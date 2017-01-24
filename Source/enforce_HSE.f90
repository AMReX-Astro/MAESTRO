module enforce_HSE_module

  use bl_types

  implicit none

  private

  public :: enforce_HSE

contains

  subroutine enforce_HSE(rho0,p0,grav_cell)

    use geometry, only: dr, r_start_coord, r_end_coord, numdisjointchunks, &
         base_cutoff_density_coord, nr, nlevs_radial, nr_fine
    use restrict_base_module
    use bl_error_module
    use make_grav_module
    use probin_module, only: do_planar_invsq_grav, do_self_grav

    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(inout) ::   p0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell(:,0:)

    integer         :: n,l,i,r
    real(kind=dp_t) :: temp,offset
    real(kind=dp_t) :: grav_edge(nlevs_radial,0:nr_fine-1)
    real(kind=dp_t) ::    p0_old(nlevs_radial,0:nr_fine-1)

    offset = 0.d0

    call make_grav_edge(grav_edge,rho0)

    ! create a copy of the input pressure to help us with initial
    ! conditions
    p0_old = p0
    
    ! zero the new pressure so we don't leave a non-zero pressure in
    ! fine radial regions that no longer have a corresponding full
    ! state
    p0 = 0.d0
    
    ! integrate all of level 1 first
    ! use the old pressure at r=0 as a reference point
    p0(1,0) = p0_old(1,0)

    ! now integrate upwards from the bottom later, we will offset the
    ! entire pressure so we have effectively integrated from the "top"
    do r=1,min(r_end_coord(1,1),base_cutoff_density_coord(1))
       p0(1,r) = p0(1,r-1) + (dr(1)/2.d0)*(rho0(1,r)+rho0(1,r-1))*grav_edge(1,r)
    end do
    do r=base_cutoff_density_coord(1)+1,r_end_coord(1,1)
       p0(1,r) = p0(1,r-1)
    end do

    do n=2,nlevs_radial
       do i=1,numdisjointchunks(n)

          ! get pressure in the bottom cell of this disjointchunk
          if (r_start_coord(n,i) .eq. 0) then
             ! if we are at the bottom of the domain, use the old
             ! pressure as reference
             p0(n,0) = p0_old(n,0)
          else if (r_start_coord(n,i) .le. base_cutoff_density_coord(n)) then
             ! we integrate upwards starting from the nearest coarse
             ! cell at a lower physical height

             if (do_planar_invsq_grav .OR. do_self_grav) then
                ! we have variable gravity
                p0(n,r_start_coord(n,i)) = p0(n-1,r_start_coord(n,i)/2-1) &
                     + (dr(n)/4.d0)* &
                      (2.d0*rho0(n,r_start_coord(n,i))/3.d0 + &
                       4.d0*rho0(n-1,r_start_coord(n,i)/2-1)/3.d0)* &
                      (grav_edge(n,r_start_coord(n,i)) + &
                       grav_cell(n-1,r_start_coord(n,i)/2-1))  &
                     + (dr(n)/8.d0)* &
                      (5.d0*rho0(n,r_start_coord(n,i))/3.d0 + &
                       1.d0*rho0(n-1,r_start_coord(n,i)/2-1)/3.d0)* &
                      (grav_edge(n,r_start_coord(n,i)) + &
                       grav_cell(n,r_start_coord(n,i)))
             else
                ! assuming constant g here
                p0(n,r_start_coord(n,i)) = p0(n-1,r_start_coord(n,i)/2-1) &
                     + (3.d0*grav_cell(1,0)*dr(n)/4.d0)* &
                     (rho0(n-1,r_start_coord(n,i)/2-1)+rho0(n,r_start_coord(n,i)))                      
             endif
          else
             ! copy pressure from below
             p0(n,r_start_coord(n,i)) = p0(n-1,r_start_coord(n,i)/2-1)
          end if

          ! integrate upwards as normal
          do r=r_start_coord(n,i)+1,min(r_end_coord(n,i),base_cutoff_density_coord(n))
             p0(n,r) = p0(n,r-1) + (dr(n)/2.d0)*(rho0(n,r)+rho0(n,r-1))*grav_edge(n,r)
          end do
          do r=base_cutoff_density_coord(n)+1,r_end_coord(n,i)
             p0(n,r) = p0(n,r-1)
          end do

          ! now we need to look at the first coarser cell above this
          ! disjoint chunk and look at the mismatch between the
          ! integration performed over the finer cells vs. the coarse
          ! cells.  Then we need to offset the coarse cell above this
          ! point to sync up.
          
          ! first, compute the value of the pressure in the coarse
          ! cell above the disjointchunk.
          if (r_end_coord(n,i) .eq. nr(n)-1) then
             ! do nothing - we are at the top of the domain
             offset = 0.d0
          else if (r_end_coord(n,i) .le. base_cutoff_density_coord(n)) then
             ! use fine -> coarse stencil in notes
             if (do_planar_invsq_grav .OR. do_self_grav) then
                ! we have variable gravity
                temp = p0(n,r_end_coord(n,i)) &
                     + (dr(n)/4.d0)* &
                      (2.d0*rho0(n,r_end_coord(n,i))/3.d0 + &
                       4.d0*rho0(n-1,(r_end_coord(n,i)+1)/2)/3.d0)* &
                      (grav_edge(n-1,(r_end_coord(n,i)+1)/2) + &
                       grav_cell(n-1,(r_end_coord(n,i)+1)/2)) &
                     + (dr(n)/8.d0)* &
                      (5.d0*rho0(n,r_end_coord(n,i))/3.d0 + &
                       1.d0*rho0(n-1,(r_end_coord(n,i)+1)/2)/3.d0)* &
                      (grav_cell(n,r_end_coord(n,i)) + &
                       grav_edge(n-1,(r_end_coord(n,i)+1)/2))
                offset = p0(n-1,(r_end_coord(n,i)+1)/2) - temp
             else
                ! assuming constant g here
                temp = p0(n,r_end_coord(n,i)) + (3.d0*grav_cell(1,0)*dr(n)/4.d0)* &
                     (rho0(n,r_end_coord(n,i))+rho0(n-1,(r_end_coord(n,i)+1)/2))
                offset = p0(n-1,(r_end_coord(n,i)+1)/2) - temp
             endif
          else
             ! copy pressure from below
             temp = p0(n,r_end_coord(n,i))
             offset = p0(n-1,(r_end_coord(n,i)+1)/2) - temp
          end if

          ! if we are not at the top of the domain, we need to
          ! subtract the offset for all values at and above this point
          if (r_end_coord(n,i) .ne. nr(n)-1) then
             do l=n-1,1,-1
                do r=(r_end_coord(n,i)+1)/(2**(n-l)),nr(l)-1
                   p0(l,r) = p0(l,r) - offset
                end do
             end do
          end if

       end do ! end loop over disjoint chunks

    end do ! end loop over levels

    ! now compare pressure in the last cell and offset to make sure we
    ! are integrating "from the top"
    ! we use the coarsest level as the reference point
    offset = p0(1,nr(1)-1) - p0_old(1,nr(1)-1)


    ! offset level 1
    p0(1,:) = p0(1,:) - offset

    ! offset remaining levels
    do n=2,nlevs_radial
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             p0(n,r) = p0(n,r) - offset
          end do
       end do
    end do

    ! zero p0 where there is no corresponding full state array
    do n=2,nlevs_radial
       do i=1,numdisjointchunks(n)
          if (i .eq. numdisjointchunks(n)) then
             do r=r_end_coord(n,i)+1,nr(n)-1
                p0(n,r) = 0.d0
             end do
          else
             do r=r_end_coord(n,i)+1,r_start_coord(n,i+1)-1
                p0(n,r) = 0.d0
             end do
          end if
       end do
    end do

    call restrict_base(p0,.true.)
    call fill_ghost_base(p0,.true.)

  end subroutine enforce_HSE

end module enforce_HSE_module
