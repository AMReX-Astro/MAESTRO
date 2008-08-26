module enforce_HSE_module

  use bl_types

  implicit none

  private

  public :: enforce_HSE

contains

  subroutine enforce_HSE(nlevs,rho0,p0,grav_cell)

    use geometry, only: dr, r_start_coord, r_end_coord, numdisjointchunks, spherical, &
         base_cutoff_density_coord, nr
    use restrict_base_module, only: fill_ghost_base
    use bl_error_module

    integer,         intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(inout) ::   p0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell(:,0:)

    integer         :: n,l,i,r
    real(kind=dp_t) :: grav,temp,offset
    real(kind=dp_t) :: temppres

    if (spherical .eq. 0) then

       ! store the pressure at the top cell
       temppres = p0(1,nr(1)-1)

       ! gravity is constant
       grav = grav_cell(1,0)

       ! integrate all of level 1 first
       ! we start at r=1 since the pressure at r=0 is assumed correct
       do r=1,min(r_end_coord(1,1),base_cutoff_density_coord(1))
          p0(1,r) = p0(1,r-1) + (dr(1)/2.d0)*(rho0(1,r)+rho0(1,r-1))*grav
       end do
       do r=base_cutoff_density_coord(1)+1,r_end_coord(1,1)
          p0(1,r) = p0(1,r-1)
       end do

       do n=2,nlevs

          do i=1,numdisjointchunks(n)

             ! use a special stencil for the first point
             if (r_start_coord(n,i) .eq. 0) then
                ! do nothing - p0(n,r_start_coord(n,i) already contains the correct data

             else if (r_start_coord(n,i) .le. base_cutoff_density_coord(n)) then
                ! use coarse -> fine stencil in notes
                p0(n,r_start_coord(n,i)) = p0(n-1,r_start_coord(n,i)/2-1) &
                     + (3.d0*grav*dr(n)/4.d0)* &
                     (rho0(n-1,r_start_coord(n,i)/2-1)+rho0(n,r_start_coord(n,i)))
             else
                ! copy pressure from below
                p0(n,r_start_coord(n,i)) = p0(n-1,r_start_coord(n,i)/2-1)
             end if

             ! iterate normally over the rest          
             do r=r_start_coord(n,i)+1,min(r_end_coord(n,i),base_cutoff_density_coord(n))
                p0(n,r) = p0(n,r-1) + (dr(n)/2.d0)*(rho0(n,r)+rho0(n,r-1))*grav
             end do
             do r=base_cutoff_density_coord(n)+1,r_end_coord(n,i)
                p0(n,r) = p0(n,r-1)
             end do

             ! use a special stencil to get the value of the coarse cell above
             if (r_end_coord(n,i) .eq. nr(n)-1) then
                ! do nothing - we are at the top of the domain

             else if (r_end_coord(n,i) .le. base_cutoff_density_coord(n)) then
                ! use fine -> coarse stencil in notes
                temp = p0(n,r_end_coord(n,i)) + (3.d0*grav*dr(n)/4.d0)* &
                     (rho0(n,r_end_coord(n,i))+rho0(n-1,(r_end_coord(n,i)+1)/2))
                offset = p0(n-1,(r_end_coord(n,i)+1)/2) - temp
             else
                ! copy pressure from below
                temp = p0(n,r_end_coord(n,i))
                offset = p0(n-1,(r_end_coord(n,i)+1)/2) - temp
             end if

             ! if we are not at the top of the domain, we need to subtract the offset 
             ! for all values at and above this point
             if (r_end_coord(n,i) .ne. nr(n)-1) then
                do l=n-1,1,-1
                   do r=(r_end_coord(n,i)+1)/(2**(n-l)),nr(l)-1
                      p0(l,r) = p0(l,r) - offset
                   end do
                end do
             end if

          end do ! end loop over disjoint chunks

       end do ! end loop over levels

       ! now compare pressure in the last cell and offset to make sure we are !
       ! integrating "from the top"
       offset = p0(1,nr(1)-1) - temppres
       
       ! offset level 1
       p0(1,:) = p0(1,:) - offset

       ! offset remaining levels
       do n=2,nlevs
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)-1
                p0(1,r) = p0(1,r) - offset
             end do
          end do
       end do

    else
       call bl_error('Have not written enforce_HSE for spherical yet')
    end if

    call fill_ghost_base(nlevs,p0,.true.)

  end subroutine enforce_HSE

end module enforce_HSE_module
