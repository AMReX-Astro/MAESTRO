module make_grav_module
  
  use bl_types

  implicit none

  private

  public :: make_grav_cell, make_grav_edge, make_grav_edge_uniform

contains

  subroutine make_grav_cell(grav_cell,rho0)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, r_cc_loc, r_edge_loc, &
         nlevs_radial, nr, numdisjointchunks, base_cutoff_density_coord, & 
         r_start_coord, r_end_coord, dr
    use probin_module, only: grav_const, base_cutoff_density, &
         do_planar_invsq_grav, planar_invsq_mass, do_2d_planar_octant, ref_ratio
    use fundamental_constants_module, only: Gconst
    use restrict_base_module

    ! compute the base state gravitational acceleration at the cell
    ! centers.  The base state uses 0-based indexing, so grav_cell 
    ! does too.
    
    real(kind=dp_t), intent(  out) :: grav_cell(:,0:)
    real(kind=dp_t), intent(in   ) ::      rho0(:,0:)

    ! Local variables
    integer                      :: r, n, i
    real(kind=dp_t), allocatable :: m(:,:)
    real(kind=dp_t)              :: term1, term2

    if (spherical .eq. 0) then

       if (do_planar_invsq_grav)  then

          ! we are doing a plane-parallel geometry with a 1/r**2
          ! gravitational acceleration.  The mass is assumed to be
          ! at the origin.  The mass in the computational domain
          ! does not contribute to the gravitational acceleration.
          do n=1,nlevs_radial
             do r = 0, nr(n)-1
                grav_cell(n,r) = -Gconst*planar_invsq_mass / r_cc_loc(n,r)**2
             enddo
          enddo

       else if (do_2d_planar_octant .eq. 1) then

          ! compute gravity as in the spherical case

          allocate(m(nlevs_radial,0:nr_fine-1))

          n = 1
          m(n,0) = FOUR3RD*M_PI*rho0(n,0)*r_cc_loc(n,0)**3
          grav_cell(n,0) = -Gconst * m(n,0) / r_cc_loc(n,0)**2
          
          do r=1,nr(n)-1
             
             ! the mass is defined at the cell-centers, so to compute
             ! the mass at the current center, we need to add the
             ! contribution of the upper half of the zone below us and
             ! the lower half of the current zone.
             
             ! don't add any contributions from outside the star --
             ! i.e.  rho < base_cutoff_density
             if (rho0(n,r-1) > base_cutoff_density) then
                term1 = FOUR3RD*M_PI*rho0(n,r-1) * &
                     (r_edge_loc(n,r) - r_cc_loc(n,r-1)) * &
                     (r_edge_loc(n,r)**2 + &
                     r_edge_loc(n,r)*r_cc_loc(n,r-1) + &
                     r_cc_loc(n,r-1)**2)
             else
                term1 = ZERO
             endif
             
             if (rho0(n,r) > base_cutoff_density) then
                term2 = FOUR3RD*M_PI*rho0(n,r  )*&
                     (r_cc_loc(n,r) - r_edge_loc(n,r  )) * &
                     (r_cc_loc(n,r)**2 + &
                     r_cc_loc(n,r)*r_edge_loc(n,r  ) + &
                     r_edge_loc(n,r  )**2)          
             else
                term2 = ZERO
             endif
          
             m(n,r) = m(n,r-1) + term1 + term2
          
             grav_cell(n,r) = -Gconst * m(n,r) / r_cc_loc(n,r)**2
             
          enddo

          do n = 2, nlevs_radial
             do i=1,numdisjointchunks(n)

                if (r_start_coord(n,i) .eq. 0) then
                   m(n,0) = FOUR3RD*M_PI*rho0(n,0)*r_cc_loc(n,0)**3
                   grav_cell(n,0) = -Gconst * m(n,0) / r_cc_loc(n,0)**2
                else 
                   r = r_start_coord(n,i)
                   m(n,r) = m(n-1,r/ref_ratio-1)

                   ! the mass is defined at the cell-centers, so to compute
                   ! the mass at the current center, we need to add the
                   ! contribution of the upper half of the zone below us and
                   ! the lower half of the current zone.

                   ! don't add any contributions from outside the star --
                   ! i.e.  rho < base_cutoff_density
                   if (rho0(n-1,r/ref_ratio-1) > base_cutoff_density) then
                      term1 = FOUR3RD*M_PI*rho0(n-1,r/ref_ratio-1) * &
                           (r_edge_loc(n-1,r/ref_ratio) - r_cc_loc(n-1,r/ref_ratio-1)) * &
                           (r_edge_loc(n-1,r/ref_ratio)**2 + &
                           r_edge_loc(n-1,r/ref_ratio)*r_cc_loc(n-1,r/ref_ratio-1) + &
                           r_cc_loc(n-1,r/ref_ratio-1)**2)
                   else
                      term1 = ZERO
                   endif

                   if (rho0(n,r) > base_cutoff_density) then
                      term2 = FOUR3RD*M_PI*rho0(n,r  )*&
                           (r_cc_loc(n,r) - r_edge_loc(n,r  )) * &
                           (r_cc_loc(n,r)**2 + &
                           r_cc_loc(n,r)*r_edge_loc(n,r  ) + &
                           r_edge_loc(n,r  )**2)          
                   else
                      term2 = ZERO
                   endif

                   m(n,r) = m(n,r) + term1 + term2

                   grav_cell(n,r) = -Gconst * m(n,r) / r_cc_loc(n,r)**2

                end if

                do r=r_start_coord(n,i)+1,r_end_coord(n,i)

                   ! the mass is defined at the cell-centers, so to compute
                   ! the mass at the current center, we need to add the
                   ! contribution of the upper half of the zone below us and
                   ! the lower half of the current zone.

                   ! don't add any contributions from outside the star --
                   ! i.e.  rho < base_cutoff_density
                   if (rho0(n,r-1) > base_cutoff_density) then
                      term1 = FOUR3RD*M_PI*rho0(n,r-1) * &
                           (r_edge_loc(n,r) - r_cc_loc(n,r-1)) * &
                           (r_edge_loc(n,r)**2 + &
                           r_edge_loc(n,r)*r_cc_loc(n,r-1) + &
                           r_cc_loc(n,r-1)**2)
                   else
                      term1 = ZERO
                   endif

                   if (rho0(n,r) > base_cutoff_density) then
                      term2 = FOUR3RD*M_PI*rho0(n,r  )*&
                           (r_cc_loc(n,r) - r_edge_loc(n,r  )) * &
                           (r_cc_loc(n,r)**2 + &
                           r_cc_loc(n,r)*r_edge_loc(n,r  ) + &
                           r_edge_loc(n,r  )**2)          
                   else
                      term2 = ZERO
                   endif

                   m(n,r) = m(n,r-1) + term1 + term2

                   grav_cell(n,r) = -Gconst * m(n,r) / r_cc_loc(n,r)**2

                end do
             enddo
          end do

          call restrict_base(grav_cell,.true.)
          call fill_ghost_base(grav_cell,.true.)  

       else

          ! constant gravity
          grav_cell = grav_const

       endif

    else  ! spherical = 1

       allocate(m(1,0:nr_fine-1))
          
       m(1,0) = FOUR3RD*M_PI*rho0(1,0)*r_cc_loc(1,0)**3
       grav_cell(1,0) = -Gconst * m(1,0) / r_cc_loc(1,0)**2
       
       do r=1,nr_fine-1

          ! the mass is defined at the cell-centers, so to compute
          ! the mass at the current center, we need to add the
          ! contribution of the upper half of the zone below us and
          ! the lower half of the current zone.
          
          ! don't add any contributions from outside the star --
          ! i.e.  rho < base_cutoff_density
          if (rho0(1,r-1) > base_cutoff_density) then
             term1 = FOUR3RD*M_PI*rho0(1,r-1) * &
                  (r_edge_loc(1,r) - r_cc_loc(1,r-1)) * &
                  (r_edge_loc(1,r)**2 + &
                   r_edge_loc(1,r)*r_cc_loc(1,r-1) + &
                   r_cc_loc(1,r-1)**2)
          else
             term1 = ZERO
          endif

          if (rho0(1,r) > base_cutoff_density) then
             term2 = FOUR3RD*M_PI*rho0(1,r  )*&
                  (r_cc_loc(1,r) - r_edge_loc(1,r  )) * &
                  (r_cc_loc(1,r)**2 + &
                   r_cc_loc(1,r)*r_edge_loc(1,r  ) + &
                   r_edge_loc(1,r  )**2)          
          else
             term2 = ZERO
          endif
          
          m(1,r) = m(1,r-1) + term1 + term2
          
          grav_cell(1,r) = -Gconst * m(1,r) / r_cc_loc(1,r)**2

       enddo

       deallocate(m)

    end if

  end subroutine make_grav_cell

  subroutine make_grav_edge(grav_edge,rho0)

    use bl_constants_module
    use geometry, only: spherical, r_edge_loc, nr_fine, nlevs_radial, nr, &
         numdisjointchunks, base_cutoff_density_coord,r_start_coord, &
         r_end_coord,dr
    use probin_module, only: grav_const, base_cutoff_density, &
         do_planar_invsq_grav, planar_invsq_mass, do_2d_planar_octant, ref_ratio
    use fundamental_constants_module, only: Gconst
    use restrict_base_module

    ! compute the base state gravity at the cell edges (grav_edge(1)
    ! is the gravitational acceleration at the left edge of zone 1).
    ! The base state uses 0-based indexing, so grav_edge does too.

    real(kind=dp_t), intent(  out) :: grav_edge(:,0:)
    real(kind=dp_t), intent(in   ) ::      rho0(:,0:)

    ! Local variables
    integer                      :: r, n, i
    real(kind=dp_t)              :: mencl
    real(kind=dp_t), allocatable :: m(:,:)
        
    if (spherical .eq. 0) then

       if (do_planar_invsq_grav)  then       

          ! we are doing a plane-parallel geometry with a 1/r**2
          ! gravitational acceleration.  The mass is assumed to be
          ! at the origin.  The mass in the computational domain
          ! does not contribute to the gravitational acceleration.
          do n=1,nlevs_radial
             do r = 0, nr(n)-1
                grav_edge(n,r) = -Gconst*planar_invsq_mass / r_edge_loc(n,r)**2
             enddo
          enddo

       else if (do_2d_planar_octant .eq. 1) then

          ! compute gravity as in spherical geometry

          allocate(m(nlevs_radial,0:nr_fine))

          grav_edge(1,0) = zero 
          m(1,0) = ZERO

          do r=1,nr(1)-1

             ! only add to the enclosed mass if the density is 
             ! > base_cutoff_density
             if (rho0(1,r-1) > base_cutoff_density) then
                m(1,r) = m(1,r-1) + FOUR3RD*M_PI * &
                     (r_edge_loc(1,r) - r_edge_loc(1,r-1)) * &
                     (r_edge_loc(1,r)**2 + &
                     r_edge_loc(1,r)*r_edge_loc(1,r-1) + &
                     r_edge_loc(1,r-1)**2) * rho0(1,r-1)
             else
                m(1,r) = m(1,r-1)
             endif

             grav_edge(1,r) = -Gconst * m(1,r) / r_edge_loc(1,r)**2

          end do

          do n = 2, nlevs_radial
             do i=1,numdisjointchunks(n)

                if (r_start_coord(n,i) .eq. 0) then

                   m(n,0) = ZERO

                else 

                   m(n,r_start_coord(n,i)) = m(n-1,r_start_coord(n,i)/ref_ratio)
                   grav_edge(n,r_start_coord(n,i)) = grav_edge(n-1,r_start_coord(n,i)/ref_ratio)

                end if

                do r=r_start_coord(n,i)+1,r_end_coord(n,i)+1

                   ! only add to the enclosed mass if the density is 
                   ! > base_cutoff_density
                   if (rho0(n,r-1) > base_cutoff_density) then
                      m(n,r) = m(n,r-1) + FOUR3RD*M_PI * &
                           (r_edge_loc(n,r) - r_edge_loc(n,r-1)) * &
                           (r_edge_loc(n,r)**2 + &
                           r_edge_loc(n,r)*r_edge_loc(n,r-1) + &
                           r_edge_loc(n,r-1)**2) * rho0(n,r-1)
                   else
                      m(n,r) = m(n,r-1)
                   endif

                   grav_edge(n,r) = -Gconst * m(n,r) / r_edge_loc(n,r)**2

                end do
             enddo
          end do

          deallocate(m)


          call restrict_base(grav_edge,.false.)
          call fill_ghost_base(grav_edge,.false.)

       
       else
          
          ! constant gravity
          grav_edge = grav_const

       endif

    else
       
       grav_edge(1,0) = zero 
       mencl = ZERO

       do r=1,nr_fine-1

          ! only add to the enclosed mass if the density is 
          ! > base_cutoff_density
          if (rho0(1,r-1) > base_cutoff_density) then
             mencl = mencl + FOUR3RD*M_PI * &
                  (r_edge_loc(1,r) - r_edge_loc(1,r-1)) * &
                  (r_edge_loc(1,r)**2 + &
                   r_edge_loc(1,r)*r_edge_loc(1,r-1) + &
                   r_edge_loc(1,r-1)**2) * rho0(1,r-1)
          endif
          
          grav_edge(1,r) = -Gconst * mencl / r_edge_loc(1,r)**2

       end do
       
    end if
    
  end subroutine make_grav_edge


  subroutine make_grav_edge_uniform(grav_edge_fine,rho0_fine)

    ! a special version of the make_grav_edge routine that takes a
    ! uniformly-gridded, single level density array at the finest base
    ! state resolution and returns the uniformly-gridded, single-level
    ! gravity at the same resolution

    use bl_constants_module
    use geometry, only: spherical, r_edge_loc, nr_fine, nlevs_radial, nr
    use probin_module, only: grav_const, base_cutoff_density, &
         do_planar_invsq_grav, planar_invsq_mass, do_2d_planar_octant
    use fundamental_constants_module, only: Gconst

    ! compute the base state gravity at the cell edges (grav_edge(1)
    ! is the gravitational acceleration at the left edge of zone 1).
    ! The base state uses 0-based indexing, so grav_edge does too.

    real(kind=dp_t), intent(  out) :: grav_edge_fine(0:)
    real(kind=dp_t), intent(in   ) ::      rho0_fine(0:)

    ! Local variables
    integer                      :: r
    real(kind=dp_t)              :: mencl
    
    if (spherical .eq. 0) then

       if (do_planar_invsq_grav)  then       

          ! we are doing a plane-parallel geometry with a 1/r**2
          ! gravitational acceleration.  The mass is assumed to be
          ! at the origin.  The mass in the computational domain
          ! does not contribute to the gravitational acceleration.
          do r = 0, nr(nlevs_radial)-1
             grav_edge_fine(r) = -Gconst*planar_invsq_mass / &
                  r_edge_loc(nlevs_radial,r)**2
          enddo
       
       else if (do_2d_planar_octant .eq. 1) then

          ! compute gravity as in spherical geometry
          grav_edge_fine(0) = ZERO
          mencl = ZERO

          do r=1,nr_fine-1

             ! only add to the enclosed mass if the density is 
             ! > base_cutoff_density
             if (rho0_fine(r-1) > base_cutoff_density) then
                mencl = mencl + FOUR3RD*M_PI * &
                     (r_edge_loc(nlevs_radial,r) - r_edge_loc(nlevs_radial,r-1)) * &
                     (r_edge_loc(nlevs_radial,r)**2 + &
                     r_edge_loc(nlevs_radial,r)*r_edge_loc(nlevs_radial,r-1) + &
                     r_edge_loc(nlevs_radial,r-1)**2) * rho0_fine(r-1)
             endif

             grav_edge_fine(r) = -Gconst * mencl / r_edge_loc(nlevs_radial,r)**2

          end do

       else

          ! constant gravity
          grav_edge_fine(:) = grav_const

       endif

    else
       
       grav_edge_fine(0) = ZERO
       mencl = ZERO

       do r=1,nr_fine-1

          ! only add to the enclosed mass if the density is 
          ! > base_cutoff_density
          if (rho0_fine(r-1) > base_cutoff_density) then
             mencl = mencl + FOUR3RD*M_PI * &
                  (r_edge_loc(nlevs_radial,r) - r_edge_loc(nlevs_radial,r-1)) * &
                  (r_edge_loc(nlevs_radial,r)**2 + &
                   r_edge_loc(nlevs_radial,r)*r_edge_loc(nlevs_radial,r-1) + &
                   r_edge_loc(nlevs_radial,r-1)**2) * rho0_fine(r-1)
          endif
          
          grav_edge_fine(r) = -Gconst * mencl / r_edge_loc(nlevs_radial,r)**2

       end do
       
    end if
    
  end subroutine make_grav_edge_uniform
  
end module make_grav_module
