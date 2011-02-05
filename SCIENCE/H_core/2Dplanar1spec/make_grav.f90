module make_grav_module
  
  use bl_types

  implicit none

  private

  public :: make_grav_cell, make_grav_edge, make_grav_edge_uniform

contains

  subroutine make_grav_cell(grav_cell,rho0)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, r_cc_loc, r_edge_loc, nlevs_radial, nr
    use probin_module, only: grav_const, base_cutoff_density, &
         do_planar_invsq_grav, planar_invsq_mass
    use fundamental_constants_module, only: Gconst

    ! compute the base state gravitational acceleration at the cell
    ! centers.  The base state uses 0-based indexing, so grav_cell 
    ! does too.
    
    real(kind=dp_t), intent(  out) :: grav_cell(:,0:)
    real(kind=dp_t), intent(in   ) ::      rho0(:,0:)

    ! Local variables
    integer                      :: r, n
    real(kind=dp_t), allocatable :: m(:)
    real(kind=dp_t)              :: term1, term2

    if (spherical .eq. 0) then

! FIXME!!!  HACK to map octant onto planar geometry
!            not made to work multilevel
       allocate(m(0:nr_fine-1))
          
       m(0) = FOUR3RD*M_PI*rho0(1,0)*r_cc_loc(1,0)**3
       grav_cell(1,0) = -Gconst * m(0) / r_cc_loc(1,0)**2
       
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
          
          m(r) = m(r-1) + term1 + term2
          
          grav_cell(1,r) = -Gconst * m(r) / r_cc_loc(1,r)**2

       enddo

       deallocate(m)

    else  ! spherical = 1

       allocate(m(0:nr_fine-1))
          
       m(0) = FOUR3RD*M_PI*rho0(1,0)*r_cc_loc(1,0)**3
       grav_cell(1,0) = -Gconst * m(0) / r_cc_loc(1,0)**2
       
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
          
          m(r) = m(r-1) + term1 + term2
          
          grav_cell(1,r) = -Gconst * m(r) / r_cc_loc(1,r)**2

       enddo

       deallocate(m)

    end if

  end subroutine make_grav_cell

  subroutine make_grav_edge(grav_edge,rho0)

    use bl_constants_module
    use geometry, only: spherical, r_edge_loc, nr_fine, nlevs_radial, nr
    use probin_module, only: grav_const, base_cutoff_density, &
         do_planar_invsq_grav, planar_invsq_mass
    use fundamental_constants_module, only: Gconst

    ! compute the base state gravity at the cell edges (grav_edge(1)
    ! is the gravitational acceleration at the left edge of zone 1).
    ! The base state uses 0-based indexing, so grav_edge does too.

    real(kind=dp_t), intent(  out) :: grav_edge(:,0:)
    real(kind=dp_t), intent(in   ) ::      rho0(:,0:)

    ! Local variables
    integer                      :: r, n
    real(kind=dp_t)              :: mencl
    
    if (spherical .eq. 0) then

! FIXME!!!  HACK to map octant onto planar geometry
!            not made to work multilevel

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
!                   * HALF*(rho0(1,r-1)+rho0(1,r))
          endif
          
          grav_edge(1,r) = -Gconst * mencl / r_edge_loc(1,r)**2

       end do
       
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
         do_planar_invsq_grav, planar_invsq_mass
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

       if (.not. do_planar_invsq_grav)  then       
          grav_edge_fine(:) = grav_const
       
       else

          ! we are doing a plane-parallel geometry with a 1/r**2
          ! gravitational acceleration.  The mass is assumed to be
          ! at the origin.  The mass in the computational domain
          ! does not contribute to the gravitational acceleration.
          do r = 0, nr(nlevs_radial)-1
             grav_edge_fine(r) = -Gconst*planar_invsq_mass / &
                  r_edge_loc(nlevs_radial,r)**2
          enddo

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
