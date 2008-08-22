module enforce_HSE_module

  use bl_types

  implicit none

  private

  public :: enforce_HSE

contains

  subroutine enforce_HSE(nlevs,rho0,p0,grav_cell)

    use geometry, only: dr, r_start_coord, r_end_coord, numdisjointchunks, spherical, &
         base_cutoff_density_coord
    use restrict_base_module, only: fill_ghost_base

    integer,         intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(inout) ::   p0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell(:,0:)

    integer         :: n,i,r
    real(kind=dp_t) :: grav

    if (spherical .eq. 0) then

       ! gravity is constant
       grav = grav_cell(1,0)

       ! do level 1 first
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
             if (r_start_coord(n,i) .le. base_cutoff_density_coord(n)) then
                p0(n,r_start_coord(n,i)) = p0(n-1,r_start_coord(n,i)/2) &
                     + (2.d0/3.d0)*(rho0(n-1,r_start_coord(n,i)/2))*grav &
                     + (1.d0/3.d0)*(rho0(n,r_start_coord(n,i)))*grav
             else
                p0(n,r_start_coord(n,i)) = p0(n-1,r_start_coord(n,i)/2)
             end if

             ! iterate normally over the rest          
             do r=r_start_coord(n,i)+1,min(r_end_coord(n,i),base_cutoff_density_coord(n))
                p0(n,r) = p0(n,r-1) + (dr(n)/2.d0)*(rho0(n,r)+rho0(n,r-1))*grav
             end do
             do r=base_cutoff_density_coord(n)+1,r_end_coord(n,i)
                p0(n,r) = p0(n,r-1)
             end do

          end do

       end do

    else

    end if

    call fill_ghost_base(nlevs,p0,.true.)

  end subroutine enforce_HSE

end module enforce_HSE_module
