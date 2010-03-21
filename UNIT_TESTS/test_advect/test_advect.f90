module test_advect_module

  use bl_types
  use bl_constants_module
  implicit none

  private
  public :: init_density_2d, init_density_3d

  ! width of the gaussian
  real (kind=dp_t), parameter :: W = 0.05

contains

  subroutine init_density_2d(rho, rhoX, ng_s, lo, hi, dx)

    ! initialize the density field to a Gaussian centered on the domain.
    ! set the first composition variable to the density and the rest to 0.

    use probin_module, only: prob_lo, prob_hi, base_cutoff_density

    integer         , intent(in   ) :: lo(:), hi(:), ng_s
    real (kind=dp_t), intent(inout) ::  rho(lo(1)-ng_s:,lo(2)-ng_s:)
    real (kind=dp_t), intent(inout) :: rhoX(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j
    real (kind=dp_t) :: x, y
    real (kind=dp_t) :: xc, yc
    real (kind=dp_t) :: dist


    xc = HALF*(prob_hi(1) - prob_lo(1))
    yc = HALF*(prob_hi(2) - prob_lo(2))

    do j = lo(2), hi(2)
       y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

       do i = lo(1), hi(1)
          x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)
    
          dist = sqrt((x - xc)**2 + (y - yc)**2)

          rho(i,j) = max(exp(-dist**2/W**2),base_cutoff_density)
          
          rhoX(i,j,:) = ZERO
          rhoX(i,j,1) = rho(i,j)

       enddo
    enddo

  end subroutine init_density_2d


  subroutine init_density_3d(rho, rhoX, ng_s, lo, hi, dx)

    ! initialize the density field to a Gaussian centered on the domain.
    ! set the first composition variable to the density and the rest to 0.

    use probin_module, only: prob_lo, prob_hi, base_cutoff_density

    integer         , intent(in   ) :: lo(:), hi(:), ng_s
    real (kind=dp_t), intent(inout) ::  rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real (kind=dp_t), intent(inout) :: rhoX(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j, k
    real (kind=dp_t) :: x, y, z
    real (kind=dp_t) :: xc, yc, zc
    real (kind=dp_t) :: dist


    xc = HALF*(prob_hi(1) - prob_lo(1))
    yc = HALF*(prob_hi(2) - prob_lo(2))
    zc = HALF*(prob_hi(3) - prob_lo(3))


    do k = lo(3), hi(3)
       z = (dble(k)+0.5d0)*dx(3) + prob_lo(3)

       do j = lo(2), hi(2)
          y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)
             x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)
    
             dist = sqrt((x - xc)**2 + (y - yc)**2 + (z - zc)**2)

             rho(i,j,k) = max(exp(-dist**2/W**2),base_cutoff_density)

             rhoX(i,j,k,:) = ZERO
             rhoX(i,j,k,1) = rho(i,j,k)

          enddo
       enddo
    enddo

  end subroutine init_density_3d

end module test_advect_module
