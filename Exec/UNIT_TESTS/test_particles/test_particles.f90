module test_particles_module

  use bl_types
  use bl_constants_module
  implicit none

  private
  public :: init_umac_2d

  ! width of the gaussian
  real (kind=dp_t), parameter :: W = 0.05

contains

  subroutine init_umac_2d(umac, vmac, ng_um, s, ng_s, lo, hi, dx)

    ! initialize the velocity field.  The s array is initialized to
    ! the magnitude of the velocity.

    use probin_module, only: prob_lo, prob_hi, vel_amp

    integer        , intent(in   ) :: lo(:), hi(:), ng_um, ng_s
    real(kind=dp_t), intent(inout) :: umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(inout) :: vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng_s:, lo(2)-ng_s: )
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j
    real (kind=dp_t) :: x, y
    real (kind=dp_t) :: xc, yc
    real (kind=dp_t) :: dist


    xc = HALF*(prob_hi(1) - prob_lo(1))
    yc = HALF*(prob_hi(2) - prob_lo(2))

    do j = lo(2), hi(2)+1
       ! coordinate of lower edge
       y = (dble(j))*dx(2) + prob_lo(2)

       do i = lo(1), hi(1)+1
          ! coordinate of lower edge
          x = (dble(i))*dx(1) + prob_lo(1)
    
          dist = sqrt((x - xc)**2 + (y - yc)**2)

          ! u = -A r sin(theta)
          ! v =  A r cos(theta)
          !
          ! sin(theta) = (y-yc)/r
          ! cos(theta) = (x-xc)/r
          umac(i,j) = -vel_amp*(y-yc)
          vmac(i,j) =  vel_amp*(x-xc)
          
       enddo
    enddo

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          
          s(i,j) = sqrt( (HALF*(umac(i,j) + umac(i+1,j  )))**2 + &
                         (HALF*(vmac(i,j) + vmac(i  ,j+1)))**2 )

       enddo
    enddo

  end subroutine init_umac_2d

end module test_particles_module
