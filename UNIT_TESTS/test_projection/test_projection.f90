module test_projection_module

  use bl_types
  use bl_constants_module
  use multifab_module

  implicit none

  private
  public :: init_velocity

  ! width of the gaussian
  real (kind=dp_t), parameter :: W = 0.05

contains

  subroutine init_velocity(U, dx)

    integer :: n, i, ng, dm, nlevs

    type(multifab) , intent(inout) :: U(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    integer :: lo(get_dim(U(1))), hi(get_dim(U(1)))

    real(kind=dp_t), pointer :: up(:,:,:,:)

    nlevs = size(U)
    dm = get_dim(U(1))

    ng = nghost(U(1))

    do n=1,nlevs
       do i = 1, nboxes(U(n))
          if ( multifab_remote(U(n),i) ) cycle
          up => dataptr(U(n), i)
          lo = lwb(get_box(U(n), i))
          hi = upb(get_box(U(n), i))

          select case (dm)
          case (2)
             call init_velocity_2d(up(:,:,1,:), ng, lo, hi, dx(n,:))

          case (3)
             call bl_error("ERROR: init_velocity not implemented in 3d")

          end select
       end do
    end do

  end subroutine init_velocity


  subroutine init_velocity_2d(U, ng, lo, hi, dx)

    ! initialize the velocity field to a divergence-free field.  This
    ! velocity field comes from Almgren, Bell, and Szymczak 1996.

    use probin_module, only: prob_lo, prob_hi

    integer         , intent(in   ) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(inout) :: U(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j
    real (kind=dp_t) :: x, y

    do j = lo(2), hi(2)
       y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

       do i = lo(1), hi(1)
          x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)
    
          U(i,j,1) = -sin(M_PI*x)**2 * sin(TWO*M_PI*y)
          U(i,j,2) =  sin(M_PI*y)**2 * sin(TWO*M_PI*x)  

       enddo
    enddo

  end subroutine init_velocity_2d

end module test_projection_module
