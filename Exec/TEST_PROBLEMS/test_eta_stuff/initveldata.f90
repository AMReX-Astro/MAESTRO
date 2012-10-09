module init_vel_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use fill_3d_module
  use eos_module
  use variables
  use network
  use geometry
  use probin_module, only: grav_const
  use ml_layout_module
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initveldata

contains

  subroutine initveldata(nlevs,u,s0_init,p0_init,dx,bc,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(dm),hi(dm),ng
    integer :: i,n
    
    ng = u(1)%ng

    do n=1,nlevs
       do i = 1, nfabs(u(n))
          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:))
          case (3) 
             call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:))
          end select
       end do
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(u(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(u(nlevs),1,1,dm,bc(nlevs))
    else
    
       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,1,dm,fill_crse_input=.false.)
       enddo
       
    end if

  end subroutine initveldata

  subroutine initveldata_2d(u,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, prob_hi

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables
    integer :: i, j, n
    real(kind=dp_t) :: x,y,r,r0,r1,r2,temp
    real(kind=dp_t) :: spx,spy,cpx,cpy,Pi

    ! initial the velocity
    u = ZERO

    Pi = 4.d0 * datan(1.d0)

    do j = lo(2), hi(2)
       if (j .ge. hi(2)/8 .and. j .le. 3*hi(2)/8) then
            y = (dble(j) - hi(2)/8 + half)*dx(2) / (prob_hi(2)-prob_lo(2))
            spy = sin(Pi*FOUR*y)
            cpy = cos(Pi*FOUR*y)
            do i = lo(1), hi(1)
               x = prob_lo(1) + (dble(i) + half)*dx(1) / (prob_hi(1)-prob_lo(1))

               spx = sin(Pi*x)
               cpx = cos(Pi*x)

               u(i,j,1) = -1.e2*two*spy*cpy*spx**2
               u(i,j,2) = +1.e2*two*spx*cpx*spy**2

            end do
       end if
    end do

  end subroutine initveldata_2d

  subroutine initveldata_3d(u,lo,hi,ng,dx,s0_init,p0_init)

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables

    ! initial the velocity
    u = ZERO

  end subroutine initveldata_3d

end module init_vel_module
