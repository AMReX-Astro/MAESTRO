module init_vel_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use variables
  use network
  use geometry
  use ml_layout_module
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initveldata

contains

  subroutine initveldata(u,s0_init,p0_init,dx,bc,mla)

    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(get_dim(u(1))),hi(get_dim(u(1))),ng,dm,nlevs
    integer :: i,n

    dm = get_dim(u(1))
    nlevs = size(u)

    ng = u(1)%ng

    do n=1,nlevs

       do i = 1, nfabs(u(n))
          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), s0_init(n,:,:))
          case (3)
             call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), s0_init(n,:,:))
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

  subroutine initveldata_2d(u,lo,hi,ng,dx,s0_init)

    use probin_module, only: prob_lo

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real(kind=dp_t), intent(in   ) ::    s0_init(0:,:)

    ! local
    integer ndum, j
    parameter (ndum = 31)

    character(len=128) :: lamsolfile
    real(kind=dp_t) :: state1d(ndum)
    real(kind=dp_t) :: loloc,hiloc,flameloc

    lamsolfile = 'flame_4.e7_screen_left.out'

    flameloc = ONE

    do j=lo(2),hi(2)

       loloc = prob_lo(2) +  dble(j)     *dx(2) - flameloc
       hiloc = prob_lo(2) + (dble(j)+ONE)*dx(2) - flameloc

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       u(lo(1):hi(1),j,1) = 0.0d0
       u(lo(1):hi(1),j,2) = state1d(2)

    enddo

  end subroutine initveldata_2d

  subroutine initveldata_3d(u,lo,hi,ng,dx,s0_init)

    use probin_module, only: prob_lo

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real(kind=dp_t), intent(in   ) ::    s0_init(0:,:)

    ! local
    integer ndum, k
    parameter (ndum = 31)

    character(len=128) :: lamsolfile
    real(kind=dp_t) :: state1d(ndum)
    real(kind=dp_t) :: loloc,hiloc

    lamsolfile = 'flame_4.e7_screen_left.out'

    do k=lo(3),hi(3)

       loloc = prob_lo(3) +  dble(k)     *dx(3)
       hiloc = prob_lo(3) + (dble(k)+ONE)*dx(3)

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       u(lo(1):hi(1),lo(2):hi(2),k,1:2) = 0.0d0
       u(lo(1):hi(1),lo(2):hi(2),k,3) = state1d(2)

    enddo

  end subroutine initveldata_3d

end module init_vel_module
