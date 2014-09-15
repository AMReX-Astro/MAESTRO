module init_vel_module

  use multifab_module
  use ml_layout_module
  use bl_constants_module
  use ml_restrict_fill_module

  implicit none

  private
  public :: initveldata

contains

  subroutine initveldata(u,s0_init,p0_init,dx,the_bc_level,mla)

    use geometry, only: spherical

    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n,dm,nlevs
    
    dm = mla%dim
    nlevs = mla%nlevel

    ng = nghost(u(1))

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
             if (spherical .eq. 1) then
                call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                    s0_init(1,:,:), p0_init(1,:))
             else
                call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                    s0_init(n,:,:), p0_init(n,:))
             end if
          end select
       end do

    enddo

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,u,mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=1, &
                              nc=dm, &
                              ng=u(1)%ng)

  end subroutine initveldata

  subroutine initveldata_2d(u,lo,hi,ng,dx,s0_init,p0_init)

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables

    ! initial the velocity
    u = ZERO

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
