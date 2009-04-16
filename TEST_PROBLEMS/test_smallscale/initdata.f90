module init_module

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
  public :: initscalardata, initscalardata_on_level, initveldata

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla)

    use geometry, only: nlevs

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(dm),hi(dm),ng,i,n

    ng = s(1)%ng

    do n=1,nlevs

       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s0_init(n,:,:))
          case (3)
             call initscalardata_3d(n,sop(:,:,:,:), lo, hi, ng, dx(n,:), s0_init(n,:,:))
          end select
       end do

    end do ! end loop over levels

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(s(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),rho_comp,dm+rho_comp,nscal,bc(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1
          
          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
          
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,dm+rho_comp,nscal, &
                                         fill_crse_input=.false.)
          
       enddo

    end if

  end subroutine initscalardata

  subroutine initscalardata_on_level(n,s,s0_init,p0_background,dx,bc)

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t), intent(in   ) :: p0_background(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    ! local
    integer                  :: ng,i
    integer                  :: lo(dm),hi(dm)
    real(kind=dp_t), pointer :: sop(:,:,:,:)

    ng = s%ng

    do i = 1, s%nboxes
       if ( multifab_remote(s,i) ) cycle
       sop => dataptr(s,i)
       lo =  lwb(get_box(s,i))
       hi =  upb(get_box(s,i))
       select case (dm)
       case (2)
          call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx, s0_init)
       case (3)
          call initscalardata_3d(n,sop(:,:,:,:), lo, hi, ng, dx, s0_init)
       end select
    end do

    call multifab_fill_boundary(s)

    call multifab_physbc(s,rho_comp,dm+rho_comp,nscal,bc)

  end subroutine initscalardata_on_level

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init)

    use probin_module, only: perturb_model

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real(kind=dp_t), intent(in   ) ::    s0_init(0:,:)

    ! Local variables
    integer :: i,j

    if (perturb_model) then
       call bl_error('perturb_model not written for initscalardata_2d')
    end if

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          s(i,j,rho_comp) = s0_init(j,rho_comp)
          s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
          s(i,j,temp_comp) = s0_init(j,temp_comp)
          s(i,j,spec_comp:spec_comp+nspec-1) = s0_init(j,spec_comp:spec_comp+nspec-1)
          s(i,j,trac_comp:trac_comp+ntrac-1) = s0_init(j,trac_comp:trac_comp+ntrac-1)
       enddo
    enddo

  end subroutine initscalardata_2d

  subroutine initscalardata_3d(n,s,lo,hi,ng,dx,s0_init)

    use probin_module, only: perturb_model

    integer, intent(in) :: n, lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real(kind=dp_t), intent(in   ) ::    s0_init(0:,:)

    !     Local variables
    integer :: i, j, k, comp

    if (perturb_model) then
       call bl_error('perturb_model not written for initscalardata_3d')
    end if

    ! initial the domain with the base state
    s = ZERO

    if (spherical .eq. 1) then

       call bl_error('Error: initdata does not handle the spherical case')

    else 
       ! initialize the scalars
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                s(i,j,k,rho_comp) = s0_init(k,rho_comp)
                s(i,j,k,rhoh_comp) = s0_init(k,rhoh_comp)
                s(i,j,k,temp_comp) = s0_init(k,temp_comp)
                s(i,j,k,spec_comp:spec_comp+nspec-1) = s0_init(k,spec_comp:spec_comp+nspec-1)
                s(i,j,k,trac_comp:trac_comp+ntrac-1) = s0_init(k,trac_comp:trac_comp+ntrac-1)
             enddo
          enddo
       enddo
    end if

  end subroutine initscalardata_3d

  subroutine initveldata(u,s0_init,p0_init,dx,bc,mla)

    use geometry, only: nlevs

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

       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
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

       loloc = prob_lo(2) +  dble(j)     *dx(dm) - flameloc
       hiloc = prob_lo(2) + (dble(j)+ONE)*dx(dm) - flameloc

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

       loloc = prob_lo(3) +  dble(k)     *dx(dm)
       hiloc = prob_lo(3) + (dble(k)+ONE)*dx(dm)

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       u(lo(1):hi(1),lo(2):hi(2),k,1:2) = 0.0d0
       u(lo(1):hi(1),lo(2):hi(2),k,3) = state1d(2)

    enddo

  end subroutine initveldata_3d

end module init_module
