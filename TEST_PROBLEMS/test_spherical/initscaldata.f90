module init_scalar_module

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
  use ml_layout_module
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initscalardata, initscalardata_on_level

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla)

    use geometry, only: spherical

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n,nlevs,dm

    nlevs = mla%nlevel
    dm = mla%dim
    
    ng = nghost(s(1))

    do n=1,nlevs
       do i = 1, nboxes(s(n))
          if ( multifab_remote(s(n),i) ) cycle
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          
          select case (dm)
          case (2)
             call bl_error('initscalardata_2d not written')
          case (3)
             if (spherical .eq. 1) then
                call initscalardata_3d_sphr(sop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                            s0_init(1,:,:), p0_init(1,:))
             else
                call bl_error('initscalardata_3d not written')
             end if
          end select
       end do
    enddo

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
                                         bc(n-1),bc(n),rho_comp,dm+rho_comp,nscal, &
                                         fill_crse_input=.false.)

       enddo

    end if
       
  end subroutine initscalardata

  subroutine initscalardata_on_level(n,s,s0_init,p0_init,dx,bc)

    use geometry, only: spherical

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    ! local
    integer                  :: ng,i,dm
    integer                  :: lo(s%dim),hi(s%dim)
    real(kind=dp_t), pointer :: sop(:,:,:,:)

    dm = s%dim

    ng = nghost(s)

    do i = 1, nboxes(s)
       if ( multifab_remote(s,i) ) cycle
       sop => dataptr(s,i)
       lo =  lwb(get_box(s,i))
       hi =  upb(get_box(s,i))
       select case (dm)
       case (2)
          call bl_error('initscalardata_2d not written')
       case (3)
          if (spherical .eq. 1) then
             call initscalardata_3d_sphr(sop(:,:,:,:),lo,hi,ng,dx,s0_init,p0_init)
          else
             call bl_error('initscalardata_3d not written')
          end if
       end select
    end do

    call multifab_fill_boundary(s)

    call multifab_physbc(s,rho_comp,dm+rho_comp,nscal,bc)

  end subroutine initscalardata_on_level

  subroutine initscalardata_3d_sphr(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, perturb_model

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t),    intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t),    intent(in   ) :: p0_init(0:)

    !     Local variables
    integer :: comp

    if (perturb_model) then
       call bl_error('perturb_model not written for initscalardata_3d_sphr')
    end if

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0_init(:,rho_comp), &
                                      s(:,:,:,rho_comp:),lo,hi,dx,ng)

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0_init(:,rhoh_comp), &
                                      s(:,:,:,rhoh_comp:),lo,hi,dx,ng)

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0_init(:,temp_comp), &
                                      s(:,:,:,temp_comp:),lo,hi,dx,ng)

    ! initialize species
    do comp = spec_comp, spec_comp+nspec-1
       call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0_init(:,comp), &
                                         s(:,:,:,comp:),lo,hi,dx,ng)
    end do

    ! initialize tracers
    do comp = trac_comp, trac_comp+ntrac-1
       call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0_init(:,comp), &
                                         s(:,:,:,comp:),lo,hi,dx,ng)
    end do

  end subroutine initscalardata_3d_sphr

end module init_scalar_module
