! setup a simple laminar flame -- fuel and ash in pressure equilibrium

module init_module

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
          case (1)
             call initscalardata_1d(sop(:,1,1,:), lo, hi, ng, s0_init(n,:,:))
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
                                         bc(n-1),bc(n),rho_comp,dm+rho_comp,nscal, &
                                         fill_crse_input=.false.)
          
       enddo

    end if

  end subroutine initscalardata

  subroutine initscalardata_on_level(n,s,s0_init,p0_init,dx,bc)

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(0:)
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
       case (1)
          call initscalardata_1d(sop(:,1,1,:), lo, hi, ng, s0_init)
       end select
    end do

    call multifab_fill_boundary(s)

    call multifab_physbc(s,rho_comp,dm+rho_comp,nscal,bc)

  end subroutine initscalardata_on_level

  subroutine initscalardata_1d(s,lo,hi,ng,s0_init)

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,:)  
    real(kind=dp_t), intent(in   ) ::    s0_init(0:,:)

    ! Local variables
    integer :: i

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do i = lo(1), hi(1)
       s(i, rho_comp) = s0_init(i,rho_comp)
       s(i,rhoh_comp) = s0_init(i,rhoh_comp)
       s(i,temp_comp) = s0_init(i,temp_comp)
       s(i,spec_comp:spec_comp+nspec-1) = s0_init(i,spec_comp:spec_comp+nspec-1)
       s(i,trac_comp:trac_comp+ntrac-1) = s0_init(i,trac_comp:trac_comp+ntrac-1)
    enddo

  end subroutine initscalardata_1d

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
          case (1)
             call initveldata_1d(uop(:,1,1,1), lo, hi, ng)
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
    
       ! the loop over nlevs must count backwards to make sure the finer 
       ! grids are done first
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

  subroutine initveldata_1d(u,lo,hi,ng)

    use probin_module, only: vel_fuel

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:)  

    ! initialize the velocity -- we'll just set it to the inflow velocity
    ! and let the projection figure out the true field
    u(:) = vel_fuel

  end subroutine initveldata_1d

end module init_module
