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
  use geometry
  use ml_layout_module
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initscalardata, initscalardata_on_level

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla)

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(dm),hi(dm),ng
    integer :: i,n
    
    ng = s(1)%ng

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
                                    p0_init(n,:))
          case (3)
             if (spherical .eq. 1) then
                call bl_error("ERROR: spherical not implemented in initdata")
             else
                call bl_error("ERROR: 3d not implemented in initdata")
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
       case (2)
          call initscalardata_2d(sop(:,:,1,:),lo,hi,ng,dx,s0_init,p0_init)
       case (3)
          if (spherical .eq. 1) then
             call bl_error("ERROR: spherical not implemented in initdata")
          else
             call bl_error("ERROR: 3d not implemented in initdata")
          end if
       end select
    end do

    call multifab_fill_boundary(s)

    call multifab_physbc(s,rho_comp,dm+rho_comp,nscal,bc)

  end subroutine initscalardata_on_level

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, prob_hi, perturb_model, rho_1, rho_2, perturb_type
    use geometry, only: nr_fine

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables
    integer         :: i,j,comp
    real(kind=dp_t) :: x,y,pertheight, L_x
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    real(kind=dp_t) :: rhoX_1,rhoX_2

    ! initial the domain with the base state
    s = ZERO

    if (perturb_type .eq. 1) then

       L_x = prob_hi(1) - prob_lo(1)
       
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
             s(i,j,temp_comp) = s0_init(j,temp_comp)
             s(i,j,trac_comp:trac_comp+ntrac-1) = s0_init(j,trac_comp:trac_comp+ntrac-1)

          end do
       end do

       do j = lo(2), hi(2)
          y = (j+HALF)*dx(2)+prob_lo(2)
          do i = lo(1), hi(1)
             x = (i+HALF)*dx(1)+prob_lo(1)

             pertheight = 0.01d0*HALF*(cos(2.d0*M_PI*x/L_x)+cos(2.d0*M_PI*(L_x-x)/L_x)) &
                  + 0.5d0

             s(i,j,rho_comp) = rho_1 + ((rho_2-rho_1)/2.d0)*(1+tanh((y-pertheight)/0.005d0))

          end do
       end do


       do comp=spec_comp,spec_comp+nspec-1
          do j = lo(2), hi(2)
             y = (j+HALF)*dx(2)+prob_lo(2)
             do i = lo(1), hi(1)
                x = (i+HALF)*dx(1)+prob_lo(1)

                pertheight = 0.01d0*HALF*(cos(2.d0*M_PI*x/L_x)+cos(2.d0*M_PI*(L_x-x)/L_x)) &
                     + 0.5d0

                rhoX_1 = s0_init(0,comp)
                rhoX_2 = s0_init(nr_fine-1,comp)
                s(i,j,comp) = rhoX_1 &
                     + ((rhoX_2-rhoX_1)/2.d0)*(1+tanh((y-pertheight)/0.005d0))

             end do
          end do
       end do

    else if (perturb_type .eq. 2) then

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             s(i,j,rho_comp)  = s0_init(j,rho_comp)
             s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
             s(i,j,temp_comp) = s0_init(j,temp_comp)
             s(i,j,spec_comp:spec_comp+nspec-1) = s0_init(j,spec_comp:spec_comp+nspec-1)
             s(i,j,trac_comp:trac_comp+ntrac-1) = s0_init(j,trac_comp:trac_comp+ntrac-1)
          enddo
       enddo
       
    end if

  end subroutine initscalardata_2d

  subroutine initscalardata_3d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, perturb_model
    
    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    s = ZERO

  end subroutine initscalardata_3d

end module init_scalar_module
