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

  subroutine initscalardata(s,s0_init,p0_background,dx,bc,mla)

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_background(:,0:)
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
                                    p0_background(n,:))
          case (3)
             if (spherical .eq. 1) then
                call bl_error("ERROR: initscalardata not implemented in spherical")
             else
                call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
                                       p0_background(n,:))
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
          call initscalardata_2d(sop(:,:,1,:),lo,hi,ng,dx,s0_init,p0_background)
       case (3)
          if (spherical .eq. 1) then
             call bl_error("ERROR: initscalardata not implemented in spherical")
          else
             call initscalardata_3d(sop(:,:,:,:),lo,hi,ng,dx,s0_init,p0_background)
          end if
       end select
    end do

    call multifab_fill_boundary(s)

    call multifab_physbc(s,rho_comp,dm+rho_comp,nscal,bc)

  end subroutine initscalardata_on_level


  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init,p0_background)

    use probin_module, only: prob_lo, prob_hi, &
         dens_fuel, temp_fuel, xc12_fuel, vel_fuel, &
         temp_ash, frac

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_background(0:)

    ! Local variables
    integer         :: i,j
    real(kind=dp_t) :: x,y, xlen, ylen
    integer         :: ic12, io16, img24
    real(kind=dp_t) :: p_ambient, dens_ash, rhoh_fuel, rhoh_ash
    real(kind=dp_t) :: xn_fuel(nspec), xn_ash(nspec)

    ! species indices
    ic12   = network_species_index("carbon-12")
    io16   = network_species_index("oxygen-16")
    img24  = network_species_index("magnesium-24")

    ! length of the domain
    xlen = (prob_hi(1) - prob_lo(1))
    ylen = (prob_hi(2) - prob_lo(2))

    ! figure out the thermodynamics of the fuel and ash state

    ! fuel
    xn_fuel(:)    = ZERO
    xn_fuel(ic12) = xc12_fuel
    xn_fuel(io16) = 1.d0 - xc12_fuel

    den_eos(1)  = dens_fuel
    temp_eos(1) = temp_fuel
    xn_eos(1,:) = xn_fuel(:)

    call eos(eos_input_rt, den_eos, temp_eos, &
             npts, nspec, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

    ! note: p_ambient should be = p0_background
    p_ambient = p_eos(1)
    rhoh_fuel = dens_fuel*h_eos(1)


    ! ash
    xn_ash(:)     = ZERO
    xn_ash(ic12)  = ZERO
    xn_ash(io16)  = 1.d0 - xc12_fuel    
    xn_ash(img24) = xc12_fuel

    den_eos(1)  = dens_fuel    ! initial guess
    temp_eos(1) = temp_ash
    xn_eos(1,:) = xn_ash(:)
    p_eos(1) = p_ambient

    call eos(eos_input_tp, den_eos, temp_eos, &
             npts, nspec, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

    dens_ash = den_eos(1)
    rhoh_ash = dens_ash*h_eos(1)


    ! initial the scalars
    s = ZERO

    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)

       do i = lo(1), hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)          

          ! the flame propagates in the -y direction.  If we are more than
          ! fuel/ash division is frac through the domain
          if (y < prob_lo(2) + frac*ylen) then

             ! fuel
             s(i,j,rho_comp)  = dens_fuel
             s(i,j,rhoh_comp) = rhoh_fuel
             s(i,j,temp_comp) = temp_fuel
             s(i,j,spec_comp:spec_comp+nspec-1) = dens_fuel*xn_fuel(:)
             s(i,j,trac_comp:trac_comp+ntrac-1) = ZERO

          else

             ! ash
             s(i,j,rho_comp)  = dens_ash
             s(i,j,rhoh_comp) = rhoh_ash
             s(i,j,temp_comp) = temp_ash
             s(i,j,spec_comp:spec_comp+nspec-1) = dens_ash*xn_ash(:)
             s(i,j,trac_comp:trac_comp+ntrac-1) = ZERO

          endif
          
       enddo
    enddo
        
  end subroutine initscalardata_2d

  subroutine initscalardata_3d(s,lo,hi,ng,dx,s0_init,p0_background)

    use probin_module, only: prob_lo
    
    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_background(0:)

    call bl_error("ERROR: initscalardata_3d not implemented")

  end subroutine initscalardata_3d


  subroutine initveldata(u,s0_init,p0_background,dx,bc,mla)

    use geometry, only: nlevs

    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_background(:,0:)
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
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_background(n,:))
          case (3) 
             if (spherical .eq. 1) then
                call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                    s0_init(1,:,:), p0_background(1,:))
             else
                call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                    s0_init(n,:,:), p0_background(n,:))
             end if
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

  subroutine initveldata_2d(u,lo,hi,ng,dx,s0_init,p0_background)

    use probin_module, only: vel_fuel

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_background(0:)

    ! Local variables

    ! initialize the velocity -- we'll just set it to the inflow velocity
    ! and let the projection figure out the true field
    u(:,:,1) = ZERO
    u(:,:,2) = vel_fuel

  end subroutine initveldata_2d

  subroutine initveldata_3d(u,lo,hi,ng,dx,s0_init,p0_background)

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_background(0:)

    ! Local variables

    ! initialize the velocity
    u = ZERO
    call bl_error("ERROR: initveldata_3d not implemented")
    
  end subroutine initveldata_3d

end module init_module
