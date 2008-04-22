module react_state_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: react_state

contains

  subroutine react_state (nlevs,mla,s_in,s_out,rho_omegadot,rho_Hext,dt,dx,the_bc_level,time)

    use variables, only: rho_comp, nscal
    use bl_prof_module
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use heating_module

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s_in(:)
    type(multifab) , intent(inout) :: s_out(:)
    type(multifab) , intent(inout) :: rho_omegadot(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:,:),time
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local
    real(kind=dp_t), pointer:: sinp(:,:,:,:)
    real(kind=dp_t), pointer:: sotp(:,:,:,:)
    real(kind=dp_t), pointer::   rp(:,:,:,:)
    real(kind=dp_t), pointer::   hp(:,:,:,:)

    integer :: lo(s_in(1)%dim),hi(s_in(1)%dim),ng,dm
    integer :: i,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "react_state")

    ng = s_in(1)%ng
    dm = s_in(1)%dim

    call get_rho_Hext(nlevs,mla,s_in,rho_Hext,dx,time)

    do n = 1, nlevs
       do i = 1, s_in(n)%nboxes
          if ( multifab_remote(s_in(n), i) ) cycle
          sinp => dataptr(s_in(n) , i)
          sotp => dataptr(s_out(n), i)
          rp => dataptr(rho_omegadot(n), i)
          hp => dataptr(rho_Hext(n), i)
          lo =  lwb(get_box(s_in(n), i))
          hi =  upb(get_box(s_in(n), i))
          select case (dm)
          case (2)
             call react_state_2d(sinp(:,:,1,:),sotp(:,:,1,:),rp(:,:,1,:), &
                                 hp(:,:,1,1),dt,lo,hi,ng)
          case (3)
             call react_state_3d(sinp(:,:,:,:),sotp(:,:,:,:),rp(:,:,:,:), &
                                 hp(:,:,:,1),dt,lo,hi,ng)
          end select
       end do
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(s_out(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s_out(nlevs),rho_comp,dm+rho_comp,nscal,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(s_out(n-1)       ,s_out(n)       ,mla%mba%rr(n-1,:))
          call ml_cc_restriction(rho_omegadot(n-1),rho_omegadot(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(rho_Hext(n-1)    ,rho_Hext(n)    ,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s_out(n),s_out(n-1), &
                                         ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         rho_comp,dm+rho_comp,nscal)
       enddo

    end if

    call destroy(bpt)

  end subroutine react_state

  subroutine react_state_2d(s_in,s_out,rho_omegadot,rho_Hext,dt,lo,hi,ng)

    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, ebin
    use probin_module, ONLY: do_burning
    use bl_constants_module, only: zero
    use eos_module

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) :: s_in (lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(  out) :: s_out(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(  out) :: rho_omegadot(lo(1):,lo(2):,:)
    real (kind = dp_t), intent(in   ) :: rho_Hext(lo(1):,lo(2):)
    real (kind = dp_t), intent(in   ) :: dt

    !     Local variables
    integer :: i, j
    real (kind = dp_t), allocatable :: x_in(:),x_out(:),rhowdot(:)
    real (kind = dp_t) :: rho,T_in,h_in,h_out

    allocate(x_in(nspec),x_out(nspec),rhowdot(nspec))

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          rho = s_in(i,j,rho_comp)
          x_in(1:nspec) = s_in(i,j,spec_comp:spec_comp+nspec-1) / rho
          h_in = s_in(i,j,rhoh_comp) / rho
          T_in = s_in(i,j,temp_comp)

          if (do_burning) then
             call burner(rho, T_in, x_in, h_in, dt, x_out, h_out, rhowdot)
          else
             x_out = x_in
             h_out = h_in
             rhowdot = ZERO
          endif

          s_out(i,j,rho_comp) = s_in(i,j,rho_comp)
          s_out(i,j,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho
          
          rho_omegadot(i,j,1:nspec) = rhowdot(1:nspec)
          
          s_out(i,j,rhoh_comp) = rho * h_out + dt * rho_Hext(i,j)

!**********************************************
! option to compute temperature and put it into s_out
! dens, enthalpy, and xmass are inputs
!
! NOTE: if you update the temperature here, then you should not add the
! reaction term explicitly to the temperature equation force returned 
! by mktempforce in scalar advance, since you will be double counting.
!
!          den_eos(1)  = s_out(i,j,rho_comp)
!          h_eos(1)    = s_out(i,j,rhoh_comp)/s_out(i,j,rho_comp)
!          xn_eos(1,:) = s_out(i,j,spec_comp:spec_comp+nspec-1)/s_out(i,j,rho_comp)
!          temp_eos(1) = T_in
!
!          call eos(eos_input_rh, den_eos, temp_eos, &
!               npts, nspec, &
!               xn_eos, &
!               p_eos, h_eos, e_eos, &
!               cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
!               dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
!               dpdX_eos, dhdX_eos, &
!               gam1_eos, cs_eos, s_eos, &
!               dsdt_eos, dsdr_eos, &
!               do_diag)
!
!          s_out(i,j,temp_comp) = temp_eos(1)
!**********************************************

          s_out(i,j,temp_comp) = s_in(i,j,temp_comp)

          s_out(i,j,trac_comp:trac_comp+ntrac-1) = &
               s_in(i,j,trac_comp:trac_comp+ntrac-1)   

       enddo
    enddo

    deallocate(x_in,x_out,rhowdot)

  end subroutine react_state_2d

  subroutine react_state_3d(s_in,s_out,rho_omegadot,rho_Hext,dt,lo,hi,ng)

    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, ebin
    use probin_module, ONLY: do_burning
    use bl_constants_module, only: zero
    use eos_module

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) :: s_in (lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) :: s_out(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) :: rho_omegadot(lo(1):,lo(2):,lo(3):,:)
    real (kind = dp_t), intent(in   ) :: rho_Hext(lo(1):,lo(2):,lo(3):)
    real (kind = dp_t), intent(in   ) :: dt

    !     Local variables
    integer :: i, j, k
    real (kind = dp_t), allocatable :: x_in(:),x_out(:),rhowdot(:)
    real (kind = dp_t) :: rho,T_in,h_in,h_out

    allocate(x_in(nspec),x_out(nspec),rhowdot(nspec))

    do k = lo(3), hi(3)
     do j = lo(2), hi(2)

       do i = lo(1), hi(1)
          rho = s_in(i,j,k,rho_comp)
          x_in = s_in(i,j,k,spec_comp:spec_comp+nspec-1) / rho
          h_in = s_in(i,j,k,rhoh_comp) / rho
          T_in = s_in(i,j,k,temp_comp)

          if (do_burning) then
             call burner(rho, T_in, x_in, h_in, dt, x_out, h_out, rhowdot)
          else
             x_out = x_in
             h_out = h_in
             rhowdot = ZERO
          endif
          
          s_out(i,j,k,rho_comp) = s_in(i,j,k,rho_comp)
          s_out(i,j,k,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho

          rho_omegadot(i,j,k,1:nspec) = rhowdot(1:nspec)

          s_out(i,j,k,rhoh_comp) = rho * h_out + dt * rho_Hext(i,j,k)

!**********************************************
! option to compute temperature and put it into s_out
! dens, enthalpy, and xmass are inputs
!
! NOTE: if you update the temperature here, then you should not add the
! reaction term explicitly to the temperature equation force returned 
! by mktempforce in scalar advance, since you will be double counting.
!
!          den_eos(1)  = s_out(i,j,k,rho_comp)
!          h_eos(1)    = s_out(i,j,k,rhoh_comp)/s_out(i,j,k,rho_comp)
!          xn_eos(1,:) = s_out(i,j,k,spec_comp:spec_comp+nspec-1)/s_out(i,j,k,rho_comp)
!          temp_eos(1) = T_in
!
!          call eos(eos_input_rh, den_eos, temp_eos, &
!               npts, nspec, &
!               xn_eos, &
!               p_eos, h_eos, e_eos, &
!               cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
!               dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
!               dpdX_eos, dhdX_eos, &
!               gam1_eos, cs_eos, s_eos, &
!               dsdt_eos, dsdr_eos, &
!               do_diag)
!
!          s_out(i,j,k,temp_comp) = temp_eos(1)
!**********************************************

          s_out(i,j,k,temp_comp) = s_in(i,j,k,temp_comp)

          s_out(i,j,k,trac_comp:trac_comp+ntrac-1) = &
               s_in(i,j,k,trac_comp:trac_comp+ntrac-1)

       enddo
     enddo
    enddo

    deallocate(x_in,x_out,rhowdot)

  end subroutine react_state_3d

end module react_state_module
