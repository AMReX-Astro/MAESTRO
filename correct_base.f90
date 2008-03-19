module correct_base_module

  use bl_types

  implicit none

  private

  public :: correct_base

contains

  subroutine correct_base(nlevs,s0_old,s0_new,eta,dz,dt)

    use bl_prof_module
    use geometry, only: spherical

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: s0_old(:,0:,:)
    real(kind=dp_t), intent(inout) :: s0_new(:,0:,:)
    real(kind=dp_t), intent(in   ) :: eta(:,0:,:)
    real(kind=dp_t), intent(in   ) :: dz(:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    integer :: n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "correct_base")
    
    do n=1,nlevs
       if (spherical .eq. 1) then
!         call correct_base_state_spherical(n,s0_old(n,0:,:),s0_new(n,0:,:),eta(n,0:,:),dz(n),dt)
       else
          call correct_base_state_planar(n,s0_old(n,0:,:),s0_new(n,0:,:),eta(n,0:,:),dz(n),dt)
       end if
    enddo

    call destroy(bpt)
       
  end subroutine correct_base

  subroutine correct_base_state_planar(n,s0_old,s0_new,eta,dz,dt)

    use bl_constants_module
    use eos_module
    use variables, only: spec_comp, rho_comp, temp_comp, rhoh_comp
    use geometry, only: nr
    use probin_module, only: grav_const, anelastic_cutoff, enthalpy_pred_type

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:)
    real(kind=dp_t), intent(inout) :: s0_new(0:,:)
    real(kind=dp_t), intent(in   ) :: eta(0:,:)
    real(kind=dp_t), intent(in   ) :: dz,dt
    
    ! Local variables
    integer :: r,comp
    integer :: r_anel
    real(kind=dp_t) :: eta_avg
   
    ! This is used to zero the eta contribution above the anelastic_cutoff
    r_anel = nr(1)-1
    do r = 0,nr(1)-1
       if (s0_old(r,rho_comp) .lt. anelastic_cutoff .and. r_anel .eq. nr(1)-1) then
          r_anel = r
          exit
       end if
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHOX0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do comp = spec_comp,spec_comp+nspec-1
       do r = 0, r_anel-1
         s0_new(r,comp) = s0_new(r,comp) - dt/dz*(eta(r+1,comp) - eta(r,comp))
       end do
    enddo
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHO0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r = 0, r_anel-1
      s0_new(r,rho_comp) = s0_new(r,rho_comp) - dt/dz*(eta(r+1,rho_comp) - eta(r,rho_comp))
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAKE TEMP0 FROM P0 AND RHO0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   do r = 0,nr(n)-1
!      
!      den_eos(1)  = s0_new(r,rho_comp)
!      temp_eos(1) = s0_old(r,temp_comp)
!      p_eos(1)    = p0_new(r)
!      xn_eos(1,:) = s0_new(r,spec_comp:spec_comp+nspec-1)/s0_new(r,rho_comp)
!      
!      ! (rho,P) --> T, h
!      call eos(eos_input_rp, den_eos, temp_eos, &
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
!      s0_new(r,temp_comp) = temp_eos(1)
!   end do
    
  end subroutine correct_base_state_planar
  
end module correct_base_module
