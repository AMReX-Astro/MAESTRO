module correct_base_module

  use bl_types

  implicit none

  private

  public :: correct_base

contains

  subroutine correct_base(nlevs,p0_old,p0_new,s0_old,s0_new, &
                         gam1,div_coeff,eta,dz,dt)

    use bl_prof_module
    use geometry, only: spherical

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), s0_old(:,0:,:)
    real(kind=dp_t), intent(  out) :: p0_new(:,0:), s0_new(:,0:,:)
    real(kind=dp_t), intent(inout) :: gam1(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: eta(:,0:,:)
    real(kind=dp_t), intent(in   ) :: dz(:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    integer :: n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "correct_base")
    
    do n=1,nlevs
       if (spherical .eq. 0) then
          call correct_base_state_planar(n,p0_old(n,0:),p0_new(n,0:),s0_old(n,0:,:), &
                                        s0_new(n,0:,:),gam1(n,0:),eta(n,0:,:),dz(n),dt)
       end if
    enddo

    call destroy(bpt)
       
  end subroutine correct_base

  subroutine correct_base_state_planar(n,p0_old,p0_new,s0_old,s0_new,gam1,eta,dz,dt)

    use bl_constants_module
    use make_edge_state_module
    use eos_module
    use variables, only: spec_comp, rho_comp, temp_comp, rhoh_comp
    use geometry, only: nr
    use probin_module, only: grav_const, anelastic_cutoff

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(in   ) :: p0_old(0:), s0_old(0:,:)
    real(kind=dp_t), intent(  out) :: p0_new(0:), s0_new(0:,:)
    real(kind=dp_t), intent(inout) :: gam1(0:)
    real(kind=dp_t), intent(in   ) :: eta(0:,:)
    real(kind=dp_t), intent(in   ) :: dz,dt
    
    ! Local variables
    integer :: r,comp
    integer :: r_anel
    real(kind=dp_t) :: eta_avg
    
    real (kind = dp_t), allocatable :: force(:)
    real (kind = dp_t), allocatable :: edge(:)
    
    ! Cell-centered
    allocate(force(0:nr(n)-1))
    
    ! Edge-centered
    allocate(edge(0:nr(n)))
   
    ! This is used to zero the eta contribution above the anelastic_cutoff
    r_anel = nr(1)-1
    do r = 0,nr(1)-1
       if (s0_old(r,rho_comp) .lt. anelastic_cutoff .and. r_anel .eq. nr(1)-1) then
          r_anel = r
          exit
       end if
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE P0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r = 0, r_anel-1
      eta_avg = HALF * (eta(r,rho_comp)+eta(r+1,rho_comp))
      p0_new(r) = p0_new(r) + dt * eta_avg * abs(grav_const)
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
! UPDATE RHOH0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r = 0, r_anel-1
      eta_avg = HALF * (eta(r,rho_comp)+eta(r+1,rho_comp))
      s0_new(r,rhoh_comp) = s0_new(r,rhoh_comp) &
         - dt/dz * (eta(r+1,rhoh_comp) - eta(r,rhoh_comp)) &
         + dt    *  eta_avg * abs(grav_const) 
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAKE TEMP0 AND GAM1 FROM P0 AND RHO0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do r = 0,nr(n)-1
       
       den_eos(1)  = s0_new(r,rho_comp)
       temp_eos(1) = s0_old(r,temp_comp)
       p_eos(1)    = p0_new(r)
       xn_eos(1,:) = s0_new(r,spec_comp:spec_comp+nspec-1)/s0_new(r,rho_comp)
       
       ! (rho,P) --> T, h
       call eos(eos_input_rp, den_eos, temp_eos, &
                npts, nspec, &
                xn_eos, &
                p_eos, h_eos, e_eos, &
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                do_diag)
       
       s0_new(r,temp_comp) = temp_eos(1)
       gam1(r) = gam1_eos(1)
       
    end do
    
    deallocate(force,edge)
    
  end subroutine correct_base_state_planar
  
end module correct_base_module
