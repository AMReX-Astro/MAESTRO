module advect_base_module

  use bl_types

  implicit none

  private

  public :: advect_base

contains

  subroutine advect_base(which_step,nlevs,vel,Sbar_in,p0_old,p0_new, &
                         s0_old,s0_new,gamma10,div_coeff,s0_predicted_edge,psi,dz,dt)

    use bl_prof_module
    use geometry, only: spherical

    integer        , intent(in   ) :: which_step,nlevs
    real(kind=dp_t), intent(in   ) :: vel(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), s0_old(:,0:,:)
    real(kind=dp_t), intent(  out) :: p0_new(:,0:), s0_new(:,0:,:)
    real(kind=dp_t), intent(inout) :: gamma10(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(  out) :: s0_predicted_edge(:,0:,:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dz(:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    integer :: n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advect_base")
    
    do n=1,nlevs
       if (spherical .eq. 0) then
          call advect_base_state_planar(which_step,n,vel(n,0:),p0_old(n,0:),p0_new(n,0:), &
                                        s0_old(n,0:,:),s0_new(n,0:,:), &
                                        s0_predicted_edge(n,0:,:),psi(n,:),dz(n),dt)
       else
          call advect_base_state_spherical(which_step,n,vel(n,:),Sbar_in(n,:,1), &
                                           p0_old(n,:),p0_new(n,:), &
                                           s0_old(n,:,:),s0_new(n,:,:), &
                                           gamma10(n,:), &
                                           s0_predicted_edge(n,0:,:),div_coeff(n,:),dt)
       end if
    enddo

    call destroy(bpt)
       
  end subroutine advect_base


  subroutine advect_base_state_planar(which_step,n,vel,p0_old,p0_new, &
                                      s0_old,s0_new,s0_predicted_edge,psi,dz,dt)

    use bl_constants_module
    use make_edge_state_module
    use eos_module
    use variables, only: spec_comp, rho_comp, temp_comp, rhoh_comp
    use geometry, only: nr
    use probin_module, only: grav_const, anelastic_cutoff, enthalpy_pred_type
    use pred_parameters

    integer        , intent(in   ) :: which_step,n
    real(kind=dp_t), intent(in   ) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), s0_old(0:,:)
    real(kind=dp_t), intent(  out) :: p0_new(0:), s0_new(0:,:)
    real(kind=dp_t), intent(  out) :: s0_predicted_edge(0:,:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    real(kind=dp_t), intent(in   ) :: dz,dt
    
    ! Local variables
    integer :: r,comp
    integer :: r_anel
    real(kind=dp_t) :: vel_avg
    
    real (kind = dp_t), allocatable :: force(:)
    real (kind = dp_t), allocatable :: edge(:)
    real (kind = dp_t), allocatable :: X0(:)
    real (kind = dp_t), allocatable :: h0(:)

    ! Cell-centered
    allocate(force(0:nr(n)-1))
    allocate(   X0(0:nr(n)-1))
    allocate(   h0(0:nr(n)-1))

    ! Edge-centered
    allocate(edge(0:nr(n)))
   
    s0_predicted_edge(:,:) = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update p_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    force = ZERO

    call make_edge_state_1d(n,p0_old,edge,vel,force,1,dz,dt)

    do r = 0, nr(n)-1
       p0_new(r) = p0_old(r) &
            - dt / dz * HALF * (vel(r) + vel(r+1)) * (edge(r+1) - edge(r)) 
            + dt * psi(j)
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Predict rho_0 to vertical edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r = 0,nr(n)-1
       force(r) = -s0_old(r,rho_comp) * (vel(r+1) - vel(r)) / dz 
    end do
    
    call make_edge_state_1d(n,s0_old(:,rho_comp),edge,vel,force,1,dz,dt)
    
    s0_predicted_edge(:,rho_comp) = edge(:)       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update (rho X)_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do comp = spec_comp,spec_comp+nspec-1

       ! here we predict X_0 on the edges
       X0(:) = s0_old(:,comp)/s0_old(:,rho_comp)
       do r = 0,nr(n)-1
          X0(r) = max(X0(r),ZERO)
       end do
       
       force = ZERO
       
       call make_edge_state_1d(n,X0,edge,vel,force,1,dz,dt)
       
       ! s0_predicted_edge will store X_0 on the vertical edges -- we need
       ! that later
       s0_predicted_edge(:,comp) = edge(:)
       
       ! our final update needs (rho X)_0 on the edges, so compute
       ! that now
       edge(:) = s0_predicted_edge(:,rho_comp)*edge(:)
       
       ! update (rho X)_0
       do r = 0,nr(n)-1
          s0_new(r,comp) = s0_old(r,comp) &
               - dt / dz * (edge(r+1) * vel(r+1) - edge(r) * vel(r)) 
       end do
       
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update rho_0 from (rho X)_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r = 0,nr(n)-1
       s0_new(r,rho_comp) =  s0_old(r,rho_comp)
       do comp = spec_comp,spec_comp+nspec-1
          s0_new(r,rho_comp) =  s0_new(r,rho_comp) + (s0_new(r,comp)-s0_old(r,comp)) 
       end do
    end do
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update (rho h)_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( (enthalpy_pred_type .eq. predict_h       ) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h) ) then

       ! here we predict h_0 on the edges
       h0(:) = s0_old(:,rhoh_comp)/s0_old(:,rho_comp)

       force = ZERO

       call make_edge_state_1d(n,h0,edge,vel,force,1,dz,dt)

       ! s0_predicted_edge will store h_0 on the vertical edges -- we need
       ! that later
       s0_predicted_edge(:,rhoh_comp) = edge(:)

       ! our final update needs (rho X)_0 on the edges, so compute
       ! that now
       edge(:) = s0_predicted_edge(:,rho_comp)*edge(:)

    else

       ! here we predict (rho h)_0 on the edges
       do r = 0,nr(n)-1
          force(r) = -s0_old(r,rhoh_comp) * (vel(r+1) - vel(r)) / dz
       end do
       
       call make_edge_state_1d(n,s0_old(:,rhoh_comp),edge,vel,force,1,dz,dt)
       
       ! s0_predicted_edge will store (rho h)_0 on the vertical edges
       s0_predicted_edge(:,rhoh_comp) = edge(:)
       
    end if

    ! update (rho h)_0
    do r = 0,nr(n)-1
       s0_new(r,rhoh_comp) = s0_old(r,rhoh_comp) &
            - dt / dz * (edge(r+1) * vel(r+1) - edge(r) * vel(r)) 
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Make t_0 from p_0 and rho_0
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

    end do
    
    deallocate(force,edge,X0,h0)
    
  end subroutine advect_base_state_planar

  
  subroutine advect_base_state_spherical(which_step,n,vel,Sbar_in, &
                                         p0_old,p0_new, &
                                         s0_old,s0_new, &
                                         gamma10, &
                                         s0_predicted_edge,div_coeff_old,dt)

    use bl_constants_module
    use make_edge_state_module
    use eos_module
    use variables, only: spec_comp, rho_comp, rhoh_comp, temp_comp
    use geometry, only: nr, base_cc_loc, base_loedge_loc, dr, nr
    use make_grav_module
    use cell_to_edge_module
    use make_div_coeff_module
    
    integer        , intent(in   ) :: which_step,n
    real(kind=dp_t), intent(in   ) :: vel(0:),Sbar_in(0:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), s0_old(0:,:)
    real(kind=dp_t), intent(  out) :: p0_new(0:), s0_new(0:,:)
    real(kind=dp_t), intent(inout) :: gamma10(0:)
    real(kind=dp_t), intent(  out) :: s0_predicted_edge(0:,:)
    real(kind=dp_t), intent(in   ) :: div_coeff_old(0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r,comp
    real(kind=dp_t) :: dtdr,divbetaw,betahalf,factor
    real(kind=dp_t) :: div_w0

    real (kind = dp_t), allocatable :: force(:)
    real (kind = dp_t), allocatable :: psi(:)
    real (kind = dp_t), allocatable :: edge(:)
    real (kind = dp_t), allocatable :: X0(:)
    real (kind = dp_t), allocatable :: div_coeff_new(:)
    real (kind = dp_t), allocatable :: beta(:),beta_new(:),beta_nh(:)
    real (kind = dp_t), allocatable :: gamma10_old(:)
    real (kind = dp_t), allocatable :: grav_cell(:)
    
    dtdr = dt / dr(n)
    
    ! Cell-centered
    allocate(force(0:nr(n)-1))
    allocate(gamma10_old(0:nr(n)-1))
    allocate(grav_cell(0:nr(n)-1))
    allocate(div_coeff_new(0:nr(n)-1))
    allocate(psi(0:nr(n)-1))
    allocate(   X0(0:nr(n)-1))
    
    ! Edge-centered
    allocate(edge(0:nr(n)))
    allocate(beta(0:nr(n)),beta_new(0:nr(n)),beta_nh(0:nr(n)))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Predict rho_0 to vertical edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r = 0,nr(n)-1
       force(r) = -s0_old(r,rho_comp) * (vel(r+1) - vel(r)) / dr(n) - &
            2.0_dp_t*s0_old(r,rho_comp)*HALF*(vel(r) + vel(r+1))/base_cc_loc(n,r)
    end do
    
    call make_edge_state_1d(n,s0_old(:,rho_comp),edge,vel,force,1,dr(n),dt)
    
    s0_predicted_edge(:,rho_comp) = edge(:)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHOX0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do comp = spec_comp,spec_comp+nspec-1
       
       ! here we predict X_0 on the edges
       X0(:) = s0_old(:,comp)/s0_old(:,rho_comp)
       do r = 0,nr(n)-1
          X0(r) = max(X0(r),ZERO)
       end do
       
       force = ZERO
       
       call make_edge_state_1d(n,X0,edge,vel,force,1,dr(n),dt)
       
       ! s0_predicted_edge will store X_0 on the vertical edges -- we need
       ! that later
       s0_predicted_edge(:,comp) = edge(:)
       
       ! our final update needs (rho X)_0 on the edges, so compute
       ! that now
       edge(:) = s0_predicted_edge(:,rho_comp)*edge(:)
       
       ! update (rho X)_0
       do r = 0,nr(n)-1
          s0_new(r,comp) = s0_old(r,comp) &
               - dtdr/base_cc_loc(n,r)**2*(base_loedge_loc(n,r+1)**2*edge(r+1)*vel(r+1) &
               - base_loedge_loc(n,r  )**2 * edge(r  ) * vel(r  ))
       end do
       
    enddo
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHO0 FROM RHOX0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do r = 0,nr(n)-1
       s0_new(r,rho_comp) =  s0_old(r,rho_comp)
       do comp = spec_comp,spec_comp+nspec-1
          s0_new(r,rho_comp) =  s0_new(r,rho_comp) + (s0_new(r,comp)-s0_old(r,comp))
       end do
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE P0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Put beta_old on edges
    call cell_to_edge(n,div_coeff_old,beta)
    
    ! Update p0 -- predictor
    do r = 0,nr(n)-1
       divbetaw = one/(base_cc_loc(n,r)**2)*(base_loedge_loc(n,r+1)**2*beta(r+1)*vel(r+1) - &
            base_loedge_loc(n,r)**2 * beta(r) * vel(r)) / dr(n)
       betahalf = div_coeff_old(r)
       factor = half * dt * gamma10(r) * (Sbar_in(r) - divbetaw / betahalf)
       p0_new(r) = p0_old(r) * (one + factor ) / (one - factor)
       
    end do
    
    do r = 0,nr(n)-1
       ! (rho, p) --> T,h, etc
       
       den_eos(1)  = s0_new(r,rho_comp)
       temp_eos(1) = s0_old(r,temp_comp) 
       p_eos(1)    = p0_new(r)
       xn_eos(1,:) = s0_new(r,spec_comp:spec_comp+nspec-1)/s0_new(r,rho_comp)
       
       gamma10_old(r) = gamma10(r)
       
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
       
       gamma10(r) = gam1_eos(1)
       s0_new(r,temp_comp) = temp_eos(1)
    end do
    
    call make_grav_cell(n,grav_cell,s0_new(:,rho_comp))
    
    ! Define beta^n+1 at cell edges using the new gravity above
    call make_div_coeff(n,div_coeff_new,s0_new(:,rho_comp),p0_new,gamma10,grav_cell)
    call cell_to_edge(n,div_coeff_new,beta_new)
    
    ! time-centered beta
    beta_nh = HALF*(beta + beta_new)
    
    ! Update p0 -- corrector
    do r = 0,nr(n)-1
       divbetaw = one / (base_cc_loc(n,r)**2) &
            * (base_loedge_loc(n,r+1)**2 * beta_nh(r+1) * vel(r+1) - &
            base_loedge_loc(n,r  )**2 * beta_nh(r  ) * vel(r  ) ) / dr(n)
       betahalf = HALF*(div_coeff_old(r) + div_coeff_new(r))
       factor = half * dt * (Sbar_in(r) - divbetaw / betahalf)
       p0_new(r) = p0_old(r) * (one + factor * gamma10_old(r)) / (one - factor * gamma10(r))
       
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHOH0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r = 0,nr(n)-1
       
       div_w0 = (vel(r+1) - vel(r)) / dr(n)
       
       force(r) = -s0_old(r,rhoh_comp) * div_w0 - &
            2.0_dp_t*s0_old(r,rhoh_comp)*HALF*(vel(r) + vel(r+1))/base_cc_loc(n,r)
       
       ! add psi at time-level n to the force for the prediction
       psi(r) = gamma10_old(r) * p0_old(r) * (Sbar_in(r) - div_w0)
       force(r) = force(r) + psi(r)
       
       ! construct a new, time-centered psi for the final update
       psi(r) = HALF*(gamma10(r)*p0_new(r) + gamma10_old(r)*p0_old(r))* &
            (Sbar_in(r) - div_w0)
    end do
    
    call make_edge_state_1d(n,s0_old(:,rhoh_comp),edge,vel,force,1,dr(n),dt)

    s0_predicted_edge(:,rhoh_comp) = edge(:)
    
    do r = 0,nr(n)-1
       
       s0_new(r,rhoh_comp) = s0_old(r,rhoh_comp) - &
            dtdr / base_cc_loc(n,r)**2 * ( base_loedge_loc(n,r+1)**2 * edge(r+1) * vel(r+1) &
            -base_loedge_loc(n,r  )**2 * edge(r  ) * vel(r  ))
       
       s0_new(r,rhoh_comp) = s0_new(r,rhoh_comp) + dt * psi(r)
       
    end do
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAKE TEMP0 FROM P0 AND RHO0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do r = 0,nr(n)-1
       
       den_eos(1)  = s0_new(r,rho_comp)
       temp_eos(1) = s0_new(r,temp_comp)
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
       
    end do
    
    deallocate(force,psi,edge,beta,beta_new,beta_nh,div_coeff_new,gamma10_old,grav_cell,X0)
    
  end subroutine advect_base_state_spherical
  
end module advect_base_module
