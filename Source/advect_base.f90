module advect_base_module

  use bl_types

  implicit none

  private

  public :: advect_base

contains

  subroutine advect_base(nlevs,vel,Sbar_in,p0_old,p0_new,rho0_old,rho0_new, &
                         rhoh0_old,rhoh0_new, &
                         gamma1bar,div_coeff,rho0_predicted_edge,psi,dz,dt)

    use bl_prof_module
    use geometry, only: spherical

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: vel(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), rho0_old(:,0:), rhoh0_old(:,0:)
    real(kind=dp_t), intent(  out) :: p0_new(:,0:), rho0_new(:,0:), rhoh0_new(:,0:)
    real(kind=dp_t), intent(inout) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dz(:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    integer :: n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advect_base")
    
    do n=1,nlevs
       if (spherical .eq. 0) then
          call advect_base_state_planar(n,vel(n,0:),p0_old(n,0:),p0_new(n,0:), &
                                        rho0_old(n,0:),rho0_new(n,0:), &
                                        rhoh0_old(n,0:),rhoh0_new(n,0:), &
                                        rho0_predicted_edge(n,0:),psi(n,:),dz(n),dt)
       else
          call advect_base_state_spherical(n,vel(n,:),Sbar_in(n,:), &
                                           p0_old(n,:),p0_new(n,:), &
                                           rho0_old(n,:),rho0_new(n,:), &
                                           rhoh0_old(n,:),rhoh0_new(n,:), &
                                           gamma1bar(n,:), &
                                           rho0_predicted_edge(n,0:),div_coeff(n,:),dt)
       end if
    enddo

    call destroy(bpt)
       
  end subroutine advect_base


  subroutine advect_base_state_planar(n,vel,p0_old,p0_new, &
                                      rho0_old,rho0_new,rhoh0_old,rhoh0_new, &
                                      rho0_predicted_edge,psi,dz,dt)

    use bl_constants_module
    use make_edge_state_module
    use variables, only: rho_comp, rhoh_comp
    use geometry, only: r_start_coord, r_end_coord
    use probin_module, only: grav_const, enthalpy_pred_type
    use pred_parameters

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(in   ) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), rho0_old(0:), rhoh0_old(0:)
    real(kind=dp_t), intent(  out) :: p0_new(0:), rho0_new(0:), rhoh0_new(0:)
    real(kind=dp_t), intent(  out) :: rho0_predicted_edge(0:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    real(kind=dp_t), intent(in   ) :: dz,dt
    
    ! Local variables
    integer :: r
    
    real (kind = dp_t), allocatable :: force(:)
    real (kind = dp_t), allocatable :: edge(:)
    real (kind = dp_t), allocatable :: X0(:)
    real (kind = dp_t), allocatable :: h0(:)

    ! Cell-centered
    allocate(force(r_start_coord(n):r_end_coord(n)))
    allocate(   X0(r_start_coord(n):r_end_coord(n)))
    allocate(   h0(r_start_coord(n):r_end_coord(n)))

    ! Edge-centered
    allocate(edge(r_start_coord(n):r_end_coord(n)+1))
   
    rho0_predicted_edge = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update p_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    force = ZERO

    call make_edge_state_1d(n,p0_old,edge,vel,force,r_start_coord(n),dz,dt)

    do r=r_start_coord(n),r_end_coord(n)
       p0_new(r) = p0_old(r) &
            - dt / dz * HALF * (vel(r) + vel(r+1)) * (edge(r+1) - edge(r))  &
            + dt * psi(r)
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Predict rho_0 to vertical edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r=r_start_coord(n),r_end_coord(n)
       force(r) = -rho0_old(r) * (vel(r+1) - vel(r)) / dz 
    end do
    
    call make_edge_state_1d(n,rho0_old(:),edge,vel,force,r_start_coord(n),dz,dt)
    
    rho0_predicted_edge = edge

    ! update rho_0
    do r=r_start_coord(n),r_end_coord(n)
       rho0_new(r) = rho0_old(r) &
            - dt / dz * (edge(r+1) * vel(r+1) - edge(r) * vel(r)) 
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update (rho h)_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( (enthalpy_pred_type .eq. predict_h       ) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h) ) then

       ! here we predict h_0 on the edges
       h0 = rhoh0_old/rho0_old

       ! The force should be Dp_0/Dt = psi, but since in plane-parallel
       ! this is equal to eta, and therefore only depends on compositional
       ! mixing, we defer this to correct_base.
       force = ZERO

       call make_edge_state_1d(n,h0,edge,vel,force,r_start_coord(n),dz,dt)

       ! our final update needs (rho h)_0 on the edges, so compute
       ! that now
       edge = rho0_predicted_edge*edge

    else

       ! here we predict (rho h)_0 on the edges
       do r=r_start_coord(n),r_end_coord(n)
          force(r) = -rhoh0_old(r) * (vel(r+1) - vel(r)) / dz
       end do
       
       call make_edge_state_1d(n,rhoh0_old(:),edge,vel,force,r_start_coord(n),dz,dt)
       
    end if

    ! update (rho h)_0
    do r=r_start_coord(n),r_end_coord(n)
       rhoh0_new(r) = rhoh0_old(r) &
            - dt / dz * (edge(r+1) * vel(r+1) - edge(r) * vel(r)) + dt*psi(r)
    end do
    
    deallocate(force,edge,X0,h0)
    
  end subroutine advect_base_state_planar

  
  subroutine advect_base_state_spherical(n,vel,Sbar_in,p0_old,p0_new, &
                                         rho0_old,rho0_new,rhoh0_old,rhoh0_new, &
                                         gamma1bar,rho0_predicted_edge,div_coeff_old,dt)

    use bl_constants_module
    use make_edge_state_module
    use variables, only: rho_comp, rhoh_comp
    use geometry, only: r_cc_loc, r_edge_loc, dr, r_start_coord, r_end_coord
    use make_grav_module
    use cell_to_edge_module
    use make_div_coeff_module
    use probin_module, only: grav_const, enthalpy_pred_type
    use pred_parameters
    
    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(in   ) :: vel(0:),Sbar_in(0:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), rho0_old(0:), rhoh0_old(0:)
    real(kind=dp_t), intent(  out) :: p0_new(0:), rho0_new(0:), rhoh0_new(0:)
    real(kind=dp_t), intent(inout) :: gamma1bar(0:)
    real(kind=dp_t), intent(  out) :: rho0_predicted_edge(0:)
    real(kind=dp_t), intent(in   ) :: div_coeff_old(0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r
    real(kind=dp_t) :: dtdr,divbetaw,betahalf,factor
    real(kind=dp_t) :: div_w0_cart, div_w0_sph

    real(kind=dp_t) :: rho0_avg,divw,p0_avg
    real(kind=dp_t) :: w0dpdr_avg,w0dpdr_avg_1,w0dpdr_avg_2

    real (kind = dp_t), allocatable :: force(:)
    real (kind = dp_t), allocatable :: psi(:)
    real (kind = dp_t), allocatable :: edge(:)
    real (kind = dp_t), allocatable :: X0(:)
    real (kind = dp_t), allocatable :: h0(:)
    real (kind = dp_t), allocatable :: div_coeff_new(:)
    real (kind = dp_t), allocatable :: beta(:),beta_new(:),beta_nh(:)
    real (kind = dp_t), allocatable :: gamma1bar_old(:)
    real (kind = dp_t), allocatable :: grav_cell(:)
    real (kind = dp_t), allocatable :: grav_edge(:)
    
    dtdr = dt / dr(n)
    
    ! Cell-centered
    allocate(        force(r_start_coord(n):r_end_coord(n)))
    allocate(gamma1bar_old(r_start_coord(n):r_end_coord(n)))
    allocate(    grav_cell(r_start_coord(n):r_end_coord(n)))
    allocate(div_coeff_new(r_start_coord(n):r_end_coord(n)))
    allocate(          psi(r_start_coord(n):r_end_coord(n)))
    allocate(           X0(r_start_coord(n):r_end_coord(n)))
    allocate(           h0(r_start_coord(n):r_end_coord(n)))
    
    ! Edge-centered
    allocate(grav_edge(r_start_coord(n):r_end_coord(n)+1))
    allocate(     edge(r_start_coord(n):r_end_coord(n)+1))
    allocate(     beta(r_start_coord(n):r_end_coord(n)+1))
    allocate( beta_new(r_start_coord(n):r_end_coord(n)+1))
    allocate(  beta_nh(r_start_coord(n):r_end_coord(n)+1))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Predict rho_0 to vertical edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r=r_start_coord(n),r_end_coord(n)
       force(r) = -rho0_old(r) * (vel(r+1) - vel(r)) / dr(n) - &
            2.0_dp_t*rho0_old(r)*HALF*(vel(r) + vel(r+1))/r_cc_loc(n,r)
    end do
    
    call make_edge_state_1d(n,rho0_old,edge,vel,force,r_start_coord(n),dr(n),dt)
    
    rho0_predicted_edge = edge

    do r=r_start_coord(n),r_end_coord(n)
       rho0_new(r) = rho0_old(r) &
            - dtdr/r_cc_loc(n,r)**2 * &
                (r_edge_loc(n,r+1)**2 * edge(r+1) * vel(r+1) - &
                 r_edge_loc(n,r  )**2 * edge(r  ) * vel(r  ))
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE P0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!   Put beta_old on edges
!   call cell_to_edge(n,div_coeff_old,beta)
!   Update p0 -- predictor
!   do r=r_start_coord(n),r_end_coord(n)
!      divbetaw = one/(r_cc_loc(n,r)**2) * &
!           (r_edge_loc(n,r+1)**2 * beta(r+1) * vel(r+1) - &
!            r_edge_loc(n,r  )**2 * beta(r  ) * vel(r  )) / dr(n)
!      betahalf = div_coeff_old(r)
!      factor = half * dt * gamma1bar(r) * (Sbar_in(r) - divbetaw / betahalf)
!      p0_new(r) = p0_old(r) * (one + factor ) / (one - factor)
!   end do

    do r=r_start_coord(n),r_end_coord(n)
       divw = one/(r_cc_loc(n,r)**2) * &
            (r_edge_loc(n,r+1)**2 * vel(r+1) - &
             r_edge_loc(n,r  )**2 * vel(r  )) / dr(n)

       if (r .eq. 0) then
          w0dpdr_avg_2 =  vel(2) * (p0_old(2)-p0_old(1)) / dr(n)
          w0dpdr_avg_1 =  vel(1) * (p0_old(1)-p0_old(0)) / dr(n)
          w0dpdr_avg =  1.5d0 * w0dpdr_avg_1 - 0.5d0 * w0dpdr_avg_2
       else if (r .eq. r_end_coord(n)) then
          w0dpdr_avg_2 =  vel(r_end_coord(n)) * &
               (p0_old(r_end_coord(n))-p0_old(r_end_coord(n)-1)) / dr(n)
          w0dpdr_avg_1 =  vel(r_end_coord(n)-1) * &
               (p0_old(r_end_coord(n)-1)-p0_old(r_end_coord(n)-2)) / dr(n)
          w0dpdr_avg =  1.5d0 * w0dpdr_avg_2 - 0.5d0 * w0dpdr_avg_1
       else
          w0dpdr_avg =  HALF * &
               ( vel(r+1)*(p0_old(r+1)-p0_old(r)) + vel(r)*(p0_old(r)-p0_old(r-1))) / dr(n)
       end if

       factor = Sbar_in(r) - divw - 1.d0 / (gamma1bar(r)*p0_old(r)) * w0dpdr_avg
       factor = half * dt * factor

       p0_new(r) = p0_old(r) * (one + gamma1bar(r)*factor ) / (one - gamma1bar(r)*factor)

    end do
    
    gamma1bar_old = gamma1bar
    
!   Define beta^n+1 at cell edges using the new gravity above
!   call make_grav_cell(n,grav_cell,rho0_new)
!   call make_div_coeff(n,div_coeff_new,rho0_new,p0_new,gamma1bar,grav_cell)
!   call cell_to_edge(n,div_coeff_new,beta_new)
!   beta_nh = HALF*(beta + beta_new)
!   Update p0 -- corrector
!   do r=r_start_coord(n),r_end_coord(n)
!      divbetaw = one / (r_cc_loc(n,r)**2) * &
!           (r_edge_loc(n,r+1)**2 * beta_nh(r+1) * vel(r+1) - &
!            r_edge_loc(n,r  )**2 * beta_nh(r  ) * vel(r  )) / dr(n)

!      betahalf = HALF*(div_coeff_old(r) + div_coeff_new(r))
!      factor = half * dt * (Sbar_in(r) - divbetaw / betahalf)
!      p0_new(r) = p0_old(r) * &
!           (one + factor * gamma1bar_old(r)) / (one - factor * gamma1bar(r))
!   end do

    do r=r_start_coord(n),r_end_coord(n)
       divw = one/(r_cc_loc(n,r)**2) * &
            (r_edge_loc(n,r+1)**2 * vel(r+1) - &
             r_edge_loc(n,r  )**2 * vel(r  )) / dr(n)

       if (r .eq. 0) then
          w0dpdr_avg_2 =  HALF * vel(2) * ( (p0_old(2)-p0_old(1)) &
                                           +(p0_new(2)-p0_new(1)) ) / dr(n)
          w0dpdr_avg_1 =  HALF * vel(1) * ( (p0_old(1)-p0_old(0)) &
                                           +(p0_new(1)-p0_new(0)) ) / dr(n)
          w0dpdr_avg =  1.5d0 * w0dpdr_avg_1 - 0.5d0 * w0dpdr_avg_2

       else if (r .eq. r_end_coord(n)) then
          w0dpdr_avg_2 = HALF * vel(r_end_coord(n)) * &
               ((p0_old(r_end_coord(n))-p0_old(r_end_coord(n)-1)) &
               +(p0_new(r_end_coord(n))-p0_new(r_end_coord(n)-1))) / dr(n)
          w0dpdr_avg_1 = HALF * vel(r_end_coord(n)-1)* &
               ((p0_old(r_end_coord(n)-1)-p0_old(r_end_coord(n)-2)) &
               +(p0_new(r_end_coord(n)-1)-p0_new(r_end_coord(n)-2))) / dr(n)
          w0dpdr_avg =  1.5d0 * w0dpdr_avg_2 - 0.5d0 * w0dpdr_avg_1

       else
          w0dpdr_avg = HALF * HALF * &
               ( vel(r+1)*(p0_old(r+1)-p0_old(r)) + vel(r)*(p0_old(r)-p0_old(r-1)) + &
               vel(r+1)*(p0_new(r+1)-p0_new(r)) + vel(r)*(p0_new(r)-p0_new(r-1)) ) / dr(n)
       end if

       p0_avg = HALF * (p0_old(r) + p0_new(r))
       factor = Sbar_in(r) - divw - 1.d0 / (gamma1bar(r)*p0_avg) * w0dpdr_avg
       factor = half * dt * factor

       p0_new(r) = p0_old(r) * (one + gamma1bar(r)*factor ) / (one - gamma1bar(r)*factor)

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHOH0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( (enthalpy_pred_type .eq. predict_h       ) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h) ) then

       ! here we predict h_0 on the edges
       h0 = rhoh0_old/rho0_old

       do r=r_start_coord(n),r_end_coord(n)

          div_w0_sph = one/(r_cc_loc(n,r)**2) * &
               (r_edge_loc(n,r+1)**2 * vel(r+1) - &
                r_edge_loc(n,r  )**2 * vel(r  )) / dr(n)
              
          ! add psi at time-level n to the force for the prediction
          psi(r) = gamma1bar_old(r) * p0_old(r) * (Sbar_in(r) - div_w0_sph)
          force(r) = psi(r)
          
          ! construct a new, time-centered psi for the final update
          psi(r) = HALF*(gamma1bar(r)*p0_new(r) + gamma1bar_old(r)*p0_old(r))* &
               (Sbar_in(r) - div_w0_sph)
       end do

       call make_edge_state_1d(n,h0,edge,vel,force,r_start_coord(n),dr(n),dt)

       ! our final update needs (rho h)_0 on the edges, so compute
       ! that now
       edge = rho0_predicted_edge*edge

    else

       ! here we predict (rho h)_0 on the edges
       do r=r_start_coord(n),r_end_coord(n)
       
          div_w0_cart = (vel(r+1) - vel(r)) / dr(n)
          div_w0_sph = one/(r_cc_loc(n,r)**2) * &
               (r_edge_loc(n,r+1)**2 * vel(r+1) - &
                r_edge_loc(n,r  )**2 * vel(r  )) / dr(n)

          force(r) = -rhoh0_old(r) * div_w0_cart - &
               2.0_dp_t*rhoh0_old(r)*HALF*(vel(r) + vel(r+1))/r_cc_loc(n,r)
       
          ! add psi at time-level n to the force for the prediction
          psi(r) = gamma1bar_old(r) * p0_old(r) * (Sbar_in(r) - div_w0_sph)
          force(r) = force(r) + psi(r)
          
          ! construct a new, time-centered psi for the final update
          psi(r) = HALF*(gamma1bar(r)*p0_new(r) + gamma1bar_old(r)*p0_old(r))* &
               (Sbar_in(r) - div_w0_sph)
       end do
    
       call make_edge_state_1d(n,rhoh0_old,edge,vel,force,r_start_coord(n),dr(n),dt)

    endif

    ! update (rho h)_0
    do r=r_start_coord(n),r_end_coord(n)
       
       rhoh0_new(r) = rhoh0_old(r) - &
            dtdr / r_cc_loc(n,r)**2 * &
            (r_edge_loc(n,r+1)**2 * edge(r+1) * vel(r+1) - &
             r_edge_loc(n,r  )**2 * edge(r  ) * vel(r  ))
       
       rhoh0_new(r) = rhoh0_new(r) + dt * psi(r)
       
    end do
    
    deallocate(force,psi,edge,beta,beta_new,beta_nh,div_coeff_new,gamma1bar_old)
    deallocate(grav_cell,grav_edge,X0,h0)
    
  end subroutine advect_base_state_spherical
  
end module advect_base_module
