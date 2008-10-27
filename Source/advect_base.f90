module advect_base_module

  use bl_types

  implicit none

  private

  public :: advect_base_dens, advect_base_pres, advect_base_enthalpy

contains

  subroutine advect_base_dens(w0,rho0_old,rho0_new,rho0_predicted_edge,dz,dt)

    use bl_prof_module
    use geometry, only: spherical, nlevs
    use restrict_base_module

    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_new(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: dz(:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    type(bl_prof_timer), save :: bpt

    call build(bpt, "advect_base")

    if (spherical .eq. 0) then
       call advect_base_dens_planar(w0,rho0_old,rho0_new,rho0_predicted_edge,dz,dt)
    else
       call advect_base_dens_spherical(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)
    end if

    call restrict_base(rho0_new,.true.)
    call fill_ghost_base(rho0_new,.true.)

    call destroy(bpt)
       
  end subroutine advect_base_dens


  subroutine advect_base_dens_planar(w0,rho0_old,rho0_new,rho0_predicted_edge,dz,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: nr_fine, r_start_coord, r_end_coord, numdisjointchunks, nlevs

    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_new(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: dz(:),dt
    
    ! Local variables
    integer :: r, n, i
    
    real (kind=dp_t) :: force(nlevs,0:nr_fine-1)
    real (kind=dp_t) ::  edge(nlevs,0:nr_fine)
   
    rho0_predicted_edge = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Predict rho_0 to vertical edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             force(n,r) = -rho0_old(n,r) * (w0(n,r+1) - w0(n,r)) / dz(n)
          end do
       end do
    end do
       
    call make_edge_state_1d(rho0_old,edge,w0,force,dz,dt)
        
    rho0_predicted_edge = edge

    ! update rho_0
    do n=1,nlevs
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             rho0_new(n,r) = rho0_old(n,r) &
                  - dt / dz(n) * (edge(n,r+1) * w0(n,r+1) - edge(n,r) * w0(n,r)) 
          end do
       end do
    end do
    
  end subroutine advect_base_dens_planar

  
  subroutine advect_base_dens_spherical(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: r_cc_loc, r_edge_loc, dr, nr_fine, nlevs
    
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_new(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer         :: r, n
    real(kind=dp_t) :: dtdr

    real(kind=dp_t) :: force(nlevs,0:nr_fine-1)
    real(kind=dp_t) ::  edge(nlevs,0:nr_fine)
    
    dtdr = dt / dr(nlevs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Predict rho_0 to vertical edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       do r=0,nr_fine-1
          force(n,r) = -rho0_old(n,r) * (w0(n,r+1) - w0(n,r)) / dr(n) - &
               2.0_dp_t*rho0_old(n,r)*HALF*(w0(n,r) + w0(n,r+1))/r_cc_loc(n,r)
       end do
    end do
    
    call make_edge_state_1d(rho0_old,edge,w0,force,dr,dt)
    
    rho0_predicted_edge = edge

    do n=1,nlevs
       do r=0,nr_fine-1
          rho0_new(n,r) = rho0_old(n,r) &
               - dtdr/r_cc_loc(n,r)**2 * &
               (r_edge_loc(n,r+1)**2 * edge(n,r+1) * w0(n,r+1) - &
               r_edge_loc(n,r  )**2 * edge(n,r  ) * w0(n,r  ))
       end do
    end do

  end subroutine advect_base_dens_spherical


  subroutine advect_base_pres(w0,Sbar_in,rho0_new,p0_old,p0_new, &
                              gamma1bar,div_coeff,psi,etarho_cc,dz,dt)

    use bl_prof_module
    use geometry, only: spherical, nlevs
    use make_psi_module
    use restrict_base_module

    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(  out) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(inout) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: etarho_cc(:,0:)
    real(kind=dp_t), intent(in   ) :: dz(:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    type(bl_prof_timer), save :: bpt

    call build(bpt, "advect_base")

    if (spherical .eq. 0) then
       call make_psi(etarho_cc,psi,w0,gamma1bar,p0_old,p0_new,Sbar_in)
       call advect_base_pres_planar(w0,p0_old,p0_new,psi,dz,dt)
       call restrict_base(p0_new,.true.)
       call fill_ghost_base(p0_new,.true.)
    else
       call advect_base_pres_spherical(w0,Sbar_in,p0_old,p0_new,rho0_new,gamma1bar, &
                                       div_coeff,dt)
       call restrict_base(p0_new,.true.)
       call fill_ghost_base(p0_new,.true.)
       call make_psi(etarho_cc,psi,w0,gamma1bar,p0_old,p0_new,Sbar_in)
    end if

    call destroy(bpt)
       
  end subroutine advect_base_pres


  subroutine advect_base_pres_planar(w0,p0_old,p0_new,psi,dz,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: nr_fine, r_start_coord, r_end_coord, numdisjointchunks, nlevs

    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(  out) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dz(:),dt
    
    ! Local variables
    integer :: r, n, i
    
    real (kind=dp_t) :: force(nlevs,0:nr_fine-1)
    real (kind=dp_t) ::  edge(nlevs,0:nr_fine)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update p_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    force = ZERO
    
    call make_edge_state_1d(p0_old,edge,w0,force,dz,dt)
    
    do n=1,nlevs
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             p0_new(n,r) = p0_old(n,r) &
                  - dt / dz(n) * HALF * (w0(n,r) + w0(n,r+1)) * (edge(n,r+1) - edge(n,r))  &
                  + dt * psi(n,r)
          end do
       end do
    end do
    
  end subroutine advect_base_pres_planar

  
  subroutine advect_base_pres_spherical(w0,Sbar_in,p0_old,p0_new,rho0_new, &
                                        gamma1bar,div_coeff_old,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: r_cc_loc, r_edge_loc, dr, nr_fine, nlevs
    use make_grav_module
    use cell_to_edge_module
    use make_div_coeff_module

    real(kind=dp_t), intent(in   ) :: w0(:,0:),Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:)
    real(kind=dp_t), intent(  out) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff_old(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r, n
    real(kind=dp_t) :: dtdr,divbetaw,betahalf,factor

    real(kind=dp_t) :: divw,p0_avg
    real(kind=dp_t) :: w0dpdr_avg,w0dpdr_avg_1,w0dpdr_avg_2

    real (kind = dp_t) :: div_coeff_new(nlevs,0:nr_fine-1)
    real (kind = dp_t) :: gamma1bar_old(nlevs,0:nr_fine-1)
    real (kind = dp_t) ::     grav_cell(nlevs,0:nr_fine-1)
    real (kind = dp_t) ::          beta(nlevs,0:nr_fine)
    real (kind = dp_t) ::      beta_new(nlevs,0:nr_fine)
    real (kind = dp_t) ::       beta_nh(nlevs,0:nr_fine)
    
    dtdr = dt / dr(nlevs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE P0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!   Put beta_old on edges
!   NOTE: Make sure ghost cells for div_coeff_old are filled before calling this
!   call cell_to_edge(n,div_coeff_old,beta)
!   Update p0 -- predictor
!   do r=0,nr_fine-1
!      divbetaw = one/(r_cc_loc(n,r)**2) * &
!           (r_edge_loc(n,r+1)**2 * beta(n,r+1) * w0(n,r+1) - &
!            r_edge_loc(n,r  )**2 * beta(n,r  ) * w0(n,r  )) / dr(n)
!      betahalf = div_coeff_old(n,r)
!      factor = half * dt * gamma1bar(n,r) * (Sbar_in(n,r) - divbetaw / betahalf)
!      p0_new(n,r) = p0_old(n,r) * (one + factor ) / (one - factor)
!   end do

    do n=1,nlevs
       do r=0,nr_fine-1
          divw = one/(r_cc_loc(n,r)**2) * &
               (r_edge_loc(n,r+1)**2 * w0(n,r+1) - &
               r_edge_loc(n,r  )**2 * w0(n,r  )) / dr(n)
          
          if (r .eq. 0) then
             w0dpdr_avg_2 =  w0(n,2) * (p0_old(n,2)-p0_old(n,1)) / dr(n)
             w0dpdr_avg_1 =  w0(n,1) * (p0_old(n,1)-p0_old(n,0)) / dr(n)
             w0dpdr_avg =  1.5d0 * w0dpdr_avg_1 - 0.5d0 * w0dpdr_avg_2
          else if (r .eq. nr_fine-1) then
             w0dpdr_avg_2 =  w0(n,nr_fine-1) * &
                  (p0_old(n,nr_fine-1)-p0_old(n,nr_fine-2)) / dr(n)
             w0dpdr_avg_1 =  w0(n,nr_fine-2) * &
                  (p0_old(n,nr_fine-2)-p0_old(n,nr_fine-3)) / dr(n)
             w0dpdr_avg =  1.5d0 * w0dpdr_avg_2 - 0.5d0 * w0dpdr_avg_1
          else
             w0dpdr_avg =  HALF * ( w0(n,r+1)*(p0_old(n,r+1)-p0_old(n,r)) &
                  + w0(n,r)*(p0_old(n,r)-p0_old(n,r-1))) / dr(n)
          end if
          
          factor = Sbar_in(n,r) - divw - 1.d0 / (gamma1bar(n,r)*p0_old(n,r)) * w0dpdr_avg
          factor = half * dt * factor
          
          p0_new(n,r) = p0_old(n,r) * (one + gamma1bar(n,r)*factor ) / &
               (one - gamma1bar(n,r)*factor)
          
       end do
    end do
    
    gamma1bar_old = gamma1bar
    
!   Define beta^n+1 at cell edges using the new gravity above
!   call make_grav_cell(n,grav_cell,rho0_new)
!   call make_div_coeff(n,div_coeff_new,rho0_new,p0_new,gamma1bar,grav_cell)
!   NOTE: Make sure ghost cells are filled for div_coeff_new before calling this
!   call cell_to_edge(n,div_coeff_new,beta_new)
!   beta_nh = HALF*(beta + beta_new)
!   Update p0 -- corrector
!   do r=0,nr_fine-1
!      divbetaw = one / (r_cc_loc(n,r)**2) * &
!           (r_edge_loc(n,r+1)**2 * beta_nh(n,r+1) * w0(n,r+1) - &
!            r_edge_loc(n,r  )**2 * beta_nh(n,r  ) * w0(n,r  )) / dr(n)

!      betahalf = HALF*(div_coeff_old(n,r) + div_coeff_new(n,r))
!      factor = half * dt * (Sbar_in(n,r) - divbetaw / betahalf)
!      p0_new(n,r) = p0_old(n,r) * &
!           (one + factor * gamma1bar_old(n,r)) / (one - factor * gamma1bar(n,r))
!   end do

    do n=1,nlevs
       do r=0,nr_fine-1
          divw = one/(r_cc_loc(n,r)**2) * &
               (r_edge_loc(n,r+1)**2 * w0(n,r+1) - &
               r_edge_loc(n,r  )**2 * w0(n,r  )) / dr(n)
          
          if (r .eq. 0) then
             w0dpdr_avg_2 =  HALF * w0(n,2) * ( (p0_old(n,2)-p0_old(n,1)) &
                  +(p0_new(n,2)-p0_new(n,1)) ) / dr(n)
             w0dpdr_avg_1 =  HALF * w0(n,1) * ( (p0_old(n,1)-p0_old(n,0)) &
                  +(p0_new(n,1)-p0_new(n,0)) ) / dr(n)
             w0dpdr_avg =  1.5d0 * w0dpdr_avg_1 - 0.5d0 * w0dpdr_avg_2
             
          else if (r .eq. nr_fine-1) then
             w0dpdr_avg_2 = HALF * w0(n,nr_fine-1) * &
                  ((p0_old(n,nr_fine-1)-p0_old(n,nr_fine-2)) &
                  +(p0_new(n,nr_fine-1)-p0_new(n,nr_fine-2))) / dr(n)
             w0dpdr_avg_1 = HALF * w0(n,nr_fine-2)* &
                  ((p0_old(n,nr_fine-2)-p0_old(n,nr_fine-3)) &
                  +(p0_new(n,nr_fine-2)-p0_new(n,nr_fine-3))) / dr(n)
             w0dpdr_avg =  1.5d0 * w0dpdr_avg_2 - 0.5d0 * w0dpdr_avg_1
             
          else
             w0dpdr_avg = HALF * HALF * ( w0(n,r+1)*(p0_old(n,r+1)-p0_old(n,r)) + &
                  w0(n,r)*(p0_old(n,r)-p0_old(n,r-1)) + &
                  w0(n,r+1)*(p0_new(n,r+1)-p0_new(n,r)) + &
                  w0(n,r)*(p0_new(n,r)-p0_new(n,r-1)) ) / dr(n)
          end if
          
          p0_avg = HALF * (p0_old(n,r) + p0_new(n,r))
          factor = Sbar_in(n,r) - divw - 1.d0 / (gamma1bar(n,r)*p0_avg) * w0dpdr_avg
          factor = half * dt * factor
          
          p0_new(n,r) = p0_old(n,r) * (one + gamma1bar(n,r)*factor ) / &
               (one - gamma1bar(n,r)*factor)
          
       end do
    end do       
       
  end subroutine advect_base_pres_spherical


  subroutine advect_base_enthalpy(w0,Sbar_in,rho0_old, &
                                  rhoh0_old,rhoh0_new,p0_old,p0_new, &
                                  gamma1bar,rho0_predicted_edge,psi,dz,dt)

    use bl_prof_module
    use geometry, only: spherical, nlevs
    use restrict_base_module

    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    real(kind=dp_t), intent(  out) :: rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dz(:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    type(bl_prof_timer), save :: bpt

    call build(bpt, "advect_base")

    if (spherical .eq. 0) then
       call advect_base_enthalpy_planar(w0,rho0_old,rhoh0_old,rhoh0_new, &
                                        rho0_predicted_edge,psi,dz,dt)
    else
       call advect_base_enthalpy_spherical(w0,Sbar_in,rho0_old,rhoh0_old,rhoh0_new, &
                                           p0_old,p0_new,gamma1bar,rho0_predicted_edge,dt)
    end if

    call restrict_base(rhoh0_new,.true.)
    call fill_ghost_base(rhoh0_new,.true.)

    call destroy(bpt)
       
  end subroutine advect_base_enthalpy


  subroutine advect_base_enthalpy_planar(w0,rho0_old,rhoh0_old,rhoh0_new, &
                                         rho0_predicted_edge,psi,dz,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: nr_fine, r_start_coord, r_end_coord, numdisjointchunks, nlevs
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:), rhoh0_old(:,0:)
    real(kind=dp_t), intent(  out) :: rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dz(:),dt
    
    ! Local variables
    integer :: r, n, i
    
    real (kind=dp_t) :: force(nlevs,0:nr_fine-1)
    real (kind=dp_t) ::    h0(nlevs,0:nr_fine-1)
    real (kind=dp_t) ::  edge(nlevs,0:nr_fine)
   
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

       call make_edge_state_1d(h0,edge,w0,force,dz,dt)

       ! our final update needs (rho h)_0 on the edges, so compute
       ! that now
       edge = rho0_predicted_edge*edge

    else

       ! here we predict (rho h)_0 on the edges
       do n=1,nlevs
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)
                force(n,r) = -rhoh0_old(n,r) * (w0(n,r+1) - w0(n,r)) / dz(n) + psi(n,r)
             end do
          end do
       end do
          
       call make_edge_state_1d(rhoh0_old,edge,w0,force,dz,dt)
              
    end if

    ! update (rho h)_0
    do n=1,nlevs
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             rhoh0_new(n,r) = rhoh0_old(n,r) &
                  - dt/dz(n) * (edge(n,r+1) * w0(n,r+1) - edge(n,r) * w0(n,r)) + dt*psi(n,r)
          end do
       end do
    end do
    
  end subroutine advect_base_enthalpy_planar

  
  subroutine advect_base_enthalpy_spherical(w0,Sbar_in,rho0_old,rhoh0_old,rhoh0_new, &
                                            p0_old,p0_new,gamma1bar,rho0_predicted_edge,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: r_cc_loc, r_edge_loc, dr, nr_fine, nlevs
    use probin_module, only: enthalpy_pred_type
    use pred_parameters
    
    real(kind=dp_t), intent(in   ) :: w0(:,0:),Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:), rhoh0_old(:,0:)
    real(kind=dp_t), intent(  out) :: rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r, n

    real(kind=dp_t) :: dtdr,divbetaw,betahalf,factor
    real(kind=dp_t) :: div_w0_cart, div_w0_sph
    real(kind=dp_t) :: divw,p0_avg
    real(kind=dp_t) :: w0dpdr_avg,w0dpdr_avg_1,w0dpdr_avg_2

    real (kind = dp_t) ::         force(nlevs,0:nr_fine-1)
    real (kind = dp_t) ::           psi(nlevs,0:nr_fine-1)
    real (kind = dp_t) ::            h0(nlevs,0:nr_fine-1)
    real (kind = dp_t) :: gamma1bar_old(nlevs,0:nr_fine-1)
    real (kind = dp_t) ::          edge(nlevs,0:nr_fine)
    real (kind = dp_t) ::          beta(nlevs,0:nr_fine)
    real (kind = dp_t) ::      beta_new(nlevs,0:nr_fine)
    real (kind = dp_t) ::       beta_nh(nlevs,0:nr_fine)
    
    dtdr = dt / dr(nlevs)

    gamma1bar_old = gamma1bar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHOH0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( (enthalpy_pred_type .eq. predict_h       ) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h) ) then

       ! here we predict h_0 on the edges
       h0 = rhoh0_old/rho0_old
       
       do n=1,nlevs
          do r=0,nr_fine-1
             
             div_w0_sph = one/(r_cc_loc(n,r)**2) * &
                  (r_edge_loc(n,r+1)**2 * w0(n,r+1) - &
                  r_edge_loc(n,r  )**2 * w0(n,r  )) / dr(n)
             
             ! add psi at time-level n to the force for the prediction
             psi(n,r) = gamma1bar_old(n,r) * p0_old(n,r) * (Sbar_in(n,r) - div_w0_sph)
             force(n,r) = psi(n,r)
             
             ! construct a new, time-centered psi for the final update
             psi(n,r) = HALF*(gamma1bar(n,r)*p0_new(n,r) + gamma1bar_old(n,r)*p0_old(n,r))* &
                  (Sbar_in(n,r) - div_w0_sph)
          end do
       end do
          
       call make_edge_state_1d(h0,edge,w0,force,dr,dt)
       
       ! our final update needs (rho h)_0 on the edges, so compute
       ! that now
       edge = rho0_predicted_edge*edge

    else

       do n=1,nlevs

          ! here we predict (rho h)_0 on the edges
          do r=0,nr_fine-1
             
             div_w0_cart = (w0(n,r+1) - w0(n,r)) / dr(n)
             div_w0_sph = one/(r_cc_loc(n,r)**2) * &
                  (r_edge_loc(n,r+1)**2 * w0(n,r+1) - &
                  r_edge_loc(n,r  )**2 * w0(n,r  )) / dr(n)
             
             force(n,r) = -rhoh0_old(n,r) * div_w0_cart - &
                  2.0_dp_t*rhoh0_old(n,r)*HALF*(w0(n,r) + w0(n,r+1))/r_cc_loc(n,r)
             
             ! add psi at time-level n to the force for the prediction
             psi(n,r) = gamma1bar_old(n,r) * p0_old(n,r) * (Sbar_in(n,r) - div_w0_sph)
             force(n,r) = force(n,r) + psi(n,r)
             
             ! construct a new, time-centered psi for the final update
             psi(n,r) = HALF*(gamma1bar(n,r)*p0_new(n,r) + gamma1bar_old(n,r)*p0_old(n,r))* &
                  (Sbar_in(n,r) - div_w0_sph)
          end do

       end do
          
       call make_edge_state_1d(rhoh0_old,edge,w0,force,dr,dt)

    endif

    do n=1,nlevs

       ! update (rho h)_0
       do r=0,nr_fine-1
          
          rhoh0_new(n,r) = rhoh0_old(n,r) - &
               dtdr / r_cc_loc(n,r)**2 * &
               (r_edge_loc(n,r+1)**2 * edge(n,r+1) * w0(n,r+1) - &
               r_edge_loc(n,r  )**2 * edge(n,r  ) * w0(n,r  ))
          
          rhoh0_new(n,r) = rhoh0_new(n,r) + dt * psi(n,r)
          
       end do

    end do
       
  end subroutine advect_base_enthalpy_spherical
  
end module advect_base_module
