module advect_base_module

  use bl_types

  implicit none

  private

  public :: advect_base_dens, advect_base_pres, advect_base_enthalpy

contains

  subroutine advect_base_dens(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)

    use bl_prof_module
    use geometry, only: spherical, nlevs
    use restrict_base_module

    real(kind=dp_t), intent(in   ) ::                  w0(:,0:)
    real(kind=dp_t), intent(in   ) ::            rho0_old(:,0:)
    real(kind=dp_t), intent(  out) ::            rho0_new(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    type(bl_prof_timer), save :: bpt

    call build(bpt, "advect_base_dens")

    if (spherical .eq. 0) then
       call advect_base_dens_planar(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)
       call restrict_base(rho0_new,.true.)
       call fill_ghost_base(rho0_new,.true.)
    else
       call advect_base_dens_spherical(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)
    end if

    call destroy(bpt)
       
  end subroutine advect_base_dens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_dens_planar(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: nr_fine, r_start_coord, r_end_coord, numdisjointchunks, nlevs, dr

    real(kind=dp_t), intent(in   ) ::                  w0(:,0:)
    real(kind=dp_t), intent(in   ) ::            rho0_old(:,0:)
    real(kind=dp_t), intent(  out) ::            rho0_new(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
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
             force(n,r) = -rho0_old(n,r) * (w0(n,r+1) - w0(n,r)) / dr(n)
          end do
       end do
    end do
       
    call make_edge_state_1d(rho0_old,edge,w0,force,dt)
        
    rho0_predicted_edge = edge

    ! update rho_0
    do n=1,nlevs
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             rho0_new(n,r) = rho0_old(n,r) &
                  - dt / dr(n) * (edge(n,r+1) * w0(n,r+1) - edge(n,r) * w0(n,r)) 
          end do
       end do
    end do
    
  end subroutine advect_base_dens_planar

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_dens_spherical(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: r_cc_loc, r_edge_loc, dr, nr_fine
    
    real(kind=dp_t), intent(in   ) ::                  w0(:,0:)
    real(kind=dp_t), intent(in   ) ::            rho0_old(:,0:)
    real(kind=dp_t), intent(  out) ::            rho0_new(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer         :: r
    real(kind=dp_t) :: dtdr
    real(kind=dp_t) :: force(1,0:nr_fine-1)
    real(kind=dp_t) ::  edge(1,0:nr_fine)
    
    dtdr = dt / dr(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Predict rho_0 to vertical edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r=0,nr_fine-1
       force(1,r) = -rho0_old(1,r) * (w0(1,r+1) - w0(1,r)) / dr(1) - &
            2.0_dp_t*rho0_old(1,r)*HALF*(w0(1,r) + w0(1,r+1))/r_cc_loc(1,r)
    end do
    
    call make_edge_state_1d(rho0_old,edge,w0,force,dt)
    
    rho0_predicted_edge = edge

    do r=0,nr_fine-1
       rho0_new(1,r) = rho0_old(1,r) - dtdr/r_cc_loc(1,r)**2 * &
            (r_edge_loc(1,r+1)**2 * edge(1,r+1) * w0(1,r+1) - &
            r_edge_loc(1,r  )**2 * edge(1,r  ) * w0(1,r  ))
    end do
    
  end subroutine advect_base_dens_spherical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_pres(w0,Sbar_in,p0_old,p0_new,gamma1bar,psi,psi_old,etarho_cc,s, &
                              dt,dx,mla)

    use bl_prof_module
    use geometry, only: spherical, nlevs, nr_fine
    use make_psi_module
    use restrict_base_module
    use multifab_module
    use ml_layout_module

    real(kind=dp_t), intent(in   ) ::        w0(:,0:)
    real(kind=dp_t), intent(in   ) ::   Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) ::    p0_old(:,0:)
    real(kind=dp_t), intent(  out) ::    p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(inout) ::       psi(:,0:)
    real(kind=dp_t), intent(inout) ::   psi_old(:,0:)
    real(kind=dp_t), intent(in   ) :: etarho_cc(:,0:)
    type(multifab) , intent(in   ) :: s(:)
    real(kind=dp_t), intent(in   ) :: dt
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(in   ) :: mla
    
    ! local
    type(bl_prof_timer), save :: bpt

    real(kind=dp_t) :: gamma1bar_nph(1,0:nr_fine-1)
    real(kind=dp_t) :: p0_nph(1,0:nr_fine-1)

    call build(bpt, "advect_base_pres")

    if (spherical .eq. 0) then

       ! make psi
       call make_psi_planar(etarho_cc,psi)
       call fill_ghost_base(psi,.true.)
       call restrict_base(psi,.true.)

       ! advect p0
       call advect_base_pres_planar(w0,p0_old,p0_new,psi,dt)
       call restrict_base(p0_new,.true.)
       call fill_ghost_base(p0_new,.true.)

    else

       ! advect p0
       call advect_base_pres_spherical(w0,Sbar_in,p0_old,p0_nph,p0_new, &
                                       gamma1bar,gamma1bar_nph,s,dt,dx,mla)

       ! make base time and time-centered psi
       call make_psi_spherical(psi_old,w0,gamma1bar    ,p0_old,Sbar_in)
       call make_psi_spherical(psi    ,w0,gamma1bar_nph,p0_nph,Sbar_in)

    end if

    call destroy(bpt)
       
  end subroutine advect_base_pres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_pres_planar(w0,p0_old,p0_new,psi,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: nr_fine, r_start_coord, r_end_coord, numdisjointchunks, nlevs, dr

    real(kind=dp_t), intent(in   ) ::     w0(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(  out) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) ::    psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r, n, i
    
    real (kind=dp_t) :: force(nlevs,0:nr_fine-1)
    real (kind=dp_t) ::  edge(nlevs,0:nr_fine)

    force = psi
    
    call make_edge_state_1d(p0_old,edge,w0,force,dt)
    
    do n=1,nlevs
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             p0_new(n,r) = p0_old(n,r) &
                  - dt / dr(n) * HALF * (w0(n,r) + w0(n,r+1)) * (edge(n,r+1) - edge(n,r))  &
                  + dt * psi(n,r)
          end do
       end do
    end do
    
  end subroutine advect_base_pres_planar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine advect_base_pres_spherical(w0,Sbar_in,p0_old,p0_nph,p0_new, &
                                        gamma1bar,gamma1bar_nph,s,dt,dx,mla)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: r_cc_loc, r_edge_loc, dr, nr_fine, nlevs
    use make_grav_module
    use cell_to_edge_module
    use make_div_coeff_module
    use multifab_module
    use make_gamma_module
    use average_module
    use ml_layout_module

    real(kind=dp_t), intent(in   ) ::        w0(:,0:)
    real(kind=dp_t), intent(in   ) ::   Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) ::    p0_old(:,0:)
    real(kind=dp_t), intent(  out) ::    p0_nph(:,0:)
    real(kind=dp_t), intent(  out) ::    p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(  out) :: gamma1bar_nph(:,0:)
    type(multifab) , intent(in   ) :: s(:)
    real(kind=dp_t), intent(in   ) :: dt
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(in   ) :: mla
    
    ! Local variables
    integer :: r,n

    real(kind=dp_t) :: factor,divw,w0dpdr_nph,w0dpdr_nph_1,w0dpdr_nph_2
    type(multifab)  :: gamma1(nlevs)

    do r=0,nr_fine-1

       divw = one/(r_cc_loc(1,r)**2) * &
            (r_edge_loc(1,r+1)**2 * w0(1,r+1) - &
            r_edge_loc(1,r  )**2 * w0(1,r  )) / dr(1)
       
       if (r .eq. 0) then
          w0dpdr_nph_2 =  w0(1,2) * (p0_old(1,2)-p0_old(1,1)) / dr(1)
          w0dpdr_nph_1 =  w0(1,1) * (p0_old(1,1)-p0_old(1,0)) / dr(1)
          w0dpdr_nph =  1.5d0 * w0dpdr_nph_1 - 0.5d0 * w0dpdr_nph_2
       else if (r .eq. nr_fine-1) then
          w0dpdr_nph_2 =  w0(1,nr_fine-1) * (p0_old(1,nr_fine-1)-p0_old(1,nr_fine-2)) / dr(1)
          w0dpdr_nph_1 =  w0(1,nr_fine-2) * (p0_old(1,nr_fine-2)-p0_old(1,nr_fine-3)) / dr(1)
          w0dpdr_nph =  1.5d0 * w0dpdr_nph_2 - 0.5d0 * w0dpdr_nph_1
       else
          w0dpdr_nph =  HALF * ( w0(1,r+1)*(p0_old(1,r+1)-p0_old(1,r  )) &
                                +w0(1,r  )*(p0_old(1,r  )-p0_old(1,r-1)) ) / dr(1)
       end if
       
       factor = Sbar_in(1,r) - divw - 1.d0 / (gamma1bar(1,r)*p0_old(1,r)) * w0dpdr_nph
       factor = half * dt * factor
       
       p0_new(1,r) = p0_old(1,r) * (one + gamma1bar(1,r)*factor ) / &
                                   (one - gamma1bar(1,r)*factor)
       
    end do

    ! compute p0_nph
    p0_nph = HALF*(p0_old + p0_new)

    ! compute gamma1bar_star and store it in gamma1bar_nph
    do n=1,nlevs
       call multifab_build(gamma1(n), mla%la(n), 1, 0)
    end do
    
    call make_gamma(mla,gamma1,s,p0_new,dx)
    call average(mla,gamma1,gamma1bar_nph,dx,1)
    
    do n=1,nlevs
       call destroy(gamma1(n))
    end do

    ! compute gamma1bar_nph
    gamma1bar_nph = HALF*(gamma1bar + gamma1bar_nph)

    do r=0,nr_fine-1

       divw = one/(r_cc_loc(1,r)**2) * &
            (r_edge_loc(1,r+1)**2 * w0(1,r+1) - &
             r_edge_loc(1,r  )**2 * w0(1,r  )) / dr(1)
       
       if (r .eq. 0) then
          w0dpdr_nph_2 =  w0(1,2) * (p0_nph(1,2)-p0_nph(1,1)) / dr(1)
          w0dpdr_nph_1 =  w0(1,1) * (p0_nph(1,1)-p0_nph(1,0)) / dr(1)
          w0dpdr_nph =  1.5d0 * w0dpdr_nph_1 - 0.5d0 * w0dpdr_nph_2
       else if (r .eq. nr_fine-1) then
          w0dpdr_nph_2 =  w0(1,nr_fine-1) * (p0_nph(1,nr_fine-1)-p0_nph(1,nr_fine-2)) / dr(1)
          w0dpdr_nph_1 =  w0(1,nr_fine-2) * (p0_nph(1,nr_fine-2)-p0_nph(1,nr_fine-3)) / dr(1)
          w0dpdr_nph =  1.5d0 * w0dpdr_nph_2 - 0.5d0 * w0dpdr_nph_1
       else
          w0dpdr_nph =  HALF * ( w0(1,r+1)*(p0_nph(1,r+1)-p0_nph(1,r  )) &
                                +w0(1,r  )*(p0_nph(1,r  )-p0_nph(1,r-1)) ) / dr(1)
       end if
       
       factor = Sbar_in(1,r) - divw - 1.d0 / (gamma1bar_nph(1,r)*p0_nph(1,r)) * w0dpdr_nph
       factor = half * dt * factor
       
       p0_new(1,r) = p0_old(1,r) * (one + gamma1bar_nph(1,r)*factor ) / &
                                   (one - gamma1bar_nph(1,r)*factor)
       
    end do

    ! compute p0_nph
    p0_nph = HALF*(p0_old + p0_new)

    ! compute gamma1bar_star and store it in gamma1bar_nph
    do n=1,nlevs
       call multifab_build(gamma1(n), mla%la(n), 1, 0)
    end do
    
    call make_gamma(mla,gamma1,s,p0_new,dx)
    call average(mla,gamma1,gamma1bar_nph,dx,1)
    
    do n=1,nlevs
       call destroy(gamma1(n))
    end do

    ! compute gamma1bar_nph
    gamma1bar_nph = HALF*(gamma1bar + gamma1bar_nph)
       
  end subroutine advect_base_pres_spherical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_enthalpy(w0,rho0_old,rhoh0_old,rhoh0_new,rho0_predicted_edge, &
                                  psi,psi_old,dt)

    use bl_prof_module
    use geometry, only: spherical, nlevs
    use restrict_base_module
    use multifab_module

    real(kind=dp_t), intent(in   ) ::                  w0(:,0:)
    real(kind=dp_t), intent(in   ) ::            rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::           rhoh0_old(:,0:)
    real(kind=dp_t), intent(  out) ::           rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) ::                 psi(:,0:)
    real(kind=dp_t), intent(in   ) ::             psi_old(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    type(bl_prof_timer), save :: bpt

    call build(bpt, "advect_base_enthalpy")

    if (spherical .eq. 0) then
       call advect_base_enthalpy_planar(w0,rho0_old,rhoh0_old,rhoh0_new, &
                                        rho0_predicted_edge,psi,dt)
       call restrict_base(rhoh0_new,.true.)
       call fill_ghost_base(rhoh0_new,.true.)
    else
       call advect_base_enthalpy_spherical(w0,rho0_old,rhoh0_old,rhoh0_new, &
                                           rho0_predicted_edge,psi,psi_old,dt)
    end if

    call destroy(bpt)
       
  end subroutine advect_base_enthalpy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_enthalpy_planar(w0,rho0_old,rhoh0_old,rhoh0_new, &
                                         rho0_predicted_edge,psi,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: nr_fine, r_start_coord, r_end_coord, numdisjointchunks, nlevs, dr
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    real(kind=dp_t), intent(in   ) ::                  w0(:,0:)
    real(kind=dp_t), intent(in   ) ::            rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::           rhoh0_old(:,0:)
    real(kind=dp_t), intent(  out) ::           rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) ::                 psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r, i, n
    
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

       call make_edge_state_1d(h0,edge,w0,force,dt)

       ! our final update needs (rho h)_0 on the edges, so compute
       ! that now
       edge = rho0_predicted_edge*edge

    else

       ! here we predict (rho h)_0 on the edges
       do n=1,nlevs
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)
                force(n,r) = -rhoh0_old(n,r) * (w0(n,r+1) - w0(n,r)) / dr(n) + psi(n,r)
             end do
          end do
       end do
          
       call make_edge_state_1d(rhoh0_old,edge,w0,force,dt)
              
    end if

    ! update (rho h)_0
    do n=1,nlevs
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             rhoh0_new(n,r) = rhoh0_old(n,r) &
                  - dt/dr(n) * (edge(n,r+1) * w0(n,r+1) - edge(n,r) * w0(n,r)) + dt*psi(n,r)
          end do
       end do
    end do
    
  end subroutine advect_base_enthalpy_planar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine advect_base_enthalpy_spherical(w0,rho0_old,rhoh0_old,rhoh0_new, &
                                            rho0_predicted_edge,psi,psi_old,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: r_cc_loc, r_edge_loc, dr, nr_fine
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:), rhoh0_old(:,0:)
    real(kind=dp_t), intent(  out) :: rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:),psi_old(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r

    real(kind=dp_t) :: dtdr
    real(kind=dp_t) :: div_w0_cart

    real (kind = dp_t) :: force(1,0:nr_fine-1)
    real (kind = dp_t) ::    h0(1,0:nr_fine-1)
    real (kind = dp_t) ::  edge(1,0:nr_fine)
    
    dtdr = dt / dr(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHOH0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if ( (enthalpy_pred_type .eq. predict_h       ) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h) ) then

       ! here we predict h_0 on the edges
       h0 = rhoh0_old/rho0_old

       do r=0,nr_fine-1
          force(1,r) = psi(1,r)
       end do

       call make_edge_state_1d(h0,edge,w0,force,dt)

       ! our final update needs (rho h)_0 on the edges, so compute
       ! that now
       edge = rho0_predicted_edge*edge

    else

       ! here we predict (rho h)_0 on the edges
       do r=0,nr_fine-1

          div_w0_cart = (w0(1,r+1) - w0(1,r)) / dr(1)

          force(1,r) = -rhoh0_old(1,r) * div_w0_cart - &
               2.0_dp_t*rhoh0_old(1,r)*HALF*(w0(1,r) + w0(1,r+1))/r_cc_loc(1,r)

          ! add psi at time-level n to the force for the prediction
          force(1,r) = force(1,r) + psi_old(1,r)
       end do

       call make_edge_state_1d(rhoh0_old,edge,w0,force,dt)

    endif

    ! update (rho h)_0
    do r=0,nr_fine-1

       rhoh0_new(1,r) = rhoh0_old(1,r) - dtdr / r_cc_loc(1,r)**2 * &
            (r_edge_loc(1,r+1)**2 * edge(1,r+1) * w0(1,r+1) - &
            r_edge_loc(1,r  )**2 * edge(1,r  ) * w0(1,r  ))

       rhoh0_new(1,r) = rhoh0_new(1,r) + dt * psi(1,r)

    end do

  end subroutine advect_base_enthalpy_spherical

end module advect_base_module
