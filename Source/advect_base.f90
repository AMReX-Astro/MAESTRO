module advect_base_module

  use bl_types, only: dp_t

  implicit none

  private

  public :: advect_base_dens, advect_base_enthalpy

contains

  subroutine advect_base_dens(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)

    use bl_prof_module, only: bl_prof_timer, build, destroy
    use geometry, only: spherical
    use restrict_base_module, only: restrict_base, fill_ghost_base

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

    use bl_constants_module, only: ZERO
    use make_edge_state_module, only: make_edge_state_1d
    use geometry, only: nr_fine, r_start_coord, r_end_coord, numdisjointchunks, dr

    real(kind=dp_t), intent(in   ) ::                  w0(:,0:)
    real(kind=dp_t), intent(in   ) ::            rho0_old(:,0:)
    real(kind=dp_t), intent(  out) ::            rho0_new(:,0:)
    real(kind=dp_t), intent(  out) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r, n, i, nlevs
    
    real (kind=dp_t) :: force(size(w0,dim=1),0:nr_fine-1)
    real (kind=dp_t) ::  edge(size(w0,dim=1),0:nr_fine)
   
    nlevs = size(w0,dim=1)

    rho0_predicted_edge = ZERO

    ! zero the new density so we don't leave a non-zero density in fine radial
    ! regions that no longer have a corresponding full state
    rho0_new = ZERO

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

    use bl_constants_module, only: HALF
    use make_edge_state_module, only: make_edge_state_1d
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

    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr_fine-1
       force(1,r) = -rho0_old(1,r) * (w0(1,r+1) - w0(1,r)) / dr(1) - &
            2.0_dp_t*rho0_old(1,r)*HALF*(w0(1,r) + w0(1,r+1))/r_cc_loc(1,r)
    end do
    !$OMP END PARALLEL DO
    
    call make_edge_state_1d(rho0_old,edge,w0,force,dt)
    
    rho0_predicted_edge = edge

    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr_fine-1
       rho0_new(1,r) = rho0_old(1,r) - dtdr/r_cc_loc(1,r)**2 * &
            (r_edge_loc(1,r+1)**2 * edge(1,r+1) * w0(1,r+1) - &
            r_edge_loc(1,r  )**2 * edge(1,r  ) * w0(1,r  ))
    end do
    !$OMP END PARALLEL DO

  end subroutine advect_base_dens_spherical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_enthalpy(w0,rho0_old,rhoh0_old,rhoh0_new,rho0_predicted_edge, &
                                  psi,dt)

    use bl_prof_module, only: bl_prof_timer, build, destroy
    use geometry, only: spherical
    use restrict_base_module, only: restrict_base, fill_ghost_base

    real(kind=dp_t), intent(in   ) ::                  w0(:,0:)
    real(kind=dp_t), intent(in   ) ::            rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::           rhoh0_old(:,0:)
    real(kind=dp_t), intent(  out) ::           rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) ::                 psi(:,0:)
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
                                           rho0_predicted_edge,psi,dt)
    end if

    call destroy(bpt)
       
  end subroutine advect_base_enthalpy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_enthalpy_planar(w0,rho0_old,rhoh0_old,rhoh0_new, &
                                         rho0_predicted_edge,psi,dt)

    use bl_constants_module, only: ZERO
    use make_edge_state_module, only: make_edge_state_1d
    use geometry, only: nr_fine, r_start_coord, r_end_coord, numdisjointchunks, dr
    use probin_module, only: enthalpy_pred_type
    use pred_parameters, only: predict_h, predict_T_then_h

    real(kind=dp_t), intent(in   ) ::                  w0(:,0:)
    real(kind=dp_t), intent(in   ) ::            rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::           rhoh0_old(:,0:)
    real(kind=dp_t), intent(  out) ::           rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) ::                 psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r, i, n, nlevs
    
    real (kind=dp_t) :: force(size(w0,dim=1),0:nr_fine-1)
    real (kind=dp_t) ::    h0(size(w0,dim=1),0:nr_fine-1)
    real (kind=dp_t) ::  edge(size(w0,dim=1),0:nr_fine)

    nlevs = size(w0,dim=1)

    ! zero the new enthalpy so we don't leave a non-zero enthalpy in fine radial
    ! regions that no longer have a corresponding full state
    rhoh0_new = ZERO
   
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
                                            rho0_predicted_edge,psi,dt)

    use bl_constants_module
    use make_edge_state_module
    use geometry, only: r_cc_loc, r_edge_loc, dr, nr_fine
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:), rhoh0_old(:,0:)
    real(kind=dp_t), intent(  out) :: rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
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

       force(1,:) = psi(1,:)

       call make_edge_state_1d(h0,edge,w0,force,dt)

       ! our final update needs (rho h)_0 on the edges, so compute
       ! that now
       edge = rho0_predicted_edge*edge

    else

       ! here we predict (rho h)_0 on the edges
       !$OMP PARALLEL DO PRIVATE(r,div_w0_cart)
       do r=0,nr_fine-1

          div_w0_cart = (w0(1,r+1) - w0(1,r)) / dr(1)

          ! add psi at time-level n to the force for the prediction
          force(1,r) = -rhoh0_old(1,r) * div_w0_cart - &
               2.0_dp_t*rhoh0_old(1,r)*HALF*(w0(1,r) + w0(1,r+1))/r_cc_loc(1,r) + psi(1,r)
          
       end do
       !$OMP END PARALLEL DO

       call make_edge_state_1d(rhoh0_old,edge,w0,force,dt)

    endif

    ! update (rho h)_0
    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr_fine-1
       rhoh0_new(1,r) = rhoh0_old(1,r) - dtdr / r_cc_loc(1,r)**2 * &
            (r_edge_loc(1,r+1)**2 * edge(1,r+1) * w0(1,r+1) - &
            r_edge_loc(1,r  )**2 * edge(1,r  ) * w0(1,r  )) + dt * psi(1,r)
    end do
    !$OMP END PARALLEL DO

  end subroutine advect_base_enthalpy_spherical

end module advect_base_module
