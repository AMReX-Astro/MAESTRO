subroutine varden()

  use BoxLib
  use omp_module
  use f2kcli
  use list_box_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use init_module
  use base_state_module
  use box_util_module
  use bl_IO_module
  use variables
  use geometry
  use network
  use eos_module
  use make_w0_module
  use advect_base_module
  use make_grav_module
  use make_edge_state_module
  use probin_module
  use bl_constants_module
  use initialize_module
  use enforce_HSE_module
  use make_psi_module
  use fundamental_constants_module, only: Gconst
  use test_basestate_module

  implicit none

  integer :: n,r,comp,comp2,iter

  real(dp_t) :: frac,delta,sumX,time,dt,dtold
  real(dp_t) :: factor,divw,w0dpdr_nph,w0dpdr_nph_1,w0dpdr_nph_2

  real(dp_t), allocatable ::                  dx(:,:)
  real(dp_t), allocatable ::           grav_cell(:,:)
  real(dp_t), allocatable ::       gamma1bar_old(:,:)
  real(dp_t), allocatable ::       gamma1bar_nph(:,:)
  real(dp_t), allocatable ::       gamma1bar_new(:,:)
  real(dp_t), allocatable ::              s0_old(:,:,:)
  real(dp_t), allocatable ::              s0_new(:,:,:)
  real(dp_t), allocatable ::                  X0(:,:)
  real(dp_t), allocatable ::              p0_old(:,:)
  real(dp_t), allocatable ::              p0_new(:,:)
  real(dp_t), allocatable ::                  w0(:,:)
  real(dp_t), allocatable ::                 psi(:,:)
  real(dp_t), allocatable ::           etarho_ec(:,:)
  real(dp_t), allocatable ::           etarho_cc(:,:)
  real(dp_t), allocatable ::            w0_force(:,:)
  real(dp_t), allocatable ::            Sbar_old(:,:)
  real(dp_t), allocatable ::            Sbar_nph(:,:)
  real(dp_t), allocatable ::            Sbar_new(:,:)
  real(dp_t), allocatable ::  p0_minus_pthermbar(:,:)
  real(dp_t), allocatable :: rho0_predicted_edge(:,:)
  real(dp_t), allocatable ::               force(:,:)
  real(dp_t), allocatable ::                edge(:,:)
  real(dp_t), allocatable ::            Hext_bar(:,:)

  real(dp_t) :: mencl, max_hse_error, starting_rad, rloc, r_r, r_l, g, dpdr, rhog

  call probin_init()
  call init_dm()
  call init_spherical()
  center(1) = ZERO

  call init_variables()

  call network_init()
  print *, 'use_eos_coulomb = ', use_eos_coulomb
  call eos_init(use_eos_coulomb=use_eos_coulomb)

  nlevs = 1
  nlevs_radial = 1

  ! we get the number of points for the base state directly from
  ! the model file and use this to set the resolution.  It should be the 
  ! case that prob_hi(1) agrees with the maximum radius for the model file.
  nr_fine = get_model_npts(model_file)
  print *, 'number of points in model file: ', nr_fine

  allocate(r_start_coord(nlevs,1))
  r_start_coord(nlevs,1) = 0

  allocate(r_end_coord(nlevs,1))
  r_end_coord(nlevs,1) = nr_fine-1

  allocate(numdisjointchunks(nlevs))
  numdisjointchunks(:) = 1

  dr_fine = (prob_hi(1) - prob_lo(1))/nr_fine

  allocate(dx(nlevs,1))
  dx(1,1) = dr_fine

  allocate(dr(nlevs_radial))
  allocate(nr(nlevs_radial))

  allocate(  r_cc_loc(nlevs_radial,0:nr_fine-1))
  allocate(r_edge_loc(nlevs_radial,0:nr_fine))

  nr(1) = nr_fine
  dr(1) = dr_fine

  do r=0,nr(1)-1
     r_cc_loc(1,r) = prob_lo(dm) + (dble(r)+HALF)*dr(1)
  end do
  do r=0,nr(1)
     r_edge_loc(1,r) = prob_lo(dm) + (dble(r))*dr(1)
  end do

  allocate(anelastic_cutoff_coord(1))
  allocate(base_cutoff_density_coord(1))
  allocate(burning_cutoff_density_coord(1))

  allocate(          grav_cell(nlevs_radial,0:nr_fine-1))
  allocate(      gamma1bar_old(nlevs_radial,0:nr_fine-1))
  allocate(      gamma1bar_nph(nlevs_radial,0:nr_fine-1))
  allocate(      gamma1bar_new(nlevs_radial,0:nr_fine-1))
  allocate(             s0_old(nlevs_radial,0:nr_fine-1,nscal))
  allocate(             s0_new(nlevs_radial,0:nr_fine-1,nscal))
  allocate(                 X0(nlevs_radial,0:nr_fine-1))
  allocate(             p0_old(nlevs_radial,0:nr_fine-1))
  allocate(             p0_new(nlevs_radial,0:nr_fine-1))
  allocate(                 w0(nlevs_radial,0:nr_fine))
  allocate(                psi(nlevs_radial,0:nr_fine-1))
  allocate(          etarho_ec(nlevs_radial,0:nr_fine))
  allocate(          etarho_cc(nlevs_radial,0:nr_fine-1))
  allocate(           w0_force(nlevs_radial,0:nr_fine))
  allocate(           Sbar_old(nlevs_radial,0:nr_fine-1))
  allocate(           Sbar_nph(nlevs_radial,0:nr_fine-1))
  allocate(           Sbar_new(nlevs_radial,0:nr_fine-1))
  allocate( p0_minus_pthermbar(nlevs_radial,0:nr_fine-1))
  allocate(rho0_predicted_edge(nlevs_radial,0:nr_fine))
  allocate(              force(nlevs_radial,0:nr_fine-1))
  allocate(               edge(nlevs_radial,0:nr_fine))
  allocate(           Hext_bar(nlevs_radial,0:nr_fine-1))

  gamma1bar_old      = ZERO
  gamma1bar_nph      = ZERO
  gamma1bar_new      = ZERO
  w0                 = ZERO
  psi                = ZERO ! note: psi not used in spherical and is 0 in planar since
                            !       etarho is zero
  etarho_ec          = ZERO
  etarho_cc          = ZERO
  p0_minus_pthermbar = ZERO
  Hext_bar           = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do n=1,nlevs
     call init_base_state(n,model_file,s0_old(n,:,:),p0_old(n,:),dx(n,:))
  enddo

  ! output
  open(unit=10,file="base.orig")
  do r=0,nr_fine-1
     write(10,1000) r_cc_loc(1,r), s0_old(1,r,rho_comp), s0_old(1,r,temp_comp), p0_old(1,r)
  enddo
  close(unit=10)

  call compute_cutoff_coords(s0_old(:,:,rho_comp))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute initial gamma1bar_old
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do r=0,nr_fine-1

     ! (rho, p) --> gamma1bar
     den_eos(1)  = s0_old(1,r,rho_comp)
     p_eos(1)    = p0_old(1,r)
     xn_eos(1,:) = s0_old(1,r,spec_comp:spec_comp-1+nspec)/s0_old(1,r,rho_comp)
     
     temp_eos(1) = s0_old(1,r,temp_comp)
     
     call eos(eos_input_rp, den_eos, temp_eos, NP, nspec, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag)
     
     gamma1bar_old(1,r) = gam1_eos(1)

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main timestepping loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
 
  time  = ZERO
  dt    = 1.d-4
  dtold = dt

  iter = 0

  do while (time < stop_time)

     print *, 'time = ', time

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the heating term and Sbar
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do n = 1, nlevs
        call get_heating(Hext_bar(n,0:))
     enddo


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! make Sbar
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call make_Sbar(Sbar_old(1,:), s0_old(1,:,:), Hext_bar(1,:))

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute w_0
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call make_w0(w0,w0,w0_force,Sbar_old,s0_old(:,:,rho_comp),s0_old(:,:,rho_comp), &
                  p0_old,p0_old,gamma1bar_old,gamma1bar_old,p0_minus_pthermbar,psi, &
                  etarho_ec,etarho_cc,dt,dtold)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update density and compute rho0_predicted_edge
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call advect_base_dens(w0,s0_old(:,:,rho_comp),s0_new(:,:,rho_comp), &
                           rho0_predicted_edge,dt)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! recompute cutoff coordinates now that rho0 has changed
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call compute_cutoff_coords(s0_new(:,:,rho_comp))
  
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute gravity
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call make_grav_cell(grav_cell,s0_new(:,:,rho_comp))

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update species
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! In the real code, this will be done
     ! for the full state.  Here we are faking a 1-d star, so
     ! we need to do this manually, since the species are not
     ! part of the base state
     do comp=spec_comp,spec_comp+nspec-1

        ! here we predict X_0 on the edges
        X0(1,:) = s0_old(1,:,comp)/s0_old(1,:,rho_comp)
        do r=0,nr_fine-1
           X0(1,r) = max(X0(1,r),ZERO)
        end do

        force = ZERO

        call make_edge_state_1d(X0,edge,w0,force,dt)

        ! our final update needs (rho X)_0 on the edges, so compute
        ! that now
        edge(1,:) = rho0_predicted_edge(1,:)*edge(1,:)

        ! update (rho X)_0
        if (spherical .eq. 0) then
           do r=0,nr_fine-1
              s0_new(1,r,comp) = s0_old(1,r,comp) &
                   - (dt/dr(1))*(edge(1,r+1) * w0(1,r+1) - edge(1,r) * w0(1,r))
           end do
        else
           do r=0,nr_fine-1
              s0_new(1,r,comp) = s0_old(1,r,comp) &
                   - (dt/dr(1))/r_cc_loc(1,r)**2* &
                   (r_edge_loc(1,r+1)**2 * edge(1,r+1) * w0(1,r+1) - &
                   r_edge_loc(1,r  )**2 * edge(1,r  ) * w0(1,r  ))
           end do           
        endif

     enddo

     ! don't let the species leave here negative
     do comp=spec_comp,spec_comp+nspec-1
        
        if (minval(s0_new(1,:,comp)) .lt. ZERO) then
           do r=0,nr_fine-1
              if (s0_new(1,r,comp) .lt. ZERO) then
                 delta = -s0_new(1,r,comp)
                 sumX = ZERO
                 do comp2=spec_comp,spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s0_new(1,r,comp2) .ge. ZERO) then
                       sumX = sumX + s0_new(1,r,comp2)
                    endif
                 enddo
                 do comp2=spec_comp,spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s0_new(1,r,comp2) .ge. ZERO) then
                       frac = s0_new(1,r,comp2) / sumX
                       s0_new(1,r,comp2) = s0_new(1,r,comp2) - frac * delta
                    endif
                 enddo
                 s0_new(1,r,comp) = ZERO
              endif
           enddo
        endif
     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update pressure
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! update pressure
     p0_new = p0_old
     call enforce_HSE(s0_new(:,:,rho_comp),p0_new,grav_cell)

     ! compute gamma1bar_new
     do r=0,nr_fine-1

        ! (rho, p) --> gamma1bar
        den_eos(1)  = s0_new(1,r,rho_comp)
        p_eos(1)    = p0_new(1,r)
        xn_eos(1,:) = s0_new(1,r,spec_comp:spec_comp-1+nspec)/s0_new(1,r,rho_comp)

        temp_eos(1) = s0_old(1,r,temp_comp)

        call eos(eos_input_rp, den_eos, temp_eos, NP, nspec, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

        gamma1bar_new(1,r) = gam1_eos(1)

     end do

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update temperature
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do r=0,nr_fine-1

        ! (rho,p) --> T,h, etc
        den_eos(1)  = s0_new(1,r,rho_comp)
        p_eos(1)    = p0_new(1,r)
        xn_eos(1,:) = s0_new(1,r,spec_comp:spec_comp-1+nspec)/s0_new(1,r,rho_comp)

        temp_eos(1) = s0_old(1,r,temp_comp)

        call eos(eos_input_rp, den_eos, temp_eos, NP, nspec, &
                 xn_eos, &
                 p_eos, h_eos, e_eos, &
                 cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                 dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                 dpdX_eos, dhdX_eos, &
                 gam1_eos, cs_eos, s_eos, &
                 dsdt_eos, dsdr_eos, &
                 do_diag)

        s0_new(1,r,temp_comp) = temp_eos(1)

     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! SECOND HALF OF ALGORITHM
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! reset cutoff coordinates
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call compute_cutoff_coords(s0_old(:,:,rho_comp))

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! make Sbar
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call make_Sbar(Sbar_new(1,:),s0_new(1,:,:),Hext_bar(1,:))
     Sbar_nph(:,:) = HALF*(Sbar_old(:,:) + Sbar_new(:,:))

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute w_0
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call make_w0(w0,w0,w0_force,Sbar_nph,s0_old(:,:,rho_comp),s0_new(:,:,rho_comp), &
                  p0_old,p0_new,gamma1bar_old,gamma1bar_new,p0_minus_pthermbar,psi, &
                  etarho_ec,etarho_cc,dt,dtold)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update density and compute rho0_predicted_edge
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call advect_base_dens(w0,s0_old(:,:,rho_comp),s0_new(:,:,rho_comp), &
                           rho0_predicted_edge,dt)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! recompute cutoff coordinates now that rho0 has changed
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call compute_cutoff_coords(s0_new(:,:,rho_comp))
  
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute gravity
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call make_grav_cell(grav_cell,s0_new(:,:,rho_comp))

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update species
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! In the real code, this will be done
     ! for the full state.  Here we are faking a 1-d star, so
     ! we need to do this manually, since the species are not
     ! part of the base state
     do comp=spec_comp,spec_comp+nspec-1

        ! here we predict X_0 on the edges
        X0(1,:) = s0_old(1,:,comp)/s0_old(1,:,rho_comp)
        do r=0,nr_fine-1
           X0(1,r) = max(X0(1,r),ZERO)
        end do

        force = ZERO

        call make_edge_state_1d(X0,edge,w0,force,dt)

        ! our final update needs (rho X)_0 on the edges, so compute
        ! that now
        edge(1,:) = rho0_predicted_edge(1,:)*edge(1,:)

        ! update (rho X)_0
        if (spherical .eq. 0) then
           do r=0,nr_fine-1
              s0_new(1,r,comp) = s0_old(1,r,comp) &
                   - (dt/dr(1))*(edge(1,r+1) * w0(1,r+1) - edge(1,r) * w0(1,r))
           end do
        else
           do r=0,nr_fine-1
              s0_new(1,r,comp) = s0_old(1,r,comp) &
                   - (dt/dr(1))/r_cc_loc(1,r)**2* &
                   (r_edge_loc(1,r+1)**2 * edge(1,r+1) * w0(1,r+1) - &
                   r_edge_loc(1,r  )**2 * edge(1,r  ) * w0(1,r  ))
           end do           
        endif

     enddo

     ! don't let the species leave here negative
     do comp=spec_comp,spec_comp+nspec-1
        
        if (minval(s0_new(1,:,comp)) .lt. ZERO) then
           do r=0,nr_fine-1
              if (s0_new(1,r,comp) .lt. ZERO) then
                 delta = -s0_new(1,r,comp)
                 sumX = ZERO
                 do comp2=spec_comp,spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s0_new(1,r,comp2) .ge. ZERO) then
                       sumX = sumX + s0_new(1,r,comp2)
                    endif
                 enddo
                 do comp2=spec_comp,spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s0_new(1,r,comp2) .ge. ZERO) then
                       frac = s0_new(1,r,comp2) / sumX
                       s0_new(1,r,comp2) = s0_new(1,r,comp2) - frac * delta
                    endif
                 enddo
                 s0_new(1,r,comp) = ZERO
              endif
           enddo
        endif
     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update pressure
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! update pressure
     p0_new = p0_old
     call enforce_HSE(s0_new(:,:,rho_comp),p0_new,grav_cell)

     ! compute gamma1bar_new
     do r=0,nr_fine-1

        ! (rho, p) --> gamma1bar
        den_eos(1)  = s0_new(1,r,rho_comp)
        p_eos(1)    = p0_new(1,r)
        xn_eos(1,:) = s0_new(1,r,spec_comp:spec_comp-1+nspec)/s0_new(1,r,rho_comp)

        temp_eos(1) = s0_old(1,r,temp_comp)

        call eos(eos_input_rp, den_eos, temp_eos, NP, nspec, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

        gamma1bar_new(1,r) = gam1_eos(1)

     end do

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update temperature
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do r=0,nr_fine-1

        ! (rho,p) --> T,h, etc
        den_eos(1)  = s0_new(1,r,rho_comp)
        p_eos(1)    = p0_new(1,r)
        xn_eos(1,:) = s0_new(1,r,spec_comp:spec_comp-1+nspec)/s0_new(1,r,rho_comp)

        temp_eos(1) = s0_old(1,r,temp_comp)

        call eos(eos_input_rp, den_eos, temp_eos, NP, nspec, &
                 xn_eos, &
                 p_eos, h_eos, e_eos, &
                 cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                 dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                 dpdX_eos, dhdX_eos, &
                 gam1_eos, cs_eos, s_eos, &
                 dsdt_eos, dsdr_eos, &
                 do_diag)

        s0_new(1,r,temp_comp) = temp_eos(1)

     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the new timestep
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     time  = time + dt
     dtold = dt

     dt = min(1.1d0*dt,cflfac*dr_fine/maxval(abs(w0)))
     if (time+dt > stop_time) dt = stop_time - time

     iter = iter + 1

     s0_old = s0_new
     p0_old = p0_new
     gamma1bar_old = gamma1bar_new

  enddo

  ! check the HSE-ness
  if (spherical .eq. 1) then
     mencl = four3rd*m_pi*dr(1)**3*s0_new(1,0,rho_comp)
  endif
  
  max_hse_error = -1.d30

  if (spherical .eq. 0) then
     starting_rad = prob_lo(1)
  else
     starting_rad = ZERO
  endif

  do r=1,nr(1)-1

     rloc = starting_rad + (dble(r) + HALF)*dr(1)

     if (r < base_cutoff_density_coord(1)) then

        r_r = starting_rad + dble(r+1)*dr(1)
        r_l = starting_rad + dble(r)*dr(1)

        if (spherical .eq. 1) then
           g = -Gconst*mencl/r_l**2
           mencl = mencl &
                + four3rd*m_pi*dr(1)*(r_l**2+r_l*r_r+r_r**2)*s0_new(1,r,rho_comp)
        else
           g = grav_const
        endif

        dpdr = (p0_new(1,r) - p0_new(1,r-1))/dr(1)
        rhog = HALF*(s0_new(1,r,rho_comp) + s0_new(1,r-1,rho_comp))*g

        max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))

     end if

  enddo

  print *, 'maximum HSE error = ', max_hse_error

  ! output
  open(unit=10,file="base.new")
  write(10,*) "ONE=r_cc, TWO=rho, THREE=temp, FOUR=p0"
  do r=0,nr_fine-1
     write(10,1000) r_cc_loc(1,r), s0_new(1,r,rho_comp), s0_new(1,r,temp_comp), &
          p0_new(1,r)
  enddo
  close(unit=10)

  open(unit=10,file="base.new_w0")
  write(10,*) "ONE=r_edge, TWO=w0"
  do r=0,nr_fine
     write(10,1000) r_edge_loc(1,r), w0(1,r)
  enddo
  close(unit=10)

1000 format(1x,6(g20.10))

  deallocate(r_start_coord,r_end_coord,numdisjointchunks)
  deallocate(dx,dr,nr,r_cc_loc,r_edge_loc)
  deallocate(anelastic_cutoff_coord,base_cutoff_density_coord,burning_cutoff_density_coord)
  deallocate(grav_cell,gamma1bar_old,gamma1bar_nph,gamma1bar_new,s0_old,s0_new)
  deallocate(X0,p0_old,p0_new,w0,psi,etarho_ec,etarho_cc,w0_force)
  deallocate(Sbar_old,Sbar_new,Sbar_nph,p0_minus_pthermbar,rho0_predicted_edge,force,edge)

end subroutine varden
