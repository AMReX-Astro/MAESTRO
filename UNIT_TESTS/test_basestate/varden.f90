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
  use make_div_coeff_module
  use make_edge_state_module
  use probin_module
  use bl_constants_module
  use initialize_module
  use enforce_HSE_module
  use make_psi_module

  implicit none

  integer :: n,r,comp,comp2,iter

  real(dp_t) :: frac,delta,sum,time,dt,dtold,y_0,Hbar
  real(dp_t) :: factor,divw,w0dpdr_nph,w0dpdr_nph_1,w0dpdr_nph_2

  real(dp_t), allocatable ::                  dx(:,:)
  real(dp_t), allocatable ::           div_coeff(:,:)
  real(dp_t), allocatable ::           grav_cell(:,:)
  real(dp_t), allocatable ::           gamma1bar(:,:)
  real(dp_t), allocatable ::       gamma1bar_nph(:,:)
  real(dp_t), allocatable ::              s0_old(:,:,:)
  real(dp_t), allocatable ::              s0_new(:,:,:)
  real(dp_t), allocatable ::                  X0(:,:)
  real(dp_t), allocatable ::              p0_old(:,:)
  real(dp_t), allocatable ::              p0_new(:,:)
  real(dp_t), allocatable ::              p0_nph(:,:)
  real(dp_t), allocatable ::                  w0(:,:)
  real(dp_t), allocatable ::                 psi(:,:)
  real(dp_t), allocatable ::           etarho_ec(:,:)
  real(dp_t), allocatable ::           etarho_cc(:,:)
  real(dp_t), allocatable ::            w0_force(:,:)
  real(dp_t), allocatable ::             Sbar_in(:,:)
  real(dp_t), allocatable ::  p0_minus_pthermbar(:,:)
  real(dp_t), allocatable :: rho0_predicted_edge(:,:)
  real(dp_t), allocatable ::               force(:,:)
  real(dp_t), allocatable ::                edge(:,:)

  call probin_init()
  call init_dm()
  call init_spherical()
  center(1) = ZERO

  call init_variables()

  call network_init()
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

  allocate(          div_coeff(nlevs_radial,0:nr_fine-1))
  allocate(          grav_cell(nlevs_radial,0:nr_fine-1))
  allocate(          gamma1bar(nlevs_radial,0:nr_fine-1))
  allocate(      gamma1bar_nph(nlevs_radial,0:nr_fine-1))
  allocate(             s0_old(nlevs_radial,0:nr_fine-1,nscal))
  allocate(             s0_new(nlevs_radial,0:nr_fine-1,nscal))
  allocate(                 X0(nlevs_radial,0:nr_fine-1))
  allocate(             p0_old(nlevs_radial,0:nr_fine-1))
  allocate(             p0_new(nlevs_radial,0:nr_fine-1))
  allocate(             p0_nph(nlevs_radial,0:nr_fine-1))
  allocate(                 w0(nlevs_radial,0:nr_fine))
  allocate(                psi(nlevs_radial,0:nr_fine-1))
  allocate(          etarho_ec(nlevs_radial,0:nr_fine))
  allocate(          etarho_cc(nlevs_radial,0:nr_fine-1))
  allocate(           w0_force(nlevs_radial,0:nr_fine))
  allocate(            Sbar_in(nlevs_radial,0:nr_fine-1))
  allocate( p0_minus_pthermbar(nlevs_radial,0:nr_fine-1))
  allocate(rho0_predicted_edge(nlevs_radial,0:nr_fine))
  allocate(              force(nlevs_radial,0:nr_fine-1))
  allocate(               edge(nlevs_radial,0:nr_fine))

  gamma1bar          = ZERO
  w0                 = ZERO
  psi                = ZERO
  etarho_ec          = ZERO
  etarho_cc          = ZERO
  p0_minus_pthermbar = ZERO

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main timestepping loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
 
  time  = ZERO
  dt    = 1.d-4
  dtold = dt

  iter = 0

  do while (time < stop_time)

     print *, 'time = ', time

     call compute_cutoff_coords(s0_old(:,:,rho_comp))

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the heating term and gamma1bar
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     print *, 'calling the eos', nr_fine
     do r=0,nr_fine-1

        if (spherical .eq. 0) then
           ! plane-parallel -- do the heating term in paper II (section 4)
           y_0 = 4.d7
           Hbar = 1.d17 * exp(-((r_cc_loc(1,r) - y_0)**2)/ 1.d14)
        else
           ! spherical -- lower amplitude heating term
           y_0 = 4.d7
           Hbar = 1.d16 * exp(-((r_cc_loc(1,r) - y_0)**2)/ 1.d14)
        endif

        ! (rho, T) --> p,h, etc
        den_eos(1)  = s0_old(1,r,rho_comp)
        temp_eos(1) = s0_old(1,r,temp_comp)
        xn_eos(1,:) = s0_old(1,r,spec_comp:spec_comp-1+nspec)/s0_old(1,r,rho_comp)

        call eos(eos_input_rt, den_eos, temp_eos, NP, nspec, &
                 xn_eos, &
                 p_eos, h_eos, e_eos, &
                 cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                 dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                 dpdX_eos, dhdX_eos, &
                 gam1_eos, cs_eos, s_eos, &
                 dsdt_eos, dsdr_eos, &
                 do_diag)

        gamma1bar(1,r) = gam1_eos(1)

        Sbar_in(1,r) = Hbar * dpdt_eos(1) / (den_eos(1) * cp_eos(1) * dpdr_eos(1))

        ! NOTE: This was here before; I don't know why.  It seems to defeat the purpose
        !       of calling advect base pressure.  Without this line, you get wild
        !       oscillations in the solution.  Even with older versions of the code you
        !       get oscillations if you remove this line.
        p0_old(1,r) = p_eos(1)

     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute w_0
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call make_w0(w0,w0,w0_force,Sbar_in,s0_old(:,:,rho_comp),s0_old(:,:,rho_comp), &
                  p0_old,p0_old,gamma1bar,gamma1bar,p0_minus_pthermbar,psi, &
                  etarho_ec,etarho_cc,dt,dtold)
  
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute gravity
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call make_grav_cell(grav_cell,s0_old(:,:,rho_comp))

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the divergence coefficient, beta_0
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call make_div_coeff(div_coeff,s0_old(:,:,rho_comp),p0_old,gamma1bar,grav_cell)     

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update density and compute rho0_predicted_edge
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call advect_base_dens(w0,s0_old(:,:,rho_comp),s0_new(:,:,rho_comp), &
                           rho0_predicted_edge,dt)

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

        if (spherical .eq. 0) then
           call make_edge_state_1d(X0,edge,w0,force,dt)
        else
           call make_edge_state_1d(X0,edge,w0,force,dt)
        endif

        ! our final update needs (rho X)_0 on the edges, so compute
        ! that now
        edge(1,:) = rho0_predicted_edge(1,:)*edge(1,:)

        ! update (rho X)_0
        do r=0,nr_fine-1
           s0_new(1,r,comp) = s0_old(1,r,comp) &
                - (dt/dr(1))/r_cc_loc(1,r)**2* &
                (r_edge_loc(1,r+1)**2 * edge(1,r+1) * w0(1,r+1) - &
                 r_edge_loc(1,r  )**2 * edge(1,r  ) * w0(1,r  ))
        end do

     enddo

     ! don't let the species leave here negative
     do comp=spec_comp,spec_comp+nspec-1
        
        if (minval(s0_new(1,:,comp)) .lt. ZERO) then
           do r=0,nr_fine-1
              if (s0_new(1,r,comp) .lt. ZERO) then
                 delta = -s0_new(1,r,comp)
                 sum = ZERO
                 do comp2=spec_comp,spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s0_new(1,r,comp2) .ge. ZERO) then
                       sum = sum + s0_new(1,r,comp2)
                    endif
                 enddo
                 do comp2=spec_comp,spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s0_new(1,r,comp2) .ge. ZERO) then
                       frac = s0_new(1,r,comp2) / sum
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

     if (p0_update_type .eq. 1) then

        if (spherical .eq. 1) then

           do r=0,nr_fine-1
              
              divw = one/(r_cc_loc(1,r)**2) * &
                   (r_edge_loc(1,r+1)**2 * w0(1,r+1) - &
                    r_edge_loc(1,r  )**2 * w0(1,r  )) / dr(1)
              
              if (r .eq. 0) then
                 w0dpdr_nph_2 =  w0(1,2) * (p0_old(1,2)-p0_old(1,1)) / dr(1)
                 w0dpdr_nph_1 =  w0(1,1) * (p0_old(1,1)-p0_old(1,0)) / dr(1)
                 w0dpdr_nph =  1.5d0 * w0dpdr_nph_1 - 0.5d0 * w0dpdr_nph_2
              else if (r .eq. nr_fine-1) then
                 w0dpdr_nph_2 =  w0(1,nr_fine-1) &
                      * (p0_old(1,nr_fine-1)-p0_old(1,nr_fine-2)) / dr(1)
                 w0dpdr_nph_1 =  w0(1,nr_fine-2) &
                      * (p0_old(1,nr_fine-2)-p0_old(1,nr_fine-3)) / dr(1)
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

           ! compute provisionally updated gamma1bar and store it in gamma1bar_nph
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
              
              gamma1bar_nph(1,r) = gam1_eos(1)

           end do

           ! compute gamma1bar_nph
           gamma1bar_nph = HALF*(gamma1bar + gamma1bar_nph)
           
           ! compute p0_nph
           p0_nph = HALF*(p0_old + p0_new)

           do r=0,nr_fine-1

              divw = one/(r_cc_loc(1,r)**2) * &
                   (r_edge_loc(1,r+1)**2 * w0(1,r+1) - &
                   r_edge_loc(1,r  )**2 * w0(1,r  )) / dr(1)

              if (r .eq. 0) then
                 w0dpdr_nph_2 =  w0(1,2) * (p0_nph(1,2)-p0_nph(1,1)) / dr(1)
                 w0dpdr_nph_1 =  w0(1,1) * (p0_nph(1,1)-p0_nph(1,0)) / dr(1)
                 w0dpdr_nph =  1.5d0 * w0dpdr_nph_1 - 0.5d0 * w0dpdr_nph_2
              else if (r .eq. nr_fine-1) then
                 w0dpdr_nph_2 =  w0(1,nr_fine-1) &
                      * (p0_nph(1,nr_fine-1)-p0_nph(1,nr_fine-2)) / dr(1)
                 w0dpdr_nph_1 =  w0(1,nr_fine-2) &
                      * (p0_nph(1,nr_fine-2)-p0_nph(1,nr_fine-3)) / dr(1)
                 w0dpdr_nph =  1.5d0 * w0dpdr_nph_2 - 0.5d0 * w0dpdr_nph_1
              else
                 w0dpdr_nph =  HALF * ( w0(1,r+1)*(p0_nph(1,r+1)-p0_nph(1,r  )) &
                      +w0(1,r  )*(p0_nph(1,r  )-p0_nph(1,r-1)) ) / dr(1)
              end if
              
!              if (r .eq. 0) then
!                 w0dpdr_nph_2 =  HALF * w0(1,2) * ( (p0_old(1,2)-p0_old(1,1)) &
!                      +(p0_new(1,2)-p0_new(1,1)) ) / dr(1)
!                 w0dpdr_nph_1 =  HALF * w0(1,1) * ( (p0_old(1,1)-p0_old(1,0)) &
!                      +(p0_new(1,1)-p0_new(1,0)) ) / dr(1)
!                 w0dpdr_nph =  1.5d0 * w0dpdr_nph_1 - 0.5d0 * w0dpdr_nph_2
!              else if (r .eq. nr_fine-1) then
!                 w0dpdr_nph_2 = HALF * w0(1,nr_fine-1) * &
!                      ((p0_old(1,nr_fine-1)-p0_old(1,nr_fine-2)) &
!                      +(p0_new(1,nr_fine-1)-p0_new(1,nr_fine-2))) / dr(1)
!                 w0dpdr_nph_1 = HALF * w0(1,nr_fine-2)* &
!                      ((p0_old(1,nr_fine-2)-p0_old(1,nr_fine-3)) &
!                      +(p0_new(1,nr_fine-2)-p0_new(1,nr_fine-3))) / dr(1)
!                 w0dpdr_nph =  1.5d0 * w0dpdr_nph_2 - 0.5d0 * w0dpdr_nph_1
!              else
!                 w0dpdr_nph = HALF * HALF * ( w0(1,r+1)*(p0_old(1,r+1)-p0_old(1,r)) + &
!                      w0(1,r)*(p0_old(1,r)-p0_old(1,r-1)) + &
!                      w0(1,r+1)*(p0_new(1,r+1)-p0_new(1,r)) + &
!                      w0(1,r)*(p0_new(1,r)-p0_new(1,r-1)) ) / dr(1)
!              end if

              factor = Sbar_in(1,r) - divw &
                   - 1.d0 / (gamma1bar_nph(1,r)*p0_nph(1,r)) * w0dpdr_nph
              factor = half * dt * factor

              p0_new(1,r) = p0_old(1,r) * (one + gamma1bar_nph(1,r)*factor ) / &
                                          (one - gamma1bar_nph(1,r)*factor)

           end do

           ! compute updated gamma1bar and store it in gamma1bar_nph
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
              
              gamma1bar_nph(1,r) = gam1_eos(1)

           end do

           ! compute gamma1bar_nph
           gamma1bar_nph = HALF*(gamma1bar + gamma1bar_nph)

           ! compute p0_nph
           p0_nph = HALF*(p0_old + p0_new)

           ! make psi
           call make_psi_spherical(psi,w0,gamma1bar_nph,p0_nph,Sbar_in)

        else
              
           ! make psi
           call make_psi_planar(etarho_cc,psi)

           ! advect p0
           force = psi
    
           call make_edge_state_1d(p0_old,edge,w0,force,dt)
    
           do n=1,nlevs
              do r=r_start_coord(n,1),r_end_coord(n,1)
                 p0_new(n,r) = p0_old(n,r) - dt / dr(n) &
                      * HALF * (w0(n,r) + w0(n,r+1)) * (edge(n,r+1) - edge(n,r)) &
                      + dt * psi(n,r)
              end do
           end do

        endif

     else

        ! set new p0 through HSE
        p0_new = p0_old
        call enforce_HSE(s0_new(:,:,rho_comp),p0_new,grav_cell)

        ! make psi
        if (spherical .eq. 0) then
           call make_psi_planar(etarho_cc,psi)
        else

           ! compute updated gamma1bar and store it in gamma1bar_nph
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
              
              gamma1bar_nph(1,r) = gam1_eos(1)

           end do

           ! compute gamma1bar_nph
           gamma1bar_nph = HALF*(gamma1bar+gamma1bar_nph)

           ! compute p0_nph
           p0_nph = HALF*(p0_old + p0_new)

           call make_psi_spherical(psi,w0,gamma1bar_nph,p0_nph,Sbar_in)

        end if

     end if

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

     dt = min(1.1*dt,cflfac*dr_fine/maxval(abs(w0)))
     if (time+dt > stop_time) dt = stop_time - time

     iter = iter + 1

     s0_old = s0_new
     p0_old = p0_new

  enddo

  ! output
  open(unit=10,file="base.new")
  do r=0,nr_fine-1
     write(10,1000) r_cc_loc(1,r), s0_new(1,r,rho_comp), s0_new(1,r,temp_comp), &
          p0_new(1,r), w0(1,r)
  enddo
  close(unit=10)

1000 format(1x,6(g20.10))

  deallocate(r_start_coord,r_end_coord,numdisjointchunks)
  deallocate(dx,dr,nr,r_cc_loc,r_edge_loc)
  deallocate(anelastic_cutoff_coord,base_cutoff_density_coord,burning_cutoff_density_coord)
  deallocate(div_coeff,grav_cell,gamma1bar,gamma1bar_nph,s0_old,s0_new)
  deallocate(X0,p0_old,p0_new,p0_nph,w0,psi,etarho_ec,etarho_cc,w0_force)
  deallocate(Sbar_in,p0_minus_pthermbar,rho0_predicted_edge,force,edge)

end subroutine varden
