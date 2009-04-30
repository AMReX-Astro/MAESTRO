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

  real(dp_t), allocatable ::                  dx(:,:)
  real(dp_t), allocatable ::       div_coeff_old(:,:)
  real(dp_t), allocatable ::           div_coeff(:,:)
  real(dp_t), allocatable ::           grav_cell(:,:)
  real(dp_t), allocatable ::           gamma1bar(:,:)
  real(dp_t), allocatable ::              s0_old(:,:,:)
  real(dp_t), allocatable ::                  s0(:,:,:)
  real(dp_t), allocatable ::              p0_old(:,:)
  real(dp_t), allocatable ::                  p0(:,:)
  real(dp_t), allocatable ::                  w0(:,:)
  real(dp_t), allocatable ::                 psi(:,:)
  real(dp_t), allocatable ::             psi_old(:,:)
  real(dp_t), allocatable ::           etarho_ec(:,:)
  real(dp_t), allocatable ::           etarho_cc(:,:)
  real(dp_t), allocatable ::          div_etarho(:,:)
  real(dp_t), allocatable ::            w0_force(:,:)
  real(dp_t), allocatable ::             Sbar_in(:,:)
  real(dp_t), allocatable ::  p0_minus_pthermbar(:,:)
  real(dp_t), allocatable :: rho0_predicted_edge(:,:)
  real(dp_t), allocatable ::               force(:,:)
  real(dp_t), allocatable ::                  X0(:,:)
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

  allocate(             s0_old(nlevs_radial,0:nr_fine-1,nscal))
  allocate(                 s0(nlevs_radial,0:nr_fine-1,nscal))
  allocate(      div_coeff_old(nlevs_radial,0:nr_fine-1))
  allocate(          div_coeff(nlevs_radial,0:nr_fine-1))
  allocate(          grav_cell(nlevs_radial,0:nr_fine-1))
  allocate(          gamma1bar(nlevs_radial,0:nr_fine-1))
  allocate(             p0_old(nlevs_radial,0:nr_fine-1))
  allocate(                 p0(nlevs_radial,0:nr_fine-1))
  allocate(                psi(nlevs_radial,0:nr_fine-1))
  allocate(            psi_old(nlevs_radial,0:nr_fine-1))
  allocate(          etarho_cc(nlevs_radial,0:nr_fine-1))
  allocate(         div_etarho(nlevs_radial,0:nr_fine-1))
  allocate(            Sbar_in(nlevs_radial,0:nr_fine-1))
  allocate( p0_minus_pthermbar(nlevs_radial,0:nr_fine-1))
  allocate(              force(nlevs_radial,0:nr_fine-1))
  allocate(                 X0(nlevs_radial,0:nr_fine-1))
  allocate(                 w0(nlevs_radial,0:nr_fine))
  allocate(          etarho_ec(nlevs_radial,0:nr_fine))
  allocate(           w0_force(nlevs_radial,0:nr_fine))
  allocate(rho0_predicted_edge(nlevs_radial,0:nr_fine))
  allocate(               edge(nlevs_radial,0:nr_fine))

  gamma1bar          = ZERO
  w0                 = ZERO
  psi                = ZERO
  psi_old            = ZERO
  etarho_ec          = ZERO
  etarho_cc          = ZERO
  div_etarho         = ZERO
  p0_minus_pthermbar = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do n=1,nlevs
     call init_base_state(n,model_file,s0(n,:,:),p0(n,:),dx(n,:))
  enddo

  ! output
  open(unit=10,file="base.orig")
  do r=0,nr_fine-1
     write(10,1000) r_cc_loc(1,r), s0(1,r,rho_comp), s0(1,r,temp_comp), p0(1,r)
  enddo
  close(unit=10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main timestepping loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
 
  time  = ZERO
  dt    = 1.d-4
  dtold = dt

  iter = 0

  allocate(anelastic_cutoff_coord(1))
  allocate(base_cutoff_density_coord(1))
  allocate(burning_cutoff_density_coord(1))

  do while (time < stop_time)

     print *, 'time = ', time

     ! here we set the old state to the current state 
     ! advect_base will update the state
     s0_old = s0
     p0_old = p0

     call compute_cutoff_coords(s0(:,:,rho_comp))

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
        den_eos(1)  = s0(1,r,rho_comp)
        temp_eos(1) = s0(1,r,temp_comp)
        xn_eos(1,:) = s0(1,r,spec_comp:spec_comp-1+nspec)/s0(1,r,rho_comp)
        
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

     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute w_0
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call make_w0(w0,w0,w0_force,Sbar_in,s0(:,:,rho_comp),s0(:,:,rho_comp),p0,p0, &
                  gamma1bar,gamma1bar,p0_minus_pthermbar,psi,etarho_ec,etarho_cc,dt,dtold)
  
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute gravity
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call make_grav_cell(grav_cell,s0(:,:,rho_comp))

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the divergence coefficient, beta_0
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call make_div_coeff(div_coeff,s0(:,:,rho_comp),p0,gamma1bar,grav_cell)     

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update density and compute rho0_predicted_edge
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call advect_base_dens(w0,s0_old(:,:,rho_comp),s0(:,:,rho_comp),rho0_predicted_edge,dt)

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
           s0(1,r,comp) = s0_old(1,r,comp) &
                - (dt/dr(1))/r_cc_loc(1,r)**2* &
                (r_edge_loc(1,r+1)**2 * edge(1,r+1) * w0(1,r+1) - &
                 r_edge_loc(1,r  )**2 * edge(1,r  ) * w0(1,r  ))
        end do

     enddo

     ! don't let the species leave here negative
     do comp=spec_comp,spec_comp+nspec-1
        
        if (minval(s0(1,:,comp)) .lt. ZERO) then
           do r=0,nr_fine-1
              if (s0(1,r,comp) .lt. ZERO) then
                 delta = -s0(1,r,comp)
                 sum = ZERO
                 do comp2=spec_comp,spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s0(1,r,comp2) .ge. ZERO) then
                       sum = sum + s0(1,r,comp2)
                    endif
                 enddo
                 do comp2=spec_comp,spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s0(1,r,comp2) .ge. ZERO) then
                       frac = s0(1,r,comp2) / sum
                       s0(1,r,comp2) = s0(1,r,comp2) - frac * delta
                    endif
                 enddo
                 s0(1,r,comp) = ZERO
              endif
           enddo
        endif
     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update pressure
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     if (p0_update_type .eq. 1) then
!
!        call advect_base_pres(w0,Sbar_in,p0_old,p0,gamma1bar,psi,psi_old,etarho_cc,dt)
!
!     else
!
!        ! set new p0 through HSE
!        p0 = p0_old
!        call enforce_HSE(s0(:,:,rho_comp),p0,grav_cell)
!
!        ! make psi
!        if (spherical .eq. 0) then
!           call make_psi_planar(etarho_cc,psi)
!        else
!           call make_psi_spherical(psi,w0,gamma1bar,p0_old,p0,Sbar_in)
!        end if
!
!     end if

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! update temperature
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do r=0,nr_fine-1

        ! (rho, T) --> p,h, etc
        den_eos(1)  = s0(1,r,rho_comp)
        temp_eos(1) = s0(1,r,temp_comp)
        p_eos(1)    = p0(1,r)
        xn_eos(1,:) = s0(1,r,spec_comp:spec_comp-1+nspec)/s0(1,r,rho_comp)

        call eos(eos_input_rp, den_eos, temp_eos, NP, nspec, &
                 xn_eos, &
                 p_eos, h_eos, e_eos, &
                 cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                 dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                 dpdX_eos, dhdX_eos, &
                 gam1_eos, cs_eos, s_eos, &
                 dsdt_eos, dsdr_eos, &
                 do_diag)

        s0(1,r,temp_comp) = temp_eos(1)

     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the new timestep
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     time  = time + dt
     dtold = dt

     dt = min(1.1*dt,cflfac*dr_fine/maxval(abs(w0)))
     if (time+dt > stop_time) dt = stop_time - time

     iter = iter + 1

  enddo

  ! output
  open(unit=10,file="base.new")
  do r=0,nr_fine-1
     write(10,1000) r_cc_loc(1,r), s0(1,r,rho_comp), s0(1,r,temp_comp), p0(1,r), w0(1,r)
  enddo
  close(unit=10)

1000 format(1x,6(g20.10))

  deallocate(r_start_coord,r_end_coord,numdisjointchunks)
  deallocate(dx,dr,nr,r_cc_loc,r_edge_loc)
  deallocate(s0_old,s0,div_coeff_old,div_coeff,grav_cell,gamma1bar,p0_old,p0,psi,psi_old)
  deallocate(etarho_cc,div_etarho,Sbar_in,p0_minus_pthermbar,force,X0,w0)
  deallocate(etarho_ec,w0_force,rho0_predicted_edge,edge)
  deallocate(anelastic_cutoff_coord,base_cutoff_density_coord,burning_cutoff_density_coord)

end subroutine varden
