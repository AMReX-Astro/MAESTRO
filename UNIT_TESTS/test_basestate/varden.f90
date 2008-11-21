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

  implicit none

  integer :: i, comp,comp2
  integer :: iter

  real(dp_t) :: frac, delta, sum

  real(dp_t) :: time,dt,dtold

  real(dp_t) :: y_0

  real(dp_t), allocatable :: dx(:,:)

  real(dp_t), allocatable :: div_coeff_old(:,:)
  real(dp_t), allocatable :: div_coeff(:,:)
  real(dp_t), allocatable :: grav_cell(:,:)
  real(dp_t), allocatable :: gam1(:,:)
  real(dp_t), allocatable :: s0_old(:,:,:)
  real(dp_t), allocatable :: s0(:,:,:)
  real(dp_t), allocatable :: p0_old(:,:)
  real(dp_t), allocatable :: p0(:,:)
  real(dp_t), allocatable :: w0(:,:)
  real(dp_t), allocatable :: w0_old(:,:)
  real(dp_t), allocatable :: psi(:,:)
  real(dp_t), allocatable :: etarho_ec(:,:)
  real(dp_t), allocatable :: etarho_cc(:,:)
  real(dp_t), allocatable :: div_etarho(:,:)
  real(dp_t), allocatable :: f(:,:)
  real(dp_t), allocatable :: Sbar_in(:,:)
  real(dp_t), allocatable :: p0_minus_pthermbar(:,:)
  real(dp_t), allocatable :: rho0_predicted_edge(:,:)
  real(dp_t), allocatable :: force(:,:)
  real(dp_t), allocatable :: X0(:,:)
  real(dp_t), allocatable :: edge(:,:)

  real(dp_t) :: coeff, Hbar

  integer :: n
  integer :: which_step

  call probin_init()
  call init_dm()
  call init_spherical()
  center(1) = ZERO

  call init_variables()

  call network_init()
  call eos_init(use_eos_coulomb=use_eos_coulomb)

  nlevs = 1


  ! we get the number of points for the base state directly from
  ! the model file and use this to set the resolution.  It should be
  ! the case that prob_hi(1) agrees with the maximum radius for the model
  ! file.
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

  allocate(dr(1))
  allocate(nr(1))

  allocate(  r_cc_loc(1,0:nr_fine-1))
  allocate(r_edge_loc(1,0:nr_fine))

  nr(1) = nr_fine
  dr(1) = dr_fine

  do i = 0,nr(1)-1
     r_cc_loc(1,i) = prob_lo(dm) + (dble(i)+HALF)*dr(1)
  end do
  do i = 0,nr(1)
     r_edge_loc(1,i) = prob_lo(dm) + (dble(i))*dr(1)
  end do

  allocate(div_coeff_old(nlevs,0:nr_fine-1))
  allocate(    div_coeff(nlevs,0:nr_fine-1))

  allocate(grav_cell(nlevs,0:nr_fine-1))

  allocate(               gam1(nlevs,0:nr_fine-1  ))
  allocate(             s0_old(nlevs,0:nr_fine-1, nscal))
  allocate(                 s0(nlevs,0:nr_fine-1, nscal))
  allocate(             p0_old(nlevs,0:nr_fine-1  ))
  allocate(                 p0(nlevs,0:nr_fine-1  ))
  allocate(             w0_old(nlevs,0:nr_fine))
  allocate(                 w0(nlevs,0:nr_fine))
  allocate(                psi(nlevs,0:nr_fine-1))
  allocate(          etarho_ec(nlevs,0:nr_fine))
  allocate(          etarho_cc(nlevs,0:nr_fine-1))
  allocate(         div_etarho(nlevs,0:nr_fine-1))
  allocate(                  f(nlevs,0:nr_fine))
  allocate(            Sbar_in(nlevs,0:nr_fine-1))
  allocate( p0_minus_pthermbar(nlevs,0:nr_fine-1))
  allocate(rho0_predicted_edge(nlevs,0:nr_fine))

  allocate(force(nlevs,0:nr_fine-1))
  allocate(   X0(nlevs,0:nr_fine-1))
  allocate( edge(nlevs,0:nr_fine))


  gam1(:,:) = ZERO
  w0(:,:) = ZERO
  psi(:,:) = ZERO
  etarho_ec(:,:) = ZERO
  etarho_cc(:,:) = ZERO
  div_etarho(:,:) = ZERO
  p0_minus_pthermbar(:,:) = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do n = 1,nlevs
     call init_base_state(n,model_file,s0(n,:,:),p0(n,:),dx(n,:))
  enddo


  ! output
  open(unit=10,file="base.orig")
  do i = 0, nr_fine-1
     write(10,1000) r_cc_loc(1,i), s0(1,i,rho_comp), s0(1,i,temp_comp), p0(1,i)
  enddo
  close(unit=10)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main timestepping loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  time = ZERO
  dt = 0.0001_dp_t
  dtold = dt

  w0_old(:,:) = ZERO

  iter = 0

  allocate(anelastic_cutoff_coord(1))
  allocate(base_cutoff_density_coord(1))

  do while (time < stop_time)

     print *, 'time = ', time

     ! compute the anelastic cutoff
     anelastic_cutoff_coord(1) = nr_fine-1
     do i = 0, nr_fine-1
        if (s0(1,i,rho_comp) .le. anelastic_cutoff .and. &
             anelastic_cutoff_coord(1) .eq. nr_fine-1) then
           anelastic_cutoff_coord(1) = i
           exit
        endif
     enddo


     ! compute the low density cutoff
     base_cutoff_density_coord(1) = nr_fine-1
     do i = 0, nr_fine-1
        if (s0(1,i,rho_comp) .le. base_cutoff_density .and. &
             base_cutoff_density_coord(1) .eq. nr_fine-1) then
           base_cutoff_density_coord(1) = i
           exit
        endif
     enddo
     

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the heating term
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     y_0 = 4.d7

     print *, 'calling the eos', nr_fine
     do i = 0, nr_fine-1

        Hbar = 1.d16 * exp(-((r_cc_loc(1,i) - y_0)**2)/ 1.d14)
     
        ! (rho, T) --> p,h, etc
        den_eos(1)  = s0(1,i,rho_comp)
        temp_eos(1) = s0(1,i,temp_comp)
        xn_eos(1,:) = s0(1,i,spec_comp:spec_comp-1+nspec)/s0(1,i,rho_comp)
        
        !print *, 'calling EOS: ', i, den_eos(1), temp_eos(1), xn_eos(1,:)
        call eos(eos_input_rt, den_eos, temp_eos, NP, nspec, &
                 xn_eos, &
                 p_eos, h_eos, e_eos, &
                 cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                 dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                 dpdX_eos, dhdX_eos, &
                 gam1_eos, cs_eos, s_eos, &
                 dsdt_eos, dsdr_eos, &
                 do_diag)

        ! in the real Maestro code, we are updating the enthalpy by differencing
        ! the enthalpy equation with the heating term, rather than using the EOS.
        s0(1,i,rhoh_comp) = den_eos(1)*h_eos(1)

        p0(1,i) = p_eos(1)
        gam1(1,i) = gam1_eos(1)

        coeff = dpdt_eos(1)/ (den_eos(1) * cp_eos(1) * dpdr_eos(1))

        Sbar_in(1,i) = coeff*Hbar
     enddo


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute w_0
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     w0(:,:) = ZERO

     call make_w0(w0,w0_old,f,Sbar_in,s0(:,:,rho_comp),s0(:,:,rho_comp),p0,p0, &
                  gam1,gam1,p0_minus_pthermbar,psi,etarho_ec,etarho_cc,dt,dtold)
  

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the divergance coefficient (beta_0)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do n=1,nlevs
        call make_grav_cell(n,grav_cell(n,:),s0(n,:,rho_comp))
     enddo

     call make_div_coeff(div_coeff,s0(:,:,rho_comp),p0,gam1,grav_cell)     


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the new base state
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! here we set the old state to the current state -- advect_base will 
     ! update the state
     s0_old(:,:,:) = s0(:,:,:)
     p0_old(:,:) = p0(:,:)

     which_step = 1

     
     call advect_base_dens(w0,s0_old(:,:,rho_comp),s0(:,:,rho_comp), &
                           rho0_predicted_edge,dx(:,1),dt)


     call advect_base_pres(w0,Sbar_in,s0(:,:,rho_comp),p0_old,p0,gam1,div_coeff, &
                             psi,etarho_cc,dx(:,1),dt)


     call advect_base_enthalpy(w0,Sbar_in,s0_old(:,:,rho_comp), &
                               s0_old(:,:,rhoh_comp),s0(:,:,rhoh_comp), &
                               p0_old,p0,gam1,rho0_predicted_edge,psi,dx(:,1),dt)


     ! update the species.  In the real code, this will be done
     ! for the full state.  Here we are faking a 1-d star, so
     ! we need to do this manually, since the species are not
     ! part of the base state
     do comp = spec_comp,spec_comp+nspec-1

        ! here we predict X_0 on the edges
        X0(1,:) = s0_old(1,:,comp)/s0_old(1,:,rho_comp)
        do i = 0,nr_fine-1
           X0(1,i) = max(X0(1,i),ZERO)
        end do

        force = ZERO

        if (spherical .eq. 0) then
           call make_edge_state_1d(X0,edge,w0,force,dx(:,1),dt)
        else
           call make_edge_state_1d(X0,edge,w0,force,dr,dt)
        endif

        ! our final update needs (rho X)_0 on the edges, so compute
        ! that now
        edge(1,:) = rho0_predicted_edge(1,:)*edge(1,:)

        ! update (rho X)_0
        do i = 0,nr_fine-1
           s0(1,i,comp) = s0_old(1,i,comp) &
                - (dt/dr(1))/r_cc_loc(1,i)**2* &
                (r_edge_loc(1,i+1)**2 * edge(1,i+1) * w0(1,i+1) - &
                 r_edge_loc(1,i  )**2 * edge(1,i  ) * w0(1,i  ))
        end do

     enddo

     ! don't let the species leave here negative
     do comp = spec_comp, spec_comp+nspec-1
        
        if (minval(s0(1,:,comp)) .lt. ZERO) then
           do i = 0, nr_fine-1
              if (s0(1,i,comp) .lt. ZERO) then
                 delta = -s0(1,i,comp)
                 sum = ZERO
                 do comp2 = spec_comp, spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s0(1,i,comp2) .ge. ZERO) then
                       sum = sum + s0(1,i,comp2)
                    endif
                 enddo
                 do comp2 = spec_comp, spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s0(1,i,comp2) .ge. ZERO) then
                       frac = s0(1,i,comp2) / sum
                       s0(1,i,comp2) = s0(1,i,comp2) - frac * delta
                    endif
                 enddo
                 s0(1,i,comp) = ZERO
              endif
           enddo
        endif
     enddo
             


     ! update the temperature -- advect base does not do this
     do i = 0, nr_fine-1

        ! (rho, T) --> p,h, etc
        den_eos(1)  = s0(1,i,rho_comp)
        temp_eos(1) = s0(1,i,temp_comp)
        p_eos(1)    = p0(1,i)
        xn_eos(1,:) = s0(1,i,spec_comp:spec_comp-1+nspec)/s0(1,i,rho_comp)
        

        !print *, i, den_eos(1), temp_eos(1), p_eos(1), xn_eos(1,:)

        call eos(eos_input_rp, den_eos, temp_eos, NP, nspec, &
                 xn_eos, &
                 p_eos, h_eos, e_eos, &
                 cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                 dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                 dpdX_eos, dhdX_eos, &
                 gam1_eos, cs_eos, s_eos, &
                 dsdt_eos, dsdr_eos, &
                 do_diag)

        s0(1,i,temp_comp) = temp_eos(1)
     enddo

     print *, 'new base pressure', p0(1,1)
     print *, 'new base density', s0(1,1,rho_comp)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the new timestep
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     time = time + dt
     dtold = dt

     ! output
     !write(unit=base_state_name,fmt='("base_",i5.5)') iter
     !open(unit=10,file=base_state_name)
     !do i = 0, nr_fine-1
     !   write(10,1000) r_cc_loc(1,i), s0(1,i,rho_comp), s0(1,i,temp_comp), p0(1,i), w0(1,i)
     !enddo
     !close(unit=10)

     dt = min(1.1*dt,cflfac*dr_fine/maxval(abs(w0)))
     if (time+dt > stop_time) dt = stop_time - time

     ! store the old velocity
     w0_old(:,:) = w0(:,:)

     iter = iter + 1
  enddo

  ! output
  open(unit=10,file="base.new")
  do i = 0, nr_fine-1
     write(10,1000) r_cc_loc(1,i), s0(1,i,rho_comp), s0(1,i,temp_comp), p0(1,i), w0(1,i)
  enddo
  close(unit=10)
1000 format(1x,6(g20.10))

  deallocate(div_coeff_old,div_coeff,grav_cell)
  deallocate(gam1,s0_old,s0,p0_old,p0,w0,f)
  deallocate(rho0_predicted_edge)

end subroutine varden
