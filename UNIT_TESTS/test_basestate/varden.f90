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

  implicit none

  integer :: i, comp,comp2
  integer :: iter

  real(dp_t) :: frac, delta, sum
  integer    :: dm

  real(dp_t) :: time,dt,half_time,dtold
  real(dp_t) :: smin,smax

  integer :: un, ierr

  real(dp_t) :: y_0

  real(dp_t) :: dx(1,1)
  real(dp_t) :: prob_lo(1), prob_hi(1)

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
  real(dp_t), allocatable :: f(:,:)
  real(dp_t), allocatable :: Sbar_in(:,:)
  real(dp_t), allocatable :: rho0_predicted_edge(:,:)
  real(dp_t), allocatable :: force(:)
  real(dp_t), allocatable :: X0(:)
  real(dp_t), allocatable :: edge(:)

  real(dp_t) :: coeff, Hbar

  integer :: nlevs,n
  integer :: nr_fine
  integer :: which_step

  real(dp_t) :: max_dist

  character (len=10) base_state_name

  call probin_init()

  nlevs = 1

  dm = 1

  call init_spherical(spherical_in)

  center(1) = ZERO


  call init_variables(dm, nspec)
  call network_init()
  call eos_init(use_eos_coulomb=use_eos_coulomb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define the grid spacing on all levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  prob_lo(1) = prob_lo_x
  prob_hi(1) = prob_hi_x



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate storage for the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (spherical .eq. 1) then
    if (dr_base .gt. 0) then
      max_dist = prob_hi_x - prob_lo_x
      nr_fine = int(max_dist / dr_base) + 1
      if ( parallel_IOProcessor() ) then
         print *,'DISTANCE FROM CENTER TO CORNER IS ',max_dist
         print *,'DR_BASE IS ',dr_base
         print *,'SETTING nr_fine TO ',nr_fine
      end if
    else
     if ( parallel_IOProcessor() ) &
       print *,'NEED TO DEFINE DR_BASE '
      stop
    endif
  else
     ! NOTE: in the basestate test, we will always use dr_base as the
     ! input
     max_dist = prob_hi_x - prob_lo_x
     nr_fine = int(max_dist / dr_base) + 1
  end if

  dx(1,1) = dr_base

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
  allocate(                psi(nlevs,0:nr_fine))
  allocate(                  f(nlevs,0:nr_fine))
  allocate(            Sbar_in(nlevs,0:nr_fine-1))
  allocate(rho0_predicted_edge(nlevs,0:nr_fine))

  allocate(force(0:nr_fine-1))
  allocate(   X0(0:nr_fine-1))
  allocate(edge(0:nr_fine))

  gam1(:,:) = ZERO
  w0(:,:) = ZERO
  psi(:,:) = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call init_geometry(center,nr_fine,dr_base)

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
  do while (time < stop_time)

     print *, 'time = ', time

     ! compute the anelastic cutoff
     anelastic_cutoff_coord(1) = nr_fine
     do i = 0, nr_fine-1
        if (s0(1,i,rho_comp) .le. anelastic_cutoff .and. &
             anelastic_cutoff_coord(1) .eq. nr_fine-1) then
           anelastic_cutoff_coord(1) = i
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

     call make_w0(nlevs,w0,w0_old,f,Sbar_in,s0(:,:,rho_comp),p0,p0, &
                  gam1,gam1,psi,dt,dtold)
  

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the divergance coefficient (beta_0)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do n=1,nlevs
        call make_grav_cell(n,grav_cell(n,:),s0(n,:,rho_comp))
        call make_div_coeff(n,div_coeff(n,:),s0(n,:,rho_comp),p0(n,:), &
                            gam1(n,:),grav_cell(n,:))     
     enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the new base state
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! here we set the old state to the current state -- advect_base will 
     ! update the state
     s0_old(:,:,:) = s0(:,:,:)
     p0_old(:,:) = p0(:,:)

     which_step = 1

     call advect_base(nlevs,w0,Sbar_in,p0_old,p0, &
                      s0_old(:,:,rho_comp),s0(:,:,rho_comp), &
                      s0_old(:,:,rhoh_comp),s0(:,:,rhoh_comp), &
                      gam1,div_coeff, &
                      rho0_predicted_edge,psi, &
                      dx(:,1),dt)


     ! update the species.  In the real code, this will be done
     ! for the full state.  Here we are faking a 1-d star, so
     ! we need to do this manually, since the species are not
     ! part of the base state
     do comp = spec_comp,spec_comp+nspec-1

        ! here we predict X_0 on the edges
        X0(:) = s0_old(1,:,comp)/s0_old(1,:,rho_comp)
        do i = 0,nr_fine-1
           X0(i) = max(X0(i),ZERO)
        end do

        force = ZERO

        call make_edge_state_1d(1,X0,edge,w0(1,:),force,1,dr(1),dt)
        
        ! our final update needs (rho X)_0 on the edges, so compute
        ! that now
        edge(:) = rho0_predicted_edge(1,:)*edge(:)

        ! update (rho X)_0
        do i = 0,nr_fine-1
           s0(1,i,comp) = s0_old(1,i,comp) &
                - (dt/dr(1))/r_cc_loc(1,i)**2* &
                (r_edge_loc(1,i+1)**2 * edge(i+1) * w0(1,i+1) - &
                 r_edge_loc(1,i  )**2 * edge(i  ) * w0(1,i  ))
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

     dt = min(1.1*dt,cflfac*dr_base/maxval(abs(w0)))
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
