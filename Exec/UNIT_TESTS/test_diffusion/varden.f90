! THIS PROBLEM DOESN"T SUPPORT MULTILEVEL
!
subroutine varden()

  use variables
  use network
  use geometry
  use base_state_module
  use fabio_module
  use probin_module
  use runtime_init_module
  use bl_constants_module
  use define_bc_module
  use eos_module
  use rhoh_vs_t_module
  use initialize_module
  use multifab_module
  use make_explicit_thermal_module
  use thermal_conduct_module
  use estdt_module
  use conductivity_module, only: conductivity_init
  use bl_IO_module

  implicit none

  integer :: istep
  integer :: i, n, r, comp, nlevs, dm

  type(ml_layout)   :: mla
  type(bc_tower)    :: the_bc_tower

  real(kind=dp_t) :: time, dt

  type(multifab),  pointer :: s_old(:), thermal(:)
  type(multifab), allocatable :: s_new(:)
  type(multifab), allocatable :: Tcoeff1(:), hcoeff1(:), Xkcoeff1(:), pcoeff1(:)
  type(multifab), allocatable :: Tcoeff2(:), hcoeff2(:), Xkcoeff2(:), pcoeff2(:)
  type(multifab), allocatable :: analytic(:), error(:)

  real(kind=dp_t) :: L1norm, L2norm

  real(kind=dp_t), pointer :: p0(:,:), s0_init(:,:,:)
  real(kind=dp_t), pointer :: dx(:,:)

  character(len=20), allocatable :: names(:)
  character(len=256) :: outdir, sstep

  real(kind=dp_t), parameter :: SMALL = 1.e-13
  
  real(kind=dp_t) :: diffusion_coefficient

  integer :: uout

  ! initialize some things
  call runtime_init()
  call init_spherical()
  call init_center()

  call init_variables()

  call network_init()
  call eos_init()

  call conductivity_init()

  time = ZERO

  ! we are hardcoded for a certain net here
  if (nspec /= 4) then
     call bl_error("ERROR: wrong network")
  endif

  ! setup some names for the data fab's
  allocate(names(nscal))
  names(1) = "density"
  names(2) = "rhoh"
  names(3) = "X(" // trim(short_spec_names(1)) // ")*rho"
  names(4) = "X(" // trim(short_spec_names(2)) // ")*rho"
  names(5) = "X(" // trim(short_spec_names(3)) // ")*rho"
  names(6) = "X(" // trim(short_spec_names(4)) // ")*rho"
  names(7) = "temp"
  names(8) = "trac"

 ! we only use fixed grids
 if (test_set /= '') then
    call initialize_with_fixed_grids(mla, time, dt, dx, pmask, the_bc_tower, &
                                     s_old, s0_init, p0, thermal, &
                                     diffusion_coefficient)
 else
    call bl_error('test_diffusion only implemented with fixed grids!')
 endif

 dm = mla%dim
 nlevs = mla%nlevel

 if (nlevs .ne. max_levs) then
    call bl_error('varden.f90: nlevs .ne. max_levs not supported yet')
 end if

 ! error checking
 if (.not. use_thermal_diffusion) then
    call bl_error('use_thermal_diffusion = F')
 endif

 if (dm .ne. get_dim(mla%mba) .or. dm .ne. 2) then
    call bl_error('dm_in not properly set in inputs file')
 endif

 if (spherical .eq. 1 ) then
    call bl_error('spherical = 1 not supported')
 endif

 if (abs(dx(1,1) - dx(1,2)) > SMALL) then
    call bl_error('zones must be square')
 endif

 if (nlevs > 1) call bl_error("code not properly set up for multilevel")

 if (parallel_IOProcessor()) then
    print *, 'number of processors = ', parallel_nprocs()
    print *, 'number of dimensions = ', dm
    do n = 1, nlevs
       print *, 'level: ', n
       print *, '   number of boxes = ', nboxes(s_old(n)%la)
       print *, '   maximum zones   = ', (extent(mla%mba%pd(n),i),i=1,dm)
    end do
    print *, ''
 end if

 ! initialize remaining arrays
 allocate(s_new(nlevs),&
          Tcoeff1(nlevs),hcoeff1(nlevs),Xkcoeff1(nlevs),pcoeff1(nlevs),&
          Tcoeff2(nlevs),hcoeff2(nlevs),Xkcoeff2(nlevs),pcoeff2(nlevs),&
          analytic(nlevs), error(nlevs))

 do n = 1, nlevs
    ! new state 
    call build(s_new(n),   mla%la(n), nscal, s_old(n)%ng)

    ! coefficients for diffusion
    call build(Tcoeff1(n),  mla%la(n),     1,           1)
    call build(hcoeff1(n),  mla%la(n),     1,           1)
    call build(Xkcoeff1(n), mla%la(n), nspec,           1)
    call build(pcoeff1(n),  mla%la(n),     1,           1)
    call build(Tcoeff2(n),  mla%la(n),     1,           1)
    call build(hcoeff2(n),  mla%la(n),     1,           1)
    call build(Xkcoeff2(n), mla%la(n), nspec,           1)
    call build(pcoeff2(n),  mla%la(n),     1,           1)

    ! variables for calculating the norms
    call build(analytic(n), mla%la(n),     1,           0)
    call build(error(n),    mla%la(n),     1,           0)


    ! for now just copy the state data to s_new
    call multifab_copy_c(s_new(n), 1, s_old(n), 1, nscal, s_old(n)%ng)
    
    ! calculate the error
    ! build the analytic solution
    call make_analytic_solution(analytic(n), dx(n,:), time)

    ! the error is (enthalpy - analytic)
    ! build the enthalpy from the state data
    call multifab_copy_c(error(n), 1, s_new(n), rhoh_comp, 1, 0)
    call multifab_div_div_c(error(n), 1, s_new(n), rho_comp, 1, 0)

    ! subtract the analytic solution from the enthalpy
    call multifab_sub_sub(error(n), analytic(n))

    ! compute the norms
!    L1norm = multifab_norm_l1(error(n))/real(nr_fine**2,kind=dp_t)
!    L2norm = multifab_norm_l2(error(n))/real(nr_fine**2,kind=dp_t)
    L1norm = multifab_norm_l1(error(n))/multifab_norm_l1(analytic(n))
    L2norm = multifab_norm_l2(error(n))/multifab_norm_l2(analytic(n))

 enddo

 if (parallel_IOProcessor()) then
    uout = unit_new()
    open(unit=uout,file=outputfile,status='replace')
    write(uout,*) time, L1norm, L2norm
    print *, time, L1norm, L2norm
 endif

 ! calculate the timestep
 ! use fixed dt
 if (fixed_dt > ZERO) then
    dt = fixed_dt
 else ! use the explicit timestep size
    call estdt(dx, diffusion_coefficient, dt)
 endif

 if (parallel_IOProcessor()) then
    print *, '... estdt gives dt =', dt    
 endif

 dt = dt * dt_mult_factor
    
 if (parallel_IOProcessor()) &
      print *, '... dt after applying mult_factor:', dt

 istep = 0

 ! dump the initial data
 write(unit=sstep,fmt='(i5.5)')istep
 outdir = trim(plot_base_name) // sstep
 ! we only give dx at coarsest level for now
 call fabio_ml_write(s_old, mla%mba%rr(:,1), trim(outdir), &
                     names=names, time=time, &
                     prob_lo=prob_lo, prob_hi=prob_hi, &
                     dx=dx(1,:))
 call fabio_ml_write(analytic, mla%mba%rr(:,1), 'analytic'//sstep, &
                     time=time, prob_lo=prob_lo, prob_hi=prob_hi, dx=dx(1,:))
 call fabio_ml_write(error, mla%mba%rr(:,1), 'error'//sstep, &
                     time=time, prob_lo=prob_lo, prob_hi=prob_hi, dx=dx(1,:))

 ! loop
 do while (time < stop_time)

    istep = istep+1

    if (parallel_IOProcessor()) print *, 'Working on step', istep


    ! build the coeffs
    if (parallel_IOProcessor()) print *, '... building thermal coefficients'
    call make_thermal_coeffs(s_old,Tcoeff1,hcoeff1,Xkcoeff1,pcoeff1)

    ! on the first step, just copy coeffs for the time centering
    if (istep==1) then
       do n=1,nlevs
          call multifab_copy_c(Tcoeff2(n),  1, Tcoeff1(n),  1, 1,     1)
          call multifab_copy_c(hcoeff2(n),  1, hcoeff1(n),  1, 1,     1)
          call multifab_copy_c(Xkcoeff2(n), 1, Xkcoeff1(n), 1, nspec, 1)
          call multifab_copy_c(pcoeff2(n),  1, pcoeff1(n),  1, 1,     1)
       enddo
    endif


    ! diffuse the enthalpy
    if (parallel_IOProcessor()) print *, '... conducting'
    call thermal_conduct(mla, dx, dt, s_old, hcoeff1, Xkcoeff1, pcoeff1, &
                         hcoeff2, Xkcoeff2, pcoeff2, s_new, p0, p0, &
                         the_bc_tower)

    ! update temperature
    call makeTfromRhoH(s_new,p0,mla,the_bc_tower%bc_tower_array,dx)


    ! solution advanced -- update the time
    time = time + dt

    ! calculate the error
    call make_analytic_solution(analytic(nlevs), dx(nlevs,:), time)

    call multifab_copy_c(error(nlevs), 1, s_new(nlevs), rhoh_comp, 1, 0)
    call multifab_div_div_c(error(nlevs), 1, s_new(nlevs), rho_comp, 1, 0)

    call multifab_sub_sub(error(nlevs), analytic(nlevs))

!    L1norm = multifab_norm_l1(error(nlevs))/real(nr_fine**2,kind=dp_t)
!    L2norm = multifab_norm_l2(error(nlevs))/real(nr_fine**2,kind=dp_t)
    L1norm = multifab_norm_l1(error(nlevs))/multifab_norm_l1(analytic(nlevs))
    L2norm = multifab_norm_l2(error(nlevs))/multifab_norm_l2(analytic(nlevs))

    if (parallel_IOProcessor()) then
       write(uout,*) time, L1norm, L2norm
       print *, time, L1norm, L2norm
    endif

    ! dump data if needed
    if (mod(istep,plot_int) .eq. 0) then

       ! dump the analytic comparison
!       if (parallel_IOProcessor()) call dump_gnuplot_analysis(istep,time)

       write(unit=sstep,fmt='(i5.5)')istep
       outdir = trim(plot_base_name) // sstep

       if (parallel_IOProcessor()) print *, '... writing to ', trim(outdir)

       ! we only give dx at the coarsest level for now
       call fabio_ml_write(s_new, mla%mba%rr(:,1), trim(outdir), &
                           names=names, time=time, &
                           prob_lo=prob_lo, prob_hi=prob_hi, &
                           dx=dx(1,:))
       call fabio_ml_write(analytic, mla%mba%rr(:,1), 'analytic'//sstep, &
                           time=time, prob_lo=prob_lo, prob_hi=prob_hi, &
                           dx=dx(1,:))
       call fabio_ml_write(error, mla%mba%rr(:,1), 'error'//sstep, &
                           time=time, prob_lo=prob_lo, prob_hi=prob_hi, &
                           dx=dx(1,:))
    endif


    ! fill the mfs for the next timestep
     do n = 1, nlevs
        call multifab_copy_c(s_old(n), 1, s_new(n), 1, nscal, s_old(n)%ng)

        call multifab_copy_c(Tcoeff2(n),  1, Tcoeff1(n),  1, 1,     1)
        call multifab_copy_c(hcoeff2(n),  1, hcoeff1(n),  1, 1,     1)
        call multifab_copy_c(Xkcoeff2(n), 1, Xkcoeff1(n), 1, nspec, 1)
        call multifab_copy_c(pcoeff2(n),  1, pcoeff1(n),  1, 1,     1)
     enddo

    if (parallel_IOProcessor()) print *, ''

    if (time + dt > stop_time) then
       dt = stop_time - time
    endif

    if (parallel_IOProcessor()) print *, '... time = ', time

 enddo

 close(uout)

contains
  
  subroutine make_analytic_solution(solution, dx, time)

    type(multifab),  intent(inout) :: solution
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: time

    real(kind=dp_t), pointer :: sp(:,:,:,:)

    integer :: lo(dm), hi(dm)
    integer :: n, i, j, ii, jj
    real(kind=dp_t) :: xx, yy, dist
    
    do i = 1, nfabs(solution)
       sp => dataptr(solution,i)
       lo = lwb(get_box(solution,i))
       hi = upb(get_box(solution,i))

       sp = ZERO

       do jj = lo(2), hi(2)

          yy = prob_lo(2) + (dble(jj)+HALF) * dx(2)

          do ii = lo(1), hi(1)

             xx = prob_lo(1) + (dble(ii)+HALF) * dx(1)

             dist = sqrt((xx-center(1))**2 + (yy-center(2))**2)

             sp(ii,jj,1,1) = f(time,dist)


          enddo
       enddo
    enddo

  end subroutine make_analytic_solution
    
  function f(t,x) result(r)
    real(kind=dp_t) :: t, x
    real(kind=dp_t) :: r

    r = (peak_h-ambient_h)*(t0/(t+t0)) * &
         dexp(-x*x/(FOUR*diffusion_coefficient*(t+t0))) + ambient_h

  end function f

end subroutine varden



  

