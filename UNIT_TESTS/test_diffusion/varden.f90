subroutine varden()

  use variables
  use network
  use geometry
  use base_state_module
  use fabio_module
  use probin_module
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

  implicit none

  integer :: istep
  integer :: i, n, r, comp

  type(ml_layout)   :: mla
  type(bc_tower)    :: the_bc_tower

  real(kind=dp_t) :: time, dt

  type(multifab),  pointer :: s_old(:), thermal(:)
  type(multifab), allocatable :: s_new(:)
  type(multifab), allocatable :: Tcoeff1(:), hcoeff1(:), Xkcoeff1(:), pcoeff1(:)
  type(multifab), allocatable :: Tcoeff2(:), hcoeff2(:), Xkcoeff2(:), pcoeff2(:)

  real(kind=dp_t), pointer :: p0(:,:), s0_init(:,:,:)
  real(kind=dp_t), pointer :: dx(:,:)

  character(len=20), allocatable :: names(:)
  character(len=256) :: outdir, sstep

  real(kind=dp_t), parameter :: SMALL = 1.e-13
  real(kind=dp_t), parameter :: FIVE3RD = FIVE/THREE 
  
  real(kind=dp_t) :: diffusion_coefficient

  ! initialize some things
  call probin_init()
  call init_dm()
  call init_spherical()
  call init_center()

  call init_variables()

  call network_init()
  call eos_init(use_eos_coulomb=use_eos_coulomb,gamma_in=FIVE3RD)

  call conductivity_init(cond_const=thermal_conductivity)

  ! setup some names for the data fab's
  allocate(names(nscal))
  names(1) = "density"
  names(2) = "rhoh"
  names(3) = "X(He4)*rho"
  names(4) = "X(C12)*rho"
  names(5) = "X(Fe56)*rho"
  names(6) = "temp"
  names(7) = "trac"

 ! we only use fixed grids
 if (test_set /= '') then
    call initialize_with_fixed_grids(mla, time, dt, dx, pmask, the_bc_tower, &
                                     s_old, s0_init, p0, thermal, &
                                     diffusion_coefficient)
 else
    call bl_error('test_diffusion only implemented with fixed grids!')
 endif

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

 if (parallel_IOProcessor()) then
    print *, 'number of processors = ', parallel_nprocs()
    print *, 'number of dimensions = ', dm
    do n = 1, nlevs
       print *, 'level: ', n
       print *, '   number of boxes = ', s_old(n)%nboxes
       print *, '   maximum zones   = ', (extent(mla%mba%pd(n),i),i=1,dm)
    end do
    print *, ''
 end if

 ! initialize remaining arrays
 allocate(s_new(nlevs),&
          Tcoeff1(nlevs),hcoeff1(nlevs),Xkcoeff1(nlevs),pcoeff1(nlevs),&
          Tcoeff2(nlevs),hcoeff2(nlevs),Xkcoeff2(nlevs),pcoeff2(nlevs))

 do n = 1, nlevs
    call build(s_new(n),   mla%la(n), nscal, s_old(n)%ng)

    call build(Tcoeff1(n),  mla%la(n),     1,           1)
    call build(hcoeff1(n),  mla%la(n),     1,           1)
    call build(Xkcoeff1(n), mla%la(n), nspec,           1)
    call build(pcoeff1(n),  mla%la(n),     1,           1)
    call build(Tcoeff2(n),  mla%la(n),     1,           1)
    call build(hcoeff2(n),  mla%la(n),     1,           1)
    call build(Xkcoeff2(n), mla%la(n), nspec,           1)
    call build(pcoeff2(n),  mla%la(n),     1,           1)


    ! for now just copy the state data to s_new
    call multifab_copy_c(s_new(n), 1, s_old(n), 1, nscal, s_old(n)%ng)
 enddo

 ! calculate the explicit timestep
 call estdt(dx, diffusion_coefficient, dt)
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
                     problo=prob_lo, probhi=prob_hi, &
                     dx=dx(1,:))

 ! dump some parameters for the analytic comparison
 if (parallel_IOProcessor()) call dump_gnuplot_analysis(istep,time)

 ! loop
 do while (istep < max_step)

    istep = istep+1

    if (parallel_IOProcessor()) print *, 'Working on step', istep

    time = time + dt
    if (parallel_IOProcessor()) print *, '... time = ', time

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
    call makeTfromRhoH(s_new,mla,the_bc_tower%bc_tower_array)

    ! dump data if needed
    if (mod(istep,plot_int) .eq. 0) then

       ! dump the analytic comparison
       if (parallel_IOProcessor()) call dump_gnuplot_analysis(istep,time)

       write(unit=sstep,fmt='(i5.5)')istep
       outdir = trim(plot_base_name) // sstep

       if (parallel_IOProcessor()) print *, '... writing to ', trim(outdir)

       ! we only give dx at the coarsest level for now
       call fabio_ml_write(s_new, mla%mba%rr(:,1), trim(outdir), &
                           names=names, time=time, &
                           problo=prob_lo, probhi=prob_hi, &
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
 enddo

contains
  
  subroutine dump_gnuplot_analysis(istep,time)

    use bl_IO_module
    
    implicit none

    integer        , intent(in   ) :: istep
    real(kind=dp_t), intent(in   ) :: time
    
!    character(len=10) :: outputfile = "params.out"

    integer :: unit

    ! dump info to the outputfile
    unit = unit_new()
    if (istep == 0) then
       open(unit=unit, file=gnuplot_outputfile, status='replace')

!       write(unit,100) "t0", t0
!       write(unit,100) "D", diffusion_coefficient
!       write(unit,100) "hp", peak_h
!       write(unit,100) "h0", ambient_h
       write(unit,*) t0
       write(unit,*) diffusion_coefficient
       write(unit,*) peak_h
       write(unit,*) ambient_h

    else
       open(unit=unit, file=gnuplot_outputfile,status='old',position='append')
    endif

!    write(unit,101) istep,time,time
    write(unit,*) time

    close(unit)

100 format (A,"=",G15.10)
101 format ("f",I5.5,"(x)=(hp-h0)*(t0/(",G15.10,"+t0))*exp(-x*x/(4.0e0*D*(",G15.10,"+t0)))+h0")

  end subroutine dump_gnuplot_analysis
    

    
 
end subroutine varden



  

