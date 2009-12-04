subroutine varden()

  use variables
  use network
  use geometry
  use init_module
  use base_state_module
  use fabio_module
  use probin_module
  use bl_constants_module
  use bl_IO_module
  use average_module
  use define_bc_module
  use fill_3d_module
  use eos_module
  use rhoh_vs_t_module
  use initialize_module
  use box_util_module
  use multifab_module
  use ml_multifab_module
  use make_explicit_thermal_module
  use make_plotfile_module
  use thermal_conduct_module

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
  type(multifab), allocatable :: plot_coeffs(:)

  real(kind=dp_t), pointer :: p0(:,:), rho0(:,:), s0_init(:,:,:)
  real(kind=dp_t), pointer :: dx(:,:)

  real(kind=dp_t) :: dist, lenx, leny, lenz, max_dist

  character(len=20), allocatable :: names(:), coeff_names(:)
  character(len=20) :: outdir, sstep

  real(kind=dp_t), parameter :: SMALL = 1.e-13
  real(kind=dp_t), parameter :: FIVE3RD = FIVE/THREE 

  ! initialize some things
  call probin_init()
  call init_dm()
  call init_spherical()
  call init_center()

  call init_variables()

  call network_init()
  call eos_init(use_eos_coulomb=use_eos_coulomb,gamma_in=FIVE3RD)

  ! setup some names for the data fab's
  allocate(names(nscal),coeff_names(3+nspec))
  names(:) = (/ "density", "rhoh", &
       "X(He4)*rho", "X(C12)*rho", "X(Fe56)*rho", &
       "temp", "trac" /)
  coeff_names(:) = (/ "Tcoeff", &
                      "hcoeff", &
                      "He4 coeff", &
                      "C12 coeff", &
                      "Fe56 coeff", &
                      "pcoeff" /)

 ! for now we only use fixed grids
 if (test_set /= '') then
    call initialize_with_fixed_grids(mla, time, dt, dx, pmask, the_bc_tower, &
                                     s_old, s0_init, p0, rho0, thermal)
 else
    call bl_error('Currently only implemented with fixed grids!')
 endif

 ! error checking
 if (.not. use_thermal_diffusion) then
    call bl_error('use_thermal_diffusion = F')
 endif

 if (dm .ne. get_dim(mla%mba) .or. dm .ne. 2) then
    call bl_error('dm_in not properly set in inputs file')
 endif

 if (spherical .eq. 1 ) then
    call bl_error('sperical = 1 not supported')
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
          Tcoeff2(nlevs),hcoeff2(nlevs),Xkcoeff2(nlevs),pcoeff2(nlevs),&
          plot_coeffs(nlevs))
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

    call build(plot_coeffs(n), mla%la(n), 3+nspec, 0)

    ! for now just copy the state data to s_new
    call multifab_copy_c(s_new(n), 1, s_old(n), 1, nscal, s_old(n)%ng)
 enddo


 istep = 0

 ! dump the initial data
 write(unit=sstep,fmt='(i5.5)')istep
 outdir = "state_" // sstep
 call fabio_ml_write(s_old, mla%mba%rr(:,1), trim(outdir), &
                     names=names, time=time)


 ! loop
 do while (istep < max_step)

    istep = istep+1

    if (parallel_IOProcessor()) print *, 'Working on step', istep

    ! get the timestep
    call make_explicit_thermal_dt(mla,the_bc_tower,dx,s_old,dt)

    if (parallel_IOProcessor()) print *, '... dt =', dt

    time = time + dt

    ! build the coeffs
    if (parallel_IOProcessor()) print *, '... building thermal coefficients'
    call make_thermal_coeffs(s_old,Tcoeff1,hcoeff1,Xkcoeff1,pcoeff1)

    ! on the first step, just copy coeffs for the time centering
    if (istep==1) then
       do n=1,nlevs
          call multifab_copy_c(Tcoeff2(n), 1, Tcoeff1(n), 1, 1)
          call multifab_copy_c(hcoeff2(n), 1, hcoeff1(n), 1, 1)
          call multifab_copy_c(Xkcoeff2(n), 1, Xkcoeff1(n), 1, nspec)
          call multifab_copy_c(pcoeff2(n), 1, pcoeff1(n), 1, 1)
       enddo
    endif

    ! throw the coeffs1 to a mfab for output 
    do n = 1, nlevs
       call multifab_copy_c(plot_coeffs(n), 1, Tcoeff1(n),1,1)
       call multifab_copy_c(plot_coeffs(n), 2, hcoeff1(n),1,1)
       call multifab_copy_c(plot_coeffs(n), 3, Xkcoeff1(n),1,nspec)
       call multifab_copy_c(plot_coeffs(n), 3+nspec-1, pcoeff1(n),1,1)
    enddo

    write(unit=sstep,fmt='(i5.5)')istep
    outdir = "coeffs1_" // sstep

    if (parallel_IOProcessor()) print *, '... writing to ', outdir

    call fabio_ml_write(plot_coeffs, mla%mba%rr(:,1), trim(outdir), &
                        names=coeff_names, &
                        time=time)

    ! throw the coeffs2 to a mfab for output
    do n = 1, nlevs
       call multifab_copy_c(plot_coeffs(n), 1, Tcoeff2(n),1,1)
       call multifab_copy_c(plot_coeffs(n), 2, hcoeff2(n),1,1)
       call multifab_copy_c(plot_coeffs(n), 3, Xkcoeff2(n),1,nspec)
       call multifab_copy_c(plot_coeffs(n), 3+nspec-1, pcoeff2(n),1,1)
    enddo

    write(unit=sstep,fmt='(i5.5)')istep
    outdir = "coeffs2_" // sstep

    if (parallel_IOProcessor()) print *, '... writing to ', outdir

    call fabio_ml_write(plot_coeffs, mla%mba%rr(:,1), trim(outdir), &
                        names=coeff_names, &
                        time=time)

    ! conduct
    if (parallel_IOProcessor()) print *, '... conducting'
    call thermal_conduct(mla, dx, dt, s_old, hcoeff1, Xkcoeff1, pcoeff1, &
                         hcoeff2, Xkcoeff2, pcoeff2, s_new, p0, p0, &
                         the_bc_tower)

    ! dump data
    write(unit=sstep,fmt='(i5.5)')istep
    outdir = "state_" // sstep

    if (parallel_IOProcessor()) print *, '... writing to ', outdir

    call fabio_ml_write(s_new, mla%mba%rr(:,1), trim(outdir), &
                        names=names, time=time)


    ! fill the mfs for the next timestep
     do n = 1, nlevs
        call multifab_copy_c(s_old(n), 1, s_new(n), 1, nscal, s_old(n)%ng)

        call multifab_copy_c(Tcoeff2(n), 1, Tcoeff1(n), 1, 1)
        call multifab_copy_c(hcoeff2(n), 1, hcoeff1(n), 1, 1)
        call multifab_copy_c(Xkcoeff2(n), 1, Xkcoeff1(n), 1, nspec)
        call multifab_copy_c(pcoeff2(n), 1, pcoeff1(n), 1, 1)
     enddo

    if (parallel_IOProcessor()) print *, ''
 enddo
 
end subroutine varden



  

