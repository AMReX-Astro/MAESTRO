! Setup a 3D grid of smoothly varying rho, T, and user-defined X.  Then
! call react_state() on the grid and output the results.

subroutine varden()

  use BoxLib
  use f2kcli            !Grants access to command line args
  use geometry
  use ml_layout_module
  use multifab_module
  use average_module
  use variables
  use probin_module, only: small_temp, &
                           drive_initial_convection,    & 
                           use_tfromp, min_time_step,   &
                           react_its, do_burning, do_heating, &
                           dens_min, base_cutoff_density
  use runtime_init_module
  use bl_constants_module
  use bl_types
  use define_bc_module
  use network
  use eos_module
  use react_state_module
  use varden_aux
  use simple_log_module
  use cputime_module, only: start_cputime_clock

  !Local variables
  implicit none
  
  ! Conventional fluid state multifabs
  type(multifab) , allocatable :: s(:), snew(:)                             

  ! react_state output
  type(multifab) , allocatable :: rho_omegadot(:), rho_Hnuc(:), rho_Hext(:) 

  real(kind=dp_t), allocatable :: tempbar(:,:), pbar(:,:)
 
  real(kind=dp_t), pointer :: dx(:,:)
  
  type(ml_layout) :: mla
  type(bc_tower)  :: bct
 
  character(len=100) :: temp_buf 

  real(kind=dp_t) :: m_in, m_out

  integer :: i, n, res
  integer :: ii, jj, kk
  integer :: nlevs

  logical :: dbo, dho
  
  !### Execution ###
  !## Initialization ##
  !General Maestro initializations
  call runtime_init()
  call init_variables()
  call start_cputime_clock()
  call simple_log_init()

  !Check for unimplemented modes
  if (drive_initial_convection .neqv. .false.) then
     call bl_error('ERROR: Driving initial convection not currently supported')
  endif
  if (use_tfromp .neqv. .false.) then
     call bl_error('ERROR: Getting temperature from pressure not currently supported')
  endif


  ! density sanity check
  if (dens_min < base_cutoff_density) then
     call bl_error('ERROR: dens_min < base_cutoff_density (= burning_cutoff_density)')
  endif

  !Grid/Geometry initializations
  call grid_init(mla, nlevs, bct, dx)

  !Microphysics
  call network_init()
  call eos_init(small_temp=small_temp)
  
  !Initialize the varden_aux reaction data
  call react_init(mla)
  
  !Initialize multifabs & the bc_level
  allocate(s(nlevs))
  allocate(snew(nlevs))
  allocate(rho_omegadot(nlevs))
  allocate(rho_Hnuc(nlevs))
  allocate(rho_Hext(nlevs))
  
  do n = 1,nlevs
    call multifab_build(s(n)            , mla%la(n), nscal , 1)
    call multifab_build(snew(n)         , mla%la(n), nscal , 1)
    call multifab_build(rho_omegadot(n) , mla%la(n), nspec , 1)
    call multifab_build(rho_Hnuc(n)     , mla%la(n), 1     , 1)
    call multifab_build(rho_Hext(n)     , mla%la(n), 1     , 1)

    call bc_tower_level_build(bct,n,mla%la(n))
  end do
  call init_multilevel(s)   !Creates numdisjointchunks, r_start_coord, r_end_coord
                            !(average() makes use of these)
  
  !Initialization of the thermodynamic grid
  call therm_init(mla, s, tempbar, pbar)

  !## react_state() testing ##
  !Calculate react_state inputs
  call average(mla,s,tempbar,dx,temp_comp)
  call setval(rho_Hext(1), ZERO, all=.true.)

  !Check the consistency of each mode
  dbo = do_burning
  dho = do_heating

  !Explore ten orders of magnitude of the time domain using user inputs.
  do_burning = dbo
  do_heating = dho
  do i=0, react_its-1
    print *, '***react it: ', i
    print *, '***dt:       ', 10**(i)*min_time_step
    call react_state(mla,tempbar,s,snew,rho_omegadot,rho_Hnuc,rho_Hext,pbar, &
                     10**(i)*min_time_step,dx,bct%bc_tower_array)
    write(temp_buf, *) i
    temp_buf = adjustl(temp_buf)
    call react_write(snew, s, rho_omegadot, rho_Hnuc, rho_Hext, mla, &
                     10**(i)*min_time_step, "dtE+" // trim(temp_buf),bct)
  enddo

  !## Clean-up ##
  call varden_close() ! -- must do this before we destroy the layout

  !If you (or a subroutine) built it, destroy it!
  do n = 1,nlevs
    call destroy(s(n))
    call destroy(snew(n))
    call destroy(rho_omegadot(n))
    call destroy(rho_Hnuc(n))
    call destroy(rho_Hext(n))
  end do
  call bc_tower_destroy(bct)
  call destroy(mla)

  !If you (or a subroutine) alloc'd it, dealloc it!
  deallocate(tempbar) 
  deallocate(pbar)
  deallocate(s)
  deallocate(snew)
  deallocate(rho_omegadot)
  deallocate(rho_Hnuc)
  deallocate(rho_Hext)

  call runtime_close()
end subroutine varden
