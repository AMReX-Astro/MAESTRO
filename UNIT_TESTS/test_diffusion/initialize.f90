module initialize_module

  use define_bc_module
  use ml_layout_module
  use multifab_module
  use bc_module
  use probin_module
  use variables, only: nscal, rho_comp, rhoh_comp, temp_comp
  use geometry
  use network, only: nspec
  use bl_constants_module
  use base_state_module
  use base_io_module

  implicit none

  private

  public :: initialize_with_fixed_grids, initialize_bc, initialize_dx

contains
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_with_fixed_grids(mla,time,dt,dx,pmask,the_bc_tower,&
                                         s,s_init,p0,thermal,&
                                         diffusion_coefficient)

    use box_util_module
    use init_scalar_module
    use average_module
    use restrict_base_module
    
    type(ml_layout),intent(  out) :: mla
    real(dp_t)    , intent(inout) :: time,dt
    real(dp_t)    , intent(  out) :: diffusion_coefficient
    logical       , intent(in   ) :: pmask(:)
    type(bc_tower), intent(  out) :: the_bc_tower
    type(multifab),  pointer      :: s(:), thermal(:)
    real(kind=dp_t), pointer      :: s_init(:,:,:),p0(:,:),dx(:,:)

    ! local
    type(ml_boxarray) :: mba

    real(dp_t) :: lenx,leny,lenz,max_dist

    integer :: n,ng_s,dm,nlevs

    ! set time and dt
    time = ZERO
    dt = 1.d20

    ! create mba
    call read_a_hgproj_grid(mba,test_set)

    ! create mla
    call ml_layout_build(mla,mba,pmask)

    dm = mla%dim
    
    ! check for proper nesting
    if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
       call bl_error('fixed_grids not properly nested')
    end if
    
    ! initialize nlevs
    nlevs = mla%nlevel
    nlevs_radial = merge(1, nlevs, spherical .eq. 1)

    ! initialize boundary conditions
    call initialize_bc(the_bc_tower,nlevs,pmask)
    do n = 1,nlevs
       call bc_tower_level_build(the_bc_tower,n,mla%la(n))
    end do

    ! allocate states
    allocate(s(nlevs))
    allocate(thermal(nlevs))

    if (ppm_type .eq. 2) then
       ng_s = 4
    else
       ng_s = 3
    end if

    ! build states
    do n = 1,nlevs
       call multifab_build(s(n)      , mla%la(n), nscal, ng_s)
       call multifab_build(thermal(n), mla%la(n),     1,    1)

       call setval(s(n)      , ZERO, all=.true.)
       call setval(thermal(n), ZERO, all=.true.)

    end do

    ! initialize dx
    call initialize_dx(dx,mba,nlevs)

    ! initialize cutoff arrays
    call init_cutoff(nlevs)

    ! now that we have dx we can initialize nr_fine and dr_fine
    if (spherical .eq. 1) then

       ! for spherical, we will now require that dr_fine = dx
       dr_fine = dx(nlevs,1) / dble(drdxfac)
       
       lenx = HALF * (prob_hi(1) - prob_lo(1))
       leny = HALF * (prob_hi(2) - prob_lo(2))
       lenz = HALF * (prob_hi(3) - prob_lo(3))
       
       max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
       nr_fine = int(max_dist / dr_fine) + 1
       
    else
       
       nr_fine = extent(mla%mba%pd(nlevs),dm)
       dr_fine = (prob_hi(dm)-prob_lo(dm)) / dble(nr_fine)
       
    end if

    ! create numdisjointchunks, r_start_coord, r_end_coord
    call init_multilevel(s)

    ! now that we have nr_fine and dr_fine we can create nr, dr, r_cc_loc, r_edge_loc
    call init_radial(nlevs,mba)

    ! now that we have nr_fine we can allocate 1d arrays
    call initialize_1d_arrays(nlevs,s_init,p0)

    ! now that we have dr and nr we can fill initial state
    if (spherical .eq. 1) then
       call init_base_state(1,model_file,s_init(1,:,:),p0(1,:),dx(nlevs,:),&
                            diffusion_coefficient)
    else
       ! init_base_state requires loop backwards over levels
       do n=nlevs,1,-1
          call init_base_state(n,model_file,s_init(n,:,:),p0(n,:),dx(n,:),&
                               diffusion_coefficient)
       end do
    end if
    
    ! fill the s multifab
    call initscalardata(s,s_init,p0,dx,the_bc_tower%bc_tower_array,mla,&
                        diffusion_coefficient)

    call destroy(mba)


  end subroutine initialize_with_fixed_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine initialize_bc(the_bc_tower,num_levs,pmask)

    use bc_module
    use probin_module, only : bcx_lo, bcx_hi, bcy_lo, bcy_hi, bcz_lo, bcz_hi

    type(bc_tower), intent(  out) :: the_bc_tower
    integer       , intent(in   ) :: num_levs
    logical       , intent(in   ) :: pmask(:)
    
    integer :: domain_phys_bc(size(pmask),2), dm

    dm = size(pmask)

    ! Define the physical boundary conditions on the domain
    ! Put the bc values from the inputs file into domain_phys_bc
    domain_phys_bc(1,1) = bcx_lo
    domain_phys_bc(1,2) = bcx_hi
    if (pmask(1)) then
       domain_phys_bc(1,:) = BC_PER
       if (bcx_lo .ne. -1 .or. bcx_hi .ne. -1) &
            call bl_error('MUST HAVE BCX = -1 if PMASK = T')
    end if
    if (dm > 1) then
       domain_phys_bc(2,1) = bcy_lo
       domain_phys_bc(2,2) = bcy_hi
       if (pmask(2)) then
          domain_phys_bc(2,:) = BC_PER
          if (bcy_lo .ne. -1 .or. bcy_hi .ne. -1) &
               call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
       end if
    end if
    if (dm > 2) then
       domain_phys_bc(3,1) = bcz_lo
       domain_phys_bc(3,2) = bcz_hi
       if (pmask(3)) then
          domain_phys_bc(3,:) = BC_PER
          if (bcz_lo .ne. -1 .or. bcz_hi .ne. -1) &
               call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
       end if
    end if
    
    ! Initialize the_bc_tower object.
    call bc_tower_init(the_bc_tower,num_levs,dm,domain_phys_bc)
    
  end subroutine initialize_bc

  subroutine initialize_dx(dx,mba,num_levs)

    real(dp_t)       , pointer     :: dx(:,:)
    type(ml_boxarray), intent(in ) :: mba
    integer          , intent(in ) :: num_levs
    
    integer :: n,d,dm

    dm = mba%dim
    
    allocate(dx(num_levs,dm))
    
    do d=1,dm
       dx(1,d) = (prob_hi(d)-prob_lo(d)) / real(extent(mba%pd(1),d),kind=dp_t)
    end do
    do n=2,num_levs
       dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
    end do

  end subroutine initialize_dx

  subroutine initialize_1d_arrays(num_levs,s_init,p0)

    integer    , intent(in)  :: num_levs    
    real(kind=dp_t), pointer :: s_init(:,:,:), p0(:,:)
    
    if (spherical .eq. 0) then
       allocate(s_init(num_levs,0:nr_fine-1,nscal))
       allocate(p0    (num_levs,0:nr_fine-1))
    else
       allocate(s_init(1,0:nr_fine-1,nscal))
       allocate(p0    (1,0:nr_fine-1))
    end if
    
    s_init = ZERO
    p0     = ZERO

  end subroutine initialize_1d_arrays
  
end module initialize_module
