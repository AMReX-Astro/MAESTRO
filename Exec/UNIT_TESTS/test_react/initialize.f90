module initialize_module
    
  use bl_constants_module
  use ml_boxarray_module

  implicit none

  private

  public :: initialize_dx, initialize_bc

contains

  subroutine initialize_dx(dx,mba,num_levs)
    use probin_module, only: prob_lo, prob_hi

    real(dp_t)       , intent(out), pointer :: dx(:,:)
    type(ml_boxarray), intent(in )          :: mba
    integer          , intent(in )          :: num_levs
    
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

  subroutine initialize_bc(the_bc_tower,num_levs,pmask)
    use bc_module
    use define_bc_module
    use bl_error_module
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
end module initialize_module
