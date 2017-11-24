module hg_multigrid_module

  use bl_types
  use mg_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bl_constants_module

  implicit none

  private

  public :: hg_multigrid

contains 

  subroutine hg_multigrid(mla,rh,unew,rhohalf,div_coeff_cart,phi,dx,the_bc_tower, &
                          stencil_type,rel_solver_eps,abs_solver_eps, nodalrhs)

    use bl_prof_module

    use nodal_stencil_fill_module , only : stencil_fill_nodal_all_mglevels
    use ml_solve_module     , only : ml_nd_solve
    use nodal_divu_module   , only : divu, subtract_divu_from_rh
    use probin_module       , only : hg_bottom_solver, max_mg_bottom_nlevels, &
                                     mg_verbose, cg_verbose, nodal, mg_bottom_nu, &
                                     hg_cycle_type

    use variables, only: press_comp
    use mg_eps_module, only: eps_hg_bottom

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: rh(:)
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: rhohalf(:)
    type(multifab ), intent(in   ) :: div_coeff_cart(:)
    type(multifab ), intent(inout) :: phi(:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: stencil_type
    real(dp_t)     , intent(in   ) :: rel_solver_eps
    real(dp_t)     , intent(in   ) :: abs_solver_eps

    type(multifab ), intent(inout), optional :: nodalrhs(:)

    ! Local variables
    type(box     )  :: pd

    type(mg_tower) :: mgt(mla%nlevel)

    type(multifab), allocatable :: coeffs(:)

    real(dp_t) :: bottom_solver_eps

    integer :: dm, nlevs
    integer :: bottom_solver, bottom_max_iter
    integer :: max_iter
    integer :: min_width
    integer :: max_nlevel
    integer :: nu1, nu2, smoother
    integer :: d,n
    integer :: max_nlevel_in
    integer :: do_diagnostics
    integer, allocatable :: lo_inflow(:),hi_inflow(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "hg_multigrid")

    dm = mla%dim
    nlevs = mla%nlevel

    !! Defaults:
    max_nlevel        = mgt(nlevs)%max_nlevel
    max_iter          = mgt(nlevs)%max_iter
    smoother          = mgt(nlevs)%smoother
    nu1               = mgt(nlevs)%nu1
    nu2               = mgt(nlevs)%nu2
    bottom_solver     = mgt(nlevs)%bottom_solver
    bottom_solver_eps = mgt(nlevs)%bottom_solver_eps
    bottom_max_iter   = mgt(nlevs)%bottom_max_iter
    min_width         = mgt(nlevs)%min_width

    if (parallel_IOProcessor()) then
       print *, 'doing hgproject with tolerance, eps = ', rel_solver_eps
    end if

    if ( hg_bottom_solver >= 0 ) then
        if (hg_bottom_solver == 4 .and. nboxes(phi(1)%la) == 1) then
           if (parallel_IOProcessor()) then
              print *,'Dont use hg_bottom_solver == 4 with only one grid -- '
              print *,'  Reverting to default bottom solver ',bottom_solver
           end if
        else if (hg_bottom_solver == 4 .and. max_mg_bottom_nlevels < 2) then
           if (parallel_IOProcessor()) then
              print *,'Dont use hg_bottom_solver == 4 with max_mg_bottom_nlevels < 2'
              print *,'  Reverting to default bottom solver ',bottom_solver
           end if
        else
           bottom_solver = hg_bottom_solver
        end if
    end if

    bottom_solver_eps = eps_hg_bottom

    ! Note: put this here for robustness
    max_iter = 100

    do n = nlevs, 1, -1

       if (dm == 1) then
          max_nlevel_in = 1
       else if (n == 1) then
          max_nlevel_in = max_nlevel
       else
          if ( all(mla%mba%rr(n-1,:) == 2) ) then
             max_nlevel_in = 1
          else if ( all(mla%mba%rr(n-1,:) == 4) ) then
             max_nlevel_in = 2
          else 
             call bl_error("HG_MULTIGRID: confused about ref_ratio")
          end if
       end if

       pd = layout_get_pd(mla%la(n))

       if (dm .eq. 1) then
          max_iter = 200
          nu1 = 200
          nu2 = 200
       end if

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,press_comp), &
                           stencil_type, &
                           dh = dx(n,:), &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           nub = mg_bottom_nu, &
                           cycle_type = hg_cycle_type, &
                           bottom_solver = bottom_solver, &
                           bottom_max_iter = bottom_max_iter, &
                           bottom_solver_eps = bottom_solver_eps, &
                           max_iter = max_iter, &
                           max_nlevel = max_nlevel_in, &
                           max_bottom_nlevel = max_mg_bottom_nlevels, &
                           min_width = min_width, &
                           eps = rel_solver_eps, &
                           abs_eps = abs_solver_eps, &
                           verbose = mg_verbose, &
                           cg_verbose = cg_verbose, &
                           nodal = nodal)
       
    end do

    !! Fill coefficient array
    do n = nlevs,1,-1

       allocate( coeffs(mgt(n)%nlevels))

       ! Build coeffs to pass into the multigrid
       ! Set to 1 in the interior, 0 in ghost cells since those are used in Anorm.
       call multifab_build( coeffs(mgt(n)%nlevels), mla%la(n), 1, 1)
       call setval(coeffs(mgt(n)%nlevels), 0.0_dp_t, 1, all=.true.)
       call multifab_setval(coeffs(mgt(n)%nlevels),1.d0)

       ! Build coeffs(i,j,1) = (rho/beta0)
       ! (and) coeffs(i,j,2) =   1./beta0 if coeff_ncomp > 1
       !
       ! Note: either rhohalf = rho/beta_0 or rho/beta_0^2
       ! (depending on use_alt_energy_fix).

       call multifab_div_div_c(coeffs(mgt(n)%nlevels),1,rhohalf(n),1,1,0)
       if (ncomp(coeffs(mgt(n)%nlevels)) .gt. 1) &
          call multifab_div_div_c(coeffs(mgt(n)%nlevels),2,div_coeff_cart(n),1,1,0)

       call multifab_fill_boundary(coeffs(mgt(n)%nlevels))

       call stencil_fill_nodal_all_mglevels(mgt(n), coeffs)

       call destroy(coeffs(mgt(n)%nlevels))
       deallocate(coeffs)

    end do

    ! ********************************************************************************
    ! Take the divergence of U:  RHS = div(U) 
    ! ********************************************************************************

    ! Set the inflow array -- 1 if inflow, otherwise 0
    allocate(lo_inflow(dm),hi_inflow(dm))
    lo_inflow(:) = 0
    hi_inflow(:) = 0
    do d = 1,dm 
       if (the_bc_tower%bc_tower_array(1)%phys_bc_level_array(0,d,1) == INLET) then
          lo_inflow(d) = 1
       end if
       if (the_bc_tower%bc_tower_array(1)%phys_bc_level_array(0,d,2) == INLET) then
          hi_inflow(d) = 1
       end if
    end do
    call divu(nlevs,mgt,unew,rh,mla%mba%rr,nodal,lo_inflow,hi_inflow)
    deallocate(lo_inflow,hi_inflow)

    ! ********************************************************************************
    ! Subtract S:  RHS = div(U) - S
    ! ********************************************************************************

    ! Note that we now set nodalrhs at outflow and at the fine nodes
    !      on coarse-fine boundaries in a call to enforce_dirichlet_rhs from ml_nd_solve.
    if (present(nodalrhs)) then
       call subtract_divu_from_rh(nlevs,mgt,rh,nodalrhs)
    end if

    ! ********************************************************************************
    ! Call the solver
    ! ********************************************************************************

    if ( mg_verbose >= 3 ) then
       do_diagnostics = 1
    else
       do_diagnostics = 0
    end if

    call ml_nd_solve(mla,mgt,rh,phi,do_diagnostics)

    ! ********************************************************************************
    ! Clean-up ...
    ! ********************************************************************************

    do n = 1, nlevs
       call mg_tower_destroy(mgt(n))
    end do

    call destroy(bpt)

  end subroutine hg_multigrid

  !   ********************************************************************************* !

end module hg_multigrid_module
