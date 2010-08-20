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

  subroutine hg_multigrid(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                          stencil_type,rel_solver_eps,abs_solver_eps,divu_rhs)

    use bl_prof_module

    use enforce_outflow_on_divu_module, only : enforce_outflow_on_divu_rhs

    use stencil_fill_module , only : stencil_fill_nodal_all_mglevels, stencil_fill_one_sided
    use ml_solve_module     , only : ml_nd_solve
    use nodal_divu_module   , only : divu, subtract_divu_from_rh
    use probin_module       , only : hg_bottom_solver, max_mg_bottom_nlevels, &
                                     verbose, mg_verbose, cg_verbose, nodal
    use geometry, only: dm, nlevs
    use variables, only: press_comp
    use mg_eps_module, only: eps_hg_bottom

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: rh(:)
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: rhohalf(:)
    type(multifab ), intent(inout) :: phi(:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: stencil_type
    real(dp_t)     , intent(in   ) :: rel_solver_eps
    real(dp_t)     , intent(in   ) :: abs_solver_eps

    type(multifab ), intent(inout), optional :: divu_rhs(:)

    ! Local variables
    type(box     )  :: pd
    type(  layout)  :: la

    type(mg_tower) :: mgt(mla%nlevel)

    type(multifab) :: one_sided_ss(2:mla%nlevel)

    type(multifab), allocatable :: coeffs(:)

    real(dp_t) :: bottom_solver_eps
    real(dp_t) :: omega

    integer :: i, ns
    integer :: bottom_solver, bottom_max_iter
    integer :: max_iter
    integer :: min_width
    integer :: max_nlevel
    integer :: nu1, nu2, nub, gamma, cycle_type, smoother
    integer :: n
    integer :: max_nlevel_in
    integer :: do_diagnostics

    type(bl_prof_timer), save :: bpt

    call build(bpt, "hg_multigrid")

    !! Defaults:
    max_nlevel        = mgt(nlevs)%max_nlevel
    max_iter          = mgt(nlevs)%max_iter
    smoother          = mgt(nlevs)%smoother
    nu1               = mgt(nlevs)%nu1
    nu2               = mgt(nlevs)%nu2
    nub               = mgt(nlevs)%nub
    gamma             = mgt(nlevs)%gamma
    omega             = mgt(nlevs)%omega
    cycle_type        = mgt(nlevs)%cycle_type
    bottom_solver     = mgt(nlevs)%bottom_solver
    bottom_solver_eps = mgt(nlevs)%bottom_solver_eps
    bottom_max_iter   = mgt(nlevs)%bottom_max_iter
    min_width         = mgt(nlevs)%min_width

    if (parallel_IOProcessor()) then
       print *, 'doing hgproject with tolerance, eps = ', rel_solver_eps
    end if

    if ( hg_bottom_solver >= 0 ) then
        if (hg_bottom_solver == 4 .and. nboxes(phi(1)) == 1) then
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

    if (stencil_type .eq. ST_DENSE) then
       if (dm .eq. 3) then
          i = mgt(nlevs)%nlevels
          if ( (dx(nlevs,1) .eq. dx(nlevs,2)) .and. &
               (dx(nlevs,1) .eq. dx(nlevs,3)) ) then
             ns = 21
          else
             ns = 27
          end if
       else if (dm .eq. 2) then
          ns = 9
       else if (dm .eq. 1) then
          ns = 3
       end if
    else
       ns = 2*dm+1
       do n = nlevs, 2, -1
          la = mla%la(n)
          call multifab_build(one_sided_ss(n), la, ns, 0, nodal)
          call setval(one_sided_ss(n), ZERO,all=.true.)
       end do
    end if

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
          omega = 1.33d0
          nu1 = 200
          nu2 = 200
       end if

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,press_comp), &
                           dh = dx(n,:), &
                           ns = ns, &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           nub = nub, &
                           gamma = gamma, &
                           cycle_type = cycle_type, &
                           omega = omega, &
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

       allocate(coeffs(mgt(n)%nlevels))

       la = mla%la(n)

       call multifab_build(coeffs(mgt(n)%nlevels), la, 1, 1)
       call setval(coeffs(mgt(n)%nlevels), 0.0_dp_t, 1, all=.true.)

       call mkcoeffs(rhohalf(n),coeffs(mgt(n)%nlevels))
       call multifab_fill_boundary(coeffs(mgt(n)%nlevels))

       call stencil_fill_nodal_all_mglevels(mgt(n), coeffs, stencil_type)

       if (stencil_type .eq. ST_CROSS .and. n .gt. 1) then
          i = mgt(n)%nlevels
          call stencil_fill_one_sided(one_sided_ss(n), coeffs(i), &
                                      mgt(n)%dh(:,i), &
                                      mgt(n)%mm(i), mgt(n)%face_type)
       end if

       call destroy(coeffs(mgt(n)%nlevels))
       deallocate(coeffs)

    end do

    ! ********************************************************************************
    ! Create the rhs
    ! ********************************************************************************

    call divu(nlevs,mgt,unew,rh,mla%mba%rr,nodal)

    ! Do rh = rh - divu_rhs (this routine preserves rh=0 on
    !  nodes which have bc_dirichlet = true.
    if (present(divu_rhs)) then
       call enforce_outflow_on_divu_rhs(divu_rhs,the_bc_tower)
       call subtract_divu_from_rh(nlevs,mgt,rh,divu_rhs)
    end if

    ! ********************************************************************************
    ! Call the solver
    ! ********************************************************************************

    if ( mg_verbose >= 3 ) then
       do_diagnostics = 1
    else
       do_diagnostics = 0
    end if

    call ml_nd_solve(mla,mgt,rh,phi,one_sided_ss,mla%mba%rr,do_diagnostics,&
                     rel_solver_eps,abs_solver_eps)

    ! ********************************************************************************
    ! Clean-up ...
    ! ********************************************************************************

    do n = nlevs,1,-1
       call multifab_fill_boundary(phi(n))
    end do

    do n = 1, nlevs
       call mg_tower_destroy(mgt(n))
    end do

    if (stencil_type .ne. ST_DENSE) then
       do n = nlevs, 2, -1
          call destroy(one_sided_ss(n))
       end do
    end if

    call destroy(bpt)

  end subroutine hg_multigrid

  !   ********************************************************************************* !

  subroutine mkcoeffs(rho,coeffs)

    use geometry, only: dm

    type(multifab) , intent(in   ) :: rho
    type(multifab) , intent(inout) :: coeffs

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    integer :: i,ng_r,ng_c

    ng_r = nghost(rho)
    ng_c = nghost(coeffs)

    do i = 1, nboxes(rho)
       if ( multifab_remote(rho, i) ) cycle
       rp => dataptr(rho   , i)
       cp => dataptr(coeffs, i)
       select case (dm)
       case (1)
          call mkcoeffs_1d(cp(:,1,1,1), ng_c, rp(:,1,1,1), ng_r)
       case (2)
          call mkcoeffs_2d(cp(:,:,1,1), ng_c, rp(:,:,1,1), ng_r)
       case (3)
          call mkcoeffs_3d(cp(:,:,:,1), ng_c, rp(:,:,:,1), ng_r)
       end select
    end do

  end subroutine mkcoeffs

  !   *********************************************************************************** !

  subroutine mkcoeffs_1d(coeffs,ng_c,rho,ng_r)

    use bl_constants_module

    integer                        :: ng_c,ng_r
    real(kind=dp_t), intent(inout) :: coeffs(1-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(1-ng_r:)

    integer :: i,nx

    nx = size(coeffs,dim=1) - 2

    do i = 1,nx
       coeffs(i) = ONE / rho(i)
    end do

  end subroutine mkcoeffs_1d

  !   *********************************************************************************** !

  subroutine mkcoeffs_2d(coeffs,ng_c,rho,ng_r)

    use bl_constants_module

    integer                        :: ng_c,ng_r
    real(kind=dp_t), intent(inout) :: coeffs(1-ng_c:,1-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(1-ng_r:,1-ng_r:)

    integer :: i,j
    integer :: nx,ny

    nx = size(coeffs,dim=1) - 2
    ny = size(coeffs,dim=2) - 2

    do j = 1,ny
       do i = 1,nx
          coeffs(i,j) = ONE / rho(i,j)
       end do
    end do

  end subroutine mkcoeffs_2d

  !   ********************************************************************************** !

  subroutine mkcoeffs_3d(coeffs,ng_c,rho,ng_r)

      use bl_constants_module

    integer                        :: ng_c,ng_r
    real(kind=dp_t), intent(inout) :: coeffs(1-ng_c:,1-ng_c:,1-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(1-ng_r:,1-ng_r:,1-ng_r:)

    integer :: i,j,k
    integer :: nx,ny,nz

    nx = size(coeffs,dim=1) - 2
    ny = size(coeffs,dim=2) - 2
    nz = size(coeffs,dim=3) - 2

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             coeffs(i,j,k) = ONE / rho(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine mkcoeffs_3d

end module hg_multigrid_module
