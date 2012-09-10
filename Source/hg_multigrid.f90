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

  subroutine hg_multigrid(mla,rh,unew,rhohalf,div_coeff_3d,phi,dx,the_bc_tower, &
                          stencil_type,rel_solver_eps,abs_solver_eps, &
                          using_alt_energy_fix,divu_rhs)

    use bl_prof_module

    use enforce_outflow_on_divu_module, only : enforce_outflow_on_divu_rhs

    use nodal_stencil_fill_module , only : stencil_fill_nodal_all_mglevels,      &
                                           stencil_fill_one_sided
    use ml_solve_module     , only : ml_nd_solve
    use nodal_divu_module   , only : divu, subtract_divu_from_rh
    use probin_module       , only : hg_bottom_solver, max_mg_bottom_nlevels, &
                                     mg_verbose, cg_verbose, nodal, mg_bottom_nu

    use variables, only: press_comp
    use mg_eps_module, only: eps_hg_bottom

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: rh(:)
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: rhohalf(:)
    type(multifab ), intent(in   ) :: div_coeff_3d(:)
    type(multifab ), intent(inout) :: phi(:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: stencil_type
    logical        , intent(in   ) :: using_alt_energy_fix
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

    integer :: i, ns, dm, nlevs
    integer :: bottom_solver, bottom_max_iter
    integer :: max_iter
    integer :: min_width
    integer :: max_nlevel
    integer :: nu1, nu2, gamma, cycle_type, smoother
    integer :: d,n,j
    integer :: max_nlevel_in
    integer :: do_diagnostics
    integer :: coeff_ncomp
    integer, allocatable :: lo_inflow(:),hi_inflow(:)

    real(dp_t), pointer :: p(:,:,:,:)

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

    if (stencil_type .eq. ND_DENSE_STENCIL) then
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
          call multifab_build(one_sided_ss(n), la, ns, 0, nodal, stencil=.true.)

          ! set the stencil to zero manually -- multifab_setval doesn't work
          ! with stencil = .true.
          do j = 1, nfabs(one_sided_ss(n))
             p => dataptr(one_sided_ss(n), j)
             p = ZERO
          enddo

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
                           stencil_type, &
                           dh = dx(n,:), &
                           ns = ns, &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           nub = mg_bottom_nu, &
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

!   if (using_alt_energy_fix) then
!      coeff_ncomp = 2
!   else
       coeff_ncomp = 1
!   end if

    do n = nlevs,1,-1

       allocate( coeffs(mgt(n)%nlevels))

       la = mla%la(n)

       ! Build coeffs to pass into the multigrid
       call multifab_build( coeffs(mgt(n)%nlevels), la, coeff_ncomp, 1)
       call setval(coeffs(mgt(n)%nlevels), 0.0_dp_t, 1, all=.true.)

       ! Build coeffs(i,j,1) = (rho/beta0)
       ! (and) coeffs(i,j,2) =   1./beta0 if coeff_ncomp > 1
       call mkcoeffs(rhohalf(n),div_coeff_3d(n),coeffs(mgt(n)%nlevels))

       call multifab_fill_boundary(coeffs(mgt(n)%nlevels))

       call stencil_fill_nodal_all_mglevels(mgt(n), coeffs, stencil_type)

       if (stencil_type .eq. ND_CROSS_STENCIL .and. n .gt. 1) then
          i = mgt(n)%nlevels
          call stencil_fill_one_sided(one_sided_ss(n), coeffs(i), &
                                      mgt(n)%dh(:,i), &
                                      mgt(n)%mm(i), mgt(n)%face_type)
       end if

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

    ! (this routine preserves rh=0 on nodes which have bc_dirichlet = true.)
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

    if (stencil_type .ne. ND_DENSE_STENCIL) then
       do n = nlevs, 2, -1
          call destroy(one_sided_ss(n))
       end do
    end if

    call destroy(bpt)

  end subroutine hg_multigrid

  !   ********************************************************************************* !

  subroutine mkcoeffs(rho,div_coeff_3d,coeffs)

    type(multifab) , intent(in   ) :: rho
    type(multifab) , intent(in   ) :: div_coeff_3d
    type(multifab) , intent(inout) :: coeffs

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    integer :: i,ng_r,ng_c,ng_d,dm
    integer :: lo(get_dim(rho)),hi(get_dim(rho))

    dm = get_dim(rho)

    ng_r = nghost(rho)
    ng_d = nghost(div_coeff_3d)
    ng_c = nghost(coeffs)

    do i = 1, nfabs(rho)
       rp => dataptr(rho   , i)
       dp => dataptr(div_coeff_3d, i)
       cp => dataptr(coeffs, i)
       lo = lwb(get_box(rho,i))
       hi = upb(get_box(rho,i))
       select case (dm)
       case (1)
          call mkcoeffs_1d(cp(:,1,1,:), ng_c, rp(:,1,1,1), ng_r, dp(:,1,1,1), ng_d, lo, hi)
       case (2)
          call mkcoeffs_2d(cp(:,:,1,:), ng_c, rp(:,:,1,1), ng_r, dp(:,:,1,1), ng_d, lo, hi)
       case (3)
          call mkcoeffs_3d(cp(:,:,:,:), ng_c, rp(:,:,:,1), ng_r, dp(:,:,:,1), ng_d, lo, hi)
       end select
    end do

  end subroutine mkcoeffs

  !   *********************************************************************************** !

  subroutine mkcoeffs_1d(coeffs,ng_c,rho,ng_r,divcoeff,ng_d,lo,hi)

    use bl_constants_module

    integer                        :: ng_c,ng_r,ng_d,lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: coeffs(lo(1)-ng_c:,:)
    real(kind=dp_t), intent(in   ) ::    rho(lo(1)-ng_r:)
    real(kind=dp_t), intent(in   ) :: divcoeff(lo(1)-ng_d:)

    integer :: i

    do i = lo(1),hi(1)
       coeffs(i,1) = ONE / rho(i)
    end do

    if (size(coeffs,dim=2) .gt. 1) then
       do i = lo(1),hi(1)
          coeffs(i,2) = ONE / divcoeff(i)
       end do
    end if

  end subroutine mkcoeffs_1d

  !   *********************************************************************************** !

  subroutine mkcoeffs_2d(coeffs,ng_c,rho,ng_r,divcoeff,ng_d,lo,hi)

    use bl_constants_module

    integer                        :: ng_c,ng_r,ng_d,lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::   coeffs(lo(1)-ng_c:,lo(2)-ng_c:,:)
    real(kind=dp_t), intent(in   ) ::      rho(lo(1)-ng_r:,lo(2)-ng_r:)
    real(kind=dp_t), intent(in   ) :: divcoeff(lo(1)-ng_d:,lo(2)-ng_d:)

    integer :: i,j

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          coeffs(i,j,1) = ONE / rho(i,j)
       end do
    end do

    if (size(coeffs,dim=3) .gt. 1) then
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             coeffs(i,j,2) = ONE / divcoeff(i,j)
          end do
       end do
    end if

  end subroutine mkcoeffs_2d

  !   ********************************************************************************** !

  subroutine mkcoeffs_3d(coeffs,ng_c,rho,ng_r,divcoeff,ng_d,lo,hi)

      use bl_constants_module

    integer                        :: ng_c,ng_r,ng_d,lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::   coeffs(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
    real(kind=dp_t), intent(in   ) ::      rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
    real(kind=dp_t), intent(in   ) :: divcoeff(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:)

    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             coeffs(i,j,k,1) = ONE / rho(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    if (size(coeffs,dim=4) .gt. 1) then
       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             coeffs(i,j,k,2) = ONE / divcoeff(i,j,k)
          end do
       end do
       end do
    end if

  end subroutine mkcoeffs_3d

end module hg_multigrid_module
