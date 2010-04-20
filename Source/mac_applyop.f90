module mac_applyop_module

  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bl_constants_module

  implicit none

  private

  public :: mac_applyop

contains

  subroutine mac_applyop(mla,res,phi,alpha,beta,dx,the_bc_tower,bc_comp,stencil_order,ref_ratio)
    use mg_module
    use stencil_fill_module
!   use coeffs_module
    use ml_cc_module, only: ml_cc_applyop
    use probin_module, only: cg_verbose, mg_verbose
    use geometry, only: dm, nlevs

    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: stencil_order
    integer        , intent(in   ) :: ref_ratio(:,:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower) , intent(in)    :: the_bc_tower
    integer        , intent(in   ) :: bc_comp
    type(multifab) , intent(in   ) :: alpha(:), beta(:,:)
    type(multifab) , intent(inout) :: res(:), phi(:)

    type(layout  ) :: la
    type(box     ) :: pd

    type(multifab), allocatable :: coeffs(:)

    type(mg_tower)  :: mgt(mla%nlevel)
    integer         :: i, ns

    ! MG solver defaults
    integer :: bottom_solver, bottom_max_iter
    integer    :: max_iter
    integer    :: min_width
    integer    :: max_nlevel
    integer    :: n, nu1, nu2, gamma, ncycle, smoother
    real(dp_t) :: eps,abs_eps,omega,bottom_solver_eps
    real(dp_t) ::  xa(dm),  xb(dm)
    real(dp_t) :: pxa(dm), pxb(dm)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mac_applyop")

    max_nlevel        = mgt(nlevs)%max_nlevel
    max_iter          = mgt(nlevs)%max_iter
    eps               = mgt(nlevs)%eps
    abs_eps           = mgt(nlevs)%abs_eps
    smoother          = mgt(nlevs)%smoother
    nu1               = mgt(nlevs)%nu1
    nu2               = mgt(nlevs)%nu2
    gamma             = mgt(nlevs)%gamma
    omega             = mgt(nlevs)%omega
    ncycle            = mgt(nlevs)%cycle_type
    bottom_solver     = mgt(nlevs)%bottom_solver
    bottom_solver_eps = mgt(nlevs)%bottom_solver_eps
    bottom_max_iter   = mgt(nlevs)%bottom_max_iter
    min_width         = mgt(nlevs)%min_width

    ns = 1 + dm*3

    do n = nlevs, 1, -1

       pd = layout_get_pd(mla%la(n))

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp),&
                           dh = dx(n,:), &
                           ns = ns, &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           gamma = gamma, &
                           cycle_type = ncycle, &
                           omega = omega, &
                           bottom_solver = bottom_solver, &
                           bottom_max_iter = bottom_max_iter, &
                           bottom_solver_eps = bottom_solver_eps, &
                           max_iter = max_iter, &
                           max_nlevel = max_nlevel, &
                           min_width = min_width, &
                           eps = eps, &
                           abs_eps = abs_eps, &
                           verbose = mg_verbose, &
                           cg_verbose = cg_verbose, &
                           nodal = res(nlevs)%nodal)

    end do

    !! Fill coefficient array

    do n = nlevs,1,-1

       allocate(coeffs(mgt(n)%nlevels))

       la = mla%la(n)
       pd = layout_get_pd(la)

       call multifab_build(coeffs(mgt(n)%nlevels), la, 1+dm, 1)
       call multifab_copy_c(coeffs(mgt(n)%nlevels),1,alpha(n),1,1,ng=alpha(n)%ng)
       call multifab_copy_c(coeffs(mgt(n)%nlevels),2,beta(n,1),1,1,ng=beta(n,1)%ng)
       if (dm > 1) &
          call multifab_copy_c(coeffs(mgt(n)%nlevels),3,beta(n,2),1,1,ng=beta(n,2)%ng)
       if (dm > 2) &
         call multifab_copy_c(coeffs(mgt(n)%nlevels),4,beta(n,3),1,1,ng=beta(n,3)%ng)

       if (n > 1) then
          xa = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
          xb = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
       else
          xa = ZERO
          xb = ZERO
       end if

       pxa = ZERO
       pxb = ZERO

       call stencil_fill_cc_all_mglevels(mgt(n), coeffs, xa, xb, pxa, pxb, pd, stencil_order, &
                                         the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp))

       call destroy(coeffs(mgt(n)%nlevels))
       deallocate(coeffs)

    end do

    call ml_cc_applyop(mla, mgt, res, phi, ref_ratio)

    do n = 1, nlevs
       call mg_tower_destroy(mgt(n))
    end do

    call destroy(bpt)

  end subroutine mac_applyop

end module mac_applyop_module
