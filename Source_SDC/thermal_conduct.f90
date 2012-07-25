module thermal_conduct_module

  use bl_types
  use multifab_module
  use bl_constants_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: thermal_conduct_predictor, thermal_conduct_corrector, &
       make_explicit_thermal_hterm

contains 

  subroutine thermal_conduct_predictor(mla,dx,dt,sold,shat,p0_old,p0_new, &
                                       hcoeff_old,Xkcoeff_old,pcoeff_old, &
                                       the_bc_tower)

    use bl_prof_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use bndry_reg_module

    use mac_multigrid_module , only : mac_multigrid
    use mac_applyop_module   , only : mac_applyop
    use fill_3d_module       , only : put_1d_array_on_cart, put_data_on_faces
    use ml_restriction_module, only : ml_cc_restriction_c

    use variables    , only : foextrap_comp, rho_comp, spec_comp, rhoh_comp
    use network      , only : nspec
    use mg_eps_module, only : eps_mac

    type(ml_layout), intent(inout) :: mla
    real(dp_t)     , intent(in   ) :: dx(:,:),dt
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: shat(:)
    type(multifab) , intent(in   ) :: hcoeff_old(:)
    type(multifab) , intent(in   ) :: Xkcoeff_old(:)
    type(multifab) , intent(in   ) :: pcoeff_old(:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! Local
    type(multifab) :: alpha(mla%nlevel)
    type(multifab) :: beta(mla%nlevel,mla%dim)
    type(multifab) :: phi(mla%nlevel)
    type(multifab) :: Lphi(mla%nlevel)
    type(multifab) :: rhs(mla%nlevel)

    integer                     :: stencil_order,dm,nlevs
    integer                     :: i,n,comp
    type(bndry_reg)             :: fine_flx(2:mla%nlevel)
    real(dp_t)                  :: abs_eps, abs_solver_eps, rel_solver_eps

    type(bl_prof_timer), save :: bpt

    ! here we will solve:
    !
    ! (alpha - div beta grad) h = RHS
    !
    ! First we will construct the RHS by adding each of the terms in turn.  
    !
    ! To actually construct each div (c grad q) term for the RHS, we will 
    ! make use of the mac_applyop routine, which constructs the quantity
    !
    !     (alpha - div beta grad) phi
    !
    ! For all RHS terms, we set alpha = 0, beta = the appropriate 
    ! coefficient, and phi = the quantity being diffused.

    call build(bpt, "therm_cond_full_alg")

    dm = mla%dim
    nlevs = mla%nlevel

    stencil_order = 2

    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    do n=1,nlevs
       call multifab_build(alpha(n), mla%la(n),  1, 0)
       call multifab_build(phi(n)  , mla%la(n),  1, 1)
       call multifab_build(Lphi(n) , mla%la(n),  1, 0)
       call multifab_build(rhs(n)  , mla%la(n),  1, 0)
       do i = 1,dm
          call multifab_build_edge(beta(n,i), mla%la(n), 1, 0, i)
       end do
    end do

    ! set alpha to zero while we evaluate the diffusion terms for the rhs
    do n=1,nlevs
       call setval(alpha(n), ZERO, all=.true.)
       call setval(rhs(n), ZERO, all=.true.)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add h^old diffusion to rhs (using old coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! put beta on faces
    call put_data_on_faces(mla,hcoeff_old,1,beta,.true.)

    ! set phi = h^old
    do n=1,nlevs
       call multifab_copy_c(phi(n),1,sold(n),rhoh_comp,1,1)
       call multifab_div_div_c(phi(n),1,sold(n),rho_comp,1,1)
    end do

    ! apply the operator
    call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                     dm+rhoh_comp,stencil_order,mla%mba%rr)

    ! add Lphi to rhs
    do n=1,nlevs
       call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,0)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add X_k^old diffusion to rhs (using old coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do comp=1,nspec

       ! put beta on faces
       call put_data_on_faces(mla,Xkcoeff_old,comp,beta,.true.)

       ! set phi = X_k^old
       do n=1,nlevs
          call multifab_copy_c(phi(n),1,sold(n),spec_comp+comp-1,1,1)
          call multifab_div_div_c(phi(n),1,sold(n),rho_comp,1,1)
       enddo

       ! apply the operator
       call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                        dm+spec_comp+comp-1,stencil_order,mla%mba%rr)

       ! add Lphi to rhs
       do n=1,nlevs
          call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,0)
       enddo

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add X_k^hat diffusion to rhs (using old coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do comp=1,nspec

       ! put beta on faces
       call put_data_on_faces(mla,Xkcoeff_old,comp,beta,.true.)

       ! set phi = X_k^hat
       do n=1,nlevs
          call multifab_copy_c(phi(n),1,shat(n),spec_comp+comp-1,1,1)
          call multifab_div_div_c(phi(n),1,shat(n),rho_comp,1,1)
       enddo

       ! apply the operator
       call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                        dm+spec_comp+comp-1,stencil_order,mla%mba%rr)

       ! add Lphi to rhs
       do n=1,nlevs
          call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,0)
       enddo

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add p0^old diffusion to rhs (using old coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! put beta on faces
    call put_data_on_faces(mla,pcoeff_old,1,beta,.true.)

    ! set phi = p0^old
    call put_1d_array_on_cart(p0_old,phi,foextrap_comp,.false.,.false., &
                              dx,the_bc_tower%bc_tower_array,mla)

    ! apply the operator
    call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                     foextrap_comp,stencil_order,mla%mba%rr)

    ! add Lphi to rhs
    do n=1,nlevs
       call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,0)
    enddo    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add p0^new diffusion to rhs (using old coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! put beta on faces
    call put_data_on_faces(mla,pcoeff_old,1,beta,.true.)

    ! set phi = p0^new
    call put_1d_array_on_cart(p0_new,phi,foextrap_comp,.false.,.false., &
                              dx,the_bc_tower%bc_tower_array,mla)

    ! apply the operator
    call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                     foextrap_comp,stencil_order,mla%mba%rr)

    ! add Lphi to rhs
    do n=1,nlevs
       call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,0)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! scale rhs and add remaining terms
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! multiply rhs by 0.5*dt
    do n=1,nlevs
       call multifab_mult_mult_s(rhs(n),0.5d0*dt,0)
    end do

    ! add (rhoh)^old + dt*(aofs + intra) to rhs
    do n=1,nlevs
       call multifab_plus_plus_c(rhs(n),1,shat(n),rhoh_comp,1,0)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setup LHS coefficients for implicit solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! copy rho^hat=rho^new into alpha
    do n=1,nlevs
       call multifab_copy_c(alpha(n),1,shat(n),rho_comp,1,0)
    enddo

    ! copy old coefficients into beta
    call put_data_on_faces(mla,hcoeff_old,1,beta,.true.)

    ! multiply beta by -dt/2
    do n=1,nlevs
       do i = 1,dm 
          call multifab_mult_mult_s_c(beta(n,i),1,-dt/2.d0,1,0)
       enddo
    enddo

    ! initialize phi to a reasonable guess, (rhoh)^hat / rho^hat
    do n=1,nlevs
       call multifab_copy_c(phi(n),1,shat(n),rhoh_comp,1,1)
       call multifab_div_div_c(phi(n),1,shat(n),rho_comp,1,1)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now do the implicit solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Compute norm(phi) to be used inside the MG solver as part of a stopping criterion
    abs_eps = -1.d0
    do n = 1,nlevs
       abs_eps = max(abs_eps, norm_inf(phi(n)) / dx(n,1))
    end do
    abs_solver_eps = eps_mac * abs_eps

    rel_solver_eps = eps_mac

    ! Call the solver to obtain h^hat (it will be stored in phi)
    ! solves (alpha - div beta grad) phi = rhs
    call mac_multigrid(mla,rhs,phi,fine_flx,alpha,beta,dx,the_bc_tower, &
                       dm+rhoh_comp,stencil_order,mla%mba%rr,rel_solver_eps,abs_solver_eps)

    ! load (rho*h)^hat into shat
    do n=1,nlevs
       call multifab_copy_c(shat(n),rhoh_comp,phi(n),1,1,0)
       call multifab_mult_mult_c(shat(n),rhoh_comp,shat(n),rho_comp,1,0)
    enddo

    do n=2,nlevs
       call destroy(fine_flx(n))
    end do

    do n=1,nlevs
       call destroy(alpha(n))
       call destroy(phi(n))
       call destroy(Lphi(n))
       call destroy(rhs(n))
       do i = 1,dm
          call destroy(beta(n,i))
       end do
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(shat(nlevs),rhoh_comp,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(shat(nlevs),rhoh_comp,dm+rhoh_comp,1, &
                            the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(shat(n-1),rhoh_comp,shat(n),rhoh_comp,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(shat(n),shat(n-1),nghost(shat(1)),mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         rhoh_comp,dm+rhoh_comp,1,fill_crse_input=.false.)
       enddo

    end if

    call destroy(bpt)

  end subroutine thermal_conduct_predictor

  subroutine thermal_conduct_corrector(mla,dx,dt,sold,shat,snew,p0_old,p0_new, &
                                       hcoeff_old,Xkcoeff_old,pcoeff_old, &
                                       hcoeff_new,Xkcoeff_new,pcoeff_new, &
                                       the_bc_tower)

    use bl_prof_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use bndry_reg_module

    use mac_multigrid_module , only : mac_multigrid
    use mac_applyop_module   , only : mac_applyop
    use fill_3d_module       , only : put_1d_array_on_cart, put_data_on_faces
    use ml_restriction_module, only : ml_cc_restriction_c

    use variables    , only : foextrap_comp, rho_comp, spec_comp, rhoh_comp
    use network      , only : nspec
    use mg_eps_module, only : eps_mac

    type(ml_layout), intent(inout) :: mla
    real(dp_t)     , intent(in   ) :: dx(:,:),dt
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: shat(:)
    type(multifab) , intent(in   ) :: snew(:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    type(multifab) , intent(in   ) :: hcoeff_old(:)
    type(multifab) , intent(in   ) :: Xkcoeff_old(:)
    type(multifab) , intent(in   ) :: pcoeff_old(:)
    type(multifab) , intent(in   ) :: hcoeff_new(:)
    type(multifab) , intent(in   ) :: Xkcoeff_new(:)
    type(multifab) , intent(in   ) :: pcoeff_new(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! Local
    type(multifab) :: alpha(mla%nlevel)
    type(multifab) :: beta(mla%nlevel,mla%dim)
    type(multifab) :: phi(mla%nlevel)
    type(multifab) :: Lphi(mla%nlevel)
    type(multifab) :: rhs(mla%nlevel)

    integer                     :: stencil_order,dm,nlevs
    integer                     :: i,n,comp
    type(bndry_reg)             :: fine_flx(2:mla%nlevel)
    real(dp_t)                  :: abs_eps, abs_solver_eps, rel_solver_eps

    type(bl_prof_timer), save :: bpt

    ! here we will solve:
    !
    ! (alpha - div beta grad) h = RHS
    !
    ! First we will construct the RHS by adding each of the terms in turn.  
    !
    ! To actually construct each div (c grad q) term for the RHS, we will 
    ! make use of the mac_applyop routine, which constructs the quantity
    !
    !     (alpha - div beta grad) phi
    !
    ! For all RHS terms, we set alpha = 0, beta = the appropriate 
    ! coefficient, and phi = the quantity being diffused.

    call build(bpt, "therm_cond_full_alg")

    dm = mla%dim
    nlevs = mla%nlevel

    stencil_order = 2

    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    do n=1,nlevs
       call multifab_build(alpha(n), mla%la(n),  1, 0)
       call multifab_build(phi(n)  , mla%la(n),  1, 1)
       call multifab_build(Lphi(n) , mla%la(n),  1, 0)
       call multifab_build(rhs(n)  , mla%la(n),  1, 0)
       do i = 1,dm
          call multifab_build_edge(beta(n,i), mla%la(n), 1, 0, i)
       end do
    end do

    ! set alpha to while we evaluate the diffusion terms for the rhs
    do n=1,nlevs
       call setval(alpha(n), ZERO, all=.true.)
       call setval(rhs(n), ZERO, all=.true.)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add h^old diffusion to rhs (using old coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! put beta on faces
    call put_data_on_faces(mla,hcoeff_old,1,beta,.true.)

    ! set phi = h^old
    do n=1,nlevs
       call multifab_copy_c(phi(n),1,sold(n),rhoh_comp,1,1)
       call multifab_div_div_c(phi(n),1,sold(n),rho_comp,1,1)
    end do

    ! apply the operator
    call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                     dm+rhoh_comp,stencil_order,mla%mba%rr)

    ! add Lphi to rhs
    do n=1,nlevs
       call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,0)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! subtract h^new enthalpy from rhs (using new coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! put beta on faces
    call put_data_on_faces(mla,hcoeff_new,1,beta,.true.)

    ! set phi = h^new
    do n=1,nlevs
       call multifab_copy_c(phi(n),1,snew(n),rhoh_comp,1,1)
       call multifab_div_div_c(phi(n),1,snew(n),rho_comp,1,1)
    end do

    ! apply the operator
    call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                     dm+rhoh_comp,stencil_order,mla%mba%rr)

    ! subtract Lphi from rhs
    do n=1,nlevs
       call multifab_sub_sub_c(rhs(n),1,Lphi(n),1,1,0)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add X_k^old diffusion to rhs (using old coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do comp=1,nspec

       ! put beta on faces
       call put_data_on_faces(mla,Xkcoeff_old,comp,beta,.true.)

       ! set phi = X_k^old
       do n=1,nlevs
          call multifab_copy_c(phi(n),1,sold(n),spec_comp+comp-1,1,1)
          call multifab_div_div_c(phi(n),1,sold(n),rho_comp,1,1)
       enddo

       ! apply the operator
       call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                        dm+spec_comp+comp-1,stencil_order,mla%mba%rr)

       ! add Lphi to rhs
       do n=1,nlevs
          call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,0)
       enddo

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add X_k^hat diffusion to rhs (using new coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do comp=1,nspec

       ! put beta on faces
       call put_data_on_faces(mla,Xkcoeff_new,comp,beta,.true.)

       ! set phi = X_k^new
       do n=1,nlevs
          call multifab_copy_c(phi(n),1,shat(n),spec_comp+comp-1,1,1)
          call multifab_div_div_c(phi(n),1,shat(n),rho_comp,1,1)
       enddo

       ! apply the operator
       call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                        dm+spec_comp+comp-1,stencil_order,mla%mba%rr)

       ! add Lphi to rhs
       do n=1,nlevs
          call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,0)
       enddo

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add p0^old diffusion to rhs (using old coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! put beta on faces
    call put_data_on_faces(mla,pcoeff_old,1,beta,.true.)

    ! set phi = p0_old
    call put_1d_array_on_cart(p0_old,phi,foextrap_comp,.false.,.false., &
                              dx,the_bc_tower%bc_tower_array,mla)

    ! apply the operator
    call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                     foextrap_comp,stencil_order,mla%mba%rr)

    ! add Lphi to rhs
    do n=1,nlevs
       call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,0)
    enddo    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add p0^new diffusion to rhs (using new coefficients)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! put beta on faces
    call put_data_on_faces(mla,pcoeff_new,1,beta,.true.)

    ! set phi = p0_new
    call put_1d_array_on_cart(p0_new,phi,foextrap_comp,.false.,.false., &
                              dx,the_bc_tower%bc_tower_array,mla)

    ! apply the operator
    call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower, &
                     foextrap_comp,stencil_order,mla%mba%rr)

    ! add Lphi to rhs
    do n=1,nlevs
       call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,0)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! scale rhs and add remaining terms
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! multiply rhs by 0.5*dt
    do n=1,nlevs
       call multifab_mult_mult_s(rhs(n),0.5d0*dt,0)
    end do

    ! add (rhoh)^old + dt*(aofs + intra) to rhs
    do n=1,nlevs
       call multifab_plus_plus_c(rhs(n),1,shat(n),rhoh_comp,1,0)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setup LHS coefficients for implicit solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! copy rho^hat\rho^new into alpha
    do n=1,nlevs
       call multifab_copy_c(alpha(n),1,shat(n),rho_comp,1,0)
    enddo

    ! copy new coefficients into beta
    call put_data_on_faces(mla,hcoeff_new,1,beta,.true.)

    ! multiply beta by -dt
    do n=1,nlevs
       do i = 1,dm 
          call multifab_mult_mult_s_c(beta(n,i),1,-dt,1,0)
       enddo
    enddo

    ! initialize phi to a reasonable guess, (rhoh)^hat / rho^hat
    do n=1,nlevs
       call multifab_copy_c(phi(n),1,shat(n),rhoh_comp,1,1)
       call multifab_div_div_c(phi(n),1,shat(n),rho_comp,1,1)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now do the implicit solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Compute norm(phi) to be used inside the MG solver as part of a stopping criterion
    abs_eps = -1.d0
    do n = 1,nlevs
       abs_eps = max(abs_eps, norm_inf(phi(n)) / dx(n,1))
    end do
    abs_solver_eps = eps_mac * abs_eps

    rel_solver_eps = eps_mac

    ! Call the solver to obtain h^hat (it will be stored in phi)
    ! solves (alpha - div beta grad) phi = rhs
    call mac_multigrid(mla,rhs,phi,fine_flx,alpha,beta,dx,the_bc_tower, &
                       dm+rhoh_comp,stencil_order,mla%mba%rr,rel_solver_eps,abs_solver_eps)

    ! load (rho*h)^hat into shat
    do n=1,nlevs
       call multifab_copy_c(shat(n),rhoh_comp,phi(n),1,1,0)
       call multifab_mult_mult_c(shat(n),rhoh_comp,shat(n),rho_comp,1,0)
    enddo

    do n=2,nlevs
       call destroy(fine_flx(n))
    end do

    do n=1,nlevs
       call destroy(alpha(n))
       call destroy(phi(n))
       call destroy(Lphi(n))
       call destroy(rhs(n))
       do i = 1,dm
          call destroy(beta(n,i))
       end do
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(shat(nlevs),rhoh_comp,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(shat(nlevs),rhoh_comp,dm+rhoh_comp,1, &
                            the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(shat(n-1),rhoh_comp,shat(n),rhoh_comp,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(shat(n),shat(n-1),nghost(shat(1)),mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         rhoh_comp,dm+rhoh_comp,1,fill_crse_input=.false.)
       enddo

    end if

    call destroy(bpt)

  end subroutine thermal_conduct_corrector

  subroutine make_explicit_thermal_hterm(mla,dx,thermal,s,hcoeff,the_bc_tower)

    use bc_module
    use bl_prof_module
    use cc_stencil_module
    use mac_applyop_module
    use network, only: nspec
    use ml_restriction_module, only : ml_cc_restriction
    use multifab_fill_ghost_module
    use bl_constants_module
    use variables, only: rho_comp, rhoh_comp
    use multifab_physbc_module
    use fill_3d_module

    type(ml_layout), intent(inout) :: mla
    real(dp_t)     , intent(in   ) :: dx(:,:)
    type(multifab) , intent(inout) :: thermal(:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: hcoeff(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! Local
    type(multifab) :: phi(mla%nlevel),alpha(mla%nlevel),beta(mla%nlevel,mla%dim)
    type(multifab) :: resid(mla%nlevel)

    integer :: comp,i,n,stencil_order,dm,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_explicit_thermal")

    dm = mla%dim
    nlevs = mla%nlevel

    stencil_order = 2

    do n=1,nlevs
       call setval(thermal(n), ZERO, all=.true.)
    end do

       ! compute thermal = del dot ( hcoeff grad h)
       
       do n=1,nlevs
          call multifab_build(phi(n),  mla%la(n), 1,  1)
          do i = 1,dm
             call multifab_build_edge(beta(n,i), mla%la(n), 1, 1, i)
          end do
       end do

       do n=1,nlevs
          ! load h into phi
          call multifab_copy_c(phi(n),1,s(n),rhoh_comp,1,1)
          call multifab_div_div_c(phi(n),1,s(n),rho_comp,1,1)
       end do

       call put_data_on_faces(mla,hcoeff,1,beta,.true.)

       do n=1,nlevs
          call multifab_build(alpha(n), mla%la(n), 1, 1)
          call multifab_build(resid(n), mla%la(n), 1, 0)
          call setval(alpha(n), ZERO, all=.true.)
       end do
       
       ! applyop to compute resid = del dot hcoeff grad h
       call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
                        stencil_order,mla%mba%rr)
       
       ! add residual to thermal
       do n=1,nlevs
          call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
       enddo
       
       do n=1,nlevs
          call destroy(resid(n))
          call destroy(phi(n))
          call destroy(alpha(n))
          do i = 1,dm 
             call destroy(beta(n,i))
          end do
       end do
    
    call destroy(bpt)
    
  end subroutine make_explicit_thermal_hterm

end module thermal_conduct_module
