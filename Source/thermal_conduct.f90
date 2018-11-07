! thermal_conduct implements thermal diffusion in the enthalpy equation.
! This is an implicit solve, using the multigrid solver.  This updates
! the enthalpy only. 

module thermal_conduct_module

  use bl_types
  use multifab_module
  use bl_constants_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: thermal_conduct

contains 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Crank-Nicholson solve for enthalpy, taking into account only the
  ! enthalpy-diffusion terms in the temperature conduction term.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine thermal_conduct(mla,dx,dt,s1,hcoeff1,Xkcoeff1,pcoeff1, &
                             hcoeff2,Xkcoeff2,pcoeff2,s2,p0_old,p0_new,the_bc_tower,&
                             snew,derivative)

    use bl_prof_module
    use bndry_reg_module
    use ml_restrict_fill_module

    use mac_multigrid_module , only : mac_multigrid
    use cc_applyop_module   ,  only : cc_applyop
    use fill_3d_module       , only : put_1d_array_on_cart, put_data_on_faces
    use ml_cc_restriction_module, only : ml_cc_restriction_c

    use variables    , only : foextrap_comp, rho_comp, spec_comp, rhoh_comp
    use network      , only : nspec
    use probin_module, only : thermal_diffusion_type
    use mg_eps_module, only : eps_mac

    type(ml_layout), intent(inout) :: mla
    real(dp_t)     , intent(in   ) :: dx(:,:),dt
    type(multifab) , intent(in   ) :: s1(:)
    type(multifab) , intent(in   ) :: hcoeff1(:)
    type(multifab) , intent(in   ) :: Xkcoeff1(:)
    type(multifab) , intent(in   ) :: pcoeff1(:)
    type(multifab) , intent(in   ) :: hcoeff2(:)
    type(multifab) , intent(in   ) :: Xkcoeff2(:)
    type(multifab) , intent(in   ) :: pcoeff2(:)
    type(multifab) , intent(inout) :: s2(:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:),p0_new(:,0:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout), optional :: snew(:)
    logical        , intent(in   ), optional :: derivative

    ! Local
    type(multifab) :: rhs_solver_alpha(mla%nlevel),lhs_solver_alpha(mla%nlevel)
    type(multifab) :: rhs_solver_beta(mla%nlevel,mla%dim),lhs_solver_beta(mla%nlevel,mla%dim)
    type(multifab) :: phi(mla%nlevel),Lphi(mla%nlevel),rhs(mla%nlevel)

    integer                     :: stencil_order,dm,nlevs
    integer                     :: i,n,comp
    type(bndry_reg)             :: fine_flx(mla%nlevel)
    real(dp_t)                  :: abs_eps, abs_solver_eps, rel_solver_eps

    type(bl_prof_timer), save :: bpt

    call build(bpt, "therm_cond_full_alg")

    dm = mla%dim
    nlevs = mla%nlevel

    stencil_order = 2

    do n = 1,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do


    ! Here we will solve equations (A2) and (A3) in 
    ! Malone et al. ApJ 728:118, 2011
    ! for thermal diffusion.  The general form looks something like:
    !
    ! (rho - dt/2 div . hcoeff grad ) h = (rho h)  + 
    !        dt/2 div .       ( hcoeff grad h^1) -
    !        dt/2 sum_k div . (Xkcoeff grad X_k + Xkcoeff grad X_k) -
    !        dt/2 div .       ( pcoeff grad p_0 +  pcoeff grad p_0)
    !
    ! or 
    !
    ! (lhs_solver_alpha - div . lhs_solver_beta grad) h = RHS
    !
    ! First we will construct the RHS by adding each of the terms in
    ! turn.  
    !
    ! To actually construct each div . (c grad q) term for the RHS, we will 
    ! make use of the cc_applyop routine, which constructs the quantity
    !
    !     (rhs_solver_alpha - div . rhs_solver_beta grad) phi = Lphi
    !
    ! For all RHS terms, we set rhs_solver_alpha = 0, rhs_solver_beta = the appropriate 
    ! coefficient, and phi = the quantity being diffused.


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add enthalpy diffusion to rhs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       do i = 1,dm
          call multifab_build_edge(rhs_solver_beta(n,i), mla%la(n), 1, 1, i)
          call setval(rhs_solver_beta(n,i), ZERO, all=.true.)
       end do
    end do

    ! put solver_beta on faces
    call put_data_on_faces(mla,hcoeff1,1,rhs_solver_beta,.true.)

    ! scale by dt/2
    if (thermal_diffusion_type .eq. 1) then
       do n=1,nlevs
          do i = 1,dm 
             call multifab_mult_mult_s_c(rhs_solver_beta(n,i),1,dt/2.d0,1,1)
          enddo
       enddo
    else
       do n=1,nlevs
          do i = 1,dm 
             call multifab_mult_mult_s_c(rhs_solver_beta(n,i),1,dt,1,1)
          end do
       enddo
    end if

    ! set rhs_solver_alpha = 0
    ! set phi = h^{(1)}
    do n=1,nlevs
       call multifab_build(phi(n), mla%la(n),  1, 1)
       call multifab_build(rhs_solver_alpha(n), mla%la(n),  1, 1)
       call multifab_build(Lphi(n), mla%la(n),  1, 0)
       call setval(rhs_solver_alpha(n), ZERO, all=.true.)
       call multifab_copy_c(phi(n),1,s1(n),rhoh_comp,1,1)
       call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,1)
    end do

    ! apply the operator
    call cc_applyop(mla,Lphi,phi,rhs_solver_alpha,rhs_solver_beta,dx,the_bc_tower, &
                    dm+rhoh_comp,stencil_order)

    ! begin construction of rhs by setting rhs = \rho^{(2)}h^{(2')}
    do n=1,nlevs
       call multifab_build(rhs(n), mla%la(n),  1, 0)
       call multifab_copy_c(rhs(n),1,s2(n),rhoh_comp,1)
    enddo

    if(thermal_diffusion_type .eq. 1) then
       ! add Lphi to rhs
       do n=1,nlevs
          call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
       enddo
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add species diffusion to rhs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! loop over species
    do comp=1,nspec

       ! do X_k^{(1)} term first
       ! put solver_beta on faces
       call put_data_on_faces(mla,Xkcoeff1,comp,rhs_solver_beta,.true.)

       ! scale by dt/2
       if (thermal_diffusion_type .eq. 1) then
          do n=1,nlevs
             do i = 1,dm 
                call multifab_mult_mult_s_c(rhs_solver_beta(n,i),1,dt/2.d0,1,1)
             enddo
          enddo
       else
          do n=1,nlevs
             do i = 1,dm 
                call multifab_mult_mult_s_c(rhs_solver_beta(n,i),1,dt,1,1)
             enddo
          enddo
       end if

       ! load phi = X_k^{(1)}
       do n=1,nlevs
          call multifab_copy_c(phi(n),1,s1(n),spec_comp+comp-1,1,1)
          call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,1)
       enddo

       ! apply the operator
       call cc_applyop(mla,Lphi,phi,rhs_solver_alpha,rhs_solver_beta,dx,the_bc_tower, &
                       dm+spec_comp+comp-1,stencil_order)

       if(thermal_diffusion_type .eq. 1) then
          ! add lphi to rhs
          do n=1,nlevs
             call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
          enddo
       end if

       ! now do X_k^{(2)} term
       ! put solver_beta on faces
       call put_data_on_faces(mla,Xkcoeff2,comp,rhs_solver_beta,.true.)

       ! scale by dt/2
       if (thermal_diffusion_type .eq. 1) then
          do n=1,nlevs
             do i = 1,dm 
                call multifab_mult_mult_s_c(rhs_solver_beta(n,i),1,dt/2.d0,1,1)
             end do
          enddo
       else
          do n=1,nlevs
             do i = 1,dm 
                call multifab_mult_mult_s_c(rhs_solver_beta(n,i),1,dt,1,1)
             enddo
          enddo
       end if

       ! load phi = X_k^{(2)}
       do n=1,nlevs
          call multifab_copy_c(phi(n),1,s2(n),spec_comp+comp-1,1,1)
          call multifab_div_div_c(phi(n),1,s2(n),rho_comp,1,1)
       enddo

       ! apply the operator
       call cc_applyop(mla,Lphi,phi,rhs_solver_alpha,rhs_solver_beta,dx,the_bc_tower, &
                       dm+spec_comp+comp-1,stencil_order)

       ! add lphi to rhs
       do n=1,nlevs
          call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
       enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add pressure diffusion to rhs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! do p0_old term first
    ! put solver_beta on faces
    call put_data_on_faces(mla,pcoeff1,1,rhs_solver_beta,.true.)

    ! scale by dt/2
    if (thermal_diffusion_type .eq. 1) then
       do n=1,nlevs
          do i = 1,dm 
             call multifab_mult_mult_s_c(rhs_solver_beta(n,i),1,dt/2.d0,1,1)
          end do
       enddo
    else
       do n=1,nlevs
          do i = 1,dm 
             call multifab_mult_mult_s_c(rhs_solver_beta(n,i),1,dt,1,1)
          end do
       enddo
    end if

    call put_1d_array_on_cart(p0_old,phi,foextrap_comp,.false.,.false., &
                              dx,the_bc_tower%bc_tower_array,mla)
    ! apply the operator
    call cc_applyop(mla,Lphi,phi,rhs_solver_alpha,rhs_solver_beta,dx,the_bc_tower, &
                    foextrap_comp,stencil_order)

    if(thermal_diffusion_type .eq. 1) then
       ! add lphi to rhs
       do n=1,nlevs
          call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
       enddo
    end if

    ! now do p0_new term
    ! put solver_beta on faces
    call put_data_on_faces(mla,pcoeff2,1,rhs_solver_beta,.true.)

    ! scale by dt/2
    if (thermal_diffusion_type .eq. 1) then
       do n=1,nlevs
          do i = 1,dm 
             call multifab_mult_mult_s_c(rhs_solver_beta(n,i),1,dt/2.d0,1,1)
          end do
       enddo
    else
       do n=1,nlevs
          do i = 1,dm 
             call multifab_mult_mult_s_c(rhs_solver_beta(n,i),1,dt,1,1)
          end do
       enddo
    end if

    call put_1d_array_on_cart(p0_new,phi,foextrap_comp,.false.,.false., &
                              dx,the_bc_tower%bc_tower_array,mla)

    ! apply the operator
    call cc_applyop(mla,Lphi,phi,rhs_solver_alpha,rhs_solver_beta,dx,the_bc_tower, &
                    foextrap_comp,stencil_order)

    do n=1,nlevs
       call destroy(rhs_solver_alpha(n))
       do i = 1,dm
          call destroy(rhs_solver_beta(n,i))
       end do
    end do

    ! add lphi to rhs
    do n=1,nlevs
       call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
    enddo

    do n=1,nlevs
       call destroy(Lphi(n))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setup LHS coefficients
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       do i = 1,dm
          call multifab_build_edge(lhs_solver_beta(n,i), mla%la(n), 1, 1, i)
          call setval(lhs_solver_beta(n,i), ZERO, all=.true.)
       end do
    end do

    ! create lhs_solver_beta = -hcoeff2 = (dt/2)k_{th}^{(2'')}/c_p^{(2'')}
    ! put solver_beta on faces (remember to scale by -dt/2 afterwards)
    call put_data_on_faces(mla,hcoeff2,1,lhs_solver_beta,.true.)

    ! scale by -dt/2
    if (thermal_diffusion_type .eq. 1) then
       do n=1,nlevs
          do i = 1,dm
             call multifab_mult_mult_s_c(lhs_solver_beta(n,i),1,-dt/2.d0,1,1)
          end do
       enddo
    else
       do n=1,nlevs
          do i = 1,dm
             call multifab_mult_mult_s_c(lhs_solver_beta(n,i),1,-dt,1,1)
          end do
       enddo
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now do the implicit solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! initialize phi to h^{(2'')} as a guess; also sets the ghost cells at inflow/outflow
    ! to a reasonable value
    do n=1,nlevs
       call multifab_copy_c(phi(n),1,s2(n),rhoh_comp,1,1)
       call multifab_div_div_c(phi(n),1,s2(n),rho_comp,1,1)
    enddo

    ! lhs_solver_alpha = \rho^{(2),*} or \rho^{(2)}
    do n=1,nlevs
       call multifab_build(lhs_solver_alpha(n), mla%la(n), 1, 1)
       call multifab_copy_c(lhs_solver_alpha(n),1,s2(n),rho_comp,1,1)
    enddo

    ! Compute norm(phi) to be used inside the MG solver as part of a stopping criterion
    abs_eps = -1.d0
    do n = 1,nlevs
       abs_eps = max(abs_eps, norm_inf(phi(n)) / dx(n,1))
    end do
    abs_solver_eps = eps_mac * abs_eps

    rel_solver_eps = eps_mac

    ! Call the solver to obtain h^(2) (it will be stored in phi)
    ! solves (solver_alpha - nabla dot solver_beta nabla)phi = rhs
    call mac_multigrid(mla,rhs,phi,fine_flx,lhs_solver_alpha,lhs_solver_beta,dx,the_bc_tower, &
                       dm+rhoh_comp,stencil_order,rel_solver_eps,abs_solver_eps)

    do n = 1,nlevs
       call destroy(fine_flx(n))
    end do

    do n=1,nlevs
       call destroy(lhs_solver_alpha(n))
       do i = 1,dm
          call destroy(lhs_solver_beta(n,i))
       end do
       call destroy(rhs(n))
    end do
    
    if (present(derivative) .and. derivative) then
        ! load new rho*h - old rho*h into snew
        do n=1,nlevs
           call multifab_mult_mult_c(phi(n),1,s2(n),rho_comp,1)
           call multifab_div_div_c(phi(n),1,s2(n),rhoh_comp,1)
           call multifab_plus_plus_c(snew(n),rhoh_comp,phi(n),1,1)
        enddo
    else    
        ! load new rho*h into s2
        do n=1,nlevs
           call multifab_copy_c(s2(n),rhoh_comp,phi(n),1,1)
           call multifab_mult_mult_c(s2(n),rhoh_comp,s2(n),rho_comp,1)
        enddo
    endif
    
    do n=1,nlevs
       call destroy(phi(n))
    end do

    ! restrict data and fill all ghost cells
    if (present(derivative) .and. derivative) then
       call ml_restrict_and_fill(nlevs,snew,mla%mba%rr,the_bc_tower%bc_tower_array, &
                              icomp=rhoh_comp, &
                              bcomp=dm+rhoh_comp, &
                              nc=1, &
                              ng=snew(1)%ng)
    else
       call ml_restrict_and_fill(nlevs,s2,mla%mba%rr,the_bc_tower%bc_tower_array, &
                              icomp=rhoh_comp, &
                              bcomp=dm+rhoh_comp, &
                              nc=1, &
                              ng=s2(1)%ng)
    endif
    
    call destroy(bpt)

  end subroutine thermal_conduct

end module thermal_conduct_module
