! The initial projection creates a first approximation to the velocity
! field by forcing the initial velocity field set by initveldata to
! satisfy the elliptic constraint equation.  Since the initial
! velocity may be zero, there is no guarantee that a well-defined
! timestep can be computed at this point, so the source term, S, used
! here only involves thermal diffusion and any external heating term,
! H_ext---no reactions are included.
!
! see paper III, section 3.3

module initial_proj_module

  use bl_constants_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private
  public :: initial_proj

contains

  subroutine initial_proj(uold,sold,pi,gpi,S_cc,normal,S_nodal,thermal, &
                          div_coeff_old,p0,gamma1bar,dx,the_bc_tower,mla)

    use variables, only: foextrap_comp
    use network, only: nspec
    use probin_module, only: use_thermal_diffusion, evolve_base_state
    use geometry, only: spherical, nr_fine, nlevs_radial
    use proj_parameters, only: initial_projection_comp
    use make_explicit_thermal_module
    use make_S_cc_module
    use average_module
    use make_S_nodal_module
    use fill_3d_module
    use hgproject_module
    use multifab_module
    use heating_module
    use mg_eps_module, only: eps_init_proj_cart, eps_init_proj_sph

    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: pi(:)
    type(multifab) , intent(inout) :: gpi(:)
    type(multifab) , intent(inout) :: S_cc(:)
    type(multifab) , intent(inout) :: normal(:)
    type(multifab) , intent(inout) :: S_nodal(:)
    type(multifab) , intent(inout) :: thermal(:)
    real(kind=dp_t), intent(in   ) :: div_coeff_old(:,0:)
    real(kind=dp_t), intent(inout) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(ml_layout), intent(inout) :: mla

    ! local
    integer    :: n,nlevs
    real(dp_t) :: dt_temp,eps_init

    type(multifab) :: delta_gamma1_term(mla%nlevel)
    type(multifab) :: delta_gamma1(mla%nlevel)
    type(multifab) :: rhohalf(mla%nlevel)
    type(multifab) :: rho_omegadot1(mla%nlevel)
    type(multifab) :: rho_Hnuc1(mla%nlevel)
    type(multifab) :: rho_Hext(mla%nlevel)
    type(multifab) :: div_coeff_cart(mla%nlevel)
    type(multifab) :: Tcoeff(mla%nlevel)
    type(multifab) :: hcoeff(mla%nlevel)
    type(multifab) :: Xkcoeff(mla%nlevel)
    type(multifab) :: pcoeff(mla%nlevel)

    real(dp_t) ::                  psi(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::                 Sbar(nlevs_radial,0:nr_fine-1)
    real(dp_t) :: delta_gamma1_termbar(nlevs_radial,0:nr_fine-1)

    nlevs = mla%nlevel

    Sbar = ZERO
    psi = ZERO

111 format(78('-'))
    if ( parallel_IOProcessor() ) then
       write (*,111)
       write (*,*) 'DOING THE INITIAL VELOCITY PROJECTION'
       write (*,111)
       write (*,*) ' '
    end if

    if(use_thermal_diffusion) then

       do n=1,nlevs
          call multifab_build(Tcoeff(n),  mla%la(n), 1,     1)
          call multifab_build(hcoeff(n),  mla%la(n), 1,     1)
          call multifab_build(Xkcoeff(n), mla%la(n), nspec, 1)
          call multifab_build(pcoeff(n),  mla%la(n), 1,     1)
       end do

       call make_thermal_coeffs(sold,Tcoeff,hcoeff,Xkcoeff,pcoeff)

       call make_explicit_thermal(mla,dx,thermal,sold,Tcoeff,hcoeff,Xkcoeff,pcoeff, &
                                  p0,the_bc_tower)

       do n=1,nlevs
          call destroy(Tcoeff(n))
          call destroy(hcoeff(n))
          call destroy(Xkcoeff(n))
          call destroy(pcoeff(n))
       end do

    else
       do n=1,nlevs
          call setval(thermal(n), ZERO, all=.true.)
       end do
    end if
    
    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1,     1)
       call multifab_build(delta_gamma1(n),      mla%la(n), 1,     1)
       call multifab_build(rho_omegadot1(n),     mla%la(n), nspec, 0)
       call multifab_build(rho_Hnuc1(n),         mla%la(n), 1,     0)
       call multifab_build(rho_Hext(n),          mla%la(n), 1,     0)
       ! we don't have a legit timestep yet, so we set rho_omegadot1, rho_Hnuc1,
       ! and rho_Hext to 0 
       call setval(     rho_omegadot1(n), ZERO, all=.true.)
       call setval(         rho_Hnuc1(n), ZERO, all=.true.)
       call setval(          rho_Hext(n), ZERO, all=.true.)
       call setval(      delta_gamma1(n), ZERO, all=.true.)
       call setval( delta_gamma1_term(n), ZERO, all=.true.)
    end do

    call make_S_cc(S_cc,delta_gamma1_term,delta_gamma1, &
                   sold,uold, &
                   normal, &
                   rho_omegadot1,rho_Hnuc1,rho_Hext,thermal, &
                   p0,gamma1bar,delta_gamma1_termbar,psi, &
                   dx,mla,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(rho_omegadot1(n))
       call destroy(rho_Hnuc1(n))
       call destroy(rho_Hext(n))
       call destroy(delta_gamma1(n))
    end do
    
    if (evolve_base_state) then
       call average(mla,S_cc,Sbar,dx,1)
    end if
    
    ! Note that we use rhohalf, filled with 1 at this point, as a temporary
    ! in order to do a constant-density initial projection.

    do n=1,nlevs
       call multifab_build(rhohalf(n), mla%la(n), 1, 1)
       call setval(rhohalf(n),ONE,1,1,all=.true.)
    end do
    
    call make_S_nodal(the_bc_tower,mla,S_nodal,S_cc,delta_gamma1_term,Sbar, &
                    div_coeff_old,dx)

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
    end do

    ! dt doesn't matter for the initial projection since we're throwing
    ! away the pi and gpi anyway
    dt_temp = ONE
    
    do n=1,nlevs
       call multifab_build(div_coeff_cart(n), mla%la(n), 1, 1)
    end do
    
    call put_1d_array_on_cart(div_coeff_old,div_coeff_cart,foextrap_comp,.false., &
                              .false.,dx,the_bc_tower%bc_tower_array,mla)

    if (spherical .eq. 1) then
       eps_init = eps_init_proj_sph
    else
       eps_init = eps_init_proj_cart
    end if

    call hgproject(initial_projection_comp,mla,uold,uold,rhohalf,pi,gpi,dx, &
                   dt_temp,the_bc_tower,div_coeff_cart,S_nodal,eps_init)
    
    do n=1,nlevs
       call destroy(div_coeff_cart(n))
       call destroy(rhohalf(n))
    end do

    do n=1,nlevs
       call setval( pi(n), 0.0_dp_t, all=.true.)
       call setval(gpi(n), 0.0_dp_t, all=.true.)
    end do
    
  end subroutine initial_proj

end module initial_proj_module
