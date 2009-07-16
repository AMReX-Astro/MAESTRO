module initial_proj_module

  implicit none

  private
  public :: initial_proj

contains

  subroutine initial_proj(uold,sold,pres,gpres,Source_old,hgrhs, &
                          div_coeff_old,p0,gamma1bar,dx,the_bc_tower,mla)

    use variables, only: press_comp, foextrap_comp, rho_comp
    use network, only: nspec
    use define_bc_module
    use bl_constants_module
    use probin_module
    use geometry, only: spherical, nr_fine, nlevs, nlevs_radial
    use proj_parameters, only: initial_projection_comp
    use make_explicit_thermal_module
    use make_S_module
    use average_module
    use hgrhs_module
    use fill_3d_module
    use hgproject_module
    use multifab_module
    use ml_layout_module
    use heating_module
    use inlet_bc_module, only: update_inlet_bcs

    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: pres(:)
    type(multifab) , intent(inout) :: gpres(:)
    type(multifab) , intent(inout) :: Source_old(:)
    type(multifab) , intent(inout) :: hgrhs(:)
    real(kind=dp_t), intent(in   ) :: div_coeff_old(:,0:)
    real(kind=dp_t), intent(inout) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(ml_layout), intent(inout) :: mla

    ! local
    integer    :: n
    real(dp_t) :: dt_temp

    type(multifab) :: delta_gamma1_term(nlevs)
    type(multifab) :: delta_gamma1(nlevs)
    type(multifab) :: thermal(nlevs)
    type(multifab) :: rhohalf(nlevs)
    type(multifab) :: rho_omegadot1(nlevs)
    type(multifab) :: rho_Hnuc1(nlevs)
    type(multifab) :: rho_Hext(nlevs)
    type(multifab) :: div_coeff_3d(nlevs)
    type(multifab) :: Tcoeff(nlevs)
    type(multifab) :: hcoeff(nlevs)
    type(multifab) :: Xkcoeff(nlevs)
    type(multifab) :: pcoeff(nlevs)

    real(dp_t) ::                  psi(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::                 Sbar(nlevs_radial,0:nr_fine-1)
    real(dp_t) :: delta_gamma1_termbar(nlevs_radial,0:nr_fine-1)

    Sbar = ZERO
    psi = ZERO

    if ( parallel_IOProcessor() ) then
       print *, 'DOING THE INITIAL VELOCITY PROJECTION'
       print *, ' '
    end if

    do n=1,nlevs
       call multifab_build(thermal(n), mla%la(n), 1, 1)
    end do
    
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

    ! set the inlet BCs, if any
    call update_inlet_bcs(ZERO,dx,mla,sold,rho_Hnuc1,rho_omegadot1,.false.)

    call make_S(Source_old,delta_gamma1_term,delta_gamma1, &
                sold,uold,rho_omegadot1,rho_Hnuc1,rho_Hext,thermal, &
                p0,gamma1bar,delta_gamma1_termbar,psi,dx, &
                mla,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(thermal(n))
       call destroy(rho_omegadot1(n))
       call destroy(rho_Hnuc1(n))
       call destroy(rho_Hext(n))
       call destroy(delta_gamma1(n))
    end do
    
    if (evolve_base_state) then
       call average(mla,Source_old,Sbar,dx,1)
    end if
    
    ! Note that we use rhohalf, filled with 1 at this point, as a temporary
    ! in order to do a constant-density initial projection.

    do n=1,nlevs
       call multifab_build(rhohalf(n), mla%la(n), 1, 1)
       call setval(rhohalf(n),ONE,1,1,all=.true.)
    end do
    
    call make_hgrhs(the_bc_tower,mla,hgrhs,Source_old,delta_gamma1_term,Sbar, &
                    div_coeff_old,dx)

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
    end do

    ! dt doesn't matter for the initial projection since we're throwing
    ! away the p and gpres anyway
    dt_temp = ONE

    if (spherical .eq. 1) then
       do n=1,nlevs
          call multifab_build(div_coeff_3d(n), mla%la(n), 1, 0)
       end do

       call put_1d_array_on_cart(div_coeff_old,div_coeff_3d,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla)

       call hgproject(initial_projection_comp,mla,uold,uold,rhohalf,pres,gpres,dx, &
                      dt_temp,the_bc_tower,press_comp, &
                      hgrhs,div_coeff_3d=div_coeff_3d,eps_in=1.d-10)
       
    else
       call hgproject(initial_projection_comp,mla,uold,uold,rhohalf,pres,gpres,dx, &
                      dt_temp,the_bc_tower,press_comp, &
                      hgrhs,div_coeff_1d=div_coeff_old)
    end if

    if(spherical .eq. 1) then
       do n=1,nlevs
          call destroy(div_coeff_3d(n))
       end do
    end if
    
    do n=1,nlevs
       call destroy(rhohalf(n))
    end do

    do n=1,nlevs
       call setval( pres(n), 0.0_dp_t, all=.true.)
       call setval(gpres(n), 0.0_dp_t, all=.true.)
    end do
    
  end subroutine initial_proj

end module initial_proj_module
