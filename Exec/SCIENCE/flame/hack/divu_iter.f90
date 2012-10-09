module divu_iter_module

  implicit none

  private

  public :: divu_iter

contains

  subroutine divu_iter(istep_divu_iter,uold,sold,pres,gpres,normal, &
                       Source_old,hgrhs,dSdt,div_coeff_old,rho0_old,p0_old,gamma1bar, &
                       w0,grav_cell,dx,dt,time,the_bc_tower,mla)

    use variables, only: press_comp, nscal, foextrap_comp, rho_comp
    use network, only: nspec
    use define_bc_module
    use bl_constants_module
    use probin_module
    use geometry, only: spherical, nr_fine, dm, nlevs, nlevs_radial
    use proj_parameters, only: divu_iters_comp
    use react_state_module
    use make_explicit_thermal_module
    use make_S_module
    use average_module
    use hgrhs_module
    use fill_3d_module
    use hgproject_module
    use estdt_module
    use multifab_module
    use ml_layout_module
    use make_w0_module
    use inlet_bc_module, only: update_inlet_bcs

    integer        , intent(in   ) :: istep_divu_iter
    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: pres(:)
    type(multifab) , intent(inout) :: gpres(:)
    type(multifab) , intent(in   ) :: normal(:)
    type(multifab) , intent(inout) :: Source_old(:)
    type(multifab) , intent(inout) :: hgrhs(:)
    type(multifab) , intent(in   ) :: dSdt(:)
    real(kind=dp_t), intent(in   ) :: div_coeff_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(inout) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(inout) :: dt
    real(kind=dp_t), intent(in   ) :: time 
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(ml_layout), intent(inout) :: mla

    ! local
    integer        :: n,ng_s
    real(dp_t)     :: halfdt,dt_temp,dt_hold

    type(multifab) :: s1(nlevs)  
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

    real(dp_t) ::            etarho_ec(nlevs_radial,0:nr_fine)
    real(dp_t) ::                 Sbar(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::             w0_force(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::                  psi(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::            etarho_cc(nlevs_radial,0:nr_fine-1)
    real(dp_t) :: delta_gamma1_termbar(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::   p0_minus_pthermbar(nlevs_radial,0:nr_fine-1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "divu_iter")

    etarho_ec = ZERO
    Sbar = ZERO
    w0_force = ZERO
    psi = ZERO
    etarho_cc = ZERO
    delta_gamma1_termbar = ZERO
    p0_minus_pthermbar = ZERO

    halfdt = HALF*dt
    ng_s = sold(1)%ng

    do n = 1,nlevs
       call multifab_build(s1(n),            mla%la(n), nscal, ng_s)      
       call multifab_build(rho_omegadot1(n), mla%la(n), nspec, 0)
       call multifab_build(rho_Hnuc1(n),     mla%la(n), 1,     0)
       call multifab_build(rho_Hext(n),      mla%la(n), 1,     0)
    end do

    ! burn to define rho_omegadot and rho_Hnuc -- needed to make S
    call react_state(mla,sold,s1,rho_omegadot1,rho_Hnuc1,rho_Hext,p0_old, &
                     halfdt,dx,the_bc_tower%bc_tower_array,time)

    do n=1,nlevs
       call destroy(s1(n))
    end do

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
                                  p0_old,the_bc_tower)

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
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(delta_gamma1(n), mla%la(n), 1, 0)
    end do

    ! set the inlet BCs, if any
    call update_inlet_bcs(time,dx,mla,sold,rho_Hnuc1,rho_omegadot1,.false.)

    call make_S(Source_old,delta_gamma1_term,delta_gamma1, &
                sold,uold,rho_omegadot1,rho_Hnuc1,rho_Hext,thermal, &
                p0_old,gamma1bar,delta_gamma1_termbar,psi,dx, &
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
       call make_w0(w0,w0,w0_force,Sbar,rho0_old,rho0_old,p0_old,p0_old, &
                    gamma1bar,gamma1bar,p0_minus_pthermbar, &
                    psi,etarho_ec,etarho_cc,dt,dt)
    end if
    
    ! This needs to be a separate loop so Sbar is fully defined before 
    ! we get here.
    call make_hgrhs(the_bc_tower,mla,hgrhs,Source_old,delta_gamma1_term, &
                    Sbar,div_coeff_old,dx)
    
    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
    end do

    ! dt doesn't matter for the div iters since we're throwing
    ! away the p and gpres anyway
    dt_temp = ONE

    do n=1,nlevs
       call multifab_build(rhohalf(n), mla%la(n), 1, 1)
       call setval(rhohalf(n),ONE,1,1,all=.true.)
    end do

    if (spherical .eq. 1) then
       do n=1,nlevs
          call multifab_build(div_coeff_3d(n), mla%la(n), 1, 0)
       end do
       
       call put_1d_array_on_cart(div_coeff_old,div_coeff_3d,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla)

       call hgproject(divu_iters_comp,mla,uold,uold,rhohalf,pres,gpres,dx,dt_temp, &
                      the_bc_tower,press_comp, &
                      hgrhs,div_coeff_3d=div_coeff_3d,eps_in=1.d-10)
       
    else
       call hgproject(divu_iters_comp,mla,uold,uold,rhohalf,pres,gpres,dx,dt_temp, &
                      the_bc_tower,press_comp, &
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

    do n = 1,nlevs
       call setval(pres(n) ,0.0_dp_t, all=.true.)
       call setval(gpres(n),0.0_dp_t, all=.true.)
    end do
    
    dt_hold = dt
    dt      = HUGE(dt)

    call estdt(mla,the_bc_tower,uold,sold,gpres,Source_old,dSdt, &
               normal,w0,rho0_old,p0_old,gamma1bar,grav_cell,div_coeff_old, &
               dx,cflfac,dt)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       print*,"Call to estdt at end of istep_divu_iter =",istep_divu_iter
       print*,"gives dt =",dt
    end if
    
    dt = dt*init_shrink
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       print*, "Multiplying dt by init_shrink; dt =",dt
    end if
    
    if(dt .gt. dt_hold) then
       if ( parallel_IOProcessor() .and. verbose .ge. 1) then
          print*, "Ignoring this new dt since it's larger than the previous dt =", &
               dt_hold
       end if
       
       dt = min(dt_hold,dt)
    end if
    
    if (fixed_dt .ne. -1.0d0) then
       dt = fixed_dt
       if (parallel_IOProcessor()) then
          print*, "Setting fixed dt =",dt
       end if
    end if
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       print*,''
    end if

    call destroy(bpt)

  end subroutine divu_iter

end module divu_iter_module
