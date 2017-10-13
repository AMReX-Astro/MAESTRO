! The ``divu'' iterations projects the velocity field from the initial
! projection to satisfy the full constraint (including reactions).  This
! is an iterative process since the reactions depend on the timestep and
! the timestep depends on the velocity field.  The number of iterations
! to take is set through the init_divu_iter runtime parameter.

module divu_iter_module

  use bl_prof_module
  use bl_constants_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: divu_iter

contains

  subroutine divu_iter(istep_divu_iter,uold,sold,pi,gpi,thermal, &
                       S_cc,normal,hgrhs,dSdt,div_coeff_old,rho0_old,p0_old,gamma1bar, &
                       tempbar_init,w0,grav_cell,dx,dt,the_bc_tower,mla)

    use variables, only: nscal, foextrap_comp
    use network, only: nspec
    use probin_module, only: use_thermal_diffusion, evolve_base_state, &
         init_divu_iter, cflfac, verbose, init_shrink, fixed_dt
    use geometry, only: spherical, nr_fine, nlevs_radial
    use proj_parameters, only: divu_iters_comp
    use react_state_module
    use make_explicit_thermal_module
    use make_S_cc_module
    use average_module
    use hgrhs_module
    use fill_3d_module
    use hgproject_module
    use estdt_module
    use multifab_module
    use make_w0_module
    use mg_eps_module, only: eps_divu_cart, eps_divu_sph, &
         divu_iter_factor, divu_level_factor

    integer        , intent(in   ) :: istep_divu_iter
    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: pi(:)
    type(multifab) , intent(inout) :: gpi(:)
    type(multifab) , intent(inout) :: thermal(:)
    type(multifab) , intent(inout) :: S_cc(:)
    type(multifab) , intent(inout) :: normal(:)
    type(multifab) , intent(inout) :: hgrhs(:)
    type(multifab) , intent(in   ) :: dSdt(:)
    real(kind=dp_t), intent(in   ) :: div_coeff_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(:,0:)
    real(kind=dp_t), intent(inout) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(inout) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(ml_layout), intent(inout) :: mla

    ! local
    integer        :: n,ng_s,nlevs
    real(dp_t)     :: halfdt,dt_temp,dt_hold
    real(dp_t)     :: eps_divu

    type(multifab) :: s1(mla%nlevel)  
    type(multifab) :: delta_gamma1_term(mla%nlevel)
    type(multifab) :: delta_gamma1(mla%nlevel)
    type(multifab) :: rhohalf(mla%nlevel)
    type(multifab) :: rho_omegadot(mla%nlevel)
    type(multifab) :: rho_Hnuc(mla%nlevel)
    type(multifab) :: rho_Hext(mla%nlevel)
    type(multifab) :: div_coeff_3d(mla%nlevel)
    type(multifab) :: Tcoeff(mla%nlevel)
    type(multifab) :: hcoeff(mla%nlevel)
    type(multifab) :: Xkcoeff(mla%nlevel)
    type(multifab) :: pcoeff(mla%nlevel)

    real(dp_t) ::            etarho_ec(nlevs_radial,0:nr_fine)
    real(dp_t) ::                 Sbar(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::             w0_force(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::                  psi(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::            etarho_cc(nlevs_radial,0:nr_fine-1)
    real(dp_t) :: delta_gamma1_termbar(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::   p0_minus_pthermbar(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::         delta_chi_w0(nlevs_radial,0:nr_fine-1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "divu_iter")

    nlevs = mla%nlevel

    etarho_ec = ZERO
    Sbar = ZERO
    w0_force = ZERO
    psi = ZERO
    etarho_cc = ZERO
    delta_gamma1_termbar = ZERO
    p0_minus_pthermbar = ZERO

    halfdt = HALF*dt
    ng_s   = nghost(sold(1))

    do n = 1,nlevs
       call multifab_build(s1(n),           mla%la(n), nscal, ng_s)
       call multifab_build(rho_Hext(n),     mla%la(n), 1,     0)
       call multifab_build(rho_omegadot(n), mla%la(n), nspec, 0)
       call multifab_build(rho_Hnuc(n),     mla%la(n), 1,     0)
    end do

    ! burn to define rho_omegadot and rho_Hnuc -- needed to make S
    call react_state(mla,tempbar_init,sold,s1,rho_omegadot,rho_Hnuc,rho_Hext,p0_old, &
                     halfdt,dx,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(s1(n))
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
    end if
    
    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(delta_gamma1(n),      mla%la(n), 1, 0)
    end do

    call make_S_cc(S_cc,delta_gamma1_term,delta_gamma1, &
                   sold,uold, &
                   normal, &
                   rho_omegadot,rho_Hnuc,rho_Hext,thermal, &
                   p0_old,gamma1bar,delta_gamma1_termbar,psi, &
                   dx,mla,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(rho_omegadot(n))
       call destroy(rho_Hnuc(n))
       call destroy(rho_Hext(n))
       call destroy(delta_gamma1(n))
    end do

    if (evolve_base_state) then
       call average(mla,S_cc,Sbar,dx,1)
       call make_w0(w0,w0,w0_force,Sbar,rho0_old,rho0_old,p0_old,p0_old, &
                    gamma1bar,gamma1bar,p0_minus_pthermbar, &
                    psi,etarho_ec,etarho_cc,dt,dt,delta_chi_w0,.true.)
    end if
    
    ! This needs to be a separate loop so Sbar is fully defined before 
    ! we get here.
    call make_hgrhs(the_bc_tower,mla,hgrhs,S_cc,delta_gamma1_term, &
                    Sbar,div_coeff_old,dx)
    
    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
    end do

    ! dt doesn't matter for the div iters since we're throwing
    ! away the pi and gpi anyway
    dt_temp = ONE

    do n=1,nlevs
       call multifab_build(rhohalf(n), mla%la(n), 1, 1)
       call setval(rhohalf(n),ONE,1,1,all=.true.)
    end do

    if (spherical .eq. 1) then
       if (istep_divu_iter .eq. init_divu_iter) then
          eps_divu = eps_divu_sph

       else if (istep_divu_iter .eq. init_divu_iter-1) then
          eps_divu = eps_divu_sph*divu_iter_factor

       else if (istep_divu_iter .le. init_divu_iter-2) then
          eps_divu = eps_divu_sph*divu_iter_factor**2
       end if

    else
       if (istep_divu_iter .eq. init_divu_iter) then
          eps_divu = min(eps_divu_cart*divu_level_factor**(nlevs-1), &
                         eps_divu_cart*divu_level_factor**2)

       else if (istep_divu_iter .eq. init_divu_iter-1) then
          eps_divu = min(eps_divu_cart*divu_iter_factor*divu_level_factor**(nlevs-1), &
                         eps_divu_cart*divu_iter_factor*divu_level_factor**2)

       else if (istep_divu_iter .le. init_divu_iter-2) then
          eps_divu = min(eps_divu_cart*divu_iter_factor**2*divu_level_factor**(nlevs-1), &
                         eps_divu_cart*divu_iter_factor**2*divu_level_factor**2)
       endif
    end if

    do n=1,nlevs
       call multifab_build(div_coeff_3d(n), mla%la(n), 1, 1)
    end do
       
    call put_1d_array_on_cart(div_coeff_old,div_coeff_3d,foextrap_comp,.false., &
                              .false.,dx,the_bc_tower%bc_tower_array,mla)


    call hgproject(divu_iters_comp,mla,uold,uold,rhohalf,pi,gpi,dx,dt_temp, &
                   the_bc_tower,div_coeff_3d,hgrhs,eps_divu)
    
    do n=1,nlevs
       call destroy(div_coeff_3d(n))
       call destroy(rhohalf(n))
    end do

    do n = 1,nlevs
       call setval(pi(n) ,0.0_dp_t, all=.true.)
       call setval(gpi(n),0.0_dp_t, all=.true.)
    end do
    
    dt_hold = dt
    dt      = HUGE(dt)

    call estdt(mla,the_bc_tower,uold,sold,gpi,S_cc,dSdt, &
               w0,rho0_old,p0_old,gamma1bar,grav_cell,dx,cflfac,dt)

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
