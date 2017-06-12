module plot_variables_module

  use bl_types
  use multifab_module
  use define_bc_module
  use probin_module, only: use_tfromp
  use bl_constants_module, only: HALF
  use variables, only: plot_t

  implicit none

  private

  public :: make_ad_excess
  public :: make_conductivity
  public :: make_tfromH, make_tfromp, make_entropypert
  public :: make_deltaT, make_divw0, make_vorticity, make_magvel, make_velrc
  public :: make_rhopert, make_rhohpert
  public :: make_processor_number, make_pidivu

contains



  !---------------------------------------------------------------------------
  ! make_ad_excess
  !---------------------------------------------------------------------------
  subroutine make_ad_excess(plotdata,comp_ad_excess,state,normal)

    use geometry, only: spherical

    type(multifab), intent(inout) :: plotdata
    integer,        intent(in   ) :: comp_ad_excess
    type(multifab), intent(in   ) :: state
    type(multifab), intent(in   ) :: normal
    
    real(kind=dp_t), pointer :: sp(:,:,:,:), cp(:,:,:,:), np(:,:,:,:)
    integer :: lo(get_dim(plotdata)), hi(get_dim(plotdata)), ng_s, ng_c, ng_n
    integer :: i, dm

    dm = get_dim(plotdata)

    ng_c = nghost(plotdata)
    ng_s = nghost(state)
    ng_n = nghost(normal)

    do i = 1, nfabs(state)
       sp => dataptr(state, i)
       cp => dataptr(plotdata, i)
       lo = lwb(get_box(state, i))
       hi = upb(get_box(state, i))
       select case (dm)
       case (1)
          call make_ad_excess_1d(cp(:,1,1,comp_ad_excess), ng_c, &
                                 sp(:,1,1,:), ng_s, &
                                 lo, hi)
       case (2)
          call make_ad_excess_2d(cp(:,:,1,comp_ad_excess), ng_c, &
                                 sp(:,:,1,:), ng_s, &
                                 lo, hi)
       case (3)
          if (spherical .eq. 1) then
             np => dataptr(normal, i)       
             call make_ad_excess_3d_sphr(cp(:,:,:,comp_ad_excess), ng_c, &
                                         sp(:,:,:,:), ng_s, &
                                         np(:,:,:,:), ng_n, lo, hi)
          else
             call make_ad_excess_3d(cp(:,:,:,comp_ad_excess), ng_c, &
                                    sp(:,:,:,:), ng_s, &
                                    lo, hi)
          endif
       end select

    enddo
  end subroutine make_ad_excess
  
  subroutine make_ad_excess_1d(ad_excess, ng_ad, state, ng_s, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: ad_excess(lo(1)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,:)

    real(kind=dp_t) :: pres(lo(1):hi(1)), nabla_ad(lo(1):hi(1))
    real(kind=dp_t) :: chi_rho, chi_t, dt, dp, nabla

    integer :: i
    
    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    do i = lo(1), hi(1)

       eos_state%rho   = state(i,rho_comp)
       eos_state%T     = state(i,temp_comp)
       eos_state%xn(:) = state(i,spec_comp:spec_comp+nspec-1)/eos_state%rho

       pt_index(:) = (/i, -1, -1/)       

       call eos(eos_input_rt, eos_state, pt_index)
       
       pres(i) = eos_state%p

       chi_rho = eos_state%rho * eos_state%dpdr / eos_state%p
       chi_t = eos_state%T * eos_state%dpdt / eos_state%p
       nabla_ad(i) = (eos_state%gam1 - chi_rho) / (chi_t * eos_state%gam1)

    enddo

    do i = lo(1), hi(1)
       if (state(i,rho_comp) <= base_cutoff_density) then
          nabla = ZERO
       else
          ! forward difference
          if (i == lo(1)) then
             dt = state(i+1,temp_comp) - state(i,temp_comp)
             dp = pres(i+1) - pres(i)
             ! backward difference
          else if (i == hi(1)) then
             dt = state(i,temp_comp) - state(i-1,temp_comp)
             dp = pres(i) - pres(i-1)
             ! centered difference
          else
             dt = state(i+1,temp_comp) - state(i-1,temp_comp)
             dp = pres(i+1) - pres(i-1)
          endif

          ! prevent Inf
          if (dp == ZERO) then
             nabla = -huge(ZERO)
          else
             nabla = pres(i) * dt / (dp * state(i,temp_comp))
          endif
       endif

       ad_excess(i) = nabla - nabla_ad(i)
    enddo

  end subroutine make_ad_excess_1d

  subroutine make_ad_excess_2d(ad_excess, ng_ad, state, ng_s, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: ad_excess(lo(1)-ng_ad:,lo(2)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,:)

    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2))
    real(kind=dp_t) :: chi_rho, chi_t, dt, dp, nabla

    integer :: i, j
    
    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          eos_state%rho   = state(i,j,rho_comp)
          eos_state%T     = state(i,j,temp_comp)
          eos_state%xn(:) = state(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho

          pt_index(:) = (/i, j, -1/)       

          call eos(eos_input_rt, eos_state, pt_index)       
       
          pres(i,j) = eos_state%p

          chi_rho = eos_state%rho * eos_state%dpdr / eos_state%p
          chi_t = eos_state%T * eos_state%dpdt / eos_state%p
          nabla_ad(i,j) = (eos_state%gam1 - chi_rho) / (chi_t * eos_state%gam1)

       enddo
    enddo

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (state(i,j,rho_comp) <= base_cutoff_density) then
             nabla = ZERO
          else
             ! forward difference
             if (j == lo(2)) then
                dt = state(i,j+1,temp_comp) - state(i,j,temp_comp)
                dp = pres(i,j+1) - pres(i,j)
                ! backward difference
             else if (j == hi(2)) then
                dt = state(i,j,temp_comp) - state(i,j-1,temp_comp)
                dp = pres(i,j) - pres(i,j-1)
                ! centered difference
             else
                dt = state(i,j+1,temp_comp) - state(i,j-1,temp_comp)
                dp = pres(i,j+1) - pres(i,j-1)
             endif

             ! prevent Inf
             if (dp == ZERO) then
                nabla = -huge(ZERO)
             else
                nabla = pres(i,j) * dt / (dp * state(i,j,temp_comp))
             endif
          endif

          ad_excess(i,j) = nabla - nabla_ad(i,j)
       enddo
    enddo

  end subroutine make_ad_excess_2d

  subroutine make_ad_excess_3d(ad_excess, ng_ad, state, ng_s, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: ad_excess(lo(1)-ng_ad:,lo(2)-ng_ad:,lo(3)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_rho, chi_t, dt, dp, nabla

    integer :: i, j, k

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    !$OMP PARALLEL DO PRIVATE(i,j,k,chi_rho,chi_t,eos_state,pt_index)    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)       

             call eos(eos_input_rt, eos_state, pt_index)       
       
             pres(i,j,k) = eos_state%p

             chi_rho = eos_state%rho * eos_state%dpdr / eos_state%p
             chi_t = eos_state%T * eos_state%dpdt / eos_state%p
             nabla_ad(i,j,k) = (eos_state%gam1 - chi_rho) / (chi_t * eos_state%gam1)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,dt,dp,nabla)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,rho_comp) <= base_cutoff_density) then
                nabla = ZERO
             else
                ! forward difference
                if (k == lo(3)) then
                   dt = state(i,j,k+1,temp_comp) - state(i,j,k,temp_comp)
                   dp = pres(i,j,k+1) - pres(i,j,k)
                   ! backward difference
                else if (k == hi(3)) then
                   dt = state(i,j,k,temp_comp) - state(i,j,k-1,temp_comp)
                   dp = pres(i,j,k) - pres(i,j,k-1)
                   ! centered difference
                else
                   dt = state(i,j,k+1,temp_comp) - state(i,j,k-1,temp_comp)
                   dp = pres(i,j,k+1) - pres(i,j,k-1)
                endif

                ! prevent Inf
                if (dp == ZERO) then
                   nabla = -huge(ZERO)
                else
                   nabla = pres(i,j,k) * dt / (dp * state(i,j,k,temp_comp))
                endif
             endif

             ad_excess(i,j,k) = nabla - nabla_ad(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_ad_excess_3d

  subroutine make_ad_excess_3d_sphr(ad_excess, ng_ad, state, ng_s, &
                                    normal, ng_n, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s, ng_n
    real(kind=dp_t), intent(  out) :: ad_excess(lo(1)-ng_ad:,lo(2)-ng_ad:,lo(3)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)  

    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_rho, chi_t, nabla
    real(kind=dp_t) :: dp(4), dt(4)

    integer :: i, j, k, c

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    !$OMP PARALLEL DO PRIVATE(i,j,k,chi_rho,chi_t,eos_state,pt_index)    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)       

             call eos(eos_input_rt, eos_state, pt_index)
       
             pres(i,j,k) = eos_state%p

             chi_rho = eos_state%rho * eos_state%dpdr / eos_state%p
             chi_t = eos_state%T * eos_state%dpdt / eos_state%p
             nabla_ad(i,j,k) = (eos_state%gam1 - chi_rho) / (chi_t * eos_state%gam1)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,dt,dp,nabla)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,rho_comp) <= base_cutoff_density) then
                nabla = ZERO
             else
                ! compute gradient

                ! forward difference
                if (k == lo(3)) then
                   dt(3) = state(i,j,k+1,temp_comp) - state(i,j,k,temp_comp)
                   dp(3) = pres(i,j,k+1) - pres(i,j,k)
                   ! backward difference
                else if (k == hi(3)) then
                   dt(3) = state(i,j,k,temp_comp) - state(i,j,k-1,temp_comp)
                   dp(3) = pres(i,j,k) - pres(i,j,k-1)
                   ! centered difference
                else
                   dt(3) = state(i,j,k+1,temp_comp) - state(i,j,k-1,temp_comp)
                   dp(3) = pres(i,j,k+1) - pres(i,j,k-1)
                endif

                if (j == lo(2)) then
                   dt(2) = state(i,j+1,k,temp_comp) - state(i,j,k,temp_comp)
                   dp(2) = pres(i,j+1,k) - pres(i,j,k)
                   ! backward difference
                else if (j == hi(2)) then
                   dt(2) = state(i,j,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp(2) = pres(i,j,k) - pres(i,j-1,k)
                   ! centered difference
                else
                   dt(2) = state(i,j+1,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp(2) = pres(i,j+1,k) - pres(i,j-1,k)
                endif

                if (i == lo(1)) then
                   dt(1) = state(i+1,j,k,temp_comp) - state(i,j,k,temp_comp)
                   dp(1) = pres(i+1,j,k) - pres(i,j,k)
                   ! backward difference
                else if (i == hi(1)) then
                   dt(1) = state(i,j,k,temp_comp) - state(i-1,j,k,temp_comp)
                   dp(1) = pres(i,j,k) - pres(i-1,j,k)
                   ! centered difference
                else
                   dt(1) = state(i+1,j,k,temp_comp) - state(i-1,j,k,temp_comp)
                   dp(1) = pres(i+1,j,k) - pres(i-1,j,k)
                endif

                ! dot into normal to get d/dr
                dp(4) = 0.d0
                dt(4) = 0.d0
                do c = 1,3
                   dp(4) = dp(4) + dp(c)*normal(i,j,k,c) 
                   dt(4) = dt(4) + dt(c)*normal(i,j,k,c) 
                enddo

                ! prevent Inf
                if (dp(4) == ZERO) then
                   nabla = -huge(ZERO)
                else
                   nabla = pres(i,j,k)*dt(4) / (dp(4)*state(i,j,k,temp_comp))
                endif
             endif

             ad_excess(i,j,k) = nabla - nabla_ad(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_ad_excess_3d_sphr


  !---------------------------------------------------------------------------
  ! make_conductivity
  !---------------------------------------------------------------------------
  subroutine make_conductivity(plotdata,comp_cond,state)

    type(multifab),  intent(inout) :: plotdata
    integer,         intent(in   ) :: comp_cond
    type(multifab),  intent(in   ) :: state

    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer :: lo(get_dim(plotdata)), hi(get_dim(plotdata)), ng_s, ng_c
    integer :: i,dm

    dm = get_dim(plotdata)

    ng_s = nghost(state)
    ng_c = nghost(plotdata)

    do i = 1, nfabs(state)
       sp => dataptr(state,i)
       cp => dataptr(plotdata,i)
       lo = lwb(get_box(state,i))
       hi = upb(get_box(state,i))
       select case (dm)
       case (1)
          call make_conductivity_1d(cp(:,1,1,comp_cond), ng_c, &
                                    sp(:,1,1,:), ng_s, lo, hi)
       case (2)
          call make_conductivity_2d(cp(:,:,1,comp_cond), ng_c, &
                                    sp(:,:,1,:), ng_s, lo, hi)
       case (3)
          call make_conductivity_3d(cp(:,:,:,comp_cond), ng_c, &
                                    sp(:,:,:,:), ng_s, lo, hi)
       end select

    enddo

  end subroutine make_conductivity

  subroutine make_conductivity_1d(cond,ng_c,state,ng_s,lo,hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use conductivity_module
    use network, only: nspec

    integer,         intent(in   ) :: lo(:), hi(:), ng_c, ng_s
    real(kind=dp_t), intent(  out) :: cond(lo(1)-ng_c:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,:)

    ! local
    integer :: i

    type (eos_t) :: eos_state
    real (kind=dp_t) :: conductivity

    do i = lo(1), hi(1)

       eos_state%rho   = state(i,rho_comp)
       eos_state%T     = state(i,temp_comp)
       eos_state%xn(:) = state(i,spec_comp:spec_comp+nspec-1) / eos_state%rho

       call conducteos(eos_input_rt, eos_state, .false., conductivity)

       cond(i) = conductivity

    enddo

  end subroutine make_conductivity_1d

  subroutine make_conductivity_2d(cond,ng_c,state,ng_s,lo,hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use conductivity_module
    use network, only: nspec

    integer,         intent(in   ) :: lo(:), hi(:), ng_c, ng_s
    real(kind=dp_t), intent(  out) :: cond(lo(1)-ng_c:,lo(2)-ng_c:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,:)

    ! local
    integer :: i, j

    type (eos_t) :: eos_state
    real (kind=dp_t) :: conductivity

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          eos_state%rho   = state(i,j,rho_comp)
          eos_state%T     = state(i,j,temp_comp)
          eos_state%xn(:) = state(i,j,spec_comp:spec_comp+nspec-1) / eos_state%rho

          call conducteos(eos_input_rt, eos_state, .false., conductivity)

          cond(i,j) = conductivity

       enddo
    enddo

  end subroutine make_conductivity_2d

  subroutine make_conductivity_3d(cond,ng_c,state,ng_s,lo,hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use conductivity_module
    use network, only: nspec

    integer,         intent(in   ) :: lo(:), hi(:), ng_c, ng_s
    real(kind=dp_t), intent(  out) :: cond(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    ! local
    integer :: i, j, k

    type (eos_t) :: eos_state
    real (kind=dp_t) :: conductivity

    !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state,conductivity)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1) / &
                           eos_state%rho

             call conducteos(eos_input_rt, eos_state, .false., conductivity)

             cond(i,j,k) = conductivity

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_conductivity_3d


  !---------------------------------------------------------------------------
  ! make_pidivu
  !---------------------------------------------------------------------------
  subroutine make_pidivu(plotdata, comp_pidivu, pi_cc, u, dx)

    use ml_layout_module
    use bc_module

    type(multifab) , intent(inout) :: plotdata
    integer        , intent(in   ) :: comp_pidivu
    type(multifab) , intent(in   ) :: pi_cc
    type(multifab) , intent(in   ) :: u
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t), pointer :: pdp(:,:,:,:)
    real(kind=dp_t), pointer :: pip(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)

    integer :: i,ng_pd,ng_pi,ng_u,dm
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata))

    dm = get_dim(plotdata)

    ng_pd = nghost(plotdata)
    ng_pi = nghost(pi_cc)
    ng_u  = nghost(u)

    do i=1,nfabs(plotdata)
       pdp => dataptr(plotdata, i)
       pip => dataptr(pi_cc, i)
       up  => dataptr(u, i)

       lo  =  lwb(get_box(plotdata, i))
       hi  =  upb(get_box(plotdata, i))
       
       select case (dm)
       case (1)
          call make_pidivu_1d(pdp(:,1,1,comp_pidivu), ng_pd, &
                              pip(:,1,1,1), ng_pi, &
                              up(:,1,1,:), ng_u, &
                              lo, hi, dx)

       case (2)
          call make_pidivu_2d(pdp(:,:,1,comp_pidivu), ng_pd, &
                              pip(:,:,1,1), ng_pi, &
                              up(:,:,1,:), ng_u, &
                              lo, hi, dx)

       case (3)
          call make_pidivu_3d(pdp(:,:,:,comp_pidivu), ng_pd, &
                              pip(:,:,:,1), ng_pi, &
                              up(:,:,:,:), ng_u, &
                              lo, hi, dx)

       end select
    end do

  end subroutine make_pidivu

  subroutine make_pidivu_1d(pidivu, ng_pd, pi_cc, ng_pi, u, ng_u, lo, hi, dx)

    integer         , intent(in   ) :: lo(:), hi(:), ng_pd, ng_pi, ng_u
    real (kind=dp_t), intent(inout) ::   pidivu(lo(1)-ng_pd:)
    real (kind=dp_t), intent(in   ) ::    pi_cc(lo(1)-ng_pi:)
    real (kind=dp_t), intent(in   ) ::        u(lo(1)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i

    do i=lo(1),hi(1)
       pidivu(i) = pi_cc(i)*( HALF*(u(i+1,1) - u(i-1,1))/dx(1) )
    end do

  end subroutine make_pidivu_1d


  subroutine make_pidivu_2d(pidivu, ng_pd, pi_cc, ng_pi, u, ng_u, lo, hi, dx)

    integer         , intent(in   ) :: lo(:), hi(:), ng_pd, ng_pi, ng_u
    real (kind=dp_t), intent(inout) ::   pidivu(lo(1)-ng_pd:,lo(2)-ng_pd:)
    real (kind=dp_t), intent(in   ) ::    pi_cc(lo(1)-ng_pi:,lo(2)-ng_pi:)
    real (kind=dp_t), intent(in   ) ::        u(lo(1)-ng_u: ,lo(2)-ng_u: ,:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i, j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          pidivu(i,j) = pi_cc(i,j)*( HALF*(u(i+1,j,1) - u(i-1,j,1))/dx(1) + &
                                     HALF*(u(i,j+1,2) - u(i,j-1,2))/dx(2) )
       enddo
    enddo

  end subroutine make_pidivu_2d


  subroutine make_pidivu_3d(pidivu, ng_pd, pi_cc, ng_pi, u, ng_u, lo, hi, dx)

    integer         , intent(in   ) :: lo(:), hi(:), ng_pd, ng_pi, ng_u
    real (kind=dp_t), intent(inout) ::   pidivu(lo(1)-ng_pd:,lo(2)-ng_pd:,lo(3)-ng_pd:)
    real (kind=dp_t), intent(in   ) ::    pi_cc(lo(1)-ng_pi:,lo(2)-ng_pi:,lo(3)-ng_pi:)
    real (kind=dp_t), intent(in   ) ::        u(lo(1)-ng_u: ,lo(2)-ng_u: ,lo(3)-ng_u: ,:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i, j, k

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             pidivu(i,j,k) = pi_cc(i,j,k)*( HALF*(u(i+1,j,k,1) - u(i-1,j,k,1))/dx(1) + & 
                                            HALF*(u(i,j+1,k,2) - u(i,j-1,k,2))/dx(2) + & 
                                            HALF*(u(i,j,k+1,3) - u(i,j-1,k,3))/dx(3) )
          enddo
       enddo
    enddo

  end subroutine make_pidivu_3d


  !---------------------------------------------------------------------------
  ! make_tfromH
  !---------------------------------------------------------------------------
  subroutine make_tfromH(plotdata,comp_t,comp_tpert,comp_dp,state,p0,tempbar,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: comp_t,comp_tpert,comp_dp
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: state
    real(kind=dp_t), intent(in   ) :: p0(0:),tempbar(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(3), hi(3)
    integer :: i,dm

    dm = get_dim(plotdata)

    lo(:) = 1; hi(:) = 1

    do i = 1, nfabs(state)
       sp => dataptr(state, i)
       tp => dataptr(plotdata, i)
       lo(1:dm) =  lwb(get_box(state, i))
       hi(1:dm) =  upb(get_box(state, i))

       if (spherical .eq. 1) then
          call make_tfromH_3d_sphr(tp, lbound(tp), ubound(tp), &
                                   comp_t, comp_tpert, comp_dp, &
                                   sp, lbound(sp), ubound(sp), &
                                   lo,hi,p0,tempbar,dx)
       else
          call make_tfromH_cart(tp, lbound(tp), ubound(tp), &
                                comp_t, comp_tpert, comp_dp, &
                                sp, lbound(sp), ubound(sp), &
                                dm,lo,hi,p0,tempbar)
       end if
    end do

  end subroutine make_tfromH

  subroutine make_tfromH_cart(pdata,dlo,dhi, &
                              it,itpert,ideltaP, &
                              state,slo,shi, &
                              dm,lo,hi,p0,tempbar)

    use variables, only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use eos_module, only: eos_input_rh, eos
    use eos_type_module
    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), dlo(4), dhi(4), slo(4), shi(4), dm
    integer, intent(in) :: it, itpert,ideltaP
    real (kind = dp_t), intent(  out) :: pdata(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),dlo(4):dhi(4))
    real (kind = dp_t), intent(in   ) :: state(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),slo(4):shi(4))
    real (kind = dp_t), intent(in   ) :: p0(0:),tempbar(0:)

    ! Local variables
    integer :: i, j, k, r

    type (eos_t) :: eos_state

    !$OMP PARALLEL DO PRIVATE(i,j,k,r,eos_state)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             select case (dm)
             case (1)
                r = i
             case (2)
                r = j
             case (3)
                r = k
             end select

             ! (rho, H) --> T, p
             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%p     = p0(r)
             eos_state%T     = state(i,j,k,temp_comp)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho
             eos_state%h     = state(i,j,k,rhoh_comp)/state(i,j,k,rho_comp)

             call eos(eos_input_rh, eos_state)

             if (it > 0) pdata(i,j,k,iT) = eos_state%T
             if (.not. use_tfromp .and. itpert > 0) pdata(i,j,k,itpert) = eos_state%T - tempbar(r)

             if (ideltaP > 0) pdata(i,j,k,ideltaP) = (eos_state%p - p0(r))/ p0(r)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO    

  end subroutine make_tfromH_cart

  subroutine make_tfromH_3d_sphr(pdata,dlo,dhi, &
                                 iT,itpert,ideltaP, &
                                 state,slo,shi, &
                                 lo,hi,p0,tempbar,dx)

    use variables, only: rho_comp, rhoh_comp, spec_comp, temp_comp
    use eos_module, only: eos_input_rh, eos
    use eos_type_module
    use network, only: nspec
    use fill_3d_module

    integer, intent(in) :: lo(:), hi(:), dlo(4), dhi(4), slo(4), shi(4)
    integer, intent(in) :: it, itpert,ideltaP
    real (kind = dp_t), intent(  out) :: pdata(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),dlo(4):dhi(4))
    real (kind = dp_t), intent(in   ) :: state(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),slo(4):shi(4))
    real (kind = dp_t), intent(in   ) :: p0(0:),tempbar(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer          :: i, j, k
    real (kind=dp_t), allocatable :: p0_cart(:,:,:,:)
    real (kind=dp_t), allocatable :: tempbar_cart(:,:,:,:)

    type (eos_t) :: eos_state

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(  tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,p0,p0_cart,lo,hi,dx,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,tempbar,tempbar_cart,lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%h     = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             ! (rho, H) --> T, p
             call eos(eos_input_rh, eos_state)

             if (iT > 0) pdata(i,j,k,iT) = eos_state%T
             if (.not. use_tfromp .and. itpert > 0) pdata(i,j,k,itpert) = eos_state%T - tempbar_cart(i,j,k,1)
             
             if(ideltaP > 0) pdata(i,j,k,ideltaP) = (eos_state%p - p0_cart(i,j,k,1))/ p0_cart(i,j,k,1)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(p0_cart,tempbar_cart)

  end subroutine make_tfromH_3d_sphr


  !---------------------------------------------------------------------------
  ! make_tfromp
  !---------------------------------------------------------------------------
  subroutine make_tfromp(p, plotdata, s, tempbar, gamma1bar, p0, dx)

    use geometry, only: spherical

    type(plot_t)   , intent(in   ) :: p
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: tempbar(0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), pointer:: sp(:,:,:,:),tp(:,:,:,:)
    integer :: lo(3), hi(3), i
    integer :: dm

    dm = get_dim(plotdata)

    lo(:) = 1; hi(:) = 1

    do i = 1, nfabs(s)
       tp => dataptr(plotdata, i)
       sp => dataptr(s, i)
       lo(1:dm) =  lwb(get_box(s, i))
       hi(1:dm) =  upb(get_box(s, i))
  
       if (spherical .eq. 1) then
          call make_tfromp_3d_sphr(tp, lbound(tp), ubound(tp), &
                                   sp, lbound(sp), ubound(sp), &
                                   p, lo, hi, tempbar, gamma1bar, p0, dx)
       else
          call make_tfromp_cart(tp, lbound(tp), ubound(tp), &
                                sp, lbound(sp), ubound(sp), &
                                p, dm, lo, hi, tempbar, gamma1bar, p0)
       endif
    end do

  end subroutine make_tfromp

  subroutine make_tfromp_cart(pdata, dlo, dhi, &
                              s, slo, shi, &
                              p, dm,lo,hi,tempbar,gamma1bar,p0)
    
    use variables, only: rho_comp, spec_comp, temp_comp, pi_comp
    use eos_module, only: eos_input_rp, eos
    use eos_type_module
    use network, only: nspec
    use probin_module, only: use_pprime_in_tfromp

    integer, intent(in) :: lo(:), hi(:), dlo(4), dhi(4), slo(4), shi(4), dm
    type (plot_t), intent(in) :: p
    real (kind=dp_t), intent(  out) :: pdata(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),dlo(4):dhi(4))
    real (kind=dp_t), intent(in   ) :: s(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),slo(4):shi(4))
    real (kind=dp_t), intent(in   ) :: tempbar(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    ! Local variables
    integer          :: i, j, k, r

    type (eos_t) :: eos_state

    !$OMP PARALLEL DO PRIVATE(i,j,k,r,eos_state)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             select case (dm)
             case (1)
                r = i
             case (2)
                r = j
             case (3)
                r = k
             end select

             eos_state%rho   = s(i,j,k,rho_comp)
             eos_state%T     = s(i,j,k,temp_comp)
             if (use_pprime_in_tfromp) then
                eos_state%p     = p0(r) + s(i,j,k,pi_comp)
             else
                eos_state%p     = p0(r)
             endif
             eos_state%xn(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             ! (rho,P) --> T,h
             call eos(eos_input_rp, eos_state)

             if (p%icomp_tfromp > 0) pdata(i,j,k,p%icomp_tfromp) = eos_state%T
             if (use_tfromp .and. p%icomp_tpert > 0) pdata(i,j,k,p%icomp_tpert) = eos_state%T - tempbar(r)

             if (p%icomp_cs > 0) pdata(i,j,k,p%icomp_cs) = eos_state%cs

             if (p%icomp_machno > 0 .and. p%icomp_magvel > 0) then
                pdata(i,j,k,p%icomp_machno) = pdata(i,j,k,p%icomp_magvel) / eos_state%cs
             endif

             if (p%icomp_dg > 0) pdata(i,j,k,p%icomp_dg) = eos_state%gam1 - gamma1bar(r)

             if (p%icomp_entropy > 0) pdata(i,j,k,p%icomp_entropy) = eos_state%s
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_tfromp_cart

  subroutine make_tfromp_3d_sphr(pdata,dlo,dhi, &
                                 s,slo,shi, &
                                 p, lo,hi,tempbar,gamma1bar,p0,dx)

    use variables, only: rho_comp, spec_comp, temp_comp, pi_comp
    use eos_module, only: eos_input_rp, eos
    use eos_type_module
    use network, only: nspec
    use fill_3d_module
    use probin_module, only: use_pprime_in_tfromp

    integer         , intent(in   ) :: lo(:),hi(:), dlo(4), dhi(4), slo(4), shi(4)
    type (plot_t), intent(in) :: p

    real (kind=dp_t), intent(  out) :: pdata(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),dlo(4):dhi(4))
    real (kind=dp_t), intent(in   ) :: s(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),slo(4):shi(4))
    real (kind=dp_t), intent(in   ) :: tempbar(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k

    type (eos_t) :: eos_state

    real (kind=dp_t), allocatable ::   tempbar_cart(:,:,:,:)
    real (kind=dp_t), allocatable ::        p0_cart(:,:,:,:)
    real (kind=dp_t), allocatable :: gamma1bar_cart(:,:,:,:)

    allocate(  tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(       p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(gamma1bar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,tempbar,tempbar_cart,lo,hi,dx,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,p0,p0_cart,lo,hi,dx,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,gamma1bar,gamma1bar_cart,lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! Then compute the perturbation and Mach number
             eos_state%rho   = s(i,j,k,rho_comp)
             eos_state%T     = s(i,j,k,temp_comp)
             if (use_pprime_in_tfromp) then
                eos_state%p     = p0_cart(i,j,k,1) + s(i,j,k,pi_comp)
             else
                eos_state%p     = p0_cart(i,j,k,1)
             endif
             eos_state%xn(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             ! (rho,P) --> T,h
             call eos(eos_input_rp, eos_state)

             if (p%icomp_tfromp > 0) pdata(i,j,k,p%icomp_tfromp) = eos_state%T
             if (use_tfromp .and. p%icomp_tpert > 0) pdata(i,j,k,p%icomp_tpert) = eos_state%T - tempbar_cart(i,j,k,1)

             if (p%icomp_cs > 0) pdata(i,j,k,p%icomp_cs) = eos_state%cs

             if (p%icomp_machno > 0 .and. p%icomp_magvel > 0) then
                pdata(i,j,k,p%icomp_machno) = pdata(i,j,k,p%icomp_magvel) / eos_state%cs
             endif

             if (p%icomp_dg > 0) pdata(i,j,k,p%icomp_dg) = eos_state%gam1 - gamma1bar_cart(i,j,k,1)

             if (p%icomp_entropy > 0) pdata(i,j,k,p%icomp_entropy) = eos_state%s
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(tempbar_cart,p0_cart,gamma1bar_cart)

  end subroutine make_tfromp_3d_sphr


  !---------------------------------------------------------------------------
  ! make_rhopert
  !---------------------------------------------------------------------------
  subroutine make_rhopert(plotdata,comp_rhopert,s,rho0,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: comp_rhopert
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), pointer:: sp(:,:,:,:),tp(:,:,:,:)
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),i
    integer :: ng_p,ng_s,dm

    dm = get_dim(plotdata)

    ng_p = nghost(plotdata)
    ng_s = nghost(s)

    do i = 1, nfabs(s)
       tp => dataptr(plotdata, i)
       sp => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (1)
          call make_rhopert_1d(tp(:,1,1,comp_rhopert ),ng_p, &
                               sp(:,1,1,:), ng_s, &
                               lo, hi, rho0)
       case (2)
          call make_rhopert_2d(tp(:,:,1,comp_rhopert ),ng_p, &
                               sp(:,:,1,:), ng_s, &
                               lo, hi, rho0)
       case (3)
          if (spherical .eq. 1) then
             call make_rhopert_3d_sphr(tp(:,:,:,comp_rhopert ),ng_p, &
                                       sp(:,:,:,:), ng_s, &
                                       lo, hi, rho0, dx)
          else
             call make_rhopert_3d_cart(tp(:,:,:,comp_rhopert ),ng_p, &
                                       sp(:,:,:,:), ng_s, &
                                       lo, hi, rho0)
          endif
       end select
    end do

  end subroutine make_rhopert

  subroutine make_rhopert_1d(rhopert,ng_p,s,ng_s,lo,hi,rho0)

    use variables, only: rho_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)

    !     Local variables
    integer          :: i

    ! Then compute the perturbation
    do i = lo(1), hi(1)
       rhopert(i)  = s(i,rho_comp)  - rho0(i)
    enddo

  end subroutine make_rhopert_1d

  subroutine make_rhopert_2d(rhopert,ng_p,s,ng_s,lo,hi,rho0)

    use variables, only: rho_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)

    !     Local variables
    integer          :: i, j

    ! Then compute the perturbation
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rhopert(i,j)  = s(i,j,rho_comp)  - rho0(j)
       enddo
    enddo

  end subroutine make_rhopert_2d

  subroutine make_rhopert_3d_cart(rhopert,ng_p,s,ng_s,lo,hi,rho0)

    use variables, only: rho_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)

    ! Local variables
    integer          :: i, j, k

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhopert(i,j,k)  = s(i,j,k,rho_comp)  - rho0(k)
          enddo
       enddo
    enddo

  end subroutine make_rhopert_3d_cart

  subroutine make_rhopert_3d_sphr(rhopert,ng_p,s,ng_s,lo,hi,rho0,dx)

    use variables, only: rho_comp
    use fill_3d_module

    integer         , intent(in   ) :: lo(:),hi(:),ng_p,ng_s
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k

    real (kind=dp_t), allocatable ::      rho0_cart(:,:,:,:)

    allocate(     rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,rho0,rho0_cart,lo,hi,dx,0)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhopert(i,j,k)  = s(i,j,k,rho_comp)  -  rho0_cart(i,j,k,1)
          enddo
       enddo
    enddo

    deallocate(rho0_cart)

  end subroutine make_rhopert_3d_sphr


  !---------------------------------------------------------------------------
  ! make_rhohpert
  !---------------------------------------------------------------------------
  subroutine make_rhohpert(plotdata,comp_rhohpert,s,rhoh0,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: comp_rhohpert
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: rhoh0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), pointer:: sp(:,:,:,:),tp(:,:,:,:)
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),i
    integer :: ng_p,ng_s,dm

    dm = get_dim(plotdata)

    ng_p = nghost(plotdata)
    ng_s = nghost(s)

    do i = 1, nfabs(s)
       tp => dataptr(plotdata, i)
       sp => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (1)
          call make_rhohpert_1d(tp(:,1,1,comp_rhohpert ),ng_p, &
                                sp(:,1,1,:), ng_s, &
                                lo, hi, rhoh0)
       case (2)
          call make_rhohpert_2d(tp(:,:,1,comp_rhohpert ),ng_p, &
                                sp(:,:,1,:), ng_s, &
                                lo, hi, rhoh0)
       case (3)
          if (spherical .eq. 1) then
             call make_rhohpert_3d_sphr(tp(:,:,:,comp_rhohpert ),ng_p, &
                                        sp(:,:,:,:), ng_s, &
                                        lo, hi, rhoh0, dx)
          else
             call make_rhohpert_3d_cart(tp(:,:,:,comp_rhohpert ),ng_p, &
                                        sp(:,:,:,:), ng_s, &
                                        lo, hi, rhoh0)
          endif
       end select
    end do

  end subroutine make_rhohpert

  subroutine make_rhohpert_1d(rhohpert,ng_p,s,ng_s,lo,hi,rhoh0)

    use variables, only: rhoh_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::   rhohpert(lo(1)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)

    !     Local variables
    integer          :: i

    ! Then compute the perturbation
    do i = lo(1), hi(1)
       rhohpert(i)  = s(i,rhoh_comp)  - rhoh0(i)
    enddo

  end subroutine make_rhohpert_1d

  subroutine make_rhohpert_2d(rhohpert,ng_p,s,ng_s,lo,hi,rhoh0)

    use variables, only: rhoh_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::   rhohpert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)

    !     Local variables
    integer          :: i, j

    ! Then compute the perturbation
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rhohpert(i,j)  = s(i,j,rhoh_comp)  - rhoh0(j)
       enddo
    enddo

  end subroutine make_rhohpert_2d

  subroutine make_rhohpert_3d_cart(rhohpert,ng_p,s,ng_s,lo,hi,rhoh0)

    use variables, only: rhoh_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::   rhohpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)

    ! Local variables
    integer          :: i, j, k

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhohpert(i,j,k)  = s(i,j,k,rhoh_comp)  - rhoh0(k)
          enddo
       enddo
    enddo

  end subroutine make_rhohpert_3d_cart

  subroutine make_rhohpert_3d_sphr(rhohpert,ng_p,s,ng_s,lo,hi,rhoh0,dx)

    use variables, only: rhoh_comp
    use fill_3d_module

    integer         , intent(in   ) :: lo(:),hi(:),ng_p,ng_s
    real (kind=dp_t), intent(  out) ::    rhohpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k

    real (kind=dp_t), allocatable ::      rhoh0_cart(:,:,:,:)

    allocate(     rhoh0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,rhoh0,rhoh0_cart,lo,hi,dx,0)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhohpert(i,j,k)  = s(i,j,k,rhoh_comp)  -  rhoh0_cart(i,j,k,1)
          enddo
       enddo
    enddo

    deallocate(rhoh0_cart)

  end subroutine make_rhohpert_3d_sphr


  !---------------------------------------------------------------------------
  ! make_entropypert
  !---------------------------------------------------------------------------
  subroutine make_entropypert(plotdata,comp_entropy,comp_entropypert,entropybar,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: comp_entropy,comp_entropypert
    type(multifab) , intent(inout) :: plotdata(:)
    real(kind=dp_t), intent(in   ) :: entropybar(:,0:)
    real(kind=dp_t), intent(in   ) ::         dx(:,:)

    ! local
    real(kind=dp_t), pointer :: tp(:,:,:,:)
    integer                  :: lo(get_dim(plotdata(1))),hi(get_dim(plotdata(1)))
    integer                  :: ng_p,n,i,dm,nlevs

    dm = get_dim(plotdata(1))
    nlevs = size(plotdata)

    ng_p = nghost(plotdata(1))

    do n=1,nlevs
       do i = 1, nfabs(plotdata(n))
          tp => dataptr(plotdata(n), i)
          lo =  lwb(get_box(plotdata(n), i))
          hi =  upb(get_box(plotdata(n), i))
          select case (dm)
          case (1)
             call make_entropypert_1d(tp(:,1,1,comp_entropy), &
                                      tp(:,1,1,comp_entropypert),ng_p, &
                                      lo, hi, entropybar(n,:))
          case (2)
             call make_entropypert_2d(tp(:,:,1,comp_entropy), &
                                      tp(:,:,1,comp_entropypert),ng_p, &
                                      lo, hi, entropybar(n,:))
          case (3)
             if (spherical .eq. 1) then
                call make_entropypert_3d_sphr(tp(:,:,:,comp_entropy), &
                                              tp(:,:,:,comp_entropypert),ng_p, &
                                              lo, hi, entropybar(1,:), dx(n,:))
             else
                call make_entropypert_3d_cart(tp(:,:,:,comp_entropy), &
                                              tp(:,:,:,comp_entropypert),ng_p, &
                                              lo, hi, entropybar(n,:))
             endif
          end select
       end do
    end do

  end subroutine make_entropypert

  subroutine make_entropypert_1d(entropy,entropypert,ng_p,lo,hi,entropybar)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind=dp_t), intent(inout) ::     entropy(lo(1)-ng_p:)
    real (kind=dp_t), intent(  out) :: entropypert(lo(1)-ng_p:)
    real (kind=dp_t), intent(in   ) :: entropybar(0:)

    !     Local variables
    integer          :: i

    ! Compute entropy - entropybar
    do i = lo(1), hi(1)
       entropypert(i) = (entropy(i) - entropybar(i))/entropybar(i)
    enddo

  end subroutine make_entropypert_1d

  subroutine make_entropypert_2d(entropy,entropypert,ng_p,lo,hi,entropybar)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind=dp_t), intent(inout) ::     entropy(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) :: entropypert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(in   ) :: entropybar(0:)

    !     Local variables
    integer          :: i, j

    ! Compute entropy - entropybar
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          entropypert(i,j) = (entropy(i,j) - entropybar(j))/entropybar(j)
       enddo
    enddo

  end subroutine make_entropypert_2d

  subroutine make_entropypert_3d_cart(entropy,entropypert,ng_p,lo,hi,entropybar)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind=dp_t), intent(inout) ::     entropy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) :: entropypert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) :: entropybar(0:)

    ! Local variables
    integer          :: i, j, k


    ! Compute entropy - entropybar

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             entropypert(i,j,k) = (entropy(i,j,k) - entropybar(k))/entropybar(k)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_entropypert_3d_cart

  subroutine make_entropypert_3d_sphr(entropy,entropypert,ng_p,lo,hi,entropybar,dx)

    use fill_3d_module

    integer         , intent(in)    :: lo(:),hi(:),ng_p
    real (kind=dp_t), intent(inout) ::     entropy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) :: entropypert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) :: entropybar(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k
    real (kind=dp_t), allocatable :: entropybar_cart(:,:,:,:)

    allocate(entropybar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,entropybar,entropybar_cart, &
                                      lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ! Compute entropy-entropybar
             entropypert(i,j,k) = &
                  (entropy(i,j,k) - entropybar_cart(i,j,k,1))/entropybar_cart(i,j,k,1)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(entropybar_cart)

  end subroutine make_entropypert_3d_sphr


  !---------------------------------------------------------------------------
  ! make_deltaT
  !---------------------------------------------------------------------------
  subroutine make_deltaT(plotdata,comp_dT,comp_tfromH,comp_tfromp)

    integer        , intent(in   ) :: comp_dT, comp_tfromH, comp_tfromp
    type(multifab) , intent(inout) :: plotdata

    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata))
    integer :: i,ng_p,dm

    dm = get_dim(plotdata)

    ng_p = nghost(plotdata)

    do i = 1, nfabs(plotdata)
       tp => dataptr(plotdata, i)
       lo =  lwb(get_box(plotdata, i))
       hi =  upb(get_box(plotdata, i))
       select case (dm)
       case (1)
          call make_deltaT_1d(tp(:,1,1,comp_dT),tp(:,1,1,comp_tfromH), &
                              tp(:,1,1,comp_tfromp), ng_p, lo, hi)
       case (2)
          call make_deltaT_2d(tp(:,:,1,comp_dT),tp(:,:,1,comp_tfromH), &
                              tp(:,:,1,comp_tfromp), ng_p, lo, hi)
       case (3)
          call make_deltaT_3d(tp(:,:,:,comp_dT),tp(:,:,:,comp_tfromH), &
                              tp(:,:,:,comp_tfromp), ng_p, lo, hi)
       end select
    end do

  end subroutine make_deltaT

  subroutine make_deltaT_1d(dT,tfromH,tfromp,ng_p,lo,hi)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind = dp_t), intent(  out) ::     dT(lo(1)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromH(lo(1)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromp(lo(1)-ng_p:)

    !     Local variables
    integer :: i

    do i = lo(1), hi(1)
       dT(i) = (tfromH(i) - tfromp(i)) / tfromH(i)
    end do

  end subroutine make_deltaT_1d

  subroutine make_deltaT_2d(dT,tfromH,tfromp,ng_p,lo,hi)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind = dp_t), intent(  out) ::     dT(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromH(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromp(lo(1)-ng_p:,lo(2)-ng_p:)

    !     Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          dT(i,j) = (tfromH(i,j) - tfromp(i,j)) / tfromH(i,j)
       end do
    end do

  end subroutine make_deltaT_2d

  subroutine make_deltaT_3d(dT,tfromH,tfromp,ng_p,lo,hi)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind = dp_t), intent(  out) ::     dT(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromH(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromp(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)

    !     Local variables
    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dT(i,j,k) = (tfromH(i,j,k) - tfromp(i,j,k)) / tfromH(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine make_deltaT_3d


  !---------------------------------------------------------------------------
  ! make_divw0
  !---------------------------------------------------------------------------
  subroutine make_divw0(divw0,comp_divw0,w0,w0mac,dx)

    use geometry, only: spherical

    type(multifab) , intent(inout) :: divw0
    integer        , intent(in   ) :: comp_divw0
    real(kind=dp_t), intent(in   ) :: w0(0:)
    type(multifab) , intent(in   ) :: w0mac(:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    real(kind=dp_t), pointer :: dwp(:,:,:,:)
    integer                  :: lo(get_dim(divw0)),hi(get_dim(divw0))
    integer                  :: i,ng_w0,ng_dw,dm

    dm = get_dim(divw0)

    ng_dw = nghost(divw0)

    do i=1,nfabs(divw0)
       dwp => dataptr(divw0, i)
       lo = lwb(get_box(divw0, i))
       hi = upb(get_box(divw0, i))
       select case (dm)
       case (1)
          call make_divw0_1d(w0, dwp(:,1,1,comp_divw0), ng_dw, lo, hi, dx)
       case (2)
          call make_divw0_2d(w0, dwp(:,:,1,comp_divw0), ng_dw, lo, hi, dx)
       case (3)
          if(spherical .eq. 1) then
             ng_w0 = nghost(w0mac(1))
             w0xp => dataptr(w0mac(1), i)
             w0yp => dataptr(w0mac(2), i)
             w0zp => dataptr(w0mac(3), i)
             call make_divw0_3d_sphr(w0xp(:,:,:,1), w0yp(:,:,:,1), w0zp(:,:,:,1), ng_w0, &
                                     dwp(:,:,:,comp_divw0), ng_dw, lo, hi, dx)
          else
             call make_divw0_3d(w0, dwp(:,:,:,comp_divw0), ng_dw, lo, hi, dx)
          end if
       end select
    end do

  end subroutine make_divw0

  subroutine make_divw0_1d(w0,divw0,ng_dw,lo,hi,dx)

    integer, intent(in)            :: lo(:), hi(:), ng_dw
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(inout) :: divw0(lo(1)-ng_dw:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer :: i

    do i = lo(1), hi(1)
       divw0(i) = ( w0(i+1) - w0(i) )/dx(1)
    end do

  end subroutine make_divw0_1d

  subroutine make_divw0_2d(w0,divw0,ng_dw,lo,hi,dx)

    integer, intent(in)            :: lo(:), hi(:), ng_dw
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(inout) :: divw0(lo(1)-ng_dw:,lo(2)-ng_dw:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer :: j

    do j = lo(2), hi(2)
       divw0(:,j) = ( w0(j+1) - w0(j) )/dx(2)
    end do

  end subroutine make_divw0_2d

  subroutine make_divw0_3d(w0,divw0,ng_dw,lo,hi,dx)

    integer, intent(in)            :: lo(:), hi(:), ng_dw
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(inout) :: divw0(lo(1)-ng_dw:,lo(2)-ng_dw:,lo(3)-ng_dw:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer :: k

    do k = lo(3), hi(3)
       divw0(:,:,k) = ( w0(k+1) - w0(k) )/dx(3)
    end do

  end subroutine make_divw0_3d

  subroutine make_divw0_3d_sphr(w0macx,w0macy,w0macz,ng_w0,divw0,ng_dw,lo,hi,dx)

    integer, intent(in)               :: lo(:), hi(:), ng_w0, ng_dw
    real (kind = dp_t), intent(inout) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real (kind = dp_t), intent(inout) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real (kind = dp_t), intent(inout) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real (kind = dp_t), intent(inout) ::  divw0(lo(1)-ng_dw:,lo(2)-ng_dw:,lo(3)-ng_dw:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             divw0(i,j,k) = (w0macx(i+1,j,k)-w0macx(i,j,k)) / dx(1) &
                           +(w0macy(i,j+1,k)-w0macy(i,j,k)) / dx(2) &
                           +(w0macz(i,j,k+1)-w0macz(i,j,k)) / dx(3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine make_divw0_3d_sphr


  !---------------------------------------------------------------------------
  ! make_vorticity
  !---------------------------------------------------------------------------
  subroutine make_vorticity(vort,comp,u,dx,bc)

    use bl_prof_module

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: vort
    type(multifab) , intent(in   ) :: u
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    real(kind=dp_t), pointer:: up(:,:,:,:)
    real(kind=dp_t), pointer:: vp(:,:,:,:)
    integer :: lo(get_dim(vort)),hi(get_dim(vort))
    integer :: i,ng_u,ng_v,dm

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_vort")

    dm = get_dim(vort)

    ng_u = nghost(u)
    ng_v = nghost(vort)

    do i = 1, nfabs(u)
       up  => dataptr(u, i)
       vp  => dataptr(vort, i)
       lo  =  lwb(get_box(u, i))
       hi  =  upb(get_box(u, i))
       select case (dm)
       case (1)
          call setval(vort,0.d0,comp,1,all=.true.)
       case (2)
          call makevort_2d(vp(:,:,1,comp),ng_v,up(:,:,1,:),ng_u,lo,hi,dx, &
                           bc%phys_bc_level_array(i,:,:))
       case (3)
          call makevort_3d(vp(:,:,:,comp),ng_v,up(:,:,:,:),ng_u,lo,hi,dx, &
                           bc%phys_bc_level_array(i,:,:))
       end select
    end do

    call destroy(bpt)

  end subroutine make_vorticity

  subroutine makevort_2d(vort,ng_v,u,ng_u,lo,hi,dx,bc)

    use bc_module
    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_v, ng_u
    real (kind = dp_t), intent(  out) :: vort(lo(1)-ng_v:,lo(2)-ng_v:)  
    real (kind = dp_t), intent(in   ) ::    u(lo(1)-ng_u:,lo(2)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    integer           , intent(in   ) :: bc(:,:)

    !     Local variables
    integer :: i, j
    real (kind = dp_t) :: vx,uy

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
          vort(i,j) = vx - uy
       enddo
    enddo

    if (bc(1,1) .eq. INLET .or. bc(1,1) .eq. SLIP_WALL .or. bc(1,1) .eq. NO_SLIP_WALL) then
       i = lo(1)
       do j = lo(2), hi(2)
          vx = (u(i+1,j,2) + 3.d0*u(i,j,2) - 4.d0*u(i-1,j,2)) / dx(1)
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
          vort(i,j) = vx - uy
       end do
    end if

    if (bc(1,2) .eq. INLET .or. bc(1,2) .eq. SLIP_WALL .or. bc(1,2) .eq. NO_SLIP_WALL) then
       i = hi(1)
       do j = lo(2), hi(2)
          vx = -(u(i-1,j,2) + 3.d0*u(i,j,2) - 4.d0*u(i+1,j,2)) / dx(1)
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
          vort(i,j) = vx - uy
       end do
    end if

    if (bc(2,1) .eq. INLET .or. bc(2,1) .eq. SLIP_WALL .or. bc(2,1) .eq. NO_SLIP_WALL) then
       j = lo(2)
       do i = lo(1), hi(1)
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
          uy = (u(i,j+1,1) + 3.d0*u(i,j,1) - 4.d0*u(i,j-1,1)) / dx(2)
          vort(i,j) = vx - uy
       end do
    end if

    if (bc(2,2) .eq. INLET .or. bc(2,2) .eq. SLIP_WALL .or. bc(2,2) .eq. NO_SLIP_WALL) then
       j = hi(2)
       do i = lo(1), hi(1)
          vx =  (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
          uy = -(u(i,j-1,1) + 3.d0*u(i,j,1) - 4.d0*u(i,j+1,1)) / dx(2)
          vort(i,j) = vx - uy
       end do
    end if

  end subroutine makevort_2d

  subroutine makevort_3d(vort,ng_v,u,ng_u,lo,hi,dx,bc)

    use bc_module
    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_v, ng_u
    real (kind = dp_t), intent(  out) :: vort(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    real (kind = dp_t), intent(in   ) ::    u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    integer           , intent(in   ) :: bc(:,:)

    !     Local variables
    integer :: i, j, k
    logical :: fix_lo_x,fix_hi_x,fix_lo_y,fix_hi_y,fix_lo_z,fix_hi_z
    real (kind = dp_t) :: wy,vz,uz,wx,vx,uy

    !$OMP PARALLEL DO PRIVATE(i,j,k,uy,uz,vx,vz,wx,wy)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uy = uycen(i,j,k)
             uz = uzcen(i,j,k)
             vx = vxcen(i,j,k)
             vz = vzcen(i,j,k)
             wx = wxcen(i,j,k)
             wy = wycen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    fix_lo_x = ( bc(1,1) .eq. INLET .or. bc(1,1) .eq. NO_SLIP_WALL )
    fix_hi_x = ( bc(1,2) .eq. INLET .or. bc(1,2) .eq. NO_SLIP_WALL )

    fix_lo_y = ( bc(2,1) .eq. INLET .or. bc(2,1) .eq. NO_SLIP_WALL )
    fix_hi_y = ( bc(2,2) .eq. INLET .or. bc(2,2) .eq. NO_SLIP_WALL )

    fix_lo_z = ( bc(3,1) .eq. INLET .or. bc(3,1) .eq. NO_SLIP_WALL )
    fix_hi_z = ( bc(3,2) .eq. INLET .or. bc(3,2) .eq. NO_SLIP_WALL )

    !
    !     First do all the faces
    !
    if (fix_lo_x) then
       i = lo(1)
       !$OMP PARALLEL DO PRIVATE(j,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             vx = vxlo(i,j,k)
             wx = wxlo(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_x) then
       i = hi(1)
       !$OMP PARALLEL DO PRIVATE(j,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             vx = vxhi(i,j,k)
             wx = wxhi(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_y) then
       j = lo(2)
       !$OMP PARALLEL DO PRIVATE(i,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(j)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uylo(i,j,k)
             wy = wylo(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_y) then
       j = hi(2)
       !$OMP PARALLEL DO PRIVATE(i,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(j)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uyhi(i,j,k)
             wy = wyhi(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_z) then
       k = lo(3)
       !$OMP PARALLEL DO PRIVATE(i,j,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(k)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzlo(i,j,k)
             vz = vzlo(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_z) then
       k = hi(3)
       !$OMP PARALLEL DO PRIVATE(i,j,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(k)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzhi(i,j,k)
             vz = vzhi(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if
    !
    !     Next do all the edges
    !
    if (fix_lo_x .and. fix_lo_y) then
       i = lo(1)
       j = lo(2)
       do k = lo(3),hi(3)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_lo_y) then
       i = hi(1)
       j = lo(2)
       do k = lo(3),hi(3)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_x .and. fix_hi_y) then
       i = lo(1)
       j = hi(2)
       do k = lo(3),hi(3)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_hi_y) then
       i = hi(1)
       j = hi(2)
       do k = lo(3),hi(3)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_x .and. fix_lo_z) then
       i = lo(1)
       k = lo(3)
       do j = lo(2),hi(2)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_lo_z) then
       i = hi(1)
       k = lo(3)
       do j = lo(2),hi(2)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_x .and. fix_hi_z) then
       i = lo(1)
       k = hi(3)
       do j = lo(2),hi(2)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_hi_z) then
       i = hi(1)
       k = hi(3)
       do j = lo(2),hi(2)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_y .and. fix_lo_z) then
       j = lo(2)
       k = lo(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_y .and. fix_lo_z) then
       j = hi(2)
       k = lo(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_y .and. fix_hi_z) then
       j = lo(2)
       k = hi(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_y .and. fix_hi_z) then
       j = hi(2)
       k = hi(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if
    !
    !     Finally do all the corners
    !
    if (fix_lo_x .and. fix_lo_y .and. fix_lo_z) then
       i = lo(1)
       j = lo(2)
       k = lo(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_lo_y .and. fix_lo_z) then
       i = hi(1)
       j = lo(2)
       k = lo(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_hi_y .and. fix_lo_z) then
       i = lo(1)
       j = hi(2)
       k = lo(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_hi_y .and. fix_lo_z) then
       i = hi(1)
       j = hi(2)
       k = lo(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_lo_y .and. fix_hi_z) then
       i = lo(1)
       j = lo(2)
       k = hi(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_lo_y .and. fix_hi_z) then
       i = hi(1)
       j = lo(2)
       k = hi(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_hi_y .and. fix_hi_z) then
       i = lo(1)
       j = hi(2)
       k = hi(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_hi_y .and. fix_hi_z) then
       i = hi(1)
       j = hi(2)
       k = hi(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

  contains

    function uycen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j+1,k,1)-u(i,j-1,k,1))/dx(2)
    end function uycen

    function uylo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j+1,k,1)+THREE*u(i,j,k,1)-FOUR*u(i,j-1,k,1))/(THREE*dx(2))
    end function uylo

    function uyhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = -(u(i,j-1,k,1)+THREE*u(i,j,k,1)-FOUR*u(i,j+1,k,1))/(THREE*dx(2))
    end function uyhi

    function uzcen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j,k+1,1)-u(i,j,k-1,1))/dx(3)
    end function uzcen

    function uzlo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j,k+1,1)+THREE*u(i,j,k,1)-FOUR*u(i,j,k-1,1))/(THREE*dx(3))
    end function uzlo

    function uzhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i,j,k-1,1)+THREE*u(i,j,k,1)-FOUR*u(i,j,k+1,1))/(THREE*dx(3))
    end function uzhi

    function vxcen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i+1,j,k,2)-u(i-1,j,k,2))/dx(1)
    end function vxcen

    function vxlo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i+1,j,k,2)+THREE*u(i,j,k,2)-FOUR*u(i-1,j,k,2))/(THREE*dx(1))
    end function vxlo

    function vxhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i-1,j,k,2)+THREE*u(i,j,k,2)-FOUR*u(i+1,j,k,2))/(THREE*dx(1))
    end function vxhi

    function vzcen(i,j,k) result(r) 
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j,k+1,2)-u(i,j,k-1,2))/dx(3)
    end function vzcen

    function vzlo(i,j,k) result(r) 
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j,k+1,2)+THREE*u(i,j,k,2)-FOUR*u(i,j,k-1,2))/(THREE*dx(3))
    end function vzlo

    function vzhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i,j,k-1,2)+THREE*u(i,j,k,2)-FOUR*u(i,j,k+1,2))/(THREE*dx(3))
    end function vzhi

    function wxcen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i+1,j,k,3)-u(i-1,j,k,3))/dx(1)
    end function wxcen

    function wxlo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i+1,j,k,3)+THREE*u(i,j,k,3)-FOUR*u(i-1,j,k,3))/(THREE*dx(1))
    end function wxlo

    function wxhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i-1,j,k,3)+THREE*u(i,j,k,3)-FOUR*u(i+1,j,k,3))/(THREE*dx(1))
    end function wxhi

    function wycen(i,j,k) result(r) 
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j+1,k,3)-u(i,j-1,k,3))/dx(2)
    end function wycen

    function wylo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j+1,k,3)+THREE*u(i,j,k,3)-FOUR*u(i,j-1,k,3))/(THREE*dx(2))
    end function wylo

    function wyhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i,j-1,k,3)+THREE*u(i,j,k,3)-FOUR*u(i,j+1,k,3))/(THREE*dx(2))
    end function wyhi

    function vorfun(uy,uz,vx,vz,wx,wy) result(r)
      real(dp_t) :: uy,uz,vx,vz,wx,wy
      real(dp_t) :: r
      r = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
    end function vorfun

  end subroutine makevort_3d


  !---------------------------------------------------------------------------
  ! make_magvel
  !---------------------------------------------------------------------------
  subroutine make_magvel(plotdata,comp_magvel,comp_mom,s,u,w0,w0mac)

    use bc_module
    use bl_constants_module
    use geometry, only : spherical
    use variables, only : rho_comp

    integer        , intent(in   ) :: comp_magvel, comp_mom
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s, u
    real(kind=dp_t), intent(in   ) :: w0(0:)
    type(multifab) , intent(in   ) :: w0mac(:)

    real(kind=dp_t), pointer:: pp(:,:,:,:)
    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: up(:,:,:,:)
    real(kind=dp_t), pointer:: wxp(:,:,:,:)
    real(kind=dp_t), pointer:: wyp(:,:,:,:)
    real(kind=dp_t), pointer:: wzp(:,:,:,:)

    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),ng_p,ng_s,ng_u,ng_w
    integer :: i,dm

    dm = get_dim(plotdata)

    ng_s = nghost(s)
    ng_u = nghost(u)
    ng_p = nghost(plotdata)

    do i = 1, nfabs(u)
       pp => dataptr(plotdata, i)
       sp => dataptr(s, i)
       up => dataptr(u, i)
       lo =  lwb(get_box(u, i))
       hi =  upb(get_box(u, i))
       select case (dm)
       case (1)
          call makemagvel_1d(pp(:,1,1,:), comp_magvel, comp_mom, ng_p, &
                             sp(:,1,1,rho_comp),ng_s,up(:,1,1,1),ng_u,w0,lo,hi)
       case (2)
          call makemagvel_2d(pp(:,:,1,:), comp_magvel, comp_mom, ng_p, &
                             sp(:,:,1,rho_comp),ng_s,up(:,:,1,:),ng_u,w0,lo,hi)
       case (3)
          if (spherical .eq. 1) then
             wxp => dataptr(w0mac(1), i)
             wyp => dataptr(w0mac(2), i)
             wzp => dataptr(w0mac(3), i)
             ng_w = nghost(w0mac(1))
             call makemagvel_3d_sphr(pp(:,:,:,:), comp_magvel, comp_mom, ng_p, &
                                     sp(:,:,:,rho_comp),ng_s,up(:,:,:,:),ng_u, &
                                     wxp(:,:,:,1),wyp(:,:,:,1),wzp(:,:,:,1),ng_w,lo,hi)
          else
             call makemagvel_3d_cart(pp(:,:,:,:), comp_magvel, comp_mom, ng_p, &
                                     sp(:,:,:,rho_comp),ng_s,up(:,:,:,:),ng_u,w0,lo,hi)
          end if
       end select
    end do

  end subroutine make_magvel

  subroutine makemagvel_1d(pdata, imagvel, imom,ng_p,rho,ng_s,u,ng_u,w0,lo,hi)

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_s
    integer           , intent(in   ) :: imagvel, imom
    real (kind = dp_t), intent(  out) :: pdata(lo(1)-ng_p:,:)
    real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng_s:)
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:)
    real (kind = dp_t), intent(in   ) :: w0(0:)

    !     Local variables
    integer :: i
    real (kind = dp_t) :: w0_cent

    ! Recall w0 is edge-centered
    do i = lo(1), hi(1)
       w0_cent = 0.5d0 * (w0(i) + w0(i+1))
       pdata(i,imagvel) = abs(u(i)+w0_cent)
       pdata(i,imom) = rho(i)*pdata(i,imagvel)
    enddo

  end subroutine makemagvel_1d

  subroutine makemagvel_2d(pdata,imagvel,imom,ng_p,rho,ng_s,u,ng_u,w0,lo,hi)

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_s
    integer           , intent(in   ) :: imagvel, imom
    real (kind = dp_t), intent(  out) :: pdata(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng_s:,lo(2)-ng_s:)  
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: w0(0:)

    !     Local variables
    integer :: i, j
    real (kind = dp_t) :: w0_cent

    ! Recall w0 is edge-centered
    do j = lo(2), hi(2)
       w0_cent = 0.5d0 * (w0(j) + w0(j+1))
       do i = lo(1), hi(1)
          pdata(i,j,imagvel) = sqrt( u(i,j,1)**2 + (u(i,j,2)+w0_cent)**2 )
          pdata(i,j,imom) = rho(i,j)*pdata(i,j,imagvel)
       enddo
    enddo

  end subroutine makemagvel_2d

  subroutine makemagvel_3d_cart(pdata,imagvel,imom,ng_p,rho,ng_s,u,ng_u,w0,lo,hi)

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_s
    integer           , intent(in   ) :: imagvel, imom
    real (kind = dp_t), intent(  out) :: pdata(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:) 
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:) 
    real (kind = dp_t), intent(in   ) :: w0(0:)

    !     Local variables
    integer :: i, j, k
    real (kind = dp_t) :: w0_cent

    ! Recall w0 is edge-centered

    !$OMP PARALLEL DO PRIVATE(i,j,k,w0_cent)
    do k = lo(3), hi(3)
       w0_cent = 0.5d0 * (w0(k) + w0(k+1))
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             pdata(i,j,k,imagvel) = sqrt(u(i,j,k,1)**2 + u(i,j,k,2)**2 + (u(i,j,k,3)+w0_cent)**2)
             pdata(i,j,k,imom) = rho(i,j,k)*pdata(i,j,k,imagvel)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine makemagvel_3d_cart

  subroutine makemagvel_3d_sphr(pdata,imagvel,imom,ng_p,rho,ng_s,u,ng_u, &
                                w0macx,w0macy,w0macz,ng_w,lo,hi)


    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_w, ng_s
    integer           , intent(in   ) :: imagvel, imom
    real (kind = dp_t), intent(  out) :: pdata(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:) 
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:) 
    real (kind = dp_t), intent(in   ) :: w0macx(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
    real (kind = dp_t), intent(in   ) :: w0macy(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
    real (kind = dp_t), intent(in   ) :: w0macz(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)

    !     Local variables
    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             pdata(i,j,k,imagvel) = sqrt( (u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)))**2 + &
                                         (u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)))**2 + &
                                         (u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)))**2)
             pdata(i,j,k,imom) = rho(i,j,k)*pdata(i,j,k,imagvel)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine makemagvel_3d_sphr


  !---------------------------------------------------------------------------
  ! make_velrc
  !---------------------------------------------------------------------------
  subroutine make_velrc(plotdata,comp_velr,comp_velc,u,w0r_cart,normal)

    use bc_module
    use bl_constants_module
    use geometry, only: spherical

    integer        , intent(in   ) :: comp_velr, comp_velc
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: u
    type(multifab) , intent(in   ) :: w0r_cart
    type(multifab) , intent(in   ) :: normal

    ! local
    real(kind=dp_t), pointer:: pp(:,:,:,:)
    real(kind=dp_t), pointer:: up(:,:,:,:)
    real(kind=dp_t), pointer:: nop(:,:,:,:)
    real(kind=dp_t), pointer:: w0rp(:,:,:,:)
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),ng_p,ng_u,ng_n,ng_w
    integer :: i,dm

    dm = get_dim(plotdata)

    ng_u = nghost(u)
    ng_p = nghost(plotdata)
    ng_n = nghost(normal)
    ng_w = nghost(w0r_cart)

    if (spherical .ne. 1) then
       call bl_error("unable to create radial and circumferential velocity -- not spherical geometry")
    endif

    do i = 1, nfabs(u)

       pp => dataptr(plotdata, i)
       up => dataptr(u, i)
       nop => dataptr(normal, i)
       w0rp => dataptr(w0r_cart, i)
       lo =  lwb(get_box(u, i))
       hi =  upb(get_box(u, i))

       call makevelrc_3d_sphr(pp(:,:,:,:), comp_velr, comp_velc,&
                              ng_p,up(:,:,:,:),ng_u, &
                              w0rp(:,:,:,1),ng_w,nop(:,:,:,:),ng_n,lo,hi)
    end do

  end subroutine make_velrc

  subroutine makevelrc_3d_sphr(pdata,ivelr,ivelc,ng_p, &
                               u,ng_u,w0r,ng_w,normal,ng_n,lo,hi)

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_n, ng_w
    integer           , intent(in   ) :: ivelr, ivelc
    real (kind = dp_t), intent(  out) :: pdata(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind = dp_t), intent(in   ) ::    w0r(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
    real (kind = dp_t), intent(in   ) :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)  

    !     Local variables
    integer :: i, j, k

    ! note, we are assuming here that both ivelr and ivelc are defined
    ! -- for mini-plotfiles, you can't have one without the other

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! radial tilde velocity (no w0_r)
             pdata(i,j,k,ivelr) = u(i,j,k,1)*normal(i,j,k,1) + &
                                  u(i,j,k,2)*normal(i,j,k,2) + &
                                  u(i,j,k,3)*normal(i,j,k,3) 

             !                                 ~     ~   ~   ^
             ! circumferential tilde velocity, U_c = U - U . e_r
             pdata(i,j,k,ivelc) = &
                  (u(i,j,k,1)-pdata(i,j,k,ivelr)*normal(i,j,k,1)) * &
                  (u(i,j,k,1)-pdata(i,j,k,ivelr)*normal(i,j,k,1))
             pdata(i,j,k,ivelc) = pdata(i,j,k,ivelc) + &
                  (u(i,j,k,2)-pdata(i,j,k,ivelr)*normal(i,j,k,2)) * &
                  (u(i,j,k,2)-pdata(i,j,k,ivelr)*normal(i,j,k,2))
             pdata(i,j,k,ivelc) = pdata(i,j,k,ivelc) + &
                  (u(i,j,k,3)-pdata(i,j,k,ivelr)*normal(i,j,k,3)) * &
                  (u(i,j,k,3)-pdata(i,j,k,ivelr)*normal(i,j,k,3))
             pdata(i,j,k,ivelc) = sqrt(pdata(i,j,k,ivelc))

             ! now add the base state expansion velocity to the 
             ! radial to get the full U_r
             pdata(i,j,k,ivelr) = pdata(i,j,k,ivelr) + w0r(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine makevelrc_3d_sphr


  !---------------------------------------------------------------------------
  ! make_processor_number
  !---------------------------------------------------------------------------
  subroutine make_processor_number(plotdata,comp_proc)

    type(multifab), intent(inout) :: plotdata
    integer,        intent(in   ) :: comp_proc
    
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    integer :: i

    do i = 1, nfabs(plotdata)
       pp => dataptr(plotdata, i)
       pp(:,:,:,comp_proc) = parallel_myproc()
    enddo
  end subroutine make_processor_number


end module plot_variables_module

