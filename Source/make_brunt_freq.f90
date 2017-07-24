module make_brunt_freq_module
  
  use bl_types
  use multifab_module
  use define_bc_module

  implicit none

  private

  public :: make_brunt_freq
  
  contains

  !---------------------------------------------------------------------------
  ! make_brunt_freq
  !---------------------------------------------------------------------------
  subroutine make_brunt_freq(brunt,state,rho0,p0,normal,dx)
    
    use make_grav_module
    use make_scale_module
    use geometry, only: spherical, nr_fine, nlevs_radial
    use probin_module, only: max_levs

    type(multifab), intent(inout) :: brunt(:)
    type(multifab), intent(in   ) :: state(:)
    type(multifab), intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in  ) :: rho0(nlevs_radial,0:nr_fine-1)
    real(kind=dp_t), intent(in  ) :: p0(nlevs_radial,0:nr_fine-1)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    
    real(kind=dp_t), pointer :: sp(:,:,:,:), cp(:,:,:,:), np(:,:,:,:)
    integer :: lo(get_dim(brunt(1))), hi(get_dim(brunt(1))), ng_s, ng_c, ng_n
    integer :: i, dm, n

    real(kind=dp_t) ::   hp(nlevs_radial,0:nr_fine-1)
    real(kind=dp_t) ::   grav(nlevs_radial,0:nr_fine-1)

    dm = get_dim(brunt(1))
    
    !grav(n,:) where n is the amr level
    call make_grav_cell(grav,rho0)
    call make_scale(hp,p0)
    
    
    do n=1,max_levs
	
	do i = 1, nfabs(state(n))
	  sp => dataptr(state(n), i)
	  cp => dataptr(brunt(n), i)
	  lo = lwb(get_box(state(n), i))
	  hi = upb(get_box(state(n), i))


	  ng_c = nghost(brunt(n))
	  ng_s = nghost(state(n))
	  ng_n = nghost(normal(n))
	  
	  select case (dm)
	  case (1)
	      call make_brunt_1d(cp(:,1,1,1), ng_c, &
				    sp(:,1,1,:), ng_s, &
				    grav(n,:), hp(n,:), &
				    lo, hi,dx(n,:))
          case (2)
	      call make_brunt_2d(cp(:,:,1,1), ng_c, &
				    sp(:,:,1,:), ng_s, &
				    grav(n,:), hp(n,:), &
				    lo, hi,dx(n,:))
	  case (3)
	      if (spherical .eq. 1) then
		np => dataptr(normal(n), i)       
		call make_brunt_3d_sphr(cp(:,:,:,1), ng_c, &
					    sp(:,:,:,:), ng_s, &
					    grav(1,:), hp(1,:), &
					    np(:,:,:,:), ng_n, lo, hi, dx(n,:))
	      else
		call make_brunt_3d(cp(:,:,:,1), ng_c, &
				    sp(:,:,:,:), ng_s, &
				    grav(n,:), hp(n,:), &
				    lo, hi,dx(n,:))
	      endif
	  end select

	enddo
    enddo
    
  end subroutine make_brunt_freq
  
  subroutine make_brunt_1d(brunt, ng_ad, state, ng_s, grav, hp, lo, hi,dx)
    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: brunt(lo(1)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: hp(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    real(kind=dp_t) :: pres(lo(1):hi(1))
    real(kind=dp_t) :: abar(lo(1):hi(1))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1))
    real(kind=dp_t) :: chi_rho(lo(1):hi(1))
    real(kind=dp_t) :: chi_t(lo(1):hi(1))
    real(kind=dp_t) :: chi_mu(lo(1):hi(1))
    real(kind=dp_t) :: nabla, nabla_mu
    real(kind=dp_t) :: dp, dt, dmu


    
    integer :: i, c

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    
    do i = lo(1), hi(1)

        eos_state%rho   = state(i,rho_comp)
        eos_state%T     = state(i,temp_comp)
        eos_state%xn(:) = state(i,spec_comp:spec_comp+nspec-1)/eos_state%rho
        
        pt_index(:) = (/i, -1, -1/) 
        
        call eos(eos_input_rt, eos_state, pt_index)

        pres(i) = eos_state%p
        abar(i) = eos_state%abar

        chi_rho(i) = eos_state%rho * eos_state%dpdr / eos_state%p
        chi_t(i) = eos_state%T * eos_state%dpdt / eos_state%p
        chi_mu(i) = eos_state%abar * eos_state%dpdA / eos_state%p
        nabla_ad(i) = (eos_state%gam1 - chi_rho(i)) / (chi_t(i) * eos_state%gam1)

    enddo

    
    do i = lo(1), hi(1)

        if (state(i,rho_comp) <= base_cutoff_density) then
        nabla = ZERO
        nabla_mu = ZERO
        else
        ! compute gradient

        ! forward difference
        if (i == lo(1)) then
            dt = state(i+1,temp_comp) - state(i,temp_comp)
            dp = pres(i+1) - pres(i)
            dmu = abar(i+1) - abar(i)
            ! backward difference
        else if (i == hi(1)) then
            dt = state(i,temp_comp) - state(i-1,temp_comp)
            dp = pres(i) - pres(i-1)
            dmu = abar(i) - abar(i-1)
            ! centered difference
        else
            dt = state(i+1,temp_comp) - state(i-1,temp_comp)
            dp = pres(i+1) - pres(i-1)
            dmu = abar(i+1) - abar(i-1)
        endif


        ! prevent Inf
        if (dp == ZERO) then
            nabla = -huge(ZERO)
            nabla_mu = -huge(ZERO)
        else
            nabla = pres(i)*dt / (dp*state(i,temp_comp))
            nabla_mu = pres(i)*dmu / (dp*abar(i))
        endif
        endif
        if (hp(i) == ZERO) then
        brunt(i) = -huge(ZERO)
        else 
        brunt(i) = - grav(i)*chi_t(i) / (chi_rho(i) * hp(i)) &
                    * (nabla_ad(i) - nabla &
                    - nabla_mu * chi_mu(i)/chi_t(i))
        endif
    enddo

  end subroutine make_brunt_1d

  subroutine make_brunt_2d(brunt, ng_ad, state, ng_s, grav, hp, lo, hi,dx)
    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: brunt(lo(1)-ng_ad:,lo(2)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: hp(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2))
    real(kind=dp_t) :: abar(lo(1):hi(1),lo(2):hi(2))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2))
    real(kind=dp_t) :: chi_rho(lo(1):hi(1),lo(2):hi(2))
    real(kind=dp_t) :: chi_t(lo(1):hi(1),lo(2):hi(2))
    real(kind=dp_t) :: chi_mu(lo(1):hi(1),lo(2):hi(2))
    real(kind=dp_t) :: nabla, nabla_mu
    real(kind=dp_t) :: dp, dt, dmu


    
    integer :: i, j, c

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    
    !$OMP PARALLEL DO PRIVATE(i,j,eos_state,pt_index)    
      do j = lo(2), hi(2)
	  do i = lo(1), hi(1)

	    eos_state%rho   = state(i,j,rho_comp)
	    eos_state%T     = state(i,j,temp_comp)
	    eos_state%xn(:) = state(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho
	    
	    pt_index(:) = (/i, j, -1/) 
	    
	    call eos(eos_input_rt, eos_state, pt_index)
      
	    pres(i,j) = eos_state%p
	    abar(i,j) = eos_state%abar

	    chi_rho(i,j) = eos_state%rho * eos_state%dpdr / eos_state%p
	    chi_t(i,j) = eos_state%T * eos_state%dpdt / eos_state%p
	    chi_mu(i,j) = eos_state%abar * eos_state%dpdA / eos_state%p
	    nabla_ad(i,j) = (eos_state%gam1 - chi_rho(i,j)) / (chi_t(i,j) * eos_state%gam1)

	  enddo
      enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,dt,dp,nabla,nabla_mu,dmu)

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,rho_comp) <= base_cutoff_density) then
                nabla = ZERO
                nabla_mu = ZERO
             else
                ! compute gradient

                ! forward difference
                if (j == lo(2)) then
                   dt = state(i,j+1,temp_comp) - state(i,j,temp_comp)
                   dp = pres(i,j+1) - pres(i,j)
                   dmu = abar(i,j+1) - abar(i,j)
                  ! backward difference
                else if (j == hi(2)) then
                   dt = state(i,j,temp_comp) - state(i,j-1,temp_comp)
                   dp = pres(i,j) - pres(i,j-1)
		   dmu = abar(i,j) - abar(i,j-1)
		  ! centered difference
                else
                   dt = state(i,j+1,temp_comp) - state(i,j-1,temp_comp)
                   dp = pres(i,j+1) - pres(i,j-1)
		   dmu = abar(i,j+1) - abar(i,j-1)
                endif


                ! prevent Inf
                if (dp == ZERO) then
                   nabla = -huge(ZERO)
                   nabla_mu = -huge(ZERO)
                else
                   nabla = pres(i,j)*dt / (dp*state(i,j,temp_comp))
                   nabla_mu = pres(i,j)*dmu / (dp*abar(i,j))
                endif
             endif
             if (hp(j) == ZERO) then
              brunt(i,j) = -huge(ZERO)
             else 
	      brunt(i,j) = - grav(j)*chi_t(i,j) / (chi_rho(i,j) * hp(j)) &
                            * (nabla_ad(i,j) - nabla &
                            - nabla_mu * chi_mu(i,j)/chi_t(i,j))
	     endif
          enddo
       enddo
    !$OMP END PARALLEL DO
    
    
    
  end subroutine make_brunt_2d

  subroutine make_brunt_3d(brunt, ng_ad, state, ng_s, grav, hp, &
                                    lo, hi, dx)

    use variables, only: rho_comp, temp_comp, spec_comp
    use geometry, only: nr_fine, nlevs_radial
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: brunt(lo(1)-ng_ad:,lo(2)-ng_ad:,lo(3)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: hp(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: abar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_rho(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_t(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_mu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: nabla, nabla_mu
    real(kind=dp_t) :: dp, dt, dmu

   integer :: i, j, k, c

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)
   
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state,pt_index)    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)       

             call eos(eos_input_rt, eos_state, pt_index)
       
             pres(i,j,k) = eos_state%p
             abar(i,j,k) = eos_state%abar

             chi_rho(i,j,k) = eos_state%rho * eos_state%dpdr / eos_state%p
             chi_t(i,j,k) = eos_state%T * eos_state%dpdt / eos_state%p
             chi_mu(i,j,k) = eos_state%abar * eos_state%dpdA / eos_state%p
             nabla_ad(i,j,k) = (eos_state%gam1 - chi_rho(i,j,k)) / (chi_t(i,j,k) * eos_state%gam1)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,dt,dp,nabla,nabla_mu,dmu)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,rho_comp) <= base_cutoff_density) then
                nabla = ZERO
                nabla_mu = ZERO
             else
                ! compute gradient

                ! forward difference
                if (k == lo(3)) then
                   dt = state(i,j,k+1,temp_comp) - state(i,j,k,temp_comp)
                   dp = pres(i,j,k+1) - pres(i,j,k)
                   dmu = abar(i,j,k+1) - abar(i,j,k)
                   ! backward difference
                else if (k == hi(3)) then
                   dt = state(i,j,k,temp_comp) - state(i,j,k-1,temp_comp)
                   dp = pres(i,j,k) - pres(i,j,k-1)
                   dmu = abar(i,j,k) - abar(i,j,k-1)
                ! centered difference
                else
                   dt = state(i,j,k+1,temp_comp) - state(i,j,k-1,temp_comp)
                   dp = pres(i,j,k+1) - pres(i,j,k-1)   
                   dmu = abar(i,j,k+1) - abar(i,j,k-1)
                endif


                ! prevent Inf
                if (dp == ZERO) then
                   nabla = -huge(ZERO)
                   nabla_mu = -huge(ZERO)
                else
                   nabla = pres(i,j,k)*dt / (dp*state(i,j,k,temp_comp))
                   nabla_mu = pres(i,j,k)*dmu / (dp*abar(i,j,k))
                endif
             endif
             if (hp(k) == ZERO) then
              brunt(i,j,k) = -huge(ZERO)
             else 
	      brunt(i,j,k) = - grav(k)*chi_t(i,j,k) / (chi_rho(i,j,k) * hp(k)) * &
                                (nabla_ad(i,j,k) - nabla &
                                - nabla_mu * chi_mu(i,j,k)/chi_t(i,j,k))
	     endif

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO    

  end subroutine make_brunt_3d

  subroutine make_brunt_3d_sphr(brunt, ng_ad, state, ng_s, grav, hp, &
                                    normal, ng_n, lo, hi, dx)

    use variables, only: rho_comp, temp_comp, spec_comp
    use geometry, only: nr_fine, nlevs_radial
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use fill_3d_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s, ng_n
    real(kind=dp_t), intent(  out) :: brunt(lo(1)-ng_ad:,lo(2)-ng_ad:,lo(3)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)  
    real(kind=dp_t), intent(in   ) ::   hp(0:)
    real(kind=dp_t), intent(in   ) ::   grav(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: abar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_rho(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_t(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_mu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: nabla, nabla_mu
    real(kind=dp_t) :: dp(4), dt(4), dmu(4)

    

    

    real(kind=dp_t), allocatable ::   hp_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: grav_cart(:,:,:,:)
    
    integer :: i, j, k, c

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    
    allocate(grav_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,grav(0:),grav_cart,lo,hi,dx,0)
    
    !hp(n,:) where n is the amr level
    allocate(hp_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,hp(0:),hp_cart,lo,hi,dx,0)
    
    
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state,pt_index)    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)       

             call eos(eos_input_rt, eos_state, pt_index)
       
             pres(i,j,k) = eos_state%p
             abar(i,j,k) = eos_state%abar

             chi_rho(i,j,k) = eos_state%rho * eos_state%dpdr / eos_state%p
             chi_t(i,j,k) = eos_state%T * eos_state%dpdt / eos_state%p
             chi_mu(i,j,k) = eos_state%abar * eos_state%dpdA / eos_state%p
             nabla_ad(i,j,k) = (eos_state%gam1 - chi_rho(i,j,k)) / (chi_t(i,j,k) * eos_state%gam1)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,dt,dp,nabla,nabla_mu,dmu)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,rho_comp) <= base_cutoff_density) then
                nabla = ZERO
                nabla_mu = ZERO
             else
                ! compute gradient

                ! forward difference
                if (k == lo(3)) then
                   dt(3) = state(i,j,k+1,temp_comp) - state(i,j,k,temp_comp)
                   dp(3) = pres(i,j,k+1) - pres(i,j,k)
                   dmu(3) = abar(i,j,k+1) - abar(i,j,k)
                   ! backward difference
                else if (k == hi(3)) then
                   dt(3) = state(i,j,k,temp_comp) - state(i,j,k-1,temp_comp)
                   dp(3) = pres(i,j,k) - pres(i,j,k-1)
                   dmu(3) = abar(i,j,k) - abar(i,j,k-1)
                ! centered difference
                else
                   dt(3) = state(i,j,k+1,temp_comp) - state(i,j,k-1,temp_comp)
                   dp(3) = pres(i,j,k+1) - pres(i,j,k-1)   
                   dmu(3) = abar(i,j,k+1) - abar(i,j,k-1)
                endif

                if (j == lo(2)) then
                   dt(2) = state(i,j+1,k,temp_comp) - state(i,j,k,temp_comp)
                   dp(2) = pres(i,j+1,k) - pres(i,j,k)
                   dmu(2) = abar(i,j+1,k) - abar(i,j,k)
                  ! backward difference
                else if (j == hi(2)) then
                   dt(2) = state(i,j,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp(2) = pres(i,j,k) - pres(i,j-1,k)
		   dmu(2) = abar(i,j,k) - abar(i,j-1,k)
		  ! centered difference
                else
                   dt(2) = state(i,j+1,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp(2) = pres(i,j+1,k) - pres(i,j-1,k)
		   dmu(2) = abar(i,j+1,k) - abar(i,j-1,k)
                endif

                if (i == lo(1)) then
                   dt(1) = state(i+1,j,k,temp_comp) - state(i,j,k,temp_comp)
                   dp(1) = pres(i+1,j,k) - pres(i,j,k)
                   dmu(1) = abar(i+1,j,k) - abar(i,j,k)
                   ! backward difference
                else if (i == hi(1)) then
                   dt(1) = state(i,j,k,temp_comp) - state(i-1,j,k,temp_comp)
                   dp(1) = pres(i,j,k) - pres(i-1,j,k)
                   dmu(1) = abar(i,j,k) - abar(i-1,j,k)
                   ! centered difference
                else
                   dt(1) = state(i+1,j,k,temp_comp) - state(i-1,j,k,temp_comp)
                   dp(1) = pres(i+1,j,k) - pres(i-1,j,k)
                   dmu(1) = abar(i+1,j,k) - abar(i-1,j,k)
                endif

                ! dot into normal to get d/dr
                dp(4) = 0.d0
                dt(4) = 0.d0
                dmu(4) = 0.d0
                do c = 1,3
                   dp(4) = dp(4) + dp(c)*normal(i,j,k,c) 
                   dt(4) = dt(4) + dt(c)*normal(i,j,k,c)
                   dmu(4) = dmu(4) + dmu(c)*normal(i,j,k,c)
                enddo

                ! prevent Inf
                if (dp(4) == ZERO) then
                   nabla = -huge(ZERO)
                   nabla_mu = -huge(ZERO)
                else
                   nabla = pres(i,j,k)*dt(4) / (dp(4)*state(i,j,k,temp_comp))
                   nabla_mu = pres(i,j,k)*dmu(4) / (dp(4)*abar(i,j,k))
                endif
             endif
             if (hp_cart(i,j,k,1) == ZERO) then
              brunt(i,j,k) = -huge(ZERO)
             else 
	      brunt(i,j,k) = - grav_cart(i,j,k,1)*chi_t(i,j,k) &
                                / (chi_rho(i,j,k) * hp_cart(i,j,k,1)) &
                                * (nabla_ad(i,j,k) - nabla &
                                - nabla_mu * chi_mu(i,j,k)/chi_t(i,j,k))
	     endif

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
    deallocate(grav_cart)
    deallocate(hp_cart)

  end subroutine make_brunt_3d_sphr




end module make_brunt_freq_module
  