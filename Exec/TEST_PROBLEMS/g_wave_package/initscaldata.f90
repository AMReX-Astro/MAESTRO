module init_scalar_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use fill_3d_module
  use eos_module
  use variables
  use network
  use geometry
  use ml_layout_module
  use ml_cc_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initscalardata, initscalardata_on_level

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla)

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n,dm,nlevs
    
    dm = mla%dim
    nlevs = mla%nlevel

    ng = s(1)%ng

    do n=1,nlevs
       do i = 1, nfabs(s(n))
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
                                    p0_init(n,:))
          case (3)
             if (spherical .eq. 1) then
                call bl_error("ERROR: spherical not implemented in initdata")
             else
                call bl_error("ERROR: 3d not implemented in initdata")
             end if
          end select
       end do
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(s(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),rho_comp,dm+rho_comp,nscal,bc(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),rho_comp,dm+rho_comp,nscal)

       enddo

    end if

  end subroutine initscalardata

  subroutine initscalardata_on_level(n,s,s0_init,p0_init,dx,bc)

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    ! local
    integer                  :: ng,i,dm
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    real(kind=dp_t), pointer :: sop(:,:,:,:)

    ng = nghost(s)
    dm = get_dim(s)

    do i = 1, nfabs(s)
       sop => dataptr(s,i)
       lo =  lwb(get_box(s,i))
       hi =  upb(get_box(s,i))
       select case (dm)
       case (2)
          call initscalardata_2d(sop(:,:,1,:),lo,hi,ng,dx,s0_init,p0_init)
       case (3)
          if (spherical .eq. 1) then
             call bl_error("ERROR: spherical not implemented in initdata")
          else
             call bl_error("ERROR: 3d not implemented in initdata")
          end if
       end select
    end do

    call multifab_fill_boundary(s)

    call multifab_physbc(s,rho_comp,dm+rho_comp,nscal,bc)

  end subroutine initscalardata_on_level

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only : prob_lo, prob_hi, wave_angle, mach_max, T_0, rho_0, grav_const
    use extern_probin_module, only: eos_gamma
    use bl_constants_module
    use fundamental_constants_module
    use eos_type_module

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables
    integer :: i, j, n, comp
    real (kind=dp_t) :: x, y
    
    real (kind=dp_t) :: wave_k, wave_kx, wave_ky 
    real (kind=dp_t) :: N0,hp,amplitude,R,ex
    
    type (eos_t) :: eos_state
 
 
    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          s(i,j,rho_comp)  = s0_init(j,rho_comp)
          s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
          s(i,j,temp_comp) = s0_init(j,temp_comp)
          s(i,j,spec_comp:spec_comp+nspec-1) = &
               s0_init(j,spec_comp:spec_comp+nspec-1)
          s(i,j,trac_comp:trac_comp+ntrac-1) = &
               s0_init(j,trac_comp:trac_comp+ntrac-1)
       enddo
    enddo
        
    wave_k = TWO * M_PI * TEN/ (prob_hi(2) - prob_lo(2))
    wave_kx = wave_k * cos(wave_angle)
    wave_ky = wave_k * sin(wave_angle)
    
    !compute the gas constant R
    R = k_B * n_A
    ! the pressure scale height is, assuming a molecular weight of 1
    hp = -R*T_0/(grav_const)
    ! N_0 is then
    N0 = sqrt(-grav_const/hp * (eos_gamma - ONE)/eos_gamma)
    
    ! initial the pressure perturbations
    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j)+HALF) * dx(2)
       
       amplitude = mach_max * sqrt(eos_gamma*R*T_0) * exp(-HALF*(y*FOUR/prob_lo(2))**TWO)
       
       do i = lo(1), hi(1)
          x = prob_lo(1) + (dble(i)+HALF) * dx(1)

          ex = cos(wave_kx * x + wave_ky * y)
          s(i,j,pi_comp) = s(i,j,pi_comp) - &
                           rho_0*cos(wave_angle)*N0*wave_ky/wave_kx**TWO*amplitude*ex
          
          !make the state thermodynamically consistent again
          eos_state%T     = 1d-15 !a temperature guess
          eos_state%rho   = s(i,j,rho_comp)
          eos_state%p     = p0_init(j) + s(i,j,pi_comp)
          eos_state%xn(:) = s(i,j,spec_comp:spec_comp+nspec-1)/s(i,j,rho_comp)

          ! (rho,p) --> T, h
          call eos(eos_input_rp, eos_state)

          s(i,j,rhoh_comp) = s(i,j,rho_comp) * eos_state%h
          s(i,j,temp_comp) = eos_state%T
       enddo
    enddo

  end subroutine initscalardata_2d

  subroutine initscalardata_3d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, perturb_model
    
    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    s = ZERO

  end subroutine initscalardata_3d

end module init_scalar_module
