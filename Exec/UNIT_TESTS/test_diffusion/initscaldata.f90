module init_scalar_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use fill_3d_module
  use variables
  use network
  use geometry
  use ml_layout_module
  use ml_cc_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initscalardata

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla,diffusion_coefficient)

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(inout) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla
    real(kind=dp_t), intent(in   ) :: diffusion_coefficient

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n, dm,nlevs
    integer :: ii,jj

    dm = mla%dim
    nlevs = mla%nlevel

    ng = s(1)%ng

    do n=1,nlevs
       do i = 1, nfabs(s(n))
          sp => dataptr(s(n),i)
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call initscalardata_2d(sp(:,:,1,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
                                    p0_init(n,:), diffusion_coefficient)
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

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init,p0_init,&
                               diffusion_coefficient)

    use probin_module, only: prob_lo, perturb_model, peak_h, t0, ambient_h
    use eos_module, only: eos_input_rt, eos
    use eos_type_module

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t), intent(inout) :: p0_init(0:)
    real(kind=dp_t), intent(in   ) :: diffusion_coefficient 

    ! Local variables
    integer         :: i,j,iter
    integer, parameter :: max_iter = 50
    real, parameter :: tol = 1.e-12
    real(kind=dp_t) :: x,y,dist2,dens_zone,temp_zone,del_dens, del_temp
    real(kind=dp_t) :: pres_zone, del_pres
    real(kind=dp_t) :: h_zone, dhdt
    logical :: converged
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)
    type (eos_t) :: eos_state

    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)

       y = prob_lo(2) + (dble(j)+HALF) * dx(2)

       do i = lo(1), hi(1)

          x = prob_lo(1) + (dble(i)+HALF) * dx(1)

          ! apply the guassian enthalpy pulse at constant density
          dist2 = (center(1) - x)**2 + (center(2) - y)**2

          h_zone = (peak_h - ambient_h) * &
                    exp(-dist2/(FOUR*diffusion_coefficient*t0)) + ambient_h

          temp_zone = s0_init(j,temp_comp)

          eos_state%xn(1:nspec) = s0_init(j,spec_comp:spec_comp+nspec-1) / &
                                  s0_init(j,rho_comp)

          eos_state%rho = s0_init(j,rho_comp)

          converged = .false.

          do iter = 1, max_iter
             eos_state%T = temp_zone

             call eos(eos_input_rt, eos_state)

             dhdt = eos_state%cv + eos_state%dpdt/eos_state%rho

             del_temp = -(eos_state%h - h_zone) / dhdt

             temp_zone = temp_zone + del_temp

             if (abs(del_temp) < tol*temp_zone) then
                converged = .true.
                exit
             endif
          enddo

          if (.not. converged) &
             call bl_error("iters did not converge in initscalars")

          ! call eos one last time
          eos_state%T = temp_zone

          call eos(eos_input_rt, eos_state)

          s(i,j,rho_comp)  = eos_state%rho
          s(i,j,rhoh_comp) = eos_state%rho * eos_state%h
          s(i,j,temp_comp) = eos_state%T
          s(i,j,spec_comp:spec_comp+nspec-1) = &
               eos_state%xn(1:nspec) * eos_state%rho
          s(i,j,trac_comp:trac_comp+ntrac-1) = &
                                  s0_init(j,trac_comp:trac_comp+ntrac-1)
          p0_init(j) = eos_state%p
                    
       enddo
    enddo
    
  end subroutine initscalardata_2d

end module init_scalar_module
