module init_module

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
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initscalardata

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla)

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(inout) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(dm),hi(dm),ng
    integer :: i,n
    integer :: ii,jj

    ng = s(1)%ng

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop => dataptr(s(n),i)
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
                                    p0_init(n,:))
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
                                         bc(n-1),bc(n),rho_comp,dm+rho_comp,nscal, &
                                         fill_crse_input=.false.)

       enddo

    end if

  end subroutine initscalardata

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, perturb_model, peak_temp, diff_coeff, t0

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)

    ! Local variables
    integer         :: i,j,iter
    integer, parameter :: max_iter = 50
    real, parameter :: tol = 1.e-12
    real(kind=dp_t) :: x,y,dist2,dens_zone,temp_zone,del_dens
    real(kind=dp_t) :: pres_zone, del_pres
    logical :: converged
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)

       y = prob_lo(2) + (dble(j)+HALF) * dx(2)

       do i = lo(1), hi(1)

          x = prob_lo(1) + (dble(i)+HALF) * dx(1)

          ! apply the guassian temperature pulse at constant density
          dist2 = (center(1) - x)**2 + (center(2) - y)**2

          temp_zone = (peak_temp-s0_init(j,temp_comp)) * &
               exp(-dist2/(FOUR*diff_coeff*t0))

          temp_eos(1) = temp_zone + s0_init(j,temp_comp)

          xn_eos(1,1:nspec) = s0_init(j,spec_comp:spec_comp+nspec-1) / &
                              s0_init(j,rho_comp)

!          den_eos(1) = s0_init(j,rho_comp)
          dens_zone = s0_init(j,rho_comp)
!          pres_zone = p0_init(j)
          p_eos(1) = p0_init(j)

          converged = .false.

          do iter = 1, max_iter
!             p_eos(1) = pres_zone
             den_eos(1) = dens_zone

             call eos( &
!                      eos_input_tp, &
                      eos_input_rt, &
                      den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)

!             del_pres = -(den_eos(1) - s0_init(j,rho_comp)) * dpdr_eos(1)
             del_dens = -(p_eos(1) - p0_init(j)) / dpdr_eos(1)

!             pres_zone = max(0.9*pres_zone, &
!                             min(pres_zone + del_pres, 1.1*pres_zone))
             dens_zone = max(0.9*dens_zone, &
                             min(dens_zone + del_dens, 1.1*dens_zone))

!             if (abs(del_pres) < tol*pres_zone) then
             if (abs(del_dens) < tol*dens_zone) then
                converged = .true.
                exit
             endif
          enddo

          if (.not. converged) &
             call bl_error("density iter did not converge in initscalars")

          ! call eos one last time
!          p_eos(1) = pres_zone
          den_eos(1) = dens_zone

          call eos( &
!                   eos_input_tp, &
                   eos_input_rt, &
                   den_eos, temp_eos, &
                   npts, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   do_diag)

          s(i,j,rho_comp)  = den_eos(1)
          s(i,j,rhoh_comp) = den_eos(1) * h_eos(1)
          s(i,j,temp_comp) = temp_eos(1)
          s(i,j,spec_comp:spec_comp+nspec-1) = &
               xn_eos(1,1:nspec) * den_eos(1)
          s(i,j,trac_comp:trac_comp+ntrac-1) = &
                                  s0_init(j,trac_comp:trac_comp+ntrac-1)
          p0_init(j) = p_eos(1)
                    
       enddo
    enddo
    
  end subroutine initscalardata_2d


end module init_module
