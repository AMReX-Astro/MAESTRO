module thermal_conduct_module

  use bl_types
  use multifab_module
  use bl_constants_module
  use define_bc_module
  use ml_layout_module
  use bndry_reg_module
  use fill_3d_module

  implicit none

  private

  public :: thermal_conduct
  public :: put_beta_on_faces_2d, put_beta_on_faces_3d

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Crank-Nicholson solve for enthalpy, taking into account only the
  ! enthalpy-diffusion terms in the temperature conduction term.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine thermal_conduct(mla,dx,dt,s1,s_for_new_coeff,s2,p0_old,p0_new,tempbar, &
                             the_bc_tower)

    use variables, only: foextrap_comp, rho_comp, spec_comp, rhoh_comp
    use macproject_module
    use network, only: nspec
    use rhoh_vs_t_module
    use probin_module, ONLY: thermal_diffusion_type, use_tfromp
    use bl_prof_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use geometry, only: dm, nlevs

    type(ml_layout), intent(inout) :: mla
    real(dp_t)     , intent(in   ) :: dx(:,:),dt
    type(multifab) , intent(in   ) :: s1(:)
    type(multifab) , intent(in   ) :: s_for_new_coeff(:)
    type(multifab) , intent(inout) :: s2(:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:),p0_new(:,0:),tempbar(:,0:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! Local
    type(multifab) :: rhsalpha(mla%nlevel),lhsalpha(mla%nlevel)
    type(multifab) :: rhsbeta(mla%nlevel),lhsbeta(mla%nlevel)
    type(multifab) :: phi(mla%nlevel),Lphi(mla%nlevel),rhs(mla%nlevel)
    type(multifab) :: hcoeff1(mla%nlevel),hcoeff2(mla%nlevel)
    type(multifab) :: Xkcoeff1(mla%nlevel),Xkcoeff2(mla%nlevel)
    type(multifab) :: pcoeff1(mla%nlevel),pcoeff2(mla%nlevel)

    real(kind=dp_t), pointer    :: s1p(:,:,:,:)
    real(kind=dp_t), pointer    :: s_for_new_coeffp(:,:,:,:)
    real(kind=dp_t), pointer    :: rhsbetap(:,:,:,:),lhsbetap(:,:,:,:)
    real(kind=dp_t), pointer    :: hcoeff1p(:,:,:,:),hcoeff2p(:,:,:,:)
    real(kind=dp_t), pointer    :: Xkcoeff1p(:,:,:,:),Xkcoeff2p(:,:,:,:)
    real(kind=dp_t), pointer    :: pcoeff1p(:,:,:,:),pcoeff2p(:,:,:,:)
    integer                     :: stencil_order
    integer                     :: i,n,comp
    integer                     :: lo(dm),hi(dm)
    integer                     :: ng_s,ng_h,ng_X,ng_p,ng_cc,ng_fc
    type(bndry_reg)             :: fine_flx(2:mla%nlevel)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "therm_cond_full_alg")

    stencil_order = 2

    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    do n = 1,nlevs
       call multifab_build( hcoeff1(n), mla%la(n), 1,     1)
       call multifab_build(Xkcoeff1(n), mla%la(n), nspec, 1)
       call multifab_build( pcoeff1(n), mla%la(n), 1,     1)
    end do

    ! compute hcoeff1, Xkcoeff1, and pcoeff1
    ! defined as:
    ! hcoeff1 = -(dt/2)k_{th}^{(1)}/c_p^{(1)}
    ! Xkcoeff1 = (dt/2)\xi_k^{(1)}k_{th}^{(1)}/c_p^{(1)}
    ! pcoeff1 =  (dt/2)h_p^{(1)}k_{th}^{(1)}/c_p^{(1)}

    ng_s = s1(1)%ng
    ng_h = hcoeff1(1)%ng
    ng_X = Xkcoeff1(1)%ng
    ng_p = pcoeff1(1)%ng

    do n=1,nlevs
       do i=1,s1(n)%nboxes
          if (multifab_remote(s1(n),i)) cycle
          s1p       => dataptr(s1(n),i)
          hcoeff1p  => dataptr(hcoeff1(n),i)
          Xkcoeff1p => dataptr(Xkcoeff1(n),i)
          pcoeff1p  => dataptr(pcoeff1(n),i)
          lo = lwb(get_box(s1(n), i))
          hi = upb(get_box(s1(n), i))
          select case (dm)
          case (2)
             call compute_thermo_quantities_2d(lo,hi,dt,s1p(:,:,1,:),ng_s, &
                                               hcoeff1p(:,:,1,1),ng_h,&
                                               Xkcoeff1p(:,:,1,:),ng_X, &
                                               pcoeff1p(:,:,1,1),ng_p)
          case (3)
             call compute_thermo_quantities_3d(lo,hi,dt,s1p(:,:,:,:),ng_s, &
                                               hcoeff1p(:,:,:,1),ng_h, &
                                               Xkcoeff1p(:,:,:,:),ng_X, &
                                               pcoeff1p(:,:,:,1),ng_p)
          end select
       end do
    enddo

    do n = 1,nlevs
       call multifab_build( hcoeff2(n), mla%la(n), 1,     1)
       call multifab_build(Xkcoeff2(n), mla%la(n), nspec, 1)
       call multifab_build( pcoeff2(n), mla%la(n), 1,     1)
    end do

    ! compute hcoeff2, Xkcoeff2, and pcoeff2
    ! defined as (if we're in the first call to this function):
    ! hcoeff2 = -(dt/2)k_{th}^{(1)}/c_p^{(1)}
    ! Xkcoeff2 = (dt/2)\xi_k^{(1)}k_{th}^{(1)}/c_p^{(1)}
    ! pcoeff2 =  (dt/2)h_p^{(1)}k_{th}^{(1)}/c_p^{(1)}
    ! or defined as (if we're in the second call to this function):
    ! hcoeff2 = -(dt/2)k_{th}^{(2),*}/c_p^{(2),*}
    ! Xkcoeff2 = (dt/2)\xi_k^{(2),*}k_{th}^{(2),*}/c_p^{(2),*}
    ! pcoeff2 =  (dt/2)h_p^{(2),*}k_{th}^{(2),*}/c_p^{(2),*}

    ng_s = s_for_new_coeff(1)%ng
    ng_h = hcoeff2(1)%ng
    ng_X = Xkcoeff2(1)%ng
    ng_p = pcoeff2(1)%ng

    do n=1,nlevs
       do i=1,s_for_new_coeff(n)%nboxes
          if (multifab_remote(s_for_new_coeff(n),i)) cycle
          s_for_new_coeffp => dataptr(s_for_new_coeff(n),i)
          hcoeff2p  => dataptr(hcoeff2(n),i)
          Xkcoeff2p => dataptr(Xkcoeff2(n),i)
          pcoeff2p  => dataptr(pcoeff2(n),i)
          lo = lwb(get_box(s_for_new_coeff(n), i))
          hi = upb(get_box(s_for_new_coeff(n), i))
          select case (dm)
          case (2)
             call compute_thermo_quantities_2d(lo,hi,dt,s_for_new_coeffp(:,:,1,:),ng_s, &
                                               hcoeff2p(:,:,1,1),ng_h, &
                                               Xkcoeff2p(:,:,1,:),ng_X, &
                                               pcoeff2p(:,:,1,1),ng_p)
          case (3)
             call compute_thermo_quantities_3d(lo,hi,dt,s_for_new_coeffp(:,:,:,:),ng_s, &
                                               hcoeff2p(:,:,:,1),ng_h, &
                                               Xkcoeff2p(:,:,:,:),ng_X, &
                                               pcoeff2p(:,:,:,1),ng_p)
          end select
       end do
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add enthalpy diffusion to rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       call multifab_build(rhsbeta(n), mla%la(n), dm, 1)
    end do

    ! put beta on faces
    ng_cc = hcoeff1(1)%ng
    ng_fc = rhsbeta(1)%ng

    do n=1,nlevs
       do i=1,rhsbeta(n)%nboxes
          if (multifab_remote(rhsbeta(n),i)) cycle
          hcoeff1p => dataptr(hcoeff1(n),i)
          rhsbetap => dataptr(rhsbeta(n),i)
          lo = lwb(get_box(rhsbeta(n), i))
          hi = upb(get_box(rhsbeta(n), i))
          select case (dm)
          case (2)
             call put_beta_on_faces_2d(lo,hcoeff1p(:,:,1,1),ng_cc,rhsbetap(:,:,1,:),ng_fc)
          case (3)
             call put_beta_on_faces_3d(lo,hcoeff1p(:,:,:,1),ng_cc,rhsbetap(:,:,:,:),ng_fc)
          end select
       end do
    enddo

    do n=1,nlevs
       call destroy(hcoeff1(n))
    end do

    ! set rhsalpha = 0
    ! set phi = h^{(1)}
    do n=1,nlevs
       call multifab_build(phi(n), mla%la(n),  1, 1)
       call multifab_build(rhsalpha(n), mla%la(n),  1, 1)
       call multifab_build(Lphi(n), mla%la(n),  1, 0)
       call setval(rhsalpha(n), ZERO, all=.true.)
       call multifab_copy_c(phi(n),1,s1(n),rhoh_comp,1,1)
       call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,1)
    end do

    ! apply the operator
    call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                     dm+rhoh_comp,stencil_order,mla%mba%rr)

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
       ! put beta on faces
       ng_cc = Xkcoeff1(1)%ng
       ng_fc = rhsbeta(1)%ng

       do n=1,nlevs
          do i=1,rhsbeta(n)%nboxes
             if (multifab_remote(rhsbeta(n),i)) cycle
             Xkcoeff1p => dataptr(Xkcoeff1(n),i)
             rhsbetap  => dataptr(rhsbeta(n),i)
             lo = lwb(get_box(rhsbeta(n), i))
             hi = upb(get_box(rhsbeta(n), i))
             select case (dm)
             case (2)
                call put_beta_on_faces_2d(lo,Xkcoeff1p(:,:,1,comp),ng_cc, &
                                          rhsbetap(:,:,1,:),ng_fc)
             case (3)
                call put_beta_on_faces_3d(lo,Xkcoeff1p(:,:,:,comp),ng_cc, &
                                          rhsbetap(:,:,:,:),ng_fc)
             end select
          end do
       enddo

       ! load phi = X_k^{(1)}
       do n=1,nlevs
          call multifab_copy_c(phi(n),1,s1(n),spec_comp+comp-1,1,1)
          call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,1)
       enddo

       ! apply the operator
       call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                        dm+spec_comp+comp-1,stencil_order,mla%mba%rr)

       if(thermal_diffusion_type .eq. 1) then
          ! add lphi to rhs
          do n=1,nlevs
             call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
          enddo
       end if

       ! now do X_k^{(2)} term
       ! put beta on faces
       ng_cc = Xkcoeff2(1)%ng
       ng_fc = rhsbeta(1)%ng

       do n=1,nlevs
          do i=1,rhsbeta(n)%nboxes
             if (multifab_remote(rhsbeta(n),i)) cycle
             Xkcoeff2p => dataptr(Xkcoeff2(n),i)
             rhsbetap  => dataptr(rhsbeta(n),i)
             lo = lwb(get_box(rhsbeta(n), i))
             hi = upb(get_box(rhsbeta(n), i))
             select case (dm)
             case (2)
                call put_beta_on_faces_2d(lo,Xkcoeff2p(:,:,1,comp),ng_cc, &
                                          rhsbetap(:,:,1,:),ng_fc)
             case (3)
                call put_beta_on_faces_3d(lo,Xkcoeff2p(:,:,:,comp),ng_cc, &
                                          rhsbetap(:,:,:,:),ng_fc)
             end select
          end do
       enddo

       ! load phi = X_k^{(2)}
       do n=1,nlevs
          call multifab_copy_c(phi(n),1,s2(n),spec_comp+comp-1,1,1)
          call multifab_div_div_c(phi(n),1,s2(n),rho_comp,1,1)
       enddo

       ! apply the operator
       call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                        dm+spec_comp+comp-1,stencil_order,mla%mba%rr)

       ! add lphi to rhs
       do n=1,nlevs
          call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
       enddo
    enddo

    do n=1,nlevs
       call destroy(Xkcoeff1(n))
       call destroy(Xkcoeff2(n))
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! add pressure diffusion to rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! do p0_old term first
    ! put beta on faces
    ng_cc = pcoeff1(1)%ng
    ng_fc = rhsbeta(1)%ng

    do n=1,nlevs
       do i=1,rhsbeta(n)%nboxes
          if (multifab_remote(rhsbeta(n),i)) cycle
          pcoeff1p => dataptr(pcoeff1(n),i)
          rhsbetap => dataptr(rhsbeta(n),i)
          lo = lwb(get_box(rhsbeta(n), i))
          hi = upb(get_box(rhsbeta(n), i))
          select case (dm)
          case (2)
             call put_beta_on_faces_2d(lo,pcoeff1p(:,:,1,1),ng_cc,rhsbetap(:,:,1,:),ng_fc)
          case (3)
             call put_beta_on_faces_3d(lo,pcoeff1p(:,:,:,1),ng_cc,rhsbetap(:,:,:,:),ng_fc)
          end select
       end do
    enddo

    do n=1,nlevs
       call destroy(pcoeff1(n))
    end do

    call put_1d_array_on_cart(p0_old,phi,foextrap_comp,.false.,.false., &
                              dx,the_bc_tower%bc_tower_array,mla)
    ! apply the operator
    call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                     foextrap_comp,stencil_order,mla%mba%rr)

    if(thermal_diffusion_type .eq. 1) then
       ! add lphi to rhs
       do n=1,nlevs
          call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
       enddo
    end if

    ! now do p0_new term
    ! put beta on faces
    ng_cc = pcoeff2(1)%ng
    ng_fc = rhsbeta(1)%ng

    do n=1,nlevs
       do i=1,rhsbeta(n)%nboxes
          if (multifab_remote(rhsbeta(n),i)) cycle
          pcoeff2p => dataptr(pcoeff2(n),i)
          rhsbetap => dataptr(rhsbeta(n),i)
          lo = lwb(get_box(rhsbeta(n), i))
          hi = upb(get_box(rhsbeta(n), i))
          select case (dm)
          case (2)
             call put_beta_on_faces_2d(lo,pcoeff2p(:,:,1,1),ng_cc,rhsbetap(:,:,1,:),ng_fc)
          case (3)
             call put_beta_on_faces_3d(lo,pcoeff2p(:,:,:,1),ng_cc,rhsbetap(:,:,:,:),ng_fc)
          end select
       end do
    enddo

    do n=1,nlevs
       call destroy(pcoeff2(n))
    end do

    call put_1d_array_on_cart(p0_new,phi,foextrap_comp,.false.,.false., &
                              dx,the_bc_tower%bc_tower_array,mla)

    ! apply the operator
    call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                     foextrap_comp,stencil_order,mla%mba%rr)

    do n=1,nlevs
       call destroy(rhsalpha(n))
       call destroy(rhsbeta(n))
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
       call multifab_build(lhsbeta(n), mla%la(n), dm, 1)
    end do

    ! create lhsbeta = -hcoeff2 = (dt/2)k_{th}^{(2'')}/c_p^{(2'')}
    ! put beta on faces (remember to scale by -1 afterwards)
    ng_cc = hcoeff2(1)%ng
    ng_fc = lhsbeta(1)%ng

    do n=1,nlevs
       do i=1,lhsbeta(n)%nboxes
          if (multifab_remote(lhsbeta(n),i)) cycle
          hcoeff2p => dataptr(hcoeff2(n),i)
          lhsbetap => dataptr(lhsbeta(n),i)
          lo = lwb(get_box(lhsbeta(n), i))
          hi = upb(get_box(lhsbeta(n), i))
          select case (dm)
          case (2)
             call put_beta_on_faces_2d(lo,hcoeff2p(:,:,1,1),ng_cc,lhsbetap(:,:,1,:),ng_fc)
          case (3)
             call put_beta_on_faces_3d(lo,hcoeff2p(:,:,:,1),ng_cc,lhsbetap(:,:,:,:),ng_fc)
          end select
       end do
    enddo

    do n=1,nlevs
       call destroy(hcoeff2(n))
    end do

    ! scale by -1
    do n=1,nlevs
       call multifab_mult_mult_s_c(lhsbeta(n),1,-1.0d0,dm,1)
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now do the implicit solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! initialize phi to h^{(2'')} as a guess; also sets the ghost cells at inflow/outflow
    ! to a reasonable value
    do n=1,nlevs
       call multifab_copy_c(phi(n),1,s2(n),rhoh_comp,1,1)
       call multifab_div_div_c(phi(n),1,s2(n),rho_comp,1,1)
    enddo

    ! lhsalpha = \rho^{(2),*} or \rho^{(2)}
    do n=1,nlevs
       call multifab_build(lhsalpha(n), mla%la(n), 1, 1)
       call multifab_copy_c(lhsalpha(n),1,s2(n),rho_comp,1,1)
    enddo

    ! Call the solver to obtain h^(2) (it will be stored in phi)
    ! solves (alpha - nabla dot beta nabla)phi = rh
    call mac_multigrid(mla,rhs,phi,fine_flx,lhsalpha,lhsbeta,dx,the_bc_tower, &
                       dm+rhoh_comp,stencil_order,mla%mba%rr)

    do n=2,nlevs
       call destroy(fine_flx(n))
    end do

    do n=1,nlevs
       call destroy(lhsalpha(n))
       call destroy(lhsbeta(n))
       call destroy(rhs(n))
    end do

    ! load new rho*h into s2
    do n=1,nlevs
       call multifab_copy_c(s2(n),rhoh_comp,phi(n),1,1)
       call multifab_mult_mult_c(s2(n),rhoh_comp,s2(n),rho_comp,1)
    enddo

    do n=1,nlevs
       call destroy(phi(n))
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s2(nlevs),rhoh_comp,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s2(nlevs),rhoh_comp,dm+rhoh_comp,1, &
                            the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(s2(n-1),rhoh_comp,s2(n),rhoh_comp,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s2(n),s2(n-1),s2(1)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         rhoh_comp,dm+rhoh_comp,1,fill_crse_input=.false.)
       enddo

    end if

    ! compute updated temperature
    if (use_tfromp) then
       call makeTfromRhoP(s2,p0_new,s1,mla,the_bc_tower%bc_tower_array,dx)
    else
       call makeTfromRhoH(s2,s1,mla,the_bc_tower%bc_tower_array)
    end if

    call destroy(bpt)

  end subroutine thermal_conduct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute hcoeff and Xkcoeff, defined as:
  ! hcoeff = -(dt/2)k_{th}/c_p
  ! Xkcoeff = (dt/2)\xi_k k_{th}/c_p
  ! pcoeff = (dt/2)h_p*k_{th}/c_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_thermo_quantities_2d(lo,hi,dt,s,ng_s,hcoeff,ng_h,Xkcoeff,ng_X, &
                                          pcoeff,ng_p)

    use variables, only: rho_comp, spec_comp, temp_comp
    use eos_module
    use probin_module, ONLY: thermal_diffusion_type

    integer        , intent(in   ) :: lo(:),hi(:)
    integer        , intent(in   ) :: ng_s, ng_h, ng_X, ng_p
    real(dp_t)    ,  intent(in   ) :: dt
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::  hcoeff(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(inout) :: Xkcoeff(lo(1)-ng_X:,lo(2)-ng_X:,:)
    real(kind=dp_t), intent(inout) ::  pcoeff(lo(1)-ng_p:,lo(2)-ng_p:)

    ! Local
    integer :: i,j,comp

    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1

          den_eos(1) = s(i,j,rho_comp)
          temp_eos(1) = s(i,j,temp_comp)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

          do_diag = .false.

          call conducteos(eos_input_rt, den_eos, temp_eos, &
                          npts, nspec, &
                          xn_eos, &
                          p_eos, h_eos, e_eos, & 
                          cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                          dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                          dpdX_eos, dhdX_eos, &
                          gam1_eos, cs_eos, s_eos, &
                          dsdt_eos, dsdr_eos, &
                          do_diag, conduct_eos)

          hcoeff(i,j) = -dt*conduct_eos(1)/cp_eos(1)
          pcoeff(i,j) = dt*(conduct_eos(1)/cp_eos(1))*((1.0d0/den_eos(1))* &
               (1.0d0-p_eos(1)/(den_eos(1)*dpdr_eos(1)))+dedr_eos(1)/dpdr_eos(1))

          do comp=1,nspec
             Xkcoeff(i,j,comp) = dt*conduct_eos(1)*dhdX_eos(1,comp)/cp_eos(1)
          enddo

       enddo
    enddo

    if(thermal_diffusion_type .eq. 1) then
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             hcoeff(i,j) = HALF*hcoeff(i,j)
             pcoeff(i,j) = HALF*pcoeff(i,j)
             do comp=1,nspec
                Xkcoeff(i,j,comp) = HALF*Xkcoeff(i,j,comp)
             end do
          end do
       end do
    end if

  end subroutine compute_thermo_quantities_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! hcoeff = -(dt/2)k_{th}/c_p
  ! Xkcoeff = (dt/2)\xi_k k_{th}/c_p
  ! pcoeff = (dt/2)h_p*k_{th}/c_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_thermo_quantities_3d(lo,hi,dt,s,ng_s,hcoeff,ng_h,Xkcoeff,ng_X, &
                                          pcoeff,ng_p)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use probin_module, ONLY: thermal_diffusion_type
    use geometry, only: spherical

    integer        , intent(in   ) :: lo(:),hi(:)
    integer        , intent(in   ) :: ng_s,ng_h,ng_X,ng_p
    real(dp_t)    ,  intent(in   ) :: dt
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::  hcoeff(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(inout) :: Xkcoeff(lo(1)-ng_X:,lo(2)-ng_X:,lo(3)-ng_X:,:)
    real(kind=dp_t), intent(inout) ::  pcoeff(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)

    ! Local
    integer :: i,j,k,comp

    if(spherical .eq. 1) then
       call bl_error("compute_thermo1_quantities_3d spherical case not written!")
    endif

    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = s(i,j,k,temp_comp)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             do_diag = .false.

             call conducteos(eos_input_rt, den_eos, temp_eos, &
                             npts, nspec, &
                             xn_eos, &
                             p_eos, h_eos, e_eos, & 
                             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                             dpdX_eos, dhdX_eos, &
                             gam1_eos, cs_eos, s_eos, &
                             dsdt_eos, dsdr_eos, &
                             do_diag, conduct_eos)

             hcoeff(i,j,k) = -dt*conduct_eos(1)/cp_eos(1)
             pcoeff(i,j,k) = dt*(conduct_eos(1)/cp_eos(1))*((1.0d0/den_eos(1))* &
                  (1.0d0-p_eos(1)/(den_eos(1)*dpdr_eos(1)))+dedr_eos(1)/dpdr_eos(1))

             do comp=1,nspec
                Xkcoeff(i,j,k,comp) = dt*conduct_eos(1)*dhdX_eos(1,comp)/cp_eos(1)
             enddo

          enddo
       enddo
    enddo

    if(thermal_diffusion_type .eq. 1) then
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                hcoeff(i,j,k) = HALF*hcoeff(i,j,k)
                pcoeff(i,j,k) = HALF*pcoeff(i,j,k)
                do comp=1,nspec
                   Xkcoeff(i,j,k,comp) = HALF*Xkcoeff(i,j,k,comp)
                end do
             end do
          end do
       end do
    end if

  end subroutine compute_thermo_quantities_3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! put beta on faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine put_beta_on_faces_2d(lo,ccbeta,ng_cc,beta,ng_fc)

    integer        , intent(in   ) :: lo(:), ng_cc, ng_fc
    real(kind=dp_t), intent(in   ) :: ccbeta(lo(1)-ng_cc:,lo(2)-ng_cc:)
    real(kind=dp_t), intent(inout) ::   beta(lo(1)-ng_fc:,lo(2)-ng_fc:,:)

    ! Local
    integer :: i,j
    integer :: nx,ny

    nx = size(beta,dim=1) - 2
    ny = size(beta,dim=2) - 2

    do j = lo(2),lo(2)+ny-1
       do i = lo(1),lo(1)+nx
          beta(i,j,1) = TWO*(ccbeta(i,j)*ccbeta(i-1,j))/(ccbeta(i,j) + ccbeta(i-1,j))
       end do
    end do

    do j = lo(2),lo(2)+ny
       do i = lo(1),lo(1)+nx-1
          beta(i,j,2) = TWO*(ccbeta(i,j)*ccbeta(i,j-1))/(ccbeta(i,j) + ccbeta(i,j-1))
       end do
    end do

  end subroutine put_beta_on_faces_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! put beta on faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine put_beta_on_faces_3d(lo,ccbeta,ng_cc,beta,ng_fc)

    integer        , intent(in   ) :: lo(:), ng_cc, ng_fc
    real(kind=dp_t), intent(in   ) :: ccbeta(lo(1)-ng_cc:,lo(2)-ng_cc:,lo(3)-ng_cc:)
    real(kind=dp_t), intent(inout) ::   beta(lo(1)-ng_fc:,lo(2)-ng_fc:,lo(3)-ng_fc:,:)

    ! Local
    integer :: i,j,k
    integer :: nx,ny,nz

    nx = size(beta,dim=1) - 2
    ny = size(beta,dim=2) - 2
    nz = size(beta,dim=3) - 2

    do k = lo(3),lo(3)+nz-1
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx
             beta(i,j,k,1) = TWO*(ccbeta(i,j,k)*ccbeta(i-1,j,k))/(ccbeta(i,j,k) &
                  + ccbeta(i-1,j,k))
          end do
       end do
    end do

    do k = lo(3),lo(3)+nz-1
       do j = lo(2),lo(2)+ny
          do i = lo(1),lo(1)+nx-1
             beta(i,j,k,2) = TWO*(ccbeta(i,j,k)*ccbeta(i,j-1,k))/(ccbeta(i,j,k) &
                  + ccbeta(i,j-1,k))
          end do
       end do
    end do

    do k = lo(3),lo(3)+nz
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             beta(i,j,k,3) = TWO*(ccbeta(i,j,k)*ccbeta(i,j,k-1))/(ccbeta(i,j,k) &
                  + ccbeta(i,j,k-1))
          end do
       end do
    end do

  end subroutine put_beta_on_faces_3d

end module thermal_conduct_module
