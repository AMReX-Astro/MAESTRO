module thermal_conduct_module

  use bl_types
  use multifab_module
  use bl_constants_module
  use define_bc_module
  use ml_layout_module
  use bndry_reg_module

  implicit none

  private

  public :: thermal_conduct_full_alg, thermal_conduct_half_alg
  public :: put_beta_on_faces_2d, put_beta_on_faces_3d
  public :: put_base_state_on_multifab_2d, put_base_state_on_multifab_3d

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Crank-Nicholson solve for enthalpy, taking into account only the
! enthalpy-diffusion terms in the temperature conduction term.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine thermal_conduct_full_alg(mla,dx,dt,s1,s_for_new_coeff,s2,p01,p02,t01,t02, &
                                    mg_verbose,cg_verbose,the_bc_tower)

  use variables, only: foextrap_comp, rho_comp, spec_comp, rhoh_comp
  use macproject_module
  use eos_module, only: nspec
  use rhoh_vs_t_module
  use probin_module, ONLY: use_big_h
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use ml_restriction_module, only: ml_cc_restriction_c

  type(ml_layout), intent(inout) :: mla
  real(dp_t)     , intent(in   ) :: dx(:,:),dt
  type(multifab) , intent(in   ) :: s1(:)
  type(multifab) , intent(in   ) :: s_for_new_coeff(:)
  type(multifab) , intent(inout) :: s2(:)
  real(kind=dp_t), intent(in   ) :: p01(:,0:),p02(:,0:),t01(:,0:),t02(:,0:)
  integer        , intent(in   ) :: mg_verbose,cg_verbose
  type(bc_tower) , intent(in   ) :: the_bc_tower

! Local
  type(multifab), allocatable :: rhsalpha(:),lhsalpha(:),rhsbeta(:),lhsbeta(:)
  type(multifab), allocatable :: ccbeta(:),phi(:),phitemp(:),Lphi(:),rhs(:)
  type(multifab), allocatable :: p01fab(:),p02fab(:)
  type(multifab), allocatable :: hcoeff1(:),hcoeff2(:),Xkcoeff1(:),Xkcoeff2(:)
  type(multifab), allocatable :: pcoeff1(:),pcoeff2(:)
  real(kind=dp_t), pointer    :: s1p(:,:,:,:),s2p(:,:,:,:),rhsalphap(:,:,:,:)
  real(kind=dp_t), pointer    :: s_for_new_coeffp(:,:,:,:)
  real(kind=dp_t), pointer    :: rhsbetap(:,:,:,:),lhsbetap(:,:,:,:)
  real(kind=dp_t), pointer    :: ccbetap(:,:,:,:),phip(:,:,:,:),rhsp(:,:,:,:)
  real(kind=dp_t), pointer    :: p01fabp(:,:,:,:),p02fabp(:,:,:,:)
  real(kind=dp_t), pointer    :: hcoeff1p(:,:,:,:),hcoeff2p(:,:,:,:)
  real(kind=dp_t), pointer    :: Xkcoeff1p(:,:,:,:),Xkcoeff2p(:,:,:,:)
  real(kind=dp_t), pointer    :: pcoeff1p(:,:,:,:),pcoeff2p(:,:,:,:)
  integer                     :: nlevs,dm,stencil_order
  integer                     :: i,n,comp,ng
  integer                     :: lo(s1(1)%dim),hi(s1(1)%dim)
  type(bndry_reg), pointer    :: fine_flx(:) => Null()

  nlevs = mla%nlevel
  dm = mla%dim
  stencil_order = 2
  ng = s2(1)%ng

  allocate(rhsalpha(nlevs),lhsalpha(nlevs))
  allocate(rhsbeta(nlevs),lhsbeta(nlevs),ccbeta(nlevs))
  allocate(phi(nlevs),phitemp(nlevs),Lphi(nlevs),rhs(nlevs))
  allocate(p01fab(nlevs),p02fab(nlevs))
  allocate(hcoeff1(nlevs),hcoeff2(nlevs))
  allocate(Xkcoeff1(nlevs),Xkcoeff2(nlevs))
  allocate(pcoeff1(nlevs),pcoeff2(nlevs))

  allocate(fine_flx(2:nlevs))
  do n = 2,nlevs
     call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
  end do

  do n = 1,nlevs
     call multifab_build(rhsalpha(n), mla%la(n),  1, 1)
     call multifab_build(lhsalpha(n), mla%la(n),  1, 1)
     call multifab_build( rhsbeta(n), mla%la(n), dm, 1)
     call multifab_build( lhsbeta(n), mla%la(n), dm, 1)
     call multifab_build(  ccbeta(n), mla%la(n),  1, 1)
     call multifab_build(     phi(n), mla%la(n),  1, 1)
     call multifab_build( phitemp(n), mla%la(n),  1, 1)
     call multifab_build(    Lphi(n), mla%la(n),  1, 0)
     call multifab_build(     rhs(n), mla%la(n),  1, 0)
     call multifab_build(  p01fab(n), mla%la(n),  1, 1)
     call multifab_build(  p02fab(n), mla%la(n),  1, 1)

     call multifab_build( hcoeff1(n), mla%la(n),  1,     1)
     call multifab_build( hcoeff2(n), mla%la(n),  1,     1)
     call multifab_build(Xkcoeff1(n), mla%la(n),  nspec, 1)
     call multifab_build(Xkcoeff2(n), mla%la(n),  nspec, 1)
     call multifab_build( pcoeff1(n), mla%la(n),  1,     1)
     call multifab_build( pcoeff2(n), mla%la(n),  1,     1)

     call setval(rhsalpha(n), ZERO, all=.true.)
     call setval(lhsalpha(n), ZERO, all=.true.)
     call setval(rhsbeta(n),  ZERO, all=.true.)
     call setval(lhsbeta(n),  ZERO, all=.true.)
     call setval(ccbeta(n),   ZERO, all=.true.)
     call setval(Lphi(n),     ZERO, all=.true.)
     call setval(phi(n),      ZERO, all=.true.)
     call setval(phitemp(n),  ZERO, all=.true.)
     call setval(rhs(n),      ZERO, all=.true.)
     call setval(p01fab(n),   ZERO, all=.true.)
     call setval(p02fab(n),   ZERO, all=.true.)

     call setval( hcoeff1(n), ZERO, all=.true.)
     call setval( hcoeff2(n), ZERO, all=.true.)
     call setval(Xkcoeff1(n), ZERO, all=.true.)
     call setval(Xkcoeff2(n), ZERO, all=.true.)
     call setval( pcoeff1(n), ZERO, all=.true.)
     call setval( pcoeff2(n), ZERO, all=.true.)
  end do

  ! create p01fab
  do n=1,nlevs
     do i=1,p01fab(n)%nboxes
        if (multifab_remote(p01fab(n),i)) cycle
        p01fabp => dataptr(p01fab(n),i)
        lo = lwb(get_box(p01fab(n), i))
        hi = upb(get_box(p01fab(n), i))
        select case (dm)
        case (2)
           call put_base_state_on_multifab_2d(lo,hi,p01(n,:),p01fabp(:,:,1,1))
        case (3)
           call put_base_state_on_multifab_3d(lo,hi,p01(n,:),p01fabp(:,:,:,1))
        end select
     end do
  enddo
  
  ! set the boundary conditions for p01
  do n=1,nlevs
     call multifab_fill_boundary(p01fab(n))
     call multifab_physbc(p01fab(n),1,foextrap_comp,1,dx(n,:),the_bc_tower%bc_tower_array(n))
  enddo

  ! create p02fab
  do n=1,nlevs
     do i=1,p02fab(n)%nboxes
        if (multifab_remote(p02fab(n),i)) cycle
        p02fabp => dataptr(p02fab(n),i)
        lo = lwb(get_box(p02fab(n), i))
        hi = upb(get_box(p02fab(n), i))
        select case (dm)
        case (2)
           call put_base_state_on_multifab_2d(lo,hi,p02(n,:),p02fabp(:,:,1,1))
        case (3)
           call put_base_state_on_multifab_3d(lo,hi,p02(n,:),p02fabp(:,:,:,1))
        end select
     end do
  enddo

  ! set the boundary conditions for p02
  do n=1,nlevs
     call multifab_fill_boundary(p02fab(n))
     call multifab_physbc(p02fab(n),1,foextrap_comp,1,dx(n,:),the_bc_tower%bc_tower_array(n))
  enddo

  ! lhsalpha = \rho^{(2),*} or \rho^{(2)}
  ! rhsalpha = 0 (already initialized above)
  ! thess will be true for this entire subroutine
  do n=1,nlevs
     call multifab_copy_c(lhsalpha(n),1,s2(n),rho_comp,1,1)
  enddo

  ! begin construction of rhs by setting rhs = \rho^{(2)}h^{(2')}
  do n=1,nlevs
     call multifab_copy_c(rhs(n),1,s2(n),rhoh_comp,1)
  enddo

  ! compute hcoeff1, Xkcoeff1, and pcoeff1
  ! defined as:
  ! hcoeff1 = -(dt/2)k_{th}^{(1)}/c_p^{(1)}
  ! Xkcoeff1 = (dt/2)\xi_k^{(1)}k_{th}^{(1)}/c_p^{(1)}
  ! pcoeff1 =  (dt/2)h_p^{(1)}k_{th}^{(1)}/c_p^{(1)}
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
           call compute_thermo_quantities_2d(lo,hi,dt, &
                                             s1p(:,:,1,:), &
                                             hcoeff1p(:,:,1,1), &
                                             Xkcoeff1p(:,:,1,:), &
                                             pcoeff1p(:,:,1,1))
        case (3)
           call compute_thermo_quantities_3d(lo,hi,dt,t01(n,:), &
                                             s1p(:,:,:,:), &
                                             hcoeff1p(:,:,:,1), &
                                             Xkcoeff1p(:,:,:,:), &
                                             pcoeff1p(:,:,:,1))
        end select
     end do
  enddo

  ! compute hcoeff^new, Xkcoeff^new, and pcoeff^new
  ! defined as (if we're in the first call to this function):
  ! hcoeff1 = -(dt/2)k_{th}^{(1)}/c_p^{(1)}
  ! Xkcoeff1 = (dt/2)\xi_k^{(1)}k_{th}^{(1)}/c_p^{(1)}
  ! pcoeff1 =  (dt/2)h_p^{(1)}k_{th}^{(1)}/c_p^{(1)}
  ! or defined as (if we're in the second call to this function):
  ! hcoeff1 = -(dt/2)k_{th}^{(2),*}/c_p^{(2),*}
  ! Xkcoeff1 = (dt/2)\xi_k^{(2),*}k_{th}^{(2),*}/c_p^{(2),*}
  ! pcoeff1 =  (dt/2)h_p^{(2),*}k_{th}^{(2),*}/c_p^{(2),*}
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
           call compute_thermo_quantities_2d(lo,hi,dt, &
                                             s_for_new_coeffp(:,:,1,:), &
                                             hcoeff2p(:,:,1,1), &
                                             Xkcoeff2p(:,:,1,:), &
                                             pcoeff2p(:,:,1,1))
        case (3)
           call compute_thermo_quantities_3d(lo,hi,dt,t02(n,:), &
                                             s_for_new_coeffp(:,:,:,:), &
                                             hcoeff2p(:,:,:,1), &
                                             Xkcoeff2p(:,:,:,:), &
                                             pcoeff2p(:,:,:,1))
        end select
     end do
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add enthalpy diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! put beta on faces
  do n=1,nlevs
     do i=1,rhsbeta(n)%nboxes
        if (multifab_remote(rhsbeta(n),i)) cycle
        hcoeff1p => dataptr(hcoeff1(n),i)
        rhsbetap => dataptr(rhsbeta(n),i)
        lo = lwb(get_box(rhsbeta(n), i))
        hi = upb(get_box(rhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,hcoeff1p(:,:,1,1), &
                                     rhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,hcoeff1p(:,:,:,1), &
                                     rhsbetap(:,:,:,:))
        end select
     end do
  enddo

  ! load phi = h^{(1)}
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,s1(n),rhoh_comp,1,1)
     call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,1)
  enddo

  ! apply the operator
  call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                   dm+rhoh_comp,stencil_order,mla%mba%rr, &
                   mg_verbose,cg_verbose)

  ! add Lphi to rhs
  do n=1,nlevs
     call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add species diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! loop over species
  do comp=1,nspec

     ! do X_k^{(1)} term first
     ! put beta on faces
     do n=1,nlevs
        do i=1,rhsbeta(n)%nboxes
           if (multifab_remote(rhsbeta(n),i)) cycle
           Xkcoeff1p => dataptr(Xkcoeff1(n),i)
           rhsbetap  => dataptr(rhsbeta(n),i)
           lo = lwb(get_box(rhsbeta(n), i))
           hi = upb(get_box(rhsbeta(n), i))
           select case (dm)
           case (2)
              call put_beta_on_faces_2d(lo,hi,Xkcoeff1p(:,:,1,comp), &
                                        rhsbetap(:,:,1,:))
           case (3)
              call put_beta_on_faces_3d(lo,hi,Xkcoeff1p(:,:,:,comp), &
                                        rhsbetap(:,:,:,:))
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
                      dm+spec_comp+comp-1,stencil_order,mla%mba%rr, &
                      mg_verbose,cg_verbose)
     
     ! add lphi to rhs
     do n=1,nlevs
        call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
     enddo

     ! now do X_k^{(2)} term
     ! put beta on faces
     do n=1,nlevs
        do i=1,rhsbeta(n)%nboxes
           if (multifab_remote(rhsbeta(n),i)) cycle
           Xkcoeff2p => dataptr(Xkcoeff2(n),i)
           rhsbetap  => dataptr(rhsbeta(n),i)
           lo = lwb(get_box(rhsbeta(n), i))
           hi = upb(get_box(rhsbeta(n), i))
           select case (dm)
           case (2)
              call put_beta_on_faces_2d(lo,hi,Xkcoeff2p(:,:,1,comp), &
                                        rhsbetap(:,:,1,:))
           case (3)
              call put_beta_on_faces_3d(lo,hi,Xkcoeff2p(:,:,:,comp), &
                                        rhsbetap(:,:,:,:))
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
                      dm+spec_comp+comp-1,stencil_order,mla%mba%rr, &
                      mg_verbose,cg_verbose)
     
     ! add lphi to rhs
     do n=1,nlevs
        call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add pressure diffusino to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! do p01 term first
  ! put beta on faces
  do n=1,nlevs
     do i=1,rhsbeta(n)%nboxes
        if (multifab_remote(rhsbeta(n),i)) cycle
        pcoeff1p => dataptr(pcoeff1(n),i)
        rhsbetap => dataptr(rhsbeta(n),i)
        lo = lwb(get_box(rhsbeta(n), i))
        hi = upb(get_box(rhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,pcoeff1p(:,:,1,1), &
                rhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,pcoeff1p(:,:,:,1), &
                rhsbetap(:,:,:,:))
        end select
     end do
  enddo

  ! load phi = p01
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,p01fab(n),1,1,1)
  enddo
  
  ! apply the operator
  call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
       foextrap_comp,stencil_order,mla%mba%rr, &
       mg_verbose,cg_verbose)
  
  ! add lphi to rhs
  do n=1,nlevs
     call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
  enddo
  
  ! now do p02 term
  ! put beta on faces
  do n=1,nlevs
     do i=1,rhsbeta(n)%nboxes
        if (multifab_remote(rhsbeta(n),i)) cycle
        pcoeff2p => dataptr(pcoeff2(n),i)
        rhsbetap => dataptr(rhsbeta(n),i)
        lo = lwb(get_box(rhsbeta(n), i))
        hi = upb(get_box(rhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,pcoeff2p(:,:,1,1), &
                rhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,pcoeff2p(:,:,:,1), &
                rhsbetap(:,:,:,:))
        end select
     end do
  enddo
  
  ! load phi = -02
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,p02fab(n),1,1,1)
  enddo
  
  ! apply the operator
  call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
       foextrap_comp,stencil_order,mla%mba%rr, &
       mg_verbose,cg_verbose)
  
  ! add lphi to rhs
  do n=1,nlevs
     call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
  enddo
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Setup LHS coefficients
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! create lhsbeta = -hcoeff2 = (dt/2)k_{th}^{(2'')}/c_p^{(2'')}
  ! put beta on faces (remember to scale by -1 afterwards)
  do n=1,nlevs
     do i=1,lhsbeta(n)%nboxes
        if (multifab_remote(lhsbeta(n),i)) cycle
        hcoeff2p => dataptr(hcoeff2(n),i)
        lhsbetap => dataptr(lhsbeta(n),i)
        lo = lwb(get_box(lhsbeta(n), i))
        hi = upb(get_box(lhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,hcoeff2p(:,:,1,1), &
                                     lhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,hcoeff2p(:,:,:,1), &
                                     lhsbetap(:,:,:,:))
        end select
     end do
  enddo

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

  ! Call the solver to obtain h^(2) (it will be stored in phi)
  ! solves (alpha - nabla dot beta nabla)phi = rh
  call mac_multigrid(mla,rhs,phi,fine_flx,lhsalpha,lhsbeta,dx,the_bc_tower, &
                     dm+rhoh_comp,stencil_order,mla%mba%rr, &
                     mg_verbose,cg_verbose)

  ! load new rho*h into s2
  do n=1,nlevs
     call multifab_copy_c(s2(n),rhoh_comp,phi(n),1,1)
     call multifab_mult_mult_c(s2(n),rhoh_comp,s2(n),rho_comp,1)

     call multifab_fill_boundary_c(s2(n),rhoh_comp,1)
     call multifab_physbc(s2(n),rhoh_comp,dm+rhoh_comp,1,dx(n,:), &
                          the_bc_tower%bc_tower_array(n))
  enddo

  do n=nlevs,2,-1
     call ml_cc_restriction_c(s2(n-1),rhoh_comp,s2(n),rhoh_comp,mla%mba%rr(n-1,:),1)
       
     call multifab_fill_ghost_cells(s2(n),s2(n-1), &
                                    ng,mla%mba%rr(n-1,:), &
                                    the_bc_tower%bc_tower_array(n-1), &
                                    the_bc_tower%bc_tower_array(n  ), &
                                    rhoh_comp,dm+rhoh_comp,1)
  enddo

  ! compute updated temperature
  call makeTfromRhoH(nlevs,s2,t02,mla,the_bc_tower%bc_tower_array,dx)

  do n = 1,nlevs
     call destroy(rhsalpha(n))
     call destroy(lhsalpha(n))
     call destroy(rhsbeta(n))
     call destroy(lhsbeta(n))
     call destroy(ccbeta(n))
     call destroy(phi(n))
     call destroy(phitemp(n))
     call destroy(Lphi(n))
     call destroy(rhs(n))
     call destroy(p01fab(n))
     call destroy(p02fab(n))
     call destroy(hcoeff1(n))
     call destroy(hcoeff2(n))
     call destroy(Xkcoeff1(n))
     call destroy(Xkcoeff2(n))
     call destroy(pcoeff1(n))
     call destroy(pcoeff2(n))
  enddo

  deallocate(rhsalpha,lhsalpha,rhsbeta,lhsbeta,ccbeta,phi,phitemp,Lphi,rhs)
  deallocate(p01fab,p02fab,fine_flx)
  deallocate(hcoeff1,hcoeff2,Xkcoeff1,Xkcoeff2,pcoeff1,pcoeff2)

end subroutine thermal_conduct_full_alg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Crank-Nicholson solve for enthalpy, taking into account only the
! enthalpy-diffusion terms in the temperature conduction term.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine thermal_conduct_half_alg(mla,dx,dt,s1,s2,p01,p02,t01,t02, &
                                    mg_verbose,cg_verbose,the_bc_tower)

  use variables, only: foextrap_comp, rho_comp, spec_comp, rhoh_comp, temp_comp
  use macproject_module
  use eos_module, only: nspec
  use rhoh_vs_t_module
  use probin_module, ONLY: use_big_h
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use ml_restriction_module, only: ml_cc_restriction_c

  type(ml_layout), intent(inout) :: mla
  real(dp_t)     , intent(in   ) :: dx(:,:),dt
  type(multifab) , intent(in   ) :: s1(:)
  type(multifab) , intent(inout) :: s2(:)
  real(kind=dp_t), intent(in   ) :: p01(:,0:),p02(:,0:),t01(:,0:),t02(:,0:)
  integer        , intent(in   ) :: mg_verbose,cg_verbose
  type(bc_tower) , intent(in   ) :: the_bc_tower

! Local
  type(multifab), allocatable :: rhsalpha(:),lhsalpha(:),rhsbeta(:),lhsbeta(:)
  type(multifab), allocatable :: ccbeta(:),phi(:),phitemp(:),Lphi(:),rhs(:)
  type(multifab), allocatable :: p01fab(:),p02fab(:)
  type(multifab), allocatable :: hcoeff1(:),hcoeff2(:),Xkcoeff1(:),Xkcoeff2(:)
  type(multifab), allocatable :: pcoeff1(:),pcoeff2(:)
  real(kind=dp_t), pointer    :: s1p(:,:,:,:),s2p(:,:,:,:),rhsalphap(:,:,:,:)
  real(kind=dp_t), pointer    :: rhsbetap(:,:,:,:),lhsbetap(:,:,:,:)
  real(kind=dp_t), pointer    :: ccbetap(:,:,:,:),phip(:,:,:,:),rhsp(:,:,:,:)
  real(kind=dp_t), pointer    :: p01fabp(:,:,:,:),p02fabp(:,:,:,:)
  real(kind=dp_t), pointer    :: hcoeff1p(:,:,:,:),hcoeff2p(:,:,:,:)
  real(kind=dp_t), pointer    :: Xkcoeff1p(:,:,:,:),Xkcoeff2p(:,:,:,:)
  real(kind=dp_t), pointer    :: pcoeff1p(:,:,:,:),pcoeff2p(:,:,:,:)
  integer                     :: nlevs,dm,stencil_order
  integer                     :: i,n,comp,ng
  integer                     :: lo(s1(1)%dim),hi(s1(1)%dim)
  type(bndry_reg), pointer    :: fine_flx(:) => Null()

  type(bc_level) ::  bc

  nlevs = mla%nlevel
  dm = mla%dim
  stencil_order = 2
  ng = s2(1)%ng

  allocate(rhsalpha(nlevs),lhsalpha(nlevs))
  allocate(rhsbeta(nlevs),lhsbeta(nlevs),ccbeta(nlevs))
  allocate(phi(nlevs),phitemp(nlevs),Lphi(nlevs),rhs(nlevs))
  allocate(p01fab(nlevs),p02fab(nlevs))
  allocate(hcoeff1(nlevs),hcoeff2(nlevs))
  allocate(Xkcoeff1(nlevs),Xkcoeff2(nlevs))
  allocate(pcoeff1(nlevs),pcoeff2(nlevs))

  allocate(fine_flx(2:nlevs))
  do n = 2,nlevs
     call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
  end do

  do n = 1,nlevs
     call multifab_build(rhsalpha(n), mla%la(n),  1, 1)
     call multifab_build(lhsalpha(n), mla%la(n),  1, 1)
     call multifab_build( rhsbeta(n), mla%la(n), dm, 1)
     call multifab_build( lhsbeta(n), mla%la(n), dm, 1)
     call multifab_build(  ccbeta(n), mla%la(n),  1, 1)
     call multifab_build(     phi(n), mla%la(n),  1, 1)
     call multifab_build( phitemp(n), mla%la(n),  1, 1)
     call multifab_build(    Lphi(n), mla%la(n),  1, 0)
     call multifab_build(     rhs(n), mla%la(n),  1, 0)
     call multifab_build(  p01fab(n), mla%la(n),  1, 1)
     call multifab_build(  p02fab(n), mla%la(n),  1, 1)

     call multifab_build( hcoeff1(n), mla%la(n),  1,     1)
     call multifab_build( hcoeff2(n), mla%la(n),  1,     1)
     call multifab_build(Xkcoeff1(n), mla%la(n),  nspec, 1)
     call multifab_build(Xkcoeff2(n), mla%la(n),  nspec, 1)
     call multifab_build( pcoeff1(n), mla%la(n),  1,     1)
     call multifab_build( pcoeff2(n), mla%la(n),  1,     1)

     call setval(rhsalpha(n), ZERO, all=.true.)
     call setval(lhsalpha(n), ZERO, all=.true.)
     call setval(rhsbeta(n),  ZERO, all=.true.)
     call setval(lhsbeta(n),  ZERO, all=.true.)
     call setval(ccbeta(n),   ZERO, all=.true.)
     call setval(Lphi(n),     ZERO, all=.true.)
     call setval(phi(n),      ZERO, all=.true.)
     call setval(phitemp(n),  ZERO, all=.true.)
     call setval(rhs(n),      ZERO, all=.true.)
     call setval(p01fab(n),   ZERO, all=.true.)
     call setval(p02fab(n),   ZERO, all=.true.)

     call setval( hcoeff1(n), ZERO, all=.true.)
     call setval( hcoeff2(n), ZERO, all=.true.)
     call setval(Xkcoeff1(n), ZERO, all=.true.)
     call setval(Xkcoeff2(n), ZERO, all=.true.)
     call setval( pcoeff1(n), ZERO, all=.true.)
     call setval( pcoeff2(n), ZERO, all=.true.)
  end do

  ! create p01fab
  do n=1,nlevs
     do i=1,p01fab(n)%nboxes
        if (multifab_remote(p01fab(n),i)) cycle
        p01fabp => dataptr(p01fab(n),i)
        lo = lwb(get_box(p01fab(n), i))
        hi = upb(get_box(p01fab(n), i))
        select case (dm)
        case (2)
           call put_base_state_on_multifab_2d(lo,hi,p01(n,:),p01fabp(:,:,1,1))
        case (3)
           call put_base_state_on_multifab_3d(lo,hi,p01(n,:),p01fabp(:,:,:,1))
        end select
     end do
  enddo

  ! set the boundary conditions for p01
  do n=1,nlevs
     call multifab_fill_boundary(p01fab(n))
     call multifab_physbc(p01fab(n),1,foextrap_comp,1,dx(n,:),the_bc_tower%bc_tower_array(n))
  enddo

  ! create p02fab
  do n=1,nlevs
     do i=1,p02fab(n)%nboxes
        if (multifab_remote(p02fab(n),i)) cycle
        p02fabp => dataptr(p02fab(n),i)
        lo = lwb(get_box(p02fab(n), i))
        hi = upb(get_box(p02fab(n), i))
        select case (dm)
        case (2)
           call put_base_state_on_multifab_2d(lo,hi,p02(n,:),p02fabp(:,:,1,1))
        case (3)
           call put_base_state_on_multifab_3d(lo,hi,p02(n,:),p02fabp(:,:,:,1))
        end select
     end do
  enddo

  ! set the boundary conditions for p02
  do n=1,nlevs
     call multifab_fill_boundary(p02fab(n))
     call multifab_physbc(p02fab(n),1,foextrap_comp,1,dx(n,:),the_bc_tower%bc_tower_array(n))
  enddo

  ! lhsalpha = \rho^{(2)}
  ! rhsalpha = 0 (already initialized above)
  ! thess will be true for this entire subroutine
  do n=1,nlevs
     call multifab_copy_c(lhsalpha(n),1,s2(n),rho_comp,1,1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!
  ! First implicit solve
  !!!!!!!!!!!!!!!!!!!!!!!

  ! compute hcoeff1 and Xkcoeff1, defined as:
  ! hcoeff1 = -(dt/2)k_{th}^{(1)}/c_p^{(1)}
  ! Xkcoeff1 = (dt/2)\xi_k^{(1)}k_{th}^{(1)}/c_p^{(1)}
  ! pcoeff1 =  (dt/2)h_p^{(1)}k_{th}^{(1)}/c_p^{(1)}
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
           call compute_thermo_quantities_2d(lo,hi,dt, &
                                             s1p(:,:,1,:), &
                                             hcoeff1p(:,:,1,1), &
                                             Xkcoeff1p(:,:,1,:), &
                                             pcoeff1p(:,:,1,1))
        case (3)
           call compute_thermo_quantities_3d(lo,hi,dt,t01(n,:), &
                                             s1p(:,:,:,:), &
                                             hcoeff1p(:,:,:,1), &
                                             Xkcoeff1p(:,:,:,:), &
                                             pcoeff1p(:,:,:,1))
        end select
     end do
  enddo

  ! begin construction of rhs by setting rhs = \rho^{(2)}h^{(2')}
  do n=1,nlevs
     call multifab_copy_c(rhs(n),1,s2(n),rhoh_comp,1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add enthalpy diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! put beta on faces
  do n=1,nlevs
     do i=1,rhsbeta(n)%nboxes
        if (multifab_remote(rhsbeta(n),i)) cycle
        hcoeff1p => dataptr(hcoeff1(n),i)
        rhsbetap => dataptr(rhsbeta(n),i)
        lo = lwb(get_box(rhsbeta(n), i))
        hi = upb(get_box(rhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,hcoeff1p(:,:,1,1), &
                                     rhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,hcoeff1p(:,:,:,1), &
                                     rhsbetap(:,:,:,:))
        end select
     end do
  enddo

  ! load phi = h^{(1)}
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,s1(n),rhoh_comp,1,1)
     call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,1)
  enddo

  ! apply the operator
  call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                   dm+rhoh_comp,stencil_order,mla%mba%rr, &
                   mg_verbose,cg_verbose)

  ! add Lphi to rhs
  do n=1,nlevs
     call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add species diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! loop over species
  do comp=1,nspec

     ! put beta on faces
     do n=1,nlevs
        do i=1,rhsbeta(n)%nboxes
           if (multifab_remote(rhsbeta(n),i)) cycle
           Xkcoeff1p => dataptr(Xkcoeff1(n),i)
           rhsbetap  => dataptr(rhsbeta(n),i)
           lo = lwb(get_box(rhsbeta(n), i))
           hi = upb(get_box(rhsbeta(n), i))
           select case (dm)
           case (2)
              call put_beta_on_faces_2d(lo,hi,Xkcoeff1p(:,:,1,comp), &
                                        rhsbetap(:,:,1,:))
           case (3)
              call put_beta_on_faces_3d(lo,hi,Xkcoeff1p(:,:,:,comp), &
                                        rhsbetap(:,:,:,:))
           end select
        end do
     enddo

     ! load phi = X_k^{(1)} + X_k^{(2)}
     do n=1,nlevs
        call multifab_copy_c(phi(n),1,s1(n),spec_comp+comp-1,1,1)
        call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,1)
        call multifab_copy_c(phitemp(n),1,s2(n),spec_comp+comp-1,1,1)
        call multifab_div_div_c(phitemp(n),1,s2(n),rho_comp,1,1)
        call multifab_plus_plus_c(phi(n),1,phitemp(n),1,1,1)
     enddo

     ! apply the operator
     call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                      dm+spec_comp+comp-1,stencil_order,mla%mba%rr, &
                      mg_verbose,cg_verbose)
     
     ! add lphi to rhs
     do n=1,nlevs
        call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add pressure diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! put beta on faces
  do n=1,nlevs
     do i=1,rhsbeta(n)%nboxes
        if (multifab_remote(rhsbeta(n),i)) cycle
        pcoeff1p => dataptr(pcoeff1(n),i)
        rhsbetap => dataptr(rhsbeta(n),i)
        lo = lwb(get_box(rhsbeta(n), i))
        hi = upb(get_box(rhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,pcoeff1p(:,:,1,1), &
                rhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,pcoeff1p(:,:,:,1), &
                rhsbetap(:,:,:,:))
        end select
     end do
  enddo
  
  ! load phi = p01 + p02
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,p01fab(n),1,1,1)
     call multifab_plus_plus_c(phi(n),1,p02fab(n),1,1,1)
  enddo
  
  ! apply the operator
  call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                   foextrap_comp,stencil_order,mla%mba%rr, &
                   mg_verbose,cg_verbose)  

  ! add lphi to rhs
  do n=1,nlevs
     call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Setup LHS coefficients
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! create lhsbeta = -hcoeff1 = (dt/2)k_{th}^{(1)}/c_p^{(1)}
  ! put beta on faces (remember to scale by -1 afterwards)
  do n=1,nlevs
     do i=1,lhsbeta(n)%nboxes
        if (multifab_remote(lhsbeta(n),i)) cycle
        hcoeff1p => dataptr(hcoeff1(n),i)
        lhsbetap => dataptr(lhsbeta(n),i)
        lo = lwb(get_box(lhsbeta(n), i))
        hi = upb(get_box(lhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,hcoeff1p(:,:,1,1), &
                                     lhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,hcoeff1p(:,:,:,1), &
                                     lhsbetap(:,:,:,:))
        end select
     end do
  enddo

  ! scale by -1
  do n=1,nlevs
     call multifab_mult_mult_s_c(lhsbeta(n),1,-1.0d0,dm,1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now do the implicit solve
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! initialize phi to h^{(2')} as a guess; also sets the ghost cells at inflow/outflow
  ! to a reasonable value
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,s2(n),rhoh_comp,1,1)
     call multifab_div_div_c(phi(n),1,s2(n),rho_comp,1,1)
  enddo

  ! Call the solver to obtain h^(2'') (it will be stored in phi)
  ! solves (alpha - nabla dot beta nabla)phi = rh
  call mac_multigrid(mla,rhs,phi,fine_flx,lhsalpha,lhsbeta,dx,the_bc_tower, &
                     dm+rhoh_comp,stencil_order,mla%mba%rr, &
                     mg_verbose,cg_verbose)

  ! need to set rhs for second implicit solve to 
  ! rho^{(2)}h^{(2')} before I overwrite s2 with rhoh^{2'}
  do n=1,nlevs
     call multifab_copy_c(rhs(n),1,s2(n),rhoh_comp,1)
  enddo

  ! load h^{2''} into s2
  do n=1,nlevs
     call multifab_copy_c(s2(n),rhoh_comp,phi(n),1,1)
     call multifab_mult_mult_c(s2(n),rhoh_comp,s2(n),rho_comp,1)


     call multifab_fill_boundary_c(s2(n),rhoh_comp,1)
     call multifab_physbc(s2(n),rhoh_comp,dm+rhoh_comp,1,dx(n,:), &
                          the_bc_tower%bc_tower_array(n))
  enddo

  do n=nlevs,2,-1
     call ml_cc_restriction_c(s2(n-1),rhoh_comp,s2(n),rhoh_comp,mla%mba%rr(n-1,:),1)
     call multifab_fill_ghost_cells(s2(n),s2(n-1), &
                                    ng,mla%mba%rr(n-1,:), &
                                    the_bc_tower%bc_tower_array(n-1), &
                                    the_bc_tower%bc_tower_array(n  ), &
                                    temp_comp,dm+temp_comp,1)
  enddo

  ! compute updated temperature
  call makeTfromRhoH(nlevs,s2,t02,mla,the_bc_tower%bc_tower_array,dx)

  !!!!!!!!!!!!!!!!!!!!!!!
  ! Second implicit solve
  !!!!!!!!!!!!!!!!!!!!!!!

  ! compute hcoeff2 and Xkcoeff2, defined as:
  ! hcoeff2 = -(dt/2)k_{th}^{(2'')}/c_p^{(2'')}
  ! Xkcoeff2 = (dt/2)\xi_k^{(2'')}k_{th}^{(2'')}/c_p^{(2'')}
  ! pcoeff2 =  (dt/2)h_p^{(2'')}k_{th}^{(2'')}/c_p^{(2'')}
  do n=1,nlevs
     do i=1,s2(n)%nboxes
        if (multifab_remote(s2(n),i)) cycle
        s2p       => dataptr(s2(n),i)
        hcoeff2p  => dataptr(hcoeff2(n),i)
        Xkcoeff2p => dataptr(Xkcoeff2(n),i)
        pcoeff2p  => dataptr(pcoeff2(n),i)
        lo = lwb(get_box(s2(n), i))
        hi = upb(get_box(s2(n), i))
        select case (dm)
        case (2)
           call compute_thermo_quantities_2d(lo,hi,dt, &
                                             s2p(:,:,1,:), &
                                             hcoeff2p(:,:,1,1), &
                                             Xkcoeff2p(:,:,1,:), &
                                             pcoeff2p(:,:,1,1))
        case (3)
           call compute_thermo_quantities_3d(lo,hi,dt,t02(n,:), &
                                             s2p(:,:,:,:), &
                                             hcoeff2p(:,:,:,1), &
                                             Xkcoeff2p(:,:,:,:), &
                                             pcoeff2p(:,:,:,1))
        end select
     end do
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add enthalpy diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! put beta on faces
  do n=1,nlevs
     do i=1,rhsbeta(n)%nboxes
        if (multifab_remote(rhsbeta(n),i)) cycle
        hcoeff1p => dataptr(hcoeff1(n),i)
        rhsbetap => dataptr(rhsbeta(n),i)
        lo = lwb(get_box(rhsbeta(n), i))
        hi = upb(get_box(rhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,hcoeff1p(:,:,1,1), &
                                     rhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,hcoeff1p(:,:,:,1), &
                                     rhsbetap(:,:,:,:))
        end select
     end do
  enddo

  ! load phi = h^{(1)}
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,s1(n),rhoh_comp,1,1)
     call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,1)
  enddo

  ! apply the operator
  call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                   dm+rhoh_comp,stencil_order,mla%mba%rr, &
                   mg_verbose,cg_verbose)

  ! add Lphi to rhs
  do n=1,nlevs
     call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add species diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! loop over species
  do comp=1,nspec

     ! do X_k^{(1)} term first
     ! put beta on faces
     do n=1,nlevs
        do i=1,rhsbeta(n)%nboxes
           if (multifab_remote(rhsbeta(n),i)) cycle
           Xkcoeff1p => dataptr(Xkcoeff1(n),i)
           rhsbetap  => dataptr(rhsbeta(n),i)
           lo = lwb(get_box(rhsbeta(n), i))
           hi = upb(get_box(rhsbeta(n), i))
           select case (dm)
           case (2)
              call put_beta_on_faces_2d(lo,hi,Xkcoeff1p(:,:,1,comp), &
                                        rhsbetap(:,:,1,:))
           case (3)
              call put_beta_on_faces_3d(lo,hi,Xkcoeff1p(:,:,:,comp), &
                                        rhsbetap(:,:,:,:))
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
                      dm+spec_comp+comp-1,stencil_order,mla%mba%rr, &
                      mg_verbose,cg_verbose)
     
     ! add lphi to rhs
     do n=1,nlevs
        call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
     enddo

     ! now do X_k^{(2)} term
     ! put beta on faces
     do n=1,nlevs
        do i=1,rhsbeta(n)%nboxes
           if (multifab_remote(rhsbeta(n),i)) cycle
           Xkcoeff2p => dataptr(Xkcoeff2(n),i)
           rhsbetap  => dataptr(rhsbeta(n),i)
           lo = lwb(get_box(rhsbeta(n), i))
           hi = upb(get_box(rhsbeta(n), i))
           select case (dm)
           case (2)
              call put_beta_on_faces_2d(lo,hi,Xkcoeff2p(:,:,1,comp), &
                                        rhsbetap(:,:,1,:))
           case (3)
              call put_beta_on_faces_3d(lo,hi,Xkcoeff2p(:,:,:,comp), &
                                        rhsbetap(:,:,:,:))
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
                      dm+spec_comp+comp-1,stencil_order,mla%mba%rr, &
                      mg_verbose,cg_verbose)
     
     ! add lphi to rhs
     do n=1,nlevs
        call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add pressure diffusino to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! do p01 term first
  ! put beta on faces
  do n=1,nlevs
     do i=1,rhsbeta(n)%nboxes
        if (multifab_remote(rhsbeta(n),i)) cycle
        pcoeff1p => dataptr(pcoeff1(n),i)
        rhsbetap => dataptr(rhsbeta(n),i)
        lo = lwb(get_box(rhsbeta(n), i))
        hi = upb(get_box(rhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,pcoeff1p(:,:,1,1), &
                rhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,pcoeff1p(:,:,:,1), &
                rhsbetap(:,:,:,:))
        end select
     end do
  enddo

  ! load phi = p01
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,p01fab(n),1,1,1)
  enddo
  
  ! apply the operator
  call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                   foextrap_comp,stencil_order,mla%mba%rr, &
                   mg_verbose,cg_verbose)
  
  ! add lphi to rhs
  do n=1,nlevs
     call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
  enddo
  
  ! now do p02 term
  ! put beta on faces
  do n=1,nlevs
     do i=1,rhsbeta(n)%nboxes
        if (multifab_remote(rhsbeta(n),i)) cycle
        pcoeff2p => dataptr(pcoeff2(n),i)
        rhsbetap => dataptr(rhsbeta(n),i)
        lo = lwb(get_box(rhsbeta(n), i))
        hi = upb(get_box(rhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,pcoeff2p(:,:,1,1), &
                rhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,pcoeff2p(:,:,:,1), &
                rhsbetap(:,:,:,:))
        end select
     end do
  enddo
  
  ! load phi = -02
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,p02fab(n),1,1,1)
  enddo
  
  ! apply the operator
  call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                   foextrap_comp,stencil_order,mla%mba%rr, &
                   mg_verbose,cg_verbose)
  
  ! add lphi to rhs
  do n=1,nlevs
     call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1)
  enddo
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Setup LHS coefficients
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! create lhsbeta = -hcoeff2 = (dt/2)k_{th}^{(2'')}/c_p^{(2'')}
  ! put beta on faces (remember to scale by -1 afterwards)
  do n=1,nlevs
     do i=1,lhsbeta(n)%nboxes
        if (multifab_remote(lhsbeta(n),i)) cycle
        hcoeff2p => dataptr(hcoeff2(n),i)
        lhsbetap => dataptr(lhsbeta(n),i)
        lo = lwb(get_box(lhsbeta(n), i))
        hi = upb(get_box(lhsbeta(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,hcoeff2p(:,:,1,1), &
                                     lhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,hcoeff2p(:,:,:,1), &
                                     lhsbetap(:,:,:,:))
        end select
     end do
  enddo

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

  ! Call the solver to obtain h^(2) (it will be stored in phi)
  ! solves (alpha - nabla dot beta nabla)phi = rh
  call mac_multigrid(mla,rhs,phi,fine_flx,lhsalpha,lhsbeta,dx,the_bc_tower, &
                     dm+rhoh_comp,stencil_order,mla%mba%rr, &
                     mg_verbose,cg_verbose)

  ! load new rho*h into s2
  do n=1,nlevs
     call multifab_copy_c(s2(n),rhoh_comp,phi(n),1,1)
     call multifab_mult_mult_c(s2(n),rhoh_comp,s2(n),rho_comp,1)

     call multifab_fill_boundary_c(s2(n),rhoh_comp,1)
     call multifab_physbc(s2(n),rhoh_comp,dm+rhoh_comp,1,dx(n,:), &
                          the_bc_tower%bc_tower_array(n))
  enddo

  do n=nlevs,2,-1
     call ml_cc_restriction_c(s2(n-1),rhoh_comp,s2(n),rhoh_comp,mla%mba%rr(n-1,:),1)
     call multifab_fill_ghost_cells(s2(n),s2(n-1), &
                                    ng,mla%mba%rr(n-1,:), &
                                    the_bc_tower%bc_tower_array(n-1), &
                                    the_bc_tower%bc_tower_array(n  ), &
                                    rhoh_comp,dm+rhoh_comp,1)
  enddo

  ! compute updated temperature
  call makeTfromRhoH(nlevs,s2,t02,mla,the_bc_tower%bc_tower_array,dx)

  do n = 1,nlevs
     call destroy(rhsalpha(n))
     call destroy(lhsalpha(n))
     call destroy(rhsbeta(n))
     call destroy(lhsbeta(n))
     call destroy(ccbeta(n))
     call destroy(phi(n))
     call destroy(phitemp(n))
     call destroy(Lphi(n))
     call destroy(rhs(n))
     call destroy(p01fab(n))
     call destroy(p02fab(n))
     call destroy(hcoeff1(n))
     call destroy(hcoeff2(n))
     call destroy(Xkcoeff1(n))
     call destroy(Xkcoeff2(n))
     call destroy(pcoeff1(n))
     call destroy(pcoeff2(n))
  enddo

  deallocate(rhsalpha,lhsalpha,rhsbeta,lhsbeta,ccbeta,phi,phitemp,Lphi,rhs)
  deallocate(p01fab,p02fab,fine_flx)
  deallocate(hcoeff1,hcoeff2,Xkcoeff1,Xkcoeff2,pcoeff1,pcoeff2)

end subroutine thermal_conduct_half_alg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute hcoeff and Xkcoeff, defined as:
! hcoeff = -(dt/2)k_{th}/c_p
! Xkcoeff = (dt/2)\xi_k k_{th}/c_p
! pcoeff = (dt/2)h_p*k_{th}/c_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_thermo_quantities_2d(lo,hi,dt,s,hcoeff,Xkcoeff,pcoeff)

  use variables, only: rho_comp, spec_comp, temp_comp
  use eos_module
  use probin_module, ONLY: use_big_h

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt
  real(kind=dp_t), intent(in   ) :: s(lo(1)-3:,lo(2)-3:,:)
  real(kind=dp_t), intent(inout) :: hcoeff(lo(1)-1:,lo(2)-1:)
  real(kind=dp_t), intent(inout) :: Xkcoeff(lo(1)-1:,lo(2)-1:,:)
  real(kind=dp_t), intent(inout) :: pcoeff(lo(1)-1:,lo(2)-1:)

  ! Local
  integer :: i,j,comp
  real(dp_t) :: qreact

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

        hcoeff(i,j) = -HALF*dt*conduct_eos(1)/cp_eos(1)
        pcoeff(i,j) = HALF*dt*(conduct_eos(1)/cp_eos(1))* &
             ((1.0d0/den_eos(1))* &
              (1.0d0-p_eos(1)/(den_eos(1)*dpdr_eos(1)))+dedr_eos(1)/dpdr_eos(1))

        if(use_big_h) then
           do comp=1,nspec
              Xkcoeff(i,j,comp) = HALF*dt*conduct_eos(1)* &
                   (dhdX_eos(1,comp)+ebin(comp))/cp_eos(1)
           enddo
        else
           do comp=1,nspec
              Xkcoeff(i,j,comp) = HALF*dt*conduct_eos(1)*dhdX_eos(1,comp)/cp_eos(1)
           enddo
        endif

     enddo
  enddo

end subroutine compute_thermo_quantities_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! hcoeff = -(dt/2)k_{th}/c_p
! Xkcoeff = (dt/2)\xi_k k_{th}/c_p
! pcoeff = (dt/2)h_p*k_{th}/c_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_thermo_quantities_3d(lo,hi,dt,t0,s,hcoeff,Xkcoeff,pcoeff)

  use variables, only: rho_comp, temp_comp, spec_comp
  use eos_module
  use probin_module, ONLY: use_big_h
  use geometry, only: spherical

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt
  real(kind=dp_t), intent(in   ) :: t0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-3:,lo(2)-3:,lo(3)-3:,:)
  real(kind=dp_t), intent(inout) :: hcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:)
  real(kind=dp_t), intent(inout) :: Xkcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
  real(kind=dp_t), intent(inout) :: pcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:)

  ! Local
  integer :: i,j,k,comp
  real(dp_t) :: qreact

  if(spherical .eq. 1) then
     print*, "compute_thermo1_quantities_3d spherical case not written!"
     stop
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

           hcoeff(i,j,k) = -HALF*dt*conduct_eos(1)/cp_eos(1)
           pcoeff(i,j,k) = HALF*dt*(conduct_eos(1)/cp_eos(1))* &
                ((1.0d0/den_eos(1))* &
                (1.0d0-p_eos(1)/(den_eos(1)*dpdr_eos(1)))+dedr_eos(1)/dpdr_eos(1))

           if(use_big_h) then
              do comp=1,nspec
                 Xkcoeff(i,j,k,comp) = HALF*dt*conduct_eos(1)* &
                      (dhdX_eos(1,comp)+ebin(comp))/cp_eos(1)

              enddo
           else
              do comp=1,nspec
                 Xkcoeff(i,j,k,comp) = HALF*dt*conduct_eos(1)*dhdX_eos(1,comp)/cp_eos(1)
              enddo
           endif
        enddo
     enddo
  enddo

end subroutine compute_thermo_quantities_3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! put beta on faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine put_beta_on_faces_2d(lo,hi,ccbeta,beta)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(kind=dp_t), intent(in   ) :: ccbeta(lo(1)-1:,lo(2)-1:)
  real(kind=dp_t), intent(inout) :: beta(lo(1)-1:,lo(2)-1:,:)

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
subroutine put_beta_on_faces_3d(lo,hi,ccbeta,beta)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(kind=dp_t), intent(in   ) :: ccbeta(lo(1)-1:,lo(2)-1:,lo(3)-1:)
  real(kind=dp_t), intent(inout) :: beta(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine put_base_state_on_multifab_2d(lo,hi,p0,phi)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(inout) :: phi(lo(1)-1:,lo(2)-1:)

  ! local
  integer :: i,j

  do j=lo(2),hi(2)
     do i=lo(1)-1,hi(1)+1
        phi(i,j) = p0(j)
     enddo
  enddo

end subroutine put_base_state_on_multifab_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine put_base_state_on_multifab_3d(lo,hi,p0,phi)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(inout) :: phi(lo(1)-1:,lo(2)-1:,lo(3)-1:)

  ! local
  integer :: i,j,k

  do k=lo(3),hi(3)
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           phi(i,j,k) = p0(k)
        enddo
     enddo
  enddo

end subroutine put_base_state_on_multifab_3d

end module thermal_conduct_module
