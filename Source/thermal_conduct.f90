module thermal_conduct_module

  use bl_types
  use bc_module
  use define_bc_module
  use multifab_module
  use boxarray_module
  use stencil_module
  use macproject_module
  use eos_module
  use fill_3d_module
  use probin_module

  implicit none

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Crank-Nicholson solve for enthalpy, taking into account only the
! enthalpy-diffusion terms in the temperature conduction term.
! See paper IV, steps 4a and 8a.
subroutine thermal_conduct_half_alg(mla,dx,dt,s1,s2,p01,p02,temp0, &
                                    mg_verbose,cg_verbose,the_bc_tower)

  type(ml_layout), intent(inout) :: mla
  real(dp_t)     , intent(in   ) :: dx(:,:),dt
  type(multifab) , intent(in   ) :: s1(:)
  type(multifab) , intent(inout) :: s2(:)
  real(kind=dp_t), intent(in   ) :: p01(0:),p02(0:),temp0(0:)
  integer        , intent(in   ) :: mg_verbose,cg_verbose
  type(bc_tower) , intent(in   ) :: the_bc_tower

! Local
  type(multifab), allocatable :: rhsalpha(:),lhsalpha(:),rhsbeta(:),lhsbeta(:)
  type(multifab), allocatable :: ccbeta(:),phi(:),phitemp(:)
  type(multifab), allocatable :: Lphi(:),rhs(:),kthovercp1(:)
  type(multifab), allocatable :: kthovercp2prime(:),xik1(:),xik2prime(:)
  real(kind=dp_t), pointer    :: s1p(:,:,:,:),s2p(:,:,:,:),rhsalphap(:,:,:,:)
  real(kind=dp_t), pointer    :: rhsbetap(:,:,:,:),lhsbetap(:,:,:,:)
  real(kind=dp_t), pointer    :: ccbetap(:,:,:,:),phip(:,:,:,:),rhsp(:,:,:,:)
  real(kind=dp_t), pointer    :: kthovercp1p(:,:,:,:),kthovercp2primep(:,:,:,:)
  real(kind=dp_t), pointer    :: xik1p(:,:,:,:),xik2primep(:,:,:,:)
  integer                     :: nlevs,dm,stencil_order
  integer                     :: i,n,spec
  integer                     :: lo(s1(1)%dim),hi(s1(1)%dim)
  type(bndry_reg), pointer    :: fine_flx(:) => Null()

  nlevs = mla%nlevel
  dm = mla%dim
  stencil_order = 2

  allocate(rhsalpha(nlevs),lhsalpha(nlevs))
  allocate(rhsbeta(nlevs),lhsbeta(nlevs),ccbeta(nlevs))
  allocate(phi(nlevs),phitemp(nlevs),Lphi(nlevs),rhs(nlevs))
  allocate(kthovercp1(nlevs),kthovercp2prime(nlevs))
  allocate(xik1(nlevs),xik2prime(nlevs))

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

     call multifab_build(     kthovercp1(n), mla%la(n),  1, 1)
     call multifab_build(kthovercp2prime(n), mla%la(n),  1, 1)
     call multifab_build(           xik1(n), mla%la(n),  nspec, 1)
     call multifab_build(      xik2prime(n), mla%la(n),  nspec, 1)

     call setval(rhsalpha(n), ZERO, all=.true.)
     call setval(lhsalpha(n), ZERO, all=.true.)
     call setval(rhsbeta(n),  ZERO, all=.true.)
     call setval(lhsbeta(n),  ZERO, all=.true.)
     call setval(ccbeta(n),   ZERO, all=.true.)
     call setval(Lphi(n),     ZERO, all=.true.)
     call setval(phi(n),      ZERO, all=.true.)
     call setval(phitemp(n),  ZERO, all=.true.)
     call setval(rhs(n),      ZERO, all=.true.)

     call setval(     kthovercp1(n), ZERO, all=.true.)
     call setval(kthovercp2prime(n), ZERO, all=.true.)
     call setval(           xik1(n), ZERO, all=.true.)
     call setval(      xik2prime(n), ZERO, all=.true.)
  end do

  ! lhsalpha = \rho^{(2')} = \rho^{(2)}
  ! rhsalpha = 0 (already initialized above)
  ! thess will be true for this entire subroutine
  do n=1,nlevs
     call multifab_copy_c(lhsalpha(n),1,s2(n),rho_comp,1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!
  ! First implicit solve
  !!!!!!!!!!!!!!!!!!!!!!!

  ! compute kthovercp1 and xik1, defined as:
  ! kthovercp1 = -(dt/2)k_{th}^{(1)}/c_p^{(1)}
  ! xik1 = (dt/2)\xi_k^{(1)}k_{th}^{(1)}/c_p^{(1)}
  do n=1,nlevs
     do i=1,s1(n)%nboxes
        if (multifab_remote(s1(n),i)) cycle
        s1p         => dataptr(s1(n),i)
        kthovercp1p => dataptr(kthovercp1(n),i)
        xik1p       => dataptr(xik1(n),i)
        lo = lwb(get_box(s1(n), i))
        hi = upb(get_box(s1(n), i))
        select case (dm)
        case (2)
           call compute_thermo_quantities_2d(lo,hi,dt,temp0, &
                                             s1p(:,:,1,:), &
                                             kthovercp1p(:,:,1,1), &
                                             xik1p(:,:,1,:))
        case (3)
           call compute_thermo_quantities_3d(lo,hi,dt,temp0, &
                                             s1p(:,:,:,:), &
                                             kthovercp1p(:,:,:,1), &
                                             xik1p(:,:,:,:))
        end select
     end do
  enddo

  ! begin construction of rhs by setting rhs = (\rho h)^{(2)}
  do n=1,nlevs
     call multifab_copy_c(rhs(n),1,s2(n),rhoh_comp,1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add enthalpy diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! put beta on faces
  do n=1,nlevs
     do i=1,s1(n)%nboxes
        if (multifab_remote(s1(n),i)) cycle
        kthovercp1p => dataptr(kthovercp1(n),i)
        rhsbetap    => dataptr(rhsbeta(n),i)
        lo = lwb(get_box(s1(n), i))
        hi = upb(get_box(s1(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,kthovercp1p(:,:,1,1), &
                                     rhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,kthovercp1p(:,:,:,1), &
                                     rhsbetap(:,:,:,:))
        end select
     end do
  enddo

  ! load phi = h^{(1)}
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,s1(n),rhoh_comp,1)
     call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,.true.)
  enddo

  ! apply the operator
  call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                   dm+rhoh_comp,stencil_order,mla%mba%rr, &
                   mg_verbose,cg_verbose)

  ! add Lphi to rhs
  do n=1,nlevs
     call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,.true.)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add species diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! loop over species
  do spec=1,nspec

     ! put beta on faces
     do n=1,nlevs
        do i=1,s1(n)%nboxes
           if (multifab_remote(s1(n),i)) cycle
           xik1p    => dataptr(xik1(n),i)
           rhsbetap => dataptr(rhsbeta(n),i)
           lo = lwb(get_box(s1(n), i))
           hi = upb(get_box(s1(n), i))
           select case (dm)
           case (2)
              call put_beta_on_faces_2d(lo,hi,xik1p(:,:,1,spec), &
                                        rhsbetap(:,:,1,:))
           case (3)
              call put_beta_on_faces_3d(lo,hi,xik1p(:,:,:,spec), &
                                        rhsbetap(:,:,:,:))
           end select
        end do
     enddo

     ! load phi = X_k^{(1)} + X_k^{(2)}
     do n=1,nlevs
        call multifab_copy_c(phi(n),1,s1(n),spec_comp+spec-1,1)
        call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,.true.)
        call multifab_copy_c(phitemp(n),1,s2(n),spec_comp+spec-1,1)
        call multifab_div_div_c(phitemp(n),1,s2(n),rho_comp,1,.true.)
        call multifab_plus_plus_c(phi(n),1,phitemp(n),1,1)
     enddo

     ! apply the operator
     call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                      dm+spec_comp+spec-1,stencil_order,mla%mba%rr, &
                      mg_verbose,cg_verbose)
     
     ! add lphi to rhs
     do n=1,nlevs
        call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,.true.)
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now do the implicit solve
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! create lhsbeta = -kthovercp1 = (dt/2)k_{th}^{(1)}/c_p^{(1)}
  ! put beta on faces (remember to scale by -1 afterwards)
  do n=1,nlevs
     do i=1,s1(n)%nboxes
        if (multifab_remote(s1(n),i)) cycle
        kthovercp1p => dataptr(kthovercp1(n),i)
        lhsbetap    => dataptr(lhsbeta(n),i)
        lo = lwb(get_box(s1(n), i))
        hi = upb(get_box(s1(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,kthovercp1p(:,:,1,1), &
                                     lhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,kthovercp1p(:,:,:,1), &
                                     lhsbetap(:,:,:,:))
        end select
     end do
  enddo

  ! scale by -1
  do n=1,nlevs
     call multifab_mult_mult_s_c(lhsbeta(n),1,-1.0d0,dm,.true.)
  enddo

  ! Call the solver to obtain h^(2') (it will be stored in phi)
  ! solves (alpha - nabla dot beta nabla)phi = rh
  call mac_multigrid(mla,rhs,phi,fine_flx,lhsalpha,lhsbeta,dx,the_bc_tower, &
                     dm+rhoh_comp,stencil_order,mla%mba%rr, &
                     mg_verbose,cg_verbose)

  ! load new h into s2
  do n=1,nlevs
     call multifab_copy_c(s2(n),rhoh_comp,phi(n),1,1,.false.)
     call multifab_mult_mult_c(s2(n),rhoh_comp,s2(n),rho_comp,1,.true.)
  enddo

  ! fill in ghost cells on s2
  do n=1,nlevs
     call multifab_fill_boundary(s2(n))

     do i = 1, s2(n)%nboxes
        if ( multifab_remote(s2(n), i) ) cycle
        s2p => dataptr(s2(n),i)
        lo =  lwb(get_box(s2(n), i))
        hi =  upb(get_box(s2(n), i))
        select case (dm)
        case (2)
           call setbc_2d(s2p(:,:,1,rhoh_comp), lo, 3, &
       the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,dm+rhoh_comp), &
       dx(n,:),dm+rhoh_comp)
        case (3)
           call setbc_3d(s2p(:,:,:,rhoh_comp), lo, 3, &
       the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,dm+rhoh_comp), &
       dx(n,:),dm+rhoh_comp)
        end select
     end do
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!
  ! Second implicit solve
  !!!!!!!!!!!!!!!!!!!!!!!

  ! compute kthovercp2prime and xik2prime, defined as:
  ! kthovercp2prime = -(dt/2)k_{th}^{(2')}/c_p^{(2')}
  ! xik2prime = (dt/2)\xi_k^{(2')}k_{th}^{(2')}/c_p^{(2')}
  do n=1,nlevs
     do i=1,s1(n)%nboxes
        if (multifab_remote(s2(n),i)) cycle
        s2p         => dataptr(s2(n),i)
        kthovercp2primep => dataptr(kthovercp2prime(n),i)
        xik2primep       => dataptr(xik2prime(n),i)
        lo = lwb(get_box(s2(n), i))
        hi = upb(get_box(s2(n), i))
        select case (dm)
        case (2)
           call compute_thermo_quantities_2d(lo,hi,dt,temp0, &
                                             s2p(:,:,1,:), &
                                             kthovercp2primep(:,:,1,1), &
                                             xik2primep(:,:,1,:))
        case (3)
           call compute_thermo_quantities_3d(lo,hi,dt,temp0, &
                                             s2p(:,:,:,:), &
                                             kthovercp2primep(:,:,:,1), &
                                             xik2primep(:,:,:,:))
        end select
     end do
  enddo

  ! begin construction of rhs by setting rhs = (\rho h)^{(2)}
  do n=1,nlevs
     call multifab_copy_c(rhs(n),1,s2(n),rhoh_comp,1)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add enthalpy diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! put beta on faces
  do n=1,nlevs
     do i=1,s1(n)%nboxes
        if (multifab_remote(s1(n),i)) cycle
        kthovercp1p => dataptr(kthovercp1(n),i)
        rhsbetap    => dataptr(rhsbeta(n),i)
        lo = lwb(get_box(s1(n), i))
        hi = upb(get_box(s1(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,kthovercp1p(:,:,1,1), &
                                     rhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,kthovercp1p(:,:,:,1), &
                                     rhsbetap(:,:,:,:))
        end select
     end do
  enddo

  ! load phi = h^{(1)}
  do n=1,nlevs
     call multifab_copy_c(phi(n),1,s1(n),rhoh_comp,1)
     call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,.true.)
  enddo

  ! apply the operator
  call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                   dm+rhoh_comp,stencil_order,mla%mba%rr, &
                   mg_verbose,cg_verbose)

  ! add Lphi to rhs
  do n=1,nlevs
     call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,.true.)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! add species diffusion to rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! loop over species
  do spec=1,nspec

     ! do X_k^{(1)} term first

     ! put beta on faces
     do n=1,nlevs
        do i=1,s1(n)%nboxes
           if (multifab_remote(s1(n),i)) cycle
           xik1p    => dataptr(xik1(n),i)
           rhsbetap => dataptr(rhsbeta(n),i)
           lo = lwb(get_box(s1(n), i))
           hi = upb(get_box(s1(n), i))
           select case (dm)
           case (2)
              call put_beta_on_faces_2d(lo,hi,xik1p(:,:,1,spec), &
                                        rhsbetap(:,:,1,:))
           case (3)
              call put_beta_on_faces_3d(lo,hi,xik1p(:,:,:,spec), &
                                        rhsbetap(:,:,:,:))
           end select
        end do
     enddo

     ! load phi = X_k^{(1)}
     do n=1,nlevs
        call multifab_copy_c(phi(n),1,s1(n),spec_comp+spec-1,1)
        call multifab_div_div_c(phi(n),1,s1(n),rho_comp,1,.true.)
     enddo

     ! apply the operator
     call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                      dm+spec_comp+spec-1,stencil_order,mla%mba%rr, &
                      mg_verbose,cg_verbose)
     
     ! add lphi to rhs
     do n=1,nlevs
        call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,.true.)
     enddo

     ! now do X_k^{(2)} term

     ! put beta on faces
     do n=1,nlevs
        do i=1,s1(n)%nboxes
           if (multifab_remote(s1(n),i)) cycle
           xik2primep    => dataptr(xik2prime(n),i)
           rhsbetap => dataptr(rhsbeta(n),i)
           lo = lwb(get_box(s1(n), i))
           hi = upb(get_box(s1(n), i))
           select case (dm)
           case (2)
              call put_beta_on_faces_2d(lo,hi,xik2primep(:,:,1,spec), &
                                        rhsbetap(:,:,1,:))
           case (3)
              call put_beta_on_faces_3d(lo,hi,xik2primep(:,:,:,spec), &
                                        rhsbetap(:,:,:,:))
           end select
        end do
     enddo

     ! load phi = X_k^{(2)}
     do n=1,nlevs
        call multifab_copy_c(phi(n),1,s2(n),spec_comp+spec-1,1)
        call multifab_div_div_c(phi(n),1,s2(n),rho_comp,1,.true.)
     enddo

     ! apply the operator
     call mac_applyop(mla,Lphi,phi,rhsalpha,rhsbeta,dx,the_bc_tower, &
                      dm+spec_comp+spec-1,stencil_order,mla%mba%rr, &
                      mg_verbose,cg_verbose)
     
     ! add lphi to rhs
     do n=1,nlevs
        call multifab_plus_plus_c(rhs(n),1,Lphi(n),1,1,.true.)
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now do the implicit solve
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! create lhsbeta = -kthovercp1 = (dt/2)k_{th}^{(2')}/c_p^{(2')}
  ! put beta on faces (remember to scale by -1 afterwards)
  do n=1,nlevs
     do i=1,s1(n)%nboxes
        if (multifab_remote(s1(n),i)) cycle
        kthovercp2primep => dataptr(kthovercp2prime(n),i)
        lhsbetap    => dataptr(lhsbeta(n),i)
        lo = lwb(get_box(s1(n), i))
        hi = upb(get_box(s1(n), i))
        select case (dm)
        case (2)
           call put_beta_on_faces_2d(lo,hi,kthovercp2primep(:,:,1,1), &
                                     lhsbetap(:,:,1,:))
        case (3)
           call put_beta_on_faces_3d(lo,hi,kthovercp2primep(:,:,:,1), &
                                     lhsbetap(:,:,:,:))
        end select
     end do
  enddo

  ! scale by -1
  do n=1,nlevs
     call multifab_mult_mult_s_c(lhsbeta(n),1,-1.0d0,dm,.true.)
  enddo

  ! Call the solver to obtain h^(3) (it will be stored in phi)
  ! solves (alpha - nabla dot beta nabla)phi = rh
  call mac_multigrid(mla,rhs,phi,fine_flx,lhsalpha,lhsbeta,dx,the_bc_tower, &
                     dm+rhoh_comp,stencil_order,mla%mba%rr, &
                     mg_verbose,cg_verbose)

  ! load new h into s2
  do n=1,nlevs
     call multifab_copy_c(s2(n),rhoh_comp,phi(n),1,1,.false.)
     call multifab_mult_mult_c(s2(n),rhoh_comp,s2(n),rho_comp,1,.true.)
  enddo

  ! fill in ghost cells on s2
  do n=1,nlevs
     call multifab_fill_boundary(s2(n))

     do i = 1, s2(n)%nboxes
        if ( multifab_remote(s2(n), i) ) cycle
        s2p => dataptr(s2(n),i)
        lo =  lwb(get_box(s2(n), i))
        hi =  upb(get_box(s2(n), i))
        select case (dm)
        case (2)
           call setbc_2d(s2p(:,:,1,rhoh_comp), lo, 3, &
       the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,dm+rhoh_comp), &
       dx(n,:),dm+rhoh_comp)
        case (3)
           call setbc_3d(s2p(:,:,:,rhoh_comp), lo, 3, &
       the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,dm+rhoh_comp), &
       dx(n,:),dm+rhoh_comp)
        end select
     end do
  enddo

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
     call destroy(kthovercp1(n))
     call destroy(kthovercp2prime(n))
     call destroy(xik1(n))
     call destroy(xik2prime(n))
  enddo

  deallocate(rhsalpha,lhsalpha,rhsbeta,lhsbeta,ccbeta,phi,phitemp,Lphi,rhs)
  deallocate(kthovercp1,kthovercp2prime,xik1,xik2prime)

end subroutine thermal_conduct_half_alg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute kthovercp and xik, defined as:
! kthovercp = -(dt/2)k_{th}/c_p
! xik = (dt/2)\xi_k k_{th}/c_p
subroutine compute_thermo_quantities_2d(lo,hi,dt,temp0,s,kthovercp,xik)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt
  real(kind=dp_t), intent(in   ) :: temp0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-3:,lo(2)-3:,:)
  real(kind=dp_t), intent(inout) :: kthovercp(lo(1)-1:,lo(2)-1:)
  real(kind=dp_t), intent(inout) :: xik(lo(1)-1:,lo(2)-1:,:)

  ! Local
  integer :: i,j,n
  real(dp_t) :: qreact

  ! density, enthalpy, and xmass are inputs with initial temperature guess
  input_flag = 2
  do_diag = .false.

  do j=lo(2)-1,hi(2)+1
     do i=lo(1)-1,hi(1)+1

        den_row(1) = s(i,j,rho_comp)
        xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

        qreact = 0.0d0
        if(use_big_h) then
           do n=1,nspec
              qreact = qreact - ebin(n)*xn_zone(n)
           enddo
           h_row(1) = s(i,j,rhoh_comp)/den_row(1) - qreact
        else
           h_row(1) = s(i,j,rhoh_comp)/den_row(1)
        endif

        if(j .lt. lo(2)) then
           temp_row(1) = temp0(lo(2))
        else if(j .gt. hi(2)) then
           temp_row(1) = temp0(hi(2))
        else
           temp_row(1) = temp0(j)
        endif

        call conducteos(input_flag, den_row, temp_row, &
                 npts, nspec, &
                 xn_zone, aion, zion, &
                 p_row, h_row, e_row, & 
                 cv_row, cp_row, xne_row, eta_row, pele_row, &
                 dpdt_row, dpdr_row, dedt_row, dedr_row, &
                 dpdX_row, dhdX_row, &
                 gam1_row, cs_row, s_row, &
                 dsdt_row, dsdr_row, &
                 do_diag, conduct_row)

        kthovercp(i,j) = -HALF*dt*conduct_row(1)/cp_row(1)

        if(use_big_h) then
           do n=1,nspec
              xik(i,j,n) = -(dhdX_row(1,n)-ebin(n))*kthovercp(i,j)
           enddo
        else
           do n=1,nspec
              xik(i,j,n) = -dhdX_row(1,n)*kthovercp(i,j)
           enddo
        endif
     enddo
  enddo

end subroutine compute_thermo_quantities_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute kthovercp1 and xik1, defined as:
! kthovercp1 = -(dt/2)k_{th}^{(1)}/c_p^{(1)}
! xik1 = (dt/2)\xi_k^{(1)}k_{th}^{(1)}/c_p^{(1)}
subroutine compute_thermo_quantities_3d(lo,hi,dt,temp0,s,kthovercp,xik)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt
  real(kind=dp_t), intent(in   ) :: temp0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-3:,lo(2)-3:,lo(3)-3:,:)
  real(kind=dp_t), intent(inout) :: kthovercp(lo(1)-1:,lo(2)-1:,lo(3)-1:)
  real(kind=dp_t), intent(inout) :: xik(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)

  ! Local
  integer :: i,j,k,n
  real(dp_t) :: qreact

  if(spherical .eq. 1) then
     print*, "compute_thermo1_quantities_3d spherical case not written!"
     stop
  endif

  ! density, enthalpy, and xmass are inputs with initial temperature guess
  input_flag = 2
  do_diag = .false.

  do k=lo(3)-1,hi(3)+1
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1

           den_row(1) = s(i,j,k,rho_comp)
           xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

           qreact = 0.0d0
           if(use_big_h) then
              do n=1,nspec
                 qreact = qreact - ebin(n)*xn_zone(n)
              enddo
              h_row(1) = s(i,j,k,rhoh_comp)/den_row(1) - qreact
           else
              h_row(1) = s(i,j,k,rhoh_comp)/den_row(1)
           endif

           if(j .lt. lo(3)) then
              temp_row(1) = temp0(lo(3))
           else if(j .gt. hi(3)) then
              temp_row(1) = temp0(hi(3))
           else
              temp_row(1) = temp0(k)
           endif

           call conducteos(input_flag, den_row, temp_row, &
                           npts, nspec, &
                           xn_zone, aion, zion, &
                           p_row, h_row, e_row, & 
                           cv_row, cp_row, xne_row, eta_row, pele_row, &
                           dpdt_row, dpdr_row, dedt_row, dedr_row, &
                           dpdX_row, dhdX_row, &
                           gam1_row, cs_row, s_row, &
                           dsdt_row, dsdr_row, &
                           do_diag, conduct_row)

           kthovercp(i,j,k) = -HALF*dt*conduct_row(1)/cp_row(1)

           if(use_big_h) then
              xik(i,j,k,n) = -(dhdX_row(1,n)-ebin(n))*kthovercp(i,j,k)
           else
              do n=1,nspec
                 xik(i,j,k,n) = -dhdX_row(1,n)*kthovercp(i,j,k)
              enddo
           endif
        enddo
     enddo
  enddo

end subroutine compute_thermo_quantities_3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! put beta on faces
subroutine put_beta_on_faces_2d(lo,hi,ccbeta,beta)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(kind=dp_t), intent(in   ) :: ccbeta(lo(1)-1:,lo(2)-1:)
  real(kind=dp_t), intent(inout) :: beta(lo(1)-1:,lo(2)-1:,:)

! Local
  integer :: i,j
  integer :: nx,ny

  nx = size(beta,dim=1) - 2
  ny = size(beta,dim=2) - 2

  do j = 0,ny-1
     do i = 0,nx
        beta(i,j,1) = TWO*(ccbeta(i,j)*ccbeta(i-1,j))/(ccbeta(i,j) + ccbeta(i-1,j))
     end do
  end do
  
  do j = 0,ny
     do i = 0,nx-1
        beta(i,j,2) = TWO*(ccbeta(i,j)*ccbeta(i,j-1))/(ccbeta(i,j) + ccbeta(i,j-1))
     end do
  end do

end subroutine put_beta_on_faces_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! put beta on faces
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

  do k = 0,nz-1
     do j = 0,ny-1
        do i = 0,nx
           beta(i,j,k,1) = (ccbeta(i,j,k) + ccbeta(i-1,j,k)) / TWO
        end do
     end do
  end do
  
  do k = 0,nz-1
     do j = 0,ny
        do i = 0,nx-1
           beta(i,j,k,2) = (ccbeta(i,j,k) + ccbeta(i,j-1,k)) / TWO
        end do
     end do
  end do
  
  do k = 0,nz
     do j = 0,ny-1
        do i = 0,nx-1
           beta(i,j,k,3) = (ccbeta(i,j,k) + ccbeta(i,j,k-1)) / TWO
        end do
     end do
  end do

end subroutine put_beta_on_faces_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output stuff in a multifab
subroutine fab_access_2d(lo,hi,fab,ng)

 integer        , intent(in   ) :: lo(:),hi(:),ng
 real(kind=dp_t), intent(in   ) :: fab(lo(1)-ng:,lo(2)-ng:,:)

! Local
  integer :: i,j



end subroutine fab_access_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output stuff in a multifab
subroutine fab_access_3d(lo,hi,fab,ng)

 integer        , intent(in   ) :: lo(:),hi(:),ng
 real(kind=dp_t), intent(in   ) :: fab(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

! Local
  integer :: i,j,k



end subroutine fab_access_3d

end module thermal_conduct_module
