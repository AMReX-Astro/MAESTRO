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

  implicit none

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Crank-Nicholson solve for enthalpy, taking into account only the
! enthalpy-diffusion terms in the temperature conduction term.
! See paper IV, steps 4a and 8a.
subroutine thermal_conduct(mla,dx,dt,sold,s2,p0old,p02, &
                           mg_verbose,cg_verbose,the_bc_tower)

  type(ml_layout), intent(inout) :: mla
  real(dp_t)     , intent(in   ) :: dx(:,:),dt
  type(multifab) , intent(in   ) :: sold(:)
  type(multifab) , intent(inout) :: s2(:)
  real(kind=dp_t), intent(in   ) :: p0old(0:),p02(0:)
  integer        , intent(in   ) :: mg_verbose,cg_verbose
  type(bc_tower) , intent(in   ) :: the_bc_tower

! Local
  type(multifab), allocatable :: rh(:),phi(:),alpha(:),beta(:),rhsbeta(:)
  type(multifab), allocatable :: const(:)
  integer                     :: i,n,nlevs,dm,ng,ng_rh,ng_s,stencil_order
  integer                     :: lo(sold(1)%dim),hi(sold(1)%dim)
  real(kind=dp_t), pointer    :: soldp(:,:,:,:),s2p(:,:,:,:)
  real(kind=dp_t), pointer    :: rhp(:,:,:,:),phip(:,:,:,:)
  real(kind=dp_t), pointer    :: betap(:,:,:,:),rhsbetap(:,:,:,:)
  real(kind=dp_t), pointer    :: constp(:,:,:,:)
  type(bndry_reg), pointer    :: fine_flx(:) => Null()

  if (parallel_IOProcessor()) print *,'... Entering thermal_conduct ...'

  nlevs = mla%nlevel
  dm    = mla%dim
  stencil_order = 2

  allocate(rh(nlevs),phi(nlevs),alpha(nlevs),beta(nlevs))
  allocate(rhsbeta(nlevs),const(nlevs))

  do n = 1,nlevs
     call multifab_build(     rh(n), mla%la(n), 1, 0)
     call multifab_build(    phi(n), mla%la(n), 1, 1)
     call multifab_build(  alpha(n), mla%la(n), 1, 1)
     call multifab_build(   beta(n), mla%la(n),dm, 1)
     call multifab_build(rhsbeta(n), mla%la(n),dm, 1)
     call multifab_build(   const(n), mla%la(n), 1, 1)  

     call setval(rh(n),ZERO,all=.true.)
     call setval(phi(n),ZERO,all=.true.)
     call setval(alpha(n),ZERO,all=.true.)
     call setval(beta(n),ZERO,all=.true.)
     call setval(rhsbeta(n),ZERO,all=.true.)
     call setval(const(n),ZERO,all=.true.)
  end do

  if (parallel_IOProcessor()) print *,'... Setting alpha = rho ...'
  if (parallel_IOProcessor()) print *,'... Computing betas and phi ...'

  do n=1,nlevs
     ng    = alpha(n)%ng
     ng_rh = rh(n)%ng
     ng_s  = sold(n)%ng
     
     ! Create beta = \frac{\Delta t k_th^(2)}{2 c_p^(2)}
     ! Create rhsbeta = dt*k_th^n/(2*c_p^n)
     ! Copy h^n into phi
     do i=1,sold(n)%nboxes
        if (multifab_remote(sold(n),i)) cycle
        soldp    => dataptr(sold(n),i)
        s2p      => dataptr(s2(n),i)
        betap    => dataptr(beta(n),i)
        rhsbetap => dataptr(rhsbeta(n),i)
        phip     => dataptr(phi(n),i)
        constp    => dataptr(const(n),i)
        lo =  lwb(get_box(sold(n), i))
        hi =  upb(get_box(sold(n), i))
        select case (dm)
        case (2)
           call make_betas_and_phi_2d(lo,hi,dt,dx(n,:),ng,ng_rh,ng_s, &
                                      p0old,p02,soldp(:,:,1,:),s2p(:,:,1,:), &
                                      betap(:,:,1,:),rhsbetap(:,:,1,:), &
                                      phip(:,:,1,1),constp(:,:,1,1))
        case (3)
           call make_betas_and_phi_3d(lo,hi,dt,dx(n,:),ng,ng_rh,ng_s, &
                                      p0old,p02,soldp(:,:,:,:),s2p(:,:,:,:), &
                                      betap(:,:,:,:),rhsbetap(:,:,:,:), &
                                      phip(:,:,:,1),constp(:,:,:,1))
        end select
     end do
  enddo

!  call fabio_ml_multifab_write_d(phi,mla%mba%rr(:,1),"h^n")
!  call fabio_ml_multifab_write_d(alpha,mla%mba%rr(:,1),"rhsalpha")
!  call fabio_ml_multifab_write_d(rhsbeta,mla%mba%rr(:,1),"rhsbeta")

  if (parallel_IOProcessor()) print *,'... Computing RHS operator residual ...'
  ! residual = -(alpha-nabla dot rhsbeta nabla)phi
  ! rh = -residual
  call mac_applyop(mla,rh,phi,alpha,rhsbeta,dx,the_bc_tower,dm+rhoh_comp, &
                   stencil_order,mla%mba%rr,mg_verbose,cg_verbose)

!  call fabio_ml_multifab_write_d(rh,mla%mba%rr(:,1),"resid")

  if (parallel_IOProcessor()) print *,'... Adding rho h to RHS ...'
  ! add (\rho h)^(2) to RHS
  ! set alpha to rho^(2)
  do n=1,nlevs
     ng    = alpha(n)%ng
     ng_rh = rh(n)%ng
     ng_s  = sold(n)%ng

     ! Copy rho^(2) directly into alpha
     call multifab_copy_c(alpha(n),1,s2(n),rho_comp,1)

     do i=1,sold(n)%nboxes
        if (multifab_remote(sold(n),i)) cycle
        s2p   => dataptr(s2(n),i)
        rhp   => dataptr(rh(n),i)
        lo =  lwb(get_box(sold(n), i))
        hi =  upb(get_box(sold(n), i))
        select case (dm)
        case (2)
           call add_rhoh_to_rh_2d(lo,hi,ng_rh,ng_s,s2p(:,:,1,:),rhp(:,:,1,1))
        case (3)
           call add_rhoh_to_rh_3d(lo,hi,ng_rh,ng_s,s2p(:,:,:,:),rhp(:,:,:,1))
        end select
     end do
  enddo

  if (parallel_IOProcessor()) print *,'... Calling solver ...'

  allocate(fine_flx(2:nlevs))
  do n = 2,nlevs
     call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
  end do

  ! Make initial guess at phi equal to h^(2)
  do n=1,nlevs
     ng    = alpha(n)%ng
     ng_rh = rh(n)%ng
     ng_s  = sold(n)%ng
     
     ! Create beta = \frac{\Delta t k_th^(2)}{2 c_p^(2)}
     ! Create rhsbeta = dt*k_th^n/(2*c_p^n)
     ! Copy h^n into phi
     do i=1,sold(n)%nboxes
        if (multifab_remote(sold(n),i)) cycle
        s2p      => dataptr(s2(n),i)
        phip     => dataptr(phi(n),i)
        lo =  lwb(get_box(sold(n), i))
        hi =  upb(get_box(sold(n), i))
        select case (dm)
        case (2)
           call make_initial_phi_2d(lo,hi,ng,ng_s,s2p(:,:,1,:),phip(:,:,1,1))
        case (3)
           call make_initial_phi_3d(lo,hi,ng,ng_s,s2p(:,:,:,:),phip(:,:,:,1))
        end select
     end do
  enddo

  ! Call the solver to obtain h^(2') (it will be stored in phi)
  ! solves (alpha - nabla dot beta nabla)phi = rh
  call mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx,the_bc_tower, &
                     dm+rhoh_comp,stencil_order,mla%mba%rr, &
                     mg_verbose,cg_verbose)

  ! Compute updated (\rho h) = \rho^(2)h^(2')
  do n=1,nlevs

     ng    = alpha(n)%ng
     ng_rh = rh(n)%ng
     ng_s  = sold(n)%ng

     do i=1,sold(n)%nboxes
        if (multifab_remote(sold(n),i)) cycle
        phip    => dataptr(phi(n),i)
        s2p     => dataptr(s2(n),i)
        lo =  lwb(get_box(sold(n), i))
        hi =  upb(get_box(sold(n), i))
        select case (dm)
        case (2)
           call compute_rhoh_2d(lo,hi,ng,ng_s,phip(:,:,1,1),s2p(:,:,1,:))
        case (3)
           call compute_rhoh_3d(lo,hi,ng,ng_s,phip(:,:,:,1),s2p(:,:,:,:))
        end select
     end do
  enddo

  ! Deallocate memory
  do n = 1,nlevs
     call destroy(rh(n))
     call destroy(phi(n))
     call destroy(alpha(n))
     call destroy(beta(n))
     call destroy(rhsbeta(n))
     call destroy(const(n))
  enddo

  deallocate(rh,phi,alpha,beta,rhsbeta,const)

end subroutine thermal_conduct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute betas and phi for 2d problems
subroutine make_betas_and_phi_2d(lo,hi,dt,dx,ng,ng_rh,ng_s, &
                                 p0old,p02,sold,s2,beta,rhsbeta,phi,const)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt,dx(:)
  integer        , intent(in   ) :: ng,ng_rh,ng_s
  real(kind=dp_t), intent(in   ) :: p0old(0:),p02(0:)
  real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_s:,lo(2)-ng_s:,:)
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,:)
  real(kind=dp_t), intent(  out) :: beta(lo(1)-ng:,lo(2)-ng:,:)
  real(kind=dp_t), intent(  out) :: rhsbeta(lo(1)-ng:,lo(2)-ng:,:)
  real(kind=dp_t), intent(  out) :: phi(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(inout) :: const(lo(1)-ng:,lo(2)-ng:)

! Local
  integer :: i,j
  integer :: nx,ny
  real(kind=dp_t) :: conductivity
  
! dens, pres, and xmass are inputs
  input_flag = 4
  do_diag = .false.

  nx = size(beta,dim=1) - 2*ng
  ny = size(beta,dim=2) - 2*ng

  ! Compute c_p^(2), k_th^2, and beta
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)

        den_row(1) = s2(i,j,rho_comp)
        p_row(1) = p02(j)
        xn_zone(:) = s2(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

        call conducteos(input_flag, den_row, temp_row, &
                 npts, nspec, &
                 xn_zone, aion, zion, &
                 p_row, h_row, e_row, & 
                 cv_row, cp_row, xne_row, eta_row, pele_row, &
                 dpdt_row, dpdr_row, dedt_row, dedr_row, &
                 dpdX_row, dhdX_row, &
                 gam1_row, cs_row, s_row, &
                 dsdt_row, dsdr_row, &
                 do_diag, conductivity)

        const(i,j) = HALF*dt*conductivity/cp_row(1)

     enddo
  enddo

  do j = 0,ny-1
     do i = 0,nx
        beta(i,j,1) = (const(i,j) + const(i-1,j)) / TWO
     end do
  end do
  
  do j = 0,ny
     do i = 0,nx-1
        beta(i,j,2) = (const(i,j) + const(i,j-1)) / TWO
     end do
  end do

 ! Compute c_p^n, k_th^n, and rhsbeta
    do j=lo(2),hi(2)
     do i=lo(1),hi(1)

        den_row(1) = sold(i,j,rho_comp)
        p_row(1) = p0old(j)
        xn_zone(:) = sold(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

        call conducteos(input_flag, den_row, temp_row, &
                 npts, nspec, &
                 xn_zone, aion, zion, &
                 p_row, h_row, e_row, & 
                 cv_row, cp_row, xne_row, eta_row, pele_row, &
                 dpdt_row, dpdr_row, dedt_row, dedr_row, &
                 dpdX_row, dhdX_row, &
                 gam1_row, cs_row, s_row, &
                 dsdt_row, dsdr_row, &
                 do_diag, conductivity)

        const(i,j) = -HALF*dt*conductivity/cp_row(1)

     enddo
  enddo

  do j = 0,ny-1
     do i = 0,nx
        rhsbeta(i,j,1) = (const(i,j) + const(i-1,j)) / TWO
     end do
  end do
  
  do j = 0,ny
     do i = 0,nx-1
        rhsbeta(i,j,2) = (const(i,j) + const(i,j-1)) / TWO
     end do
  end do

  ! set phi = h^n for applyop on RHS
    do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        phi(i,j) = sold(i,j,rhoh_comp)/sold(i,j,rho_comp)
     enddo
  enddo

end subroutine make_betas_and_phi_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute betas and phi for 3d problems
subroutine make_betas_and_phi_3d(lo,hi,dt,dx,ng,ng_rh,ng_s, &
                                 p0old,p02,sold,s2,beta,rhsbeta,phi,const)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt,dx(:)
  integer        , intent(in   ) :: ng,ng_rh,ng_s
  real(kind=dp_t), intent(in   ) :: p0old(0:),p02(0:)
  real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
  real(kind=dp_t), intent(  out) :: beta(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
  real(kind=dp_t), intent(  out) :: rhsbeta(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
  real(kind=dp_t), intent(  out) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(inout) :: const(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

! Local
  integer :: i,j,k
  integer :: nx,ny,nz
  real(kind=dp_t), allocatable :: p0_cart(:,:,:)
  real(kind=dp_t)              :: conductivity

  if (spherical .eq. 1) then
     allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
     call fill_3d_data(p0_cart,p02,lo,hi,dx,0)
  end if

! dens, pres, and xmass are inputs
  input_flag = 4
  do_diag = .false.

      nx = size(beta,dim=1) - 2*ng
      ny = size(beta,dim=2) - 2*ng
      nz = size(beta,dim=3) - 2*ng

  ! Compute c_p^(2), k_th^2, and beta
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           
           den_row(1) = s2(i,j,k,rho_comp)
           xn_zone(:) = s2(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

           if(spherical .eq. 0) then
              p_row(1) = p02(k)
           else
              p_row(1) = p0_cart(i,j,k)
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
                    do_diag, conductivity)

           const(i,j,k) = HALF*dt*conductivity/cp_row(1)

        enddo
     enddo
  enddo

  do k = 0,nz-1
     do j = 0,ny-1
        do i = 0,nx
           beta(i,j,k,1) = (const(i,j,k) + const(i-1,j,k)) / TWO
        end do
     end do
  end do
  
  do k = 0,nz-1
     do j = 0,ny
        do i = 0,nx-1
           beta(i,j,k,2) = (const(i,j,k) + const(i,j-1,k)) / TWO
        end do
     end do
  end do
  
  do k = 0,nz
     do j = 0,ny-1
        do i = 0,nx-1
           beta(i,j,k,3) = (const(i,j,k) + const(i,j,k-1)) / TWO
        end do
     end do
  end do
  
  if (spherical .eq. 1) then
     call fill_3d_data(p0_cart,p0old,lo,hi,dx,0)
  end if

 ! Compute c_p^n, k_th^n, and rhsbeta
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           
           den_row(1) = sold(i,j,k,rho_comp)
           xn_zone(:) = sold(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)
           
           if(spherical .eq. 0) then
              p_row(1) = p0old(k)
           else
              p_row(1) = p0_cart(i,j,k)
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
                    do_diag, conductivity)
           
           const(i,j,k) = -HALF*dt*conductivity/cp_row(1)

        enddo
     enddo
  enddo

  do k = 0,nz-1
     do j = 0,ny-1
        do i = 0,nx
           rhsbeta(i,j,k,1) = (const(i,j,k) + const(i-1,j,k)) / TWO
        end do
     end do
  end do
  
  do k = 0,nz-1
     do j = 0,ny
        do i = 0,nx-1
           rhsbeta(i,j,k,2) = (const(i,j,k) + const(i,j-1,k)) / TWO
        end do
     end do
  end do
  
  do k = 0,nz
     do j = 0,ny-1
        do i = 0,nx-1
           rhsbeta(i,j,k,3) = (const(i,j,k) + const(i,j,k-1)) / TWO
        end do
     end do
  end do

  ! set phi = h^n for applyop on RHS
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           phi(i,j,k) = sold(i,j,k,rhoh_comp)/sold(i,j,k,rho_comp)
        enddo
     enddo
  enddo

end subroutine make_betas_and_phi_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute betas and phi for 2d problems
subroutine make_initial_phi_2d(lo,hi,ng,ng_s,s2,phi)

  integer        , intent(in   ) :: lo(:),hi(:)
  integer        , intent(in   ) :: ng,ng_s
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,:)
  real(kind=dp_t), intent(  out) :: phi(lo(1)-ng:,lo(2)-ng:)

! Local
  integer :: i,j

  ! set phi = h^n for applyop on RHS
    do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        phi(i,j) = s2(i,j,rhoh_comp)/s2(i,j,rho_comp)
     enddo
  enddo

end subroutine make_initial_phi_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute betas and phi for 3d problems
subroutine make_initial_phi_3d(lo,hi,ng,ng_s,s2,phi)

  integer        , intent(in   ) :: lo(:),hi(:)
  integer        , intent(in   ) :: ng,ng_s
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
  real(kind=dp_t), intent(  out) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

! Local
  integer :: i,j,k

  ! set phi = h^n for applyop on RHS
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           phi(i,j,k) = s2(i,j,k,rhoh_comp)/s2(i,j,k,rho_comp)
        enddo
     enddo
  enddo

end subroutine make_initial_phi_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add rho h to RHS for 2d problems
subroutine add_rhoh_to_rh_2d(lo,hi,ng_rh,ng_s,s2,rh)

  integer        , intent(in   ) :: lo(:),hi(:),ng_rh,ng_s
  real(kind=dp_t), intent(inout) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,:)
  real(kind=dp_t), intent(inout) :: rh(lo(1)-ng_rh:,lo(2)-ng_rh:)

! Local
  integer :: i,j

  ! rh += rho h
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        rh(i,j) = rh(i,j) + s2(i,j,rhoh_comp)
     enddo
  enddo

end subroutine add_rhoh_to_rh_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add rho h to RHS for 3d problems
subroutine add_rhoh_to_rh_3d(lo,hi,ng_rh,ng_s,s2,rh)

  integer        , intent(in   ) :: lo(:),hi(:),ng_rh,ng_s
  real(kind=dp_t), intent(inout) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
  real(kind=dp_t), intent(inout) :: rh(lo(1)-ng_rh:,lo(2)-ng_rh:,lo(3)-ng_rh:)

! Local
  integer :: i,j,k

  ! rh += rho h
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           rh(i,j,k) = rh(i,j,k) + s2(i,j,k,rhoh_comp)
        enddo
     enddo
  enddo

end subroutine add_rhoh_to_rh_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute \rho h for 2d problems
subroutine compute_rhoh_2d(lo,hi,ng,ng_s,phi,s2)

  integer        , intent(in   ) :: lo(:),hi(:),ng,ng_s
  real(kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(inout) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,:)

! Local
  integer :: i,j

  ! Compute updated (\rho h) = \rho^(2)h^(2')
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        s2(i,j,rhoh_comp) = s2(i,j,rho_comp)*phi(i,j)
     enddo
  enddo

end subroutine compute_rhoh_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute \rho h for 3d problems
subroutine compute_rhoh_3d(lo,hi,ng,ng_s,phi,s2)

  integer        , intent(in   ) :: lo(:),hi(:),ng,ng_s
  real(kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(inout) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

! Local
  integer :: i,j,k

  ! Compute updated (\rho h) = \rho^(2)h^(2')
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           s2(i,j,k,rhoh_comp) = s2(i,j,k,rho_comp)*phi(i,j,k)
        enddo
     enddo
  enddo

end subroutine compute_rhoh_3d

end module thermal_conduct_module
