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
                           mg_verbose,cg_verbose)

  type(ml_layout), intent(inout) :: mla
  real(dp_t)     , intent(in   ) :: dx(:,:)
  real(dp_t)     , intent(in   ) :: dt
  type(multifab) , intent(in   ) :: sold(:)
  type(multifab) , intent(inout) :: s2(:)
  real(kind=dp_t), intent(in   ) :: p0old(0:)
  real(kind=dp_t), intent(in   ) :: p02(0:)
  integer        , intent(in   ) :: mg_verbose,cg_verbose

! Local
  type(multifab), allocatable :: rh(:),phi(:),alpha(:),beta(:)
  type(multifab), allocatable :: kthold(:),kth2(:),cpold(:),cp2(:)
  type(multifab), allocatable :: rhsbeta(:),res(:)
  integer                     :: i,n,nlevs,dm,ng,ng_rh,ng_s
  integer                     :: lo(sold(1)%dim),hi(sold(1)%dim)
  real(kind=dp_t), pointer    :: soldp(:,:,:,:),s2p(:,:,:,:)
  real(kind=dp_t), pointer    :: rhp(:,:,:,:),phip(:,:,:,:)
  real(kind=dp_t), pointer    :: betap(:,:,:,:),rhsbetap(:,:,:,:)
  real(kind=dp_t), pointer    :: ktholdp(:,:,:,:),kth2p(:,:,:,:)
  real(kind=dp_t), pointer    :: cpoldp(:,:,:,:),cp2p(:,:,:,:)

  integer                     :: test1,test2,nscal
  integer                     :: stencil_order,bc_comp
  integer       , allocatable :: domain_phys_bc(:,:)
  type(bc_tower)              :: the_bc_tower
  type(box)     , allocatable :: domain_box(:)
  integer                     :: bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi

  bcx_lo = SLIP_WALL
  bcy_lo = SLIP_WALL
  bcz_lo = SLIP_WALL
  bcx_hi = SLIP_WALL
  bcy_hi = SLIP_WALL
  bcz_hi = SLIP_WALL

  if (parallel_IOProcessor()) print *,'... Entering thermal_conduct ...'

  nlevs = mla%nlevel
  dm    = mla%dim

  allocate(rh(nlevs),phi(nlevs),alpha(nlevs),beta(nlevs))
  allocate(kthold(nlevs),kth2(nlevs),cpold(nlevs),cp2(nlevs))
  allocate(rhsbeta(nlevs),res(nlevs))

  do n = 1,nlevs
     call multifab_build(   rh(n), mla%la(n), 1, 0)
     call multifab_build(  phi(n), mla%la(n), 1, 1)
     call multifab_build(alpha(n), mla%la(n), 1, 1)
     call multifab_build( beta(n), mla%la(n), 1, 1)

     call multifab_build(kthold(n), mla%la(n), 1, 1)
     call multifab_build(  kth2(n), mla%la(n), 1, 1)
     call multifab_build( cpold(n), mla%la(n), 1, 1)
     call multifab_build(   cp2(n), mla%la(n), 1, 1)

     call multifab_build(rhsbeta(n),mla%la(n), 1, 1)
     call multifab_build(    res(n),mla%la(n), 1, 0)
  end do

  do n=1,nlevs

     ng    = alpha(n)%ng
     ng_rh = rh(n)%ng
     ng_s  = sold(n)%ng

     if (parallel_IOProcessor()) print *,'... Setting alpha = rho ...'

     ! Copy rho^(2) directly into alpha
     call multifab_copy_c(alpha(n),1,s2(n),rho_comp,1)
     
     ! Create beta = \frac{\Delta t k_th^(2)}{2 c_p^(2)}
     ! Create rhsbeta = dt*k_th^n/(2*c_p^n)
     ! Copy h^n into phi
     if (parallel_IOProcessor()) print *,'... Computing betas and phi ...'

     do i=1,sold(n)%nboxes
        if (multifab_remote(sold(n),i)) cycle
        soldp    => dataptr(sold(n),i)
        s2p      => dataptr(s2(n),i)
        betap    => dataptr(beta(n),i)
        rhsbetap => dataptr(rhsbeta(n),i)
        phip     => dataptr(phi(n),i)
        lo =  lwb(get_box(sold(n), i))
        hi =  upb(get_box(sold(n), i))
        select case (dm)
        case (2)
           call make_betas_and_phi_2d(lo,hi,dt,dx(n,:),ng,ng_rh,ng_s, &
                                      p0old,p02,soldp(:,:,1,:),s2p(:,:,1,:), &
                                      betap(:,:,1,1),rhsbetap(:,:,1,1), &
                                      phip(:,:,1,1))
        case (3)
           call make_betas_and_phi_3d(lo,hi,dt,dx(n,:),ng,ng_rh,ng_s, &
                                      p0old,p02,soldp(:,:,:,:),s2p(:,:,:,:), &
                                      betap(:,:,:,1),rhsbetap(:,:,:,1), &
                                      phip(:,:,:,1))
        end select
     end do

  enddo

  ! define a solver
  allocate(domain_box(nlevs))
  do n = 1,nlevs
     domain_box(n) = layout_get_pd(mla%la(n))
  end do

  allocate(domain_phys_bc(dm,2))
  domain_phys_bc(1,1) = bcx_lo
  domain_phys_bc(1,2) = bcx_hi
  if (dm > 1) then
     domain_phys_bc(2,1) = bcy_lo
     domain_phys_bc(2,2) = bcy_hi
  end if
  if (dm > 2) then
     domain_phys_bc(3,1) = bcz_lo
     domain_phys_bc(3,2) = bcz_hi
  end if

  nscal = 1

  ! call bc_tower_build( the_bc_tower,mla,domain_phys_bc,domain_box,nscal)

  bc_comp = rhoh_comp
  stencil_order = 2

  ! compute residual to get del dot rhsbeta grad h term in RHS
  !call mac_applyop(mla,res,phi,alpha,beta,dx,&
  !                 the_bc_tower,bc_comp,stencil_order,mla%mba%rr, &
  !                 mg_verbose,cg_verbose)

  ! multiply residual by -1



  ! add (\rho h)^(2) to RHS
  do n=1,nlevs

     ng    = alpha(n)%ng
     ng_rh = rh(n)%ng
     ng_s  = sold(n)%ng

     if (parallel_IOProcessor()) print *,'... Adding rho h to RHS ...'

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

  ! Compute solution to (alpha - \nabla\cdot\beta\nabla)\phi = RHS
  ! Then, h^(2') = phi
  if (parallel_IOProcessor()) print *,'... Calling solver ...'

  ! Define a new solver with different alpha and beta



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

     call destroy(kthold(n))
     call destroy(kth2(n))
     call destroy(cpold(n))
     call destroy(cp2(n))

     call destroy(rhsbeta(n))
     call destroy(res(n))
  enddo

  deallocate(rh,phi,alpha,beta)
  deallocate(kthold,kth2,cpold,cp2)
  deallocate(rhsbeta,res)

end subroutine thermal_conduct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute betas and phi for 2d problems
subroutine make_betas_and_phi_2d(lo,hi,dt,dx,ng,ng_rh,ng_s, &
                                 p0old,p02,sold,s2,beta,rhsbeta,phi)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt,dx(:)
  integer        , intent(in   ) :: ng,ng_rh,ng_s
  real(kind=dp_t), intent(in   ) :: p0old(0:),p02(0:)
  real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_s:,lo(2)-ng_s:,:)
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,:)
  real(kind=dp_t), intent(  out) :: beta(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(  out) :: rhsbeta(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(  out) :: phi(lo(1)-ng:,lo(2)-ng:)

! Local
  integer :: i,j
  
! dens, pres, and xmass are inputs
  input_flag = 4
  do_diag = .false.

  ! Compute c_p^(2), k_th^2, and beta
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)

        den_row(1) = s2(i,j,rho_comp)
        p_row(1) = p02(j)
        xn_zone(:) = s2(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

        call eos(input_flag, den_row, temp_row, &
                 npts, nspec, &
                 xn_zone, aion, zion, &
                 p_row, h_row, e_row, & 
                 cv_row, cp_row, xne_row, eta_row, pele_row, &
                 dpdt_row, dpdr_row, dedt_row, dedr_row, &
                 dpdX_row, dhdX_row, &
                 gam1_row, cs_row, s_row, &
                 dsdt_row, dsdr_row, &
                 do_diag)

        beta(i,j) = HALF*dt*ONE/cp_row(1) ! k_th^(2) = 1 for now
     enddo
  enddo

 ! Compute c_p^n, k_th^n, and rhsbeta
    do j=lo(2),hi(2)
     do i=lo(1),hi(1)

        den_row(1) = sold(i,j,rho_comp)
        p_row(1) = p0old(j)
        xn_zone(:) = sold(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

        call eos(input_flag, den_row, temp_row, &
                 npts, nspec, &
                 xn_zone, aion, zion, &
                 p_row, h_row, e_row, & 
                 cv_row, cp_row, xne_row, eta_row, pele_row, &
                 dpdt_row, dpdr_row, dedt_row, dedr_row, &
                 dpdX_row, dhdX_row, &
                 gam1_row, cs_row, s_row, &
                 dsdt_row, dsdr_row, &
                 do_diag)

        rhsbeta(i,j) = (dt*ONE)/(TWO*cp_row(1)) ! k_th^n = 1 for now
     enddo
  enddo

  ! set phi = h^n for applyop on RHS
    do j=lo(2),hi(2)
     do i=lo(1),hi(1)

        phi(i,j) = sold(i,j,rhoh_comp)/sold(i,j,rho_comp)

        den_row(1) = sold(i,j,rho_comp)
        p_row(1) = p0old(j)
        xn_zone(:) = sold(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

        call eos(input_flag, den_row, temp_row, &
                 npts, nspec, &
                 xn_zone, aion, zion, &
                 p_row, h_row, e_row, & 
                 cv_row, cp_row, xne_row, eta_row, pele_row, &
                 dpdt_row, dpdr_row, dedt_row, dedr_row, &
                 dpdX_row, dhdX_row, &
                 gam1_row, cs_row, s_row, &
                 dsdt_row, dsdr_row, &
                 do_diag)

        rhsbeta(i,j) = (dt*ONE)/(TWO*cp_row(1)) ! k_th^n = 1 for now
     enddo
  enddo

end subroutine make_betas_and_phi_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute betas and phi for 3d problems
subroutine make_betas_and_phi_3d(lo,hi,dt,dx,ng,ng_rh,ng_s, &
                                 p0old,p02,sold,s2,beta,rhsbeta,phi)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt,dx(:)
  integer        , intent(in   ) :: ng,ng_rh,ng_s
  real(kind=dp_t), intent(in   ) :: p0old(0:),p02(0:)
  real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
  real(kind=dp_t), intent(  out) :: beta(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(  out) :: rhsbeta(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(  out) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

! Local
  integer :: i,j,k
  real(kind=dp_t), allocatable :: p0_cart(:,:,:)

  if (spherical .eq. 1) then
     allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
     call fill_3d_data(p0_cart,p02,lo,hi,dx,0)
  end if

! dens, pres, and xmass are inputs
  input_flag = 4
  do_diag = .false.

  ! Compute c_p^(2), k_th^2, and beta
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           
           den_row(1) = s2(i,j,k,rho_comp)
           xn_zone(:) = s2(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

           if(spherical .eq. 0) then
              p_row(1) = p02(j)
           else
              p_row(1) = p0_cart(i,j,k)
           endif
           
           call eos(input_flag, den_row, temp_row, &
                    npts, nspec, &
                    xn_zone, aion, zion, &
                    p_row, h_row, e_row, & 
                    cv_row, cp_row, xne_row, eta_row, pele_row, &
                    dpdt_row, dpdr_row, dedt_row, dedr_row, &
                    dpdX_row, dhdX_row, &
                    gam1_row, cs_row, s_row, &
                    dsdt_row, dsdr_row, &
                    do_diag)

           beta(i,j,k) = HALF*dt*ONE/cp_row(1) ! k_th^(2) = 1 for now
        enddo
     enddo
  enddo

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
              p_row(1) = p0old(j)
           else
              p_row(1) = p0_cart(i,j,k)
           endif
           
           call eos(input_flag, den_row, temp_row, &
                    npts, nspec, &
                    xn_zone, aion, zion, &
                    p_row, h_row, e_row, & 
                    cv_row, cp_row, xne_row, eta_row, pele_row, &
                    dpdt_row, dpdr_row, dedt_row, dedr_row, &
                    dpdX_row, dhdX_row, &
                    gam1_row, cs_row, s_row, &
                    dsdt_row, dsdr_row, &
                    do_diag)
           
           rhsbeta(i,j,k) = (dt*ONE)/(TWO*cp_row(1)) ! k_th^n = 1 for now
        enddo
     enddo
  enddo

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
