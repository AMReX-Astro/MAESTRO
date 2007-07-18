module thermal_conduct_module

  use bl_types
  use bc_module
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
subroutine thermal_conduct(mla,dx,dt,sold,s2,p0old,p02)

  type(ml_layout), intent(inout) :: mla
  real(dp_t)     , intent(in   ) :: dx(:,:)
  real(dp_t)     , intent(in   ) :: dt
  type(multifab) , intent(in   ) :: sold(:)
  type(multifab) , intent(inout) :: s2(:)
  real(kind=dp_t), intent(in   ) :: p0old(0:)
  real(kind=dp_t), intent(in   ) :: p02(0:)

! Local
  type(multifab), allocatable :: rh(:),phi(:),alpha(:),beta(:)
  type(multifab), allocatable :: kthold(:),kth2(:),cpold(:),cp2(:)
  integer                     :: i,n,nlevs,dm,ng,ng_rh,ng_s
  integer                     :: lo(sold(1)%dim),hi(sold(1)%dim)
  real(kind=dp_t), pointer    :: soldp(:,:,:,:),s2p(:,:,:,:)
  real(kind=dp_t), pointer    :: rhp(:,:,:,:),phip(:,:,:,:)
  real(kind=dp_t), pointer    :: betap(:,:,:,:)
  real(kind=dp_t), pointer    :: ktholdp(:,:,:,:),kth2p(:,:,:,:)
  real(kind=dp_t), pointer    :: cpoldp(:,:,:,:),cp2p(:,:,:,:)

  integer                     :: test1,test2

  if (parallel_IOProcessor()) print *,'... Entering thermal_conduct ...'

  nlevs = mla%nlevel
  dm    = mla%dim

  allocate(rh(nlevs),phi(nlevs),alpha(nlevs),beta(nlevs))
  allocate(kthold(nlevs),kth2(nlevs),cpold(nlevs),cp2(nlevs))

  do n = 1,nlevs
     call multifab_build(   rh(n), mla%la(n), 1, 0)
     call multifab_build(  phi(n), mla%la(n), 1, 1)
     call multifab_build(alpha(n), mla%la(n), 1, 1)
     call multifab_build( beta(n), mla%la(n), 1, 1)

     call multifab_build(kthold(n), mla%la(n), 1, 1)
     call multifab_build(  kth2(n), mla%la(n), 1, 1)
     call multifab_build( cpold(n), mla%la(n), 1, 1)
     call multifab_build(   cp2(n), mla%la(n), 1, 1)
  end do

  do n=1,nlevs

     ng    = alpha(n)%ng
     ng_rh = rh(n)%ng
     ng_s  = sold(n)%ng

     if (parallel_IOProcessor()) print *,'... Setting alpha = rho ...'

     ! Copy rho directly into alpha
     call multifab_copy_c(alpha(n),1,s2(n),rho_comp,1)
     
     ! Create beta = \frac{\Delta t k_th^(2)}{2 c_p^(2)}
     if (parallel_IOProcessor()) print *,'... Setting beta ...'

     do i=1,sold(n)%nboxes
        if (multifab_remote(sold(n),i)) cycle
        betap => dataptr(beta(n),i)
        s2p   => dataptr(s2(n),i)
        kth2p => dataptr(kth2(n),i)
        cp2p  => dataptr(cp2(n),i)
        lo =  lwb(get_box(sold(n), i))
        hi =  upb(get_box(sold(n), i))
        select case (dm)
        case (2)
           call make_thermal_beta_2d(lo,hi,ng,ng_s,dt,betap(:,:,1,1),p02, &
                                     s2p(:,:,1,:),kth2p(:,:,1,1), &
                                     cp2p(:,:,1,1),dx(n,:))
        case (3)
           call make_thermal_beta_3d(lo,hi,ng,ng_s,dt,betap(:,:,:,1),p02, &
                                     s2p(:,:,:,:),kth2p(:,:,:,1), &
                                     cp2p(:,:,:,1),dx(n,:))
        end select
     end do
     
     ! Make RHS 
     !    = (\rho h)^(2) + \nabla\cdot(\frac{\Delta t k_th^n}{2 c_p^n}\nabla h)
     if (parallel_IOProcessor()) print *,'... Making RHS ...'
 
     do i=1,sold(n)%nboxes
        if (multifab_remote(sold(n),i)) cycle
        rhp     => dataptr(rh(n),i)
        ktholdp => dataptr(kthold(n),i)
        cpoldp  => dataptr(cpold(n),i)
        soldp   => dataptr(sold(n),i)
        s2p     => dataptr(s2(n),i)
        lo =  lwb(get_box(sold(n), i))
        hi =  upb(get_box(sold(n), i))
        select case (dm)
        case (2)
           call make_thermal_rhs_2d(lo,hi,ng,ng_rh,ng_s,dt,rhp(:,:,1,1), &
                                    p0old,ktholdp(:,:,1,1),cpoldp(:,:,1,1), &
                                    soldp(:,:,1,:),s2p(:,:,1,:),dx(n,:))
        case (3)
           call make_thermal_rhs_3d(lo,hi,ng,ng_rh,ng_s,dt,rhp(:,:,:,1), &
                                    p0old,ktholdp(:,:,:,1),cpoldp(:,:,:,1), &
                                    soldp(:,:,:,:),s2p(:,:,:,:),dx(n,:))
        end select
     end do
   
  enddo

  ! Compute solution to (alpha - \nabla\cdot\beta\nabla)\phi = RHS
  if (parallel_IOProcessor()) print *,'... Calling solver ...'






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
           call make_thermal_rhoh_2d(lo,hi,ng,ng_s,phip(:,:,1,1),s2p(:,:,1,:))
        case (3)
           call make_thermal_rhoh_3d(lo,hi,ng,ng_s,phip(:,:,:,1),s2p(:,:,:,:))
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
  enddo

  deallocate(rh,phi,alpha,beta)
  deallocate(kthold,kth2,cpold,cp2)

end subroutine thermal_conduct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute beta for 2d problems
subroutine make_thermal_beta_2d(lo,hi,ng,ng_s,dt,beta,p02,s2,kth2,cp2,dx)

  integer        , intent(in   ) :: lo(:),hi(:),ng,ng_s
  real(dp_t)     , intent(in   ) :: dt
  real(kind=dp_t), intent(  out) :: beta(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(in   ) :: p02(0:)
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,:)
  real(kind=dp_t), intent(inout) :: kth2(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(inout) :: cp2(lo(1)-ng:,lo(2)-ng:)
  real(dp_t)    ,  intent(in   ) :: dx(:)

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

        cp2(i,j) = cp_row(1)
        kth2(i,j) = ONE ! Temporarily set to 1
        beta(i,j) = HALF*dt*kth2(i,j)/cp2(i,j)
     enddo
  enddo

end subroutine make_thermal_beta_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute beta for 3d problems
subroutine make_thermal_beta_3d(lo,hi,ng,ng_s,dt,beta,p02,s2,kth2,cp2,dx)

  integer        , intent(in   ) :: lo(:),hi(:),ng,ng_s
  real(dp_t)     , intent(in   ) :: dt
  real(kind=dp_t), intent(  out) :: beta(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(in   ) :: p02(0:)
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
  real(kind=dp_t), intent(inout) :: kth2(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(inout) :: cp2(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(dp_t)    ,  intent(in   ) :: dx(:)

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

  ! Compute c_p^(2), k_th^(2), and beta
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           den_row(1) = s2(i,j,k,rho_comp)
           xn_zone(:) = s2(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

           if (spherical .eq. 0) then
              p_row(1) = p02(k)
           else
              p_row(1) = p0_cart(i,j,k)
           end if
           
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

           cp2(i,j,k) = cp_row(1)
           kth2(i,j,k) = ONE ! Temporarily set to 1
           beta(i,j,k) = HALF*dt*kth2(i,j,k)/cp2(i,j,k)
        enddo
     enddo
  enddo

end subroutine make_thermal_beta_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute RHS for 2d problems
subroutine make_thermal_rhs_2d(lo,hi,ng,ng_rh,ng_s,dt,rh,p0old,kthold,cpold, &
                               sold,s2,dx)

  integer        , intent(in   ) :: lo(:),hi(:),ng,ng_rh,ng_s
  real(dp_t)     , intent(in   ) :: dt
  real(kind=dp_t), intent(  out) :: rh(lo(1)-ng_rh:,lo(2)-ng_rh:)
  real(kind=dp_t), intent(in   ) :: p0old(0:)
  real(kind=dp_t), intent(inout) :: kthold(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(inout) :: cpold(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_s:,lo(2)-ng_s:,:)
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,:)
  real(dp_t)     , intent(in   ) :: dx(:)

! Local
  integer :: i,j

! dens, pres, and xmass are inputs
  input_flag = 4
  do_diag = .false.

  ! Compute c_p^n and k_th^n
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

        cpold(i,j) = cp_row(1)
        kthold(i,j) = ONE ! Temporarily set to 1
     enddo
  enddo

  ! Compute residual = del dot (dt*kthold/(2*cpold)) nabla h => store in rh




  ! Compute rh += (\rho h)^{(2)}
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        rh(i,j) = rh(i,j) + s2(i,j,rhoh_comp)
     enddo
  enddo

end subroutine make_thermal_rhs_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute RHS for 3d problems
subroutine make_thermal_rhs_3d(lo,hi,ng,ng_rh,ng_s,dt,rh,p0old,kthold,cpold, &
                               sold,s2,dx)

  integer        , intent(in   ) :: lo(:),hi(:),ng,ng_rh,ng_s
  real(dp_t)     , intent(in   ) :: dt
  real(kind=dp_t), intent(  out) :: rh(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
  real(kind=dp_t), intent(in   ) :: p0old(0:)
  real(kind=dp_t), intent(inout) :: kthold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(inout) :: cpold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
  real(dp_t)     , intent(in   ) :: dx(:)

! Local
  integer :: i,j,k
  real(kind=dp_t), allocatable :: p0_cart(:,:,:)

  if (spherical .eq. 1) then
     allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
     call fill_3d_data(p0_cart,p0old,lo,hi,dx,0)
  end if

! dens, pres, and xmass are inputs
  input_flag = 4
  do_diag = .false.

  ! Compute c_p^n and k_th^n
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           den_row(1) = sold(i,j,k,rho_comp)
           xn_zone(:) = sold(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

           if (spherical .eq. 0) then
              p_row(1) = p0old(k)
           else
              p_row(1) = p0_cart(i,j,k)
           end if
           
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

           cpold(i,j,k) = cp_row(1)
           kthold(i,j,k) = ONE ! Temporarily set to 1
        enddo
     enddo
  enddo

  ! Compute residual = del dot (dt*kthold/(2*cpold)) nabla h; store in rh




  ! Compute rh += (\rho h)^{(2)}
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           rh(i,j,k) = rh(i,j,k) + s2(i,j,k,rhoh_comp)
        enddo
     enddo
  enddo

end subroutine make_thermal_rhs_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute \rho h for 2d problems
subroutine make_thermal_rhoh_2d(lo,hi,ng,ng_s,phi,s2)

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

end subroutine make_thermal_rhoh_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute \rho h for 3d problems
subroutine make_thermal_rhoh_3d(lo,hi,ng,ng_s,phi,s2)

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

end subroutine make_thermal_rhoh_3d

end module thermal_conduct_module
