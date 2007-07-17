module thermal_conduct_module

  use bl_types
  use bc_module
  use multifab_module
  use boxarray_module
  use stencil_module
  use macproject_module

  implicit none

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Crank-Nicholson solve for enthalpy, taking into account only the
! enthalpy-diffusion terms in the temperature conduction term.
! See paper IV, steps 4a and 8a.
subroutine thermal_conduct(mla,dx,dt,sold,s2)

  type(ml_layout), intent(inout) :: mla
  real(dp_t)     , intent(in   ) :: dx(:,:)
  real(dp_t)     , intent(in   ) :: dt
  type(multifab) , intent(in   ) :: sold(:)
  type(multifab) , intent(inout) :: s2(:)

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

     if (parallel_IOProcessor()) print *,'... Setting alpha = rho ...'
     ! Copy rho directly into alpha
     call multifab_copy_c(alpha(n),1,s2(n),rho_comp,1)
     
     if (parallel_IOProcessor()) print *,'... Setting beta ...'
     ! Compute k_th^(2) - (temporarily set to 1 until I hook it into the eos)
     call setval(kth2(n),ONE,all=.true.)
          
     ! Compute c_p^(2) - (temporarily set to 1 until I hook it into the eos)
     call setval(cp2(n),ONE,all=.true.)
     
     ! Create beta = \frac{\Delta t k_th^(2)}{2 c_p^(2)}
     ng    = alpha(n)%ng
     ng_rh = rh(n)%ng
     ng_s  = sold(n)%ng

     do i=1,sold(n)%nboxes
        if (multifab_remote(sold(n),i)) cycle
        betap => dataptr(beta(n),i)
        kth2p => dataptr(kth2(n),i)
        cp2p  => dataptr(cp2(n),i)
        lo =  lwb(get_box(sold(n), i))
        hi =  upb(get_box(sold(n), i))
        select case (dm)
        case (2)
           call make_thermal_beta_2d(lo,hi,ng,dt,betap(:,:,1,1), &
                                     kth2p(:,:,1,1),cp2p(:,:,1,1))
        case (3)
           call make_thermal_beta_3d(lo,hi,ng,dt,betap(:,:,:,1), &
                                     kth2p(:,:,:,1),cp2p(:,:,:,1))
        end select
     end do
     
     if (parallel_IOProcessor()) print *,'... Making RHS ...'
     ! Compute k_th^n
     ! Temporarily set to 1 until I hook into eos
     call setval(kthold(n),ONE,all=.true.)
     
     ! Compute c_p^n - (temporarily set to 1 until I hook it into the eos)
     call setval(cpold(n),ONE,all=.true.)
     
     ! RHS = (\rho h)^(2)+\nabla\cdot(\frac{\Delta t k_th^n}{2 c_p^n}\nabla h)
     ng    = alpha(n)%ng
     ng_rh = rh(n)%ng
     ng_s  = sold(n)%ng

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
                                    ktholdp(:,:,1,1),cpoldp(:,:,1,1), &
                                    soldp(:,:,1,:),s2p(:,:,1,:))
        case (3)
           call make_thermal_rhs_3d(lo,hi,ng,ng_rh,ng_s,dt,rhp(:,:,:,1), &
                                    ktholdp(:,:,:,1),cpoldp(:,:,:,1), &
                                    soldp(:,:,:,:),s2p(:,:,:,:))
        end select
     end do
   
  enddo

  if (parallel_IOProcessor()) print *,'... Calling solver ...'
  ! Compute solution to (alpha - \nabla\cdot\beta\nabla)\phi = RHS




  do n=1,nlevs
     ! Compute updated (\rho h) = \rho^(2)h^(2')
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
  enddo

  deallocate(rh,phi,alpha,beta)
  deallocate(kthold,kth2,cpold,cp2)

end subroutine thermal_conduct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute beta for 2d problems
subroutine make_thermal_beta_2d(lo,hi,ng,dt,beta,kth,cp)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(dp_t)     , intent(in   ) :: dt
  real(kind=dp_t), intent(  out) :: beta(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(in   ) :: kth(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(in   ) :: cp(lo(1)-ng:,lo(2)-ng:)

! Local
  integer :: i,j

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        beta(i,j) = HALF*dt*kth(i,j)/cp(i,j)
     enddo
  enddo

end subroutine make_thermal_beta_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute beta for 3d problems
subroutine make_thermal_beta_3d(lo,hi,ng,dt,beta,kth,cp)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(dp_t)     , intent(in   ) :: dt
  real(kind=dp_t), intent(  out) :: beta(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(in   ) :: kth(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(in   ) :: cp(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

! Local
  integer :: i,j,k

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           beta(i,j,k) = HALF*dt*kth(i,j,k)/cp(i,j,k)
        enddo
     enddo
  enddo

end subroutine make_thermal_beta_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute RHS for 2d problems
subroutine make_thermal_rhs_2d(lo,hi,ng,ng_rh,ng_s,dt,rh,kphold,cpold,sold,s2)

  integer        , intent(in   ) :: lo(:),hi(:),ng,ng_rh,ng_s
  real(dp_t)     , intent(in   ) :: dt
  real(kind=dp_t), intent(  out) :: rh(lo(1)-ng_rh:,lo(2)-ng_rh:)
  real(kind=dp_t), intent(in   ) :: kphold(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(in   ) :: cpold(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_s:,lo(2)-ng_s:,:)
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,:)

! Local
  integer :: i,j

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        ! Temporary - need to fix
        rh(i,j) = ZERO
     enddo
  enddo

end subroutine make_thermal_rhs_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute RHS for 3d problems
subroutine make_thermal_rhs_3d(lo,hi,ng,ng_rh,ng_s,dt,rh,kphold,cpold,sold,s2)

  integer        , intent(in   ) :: lo(:),hi(:),ng,ng_rh,ng_s
  real(dp_t)     , intent(in   ) :: dt
  real(kind=dp_t), intent(  out) :: rh(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
  real(kind=dp_t), intent(in   ) :: kphold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(in   ) :: cpold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
  real(kind=dp_t), intent(in   ) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

! Local
  integer :: i,j,k

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           ! Temporary - need to fix
           rh(i,j,k) = ZERO
        enddo
     enddo
  enddo

end subroutine make_thermal_rhs_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute \rho h for 2d problems
subroutine compute_rhoh_2d(lo,hi,ng,ng_s,phi,s2)

  integer        , intent(in   ) :: lo(:),hi(:),ng,ng_s
  real(kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(inout) :: s2(lo(1)-ng_s:,lo(2)-ng_s:,:)

! Local
  integer :: i,j

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

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           s2(i,j,k,rhoh_comp) = s2(i,j,k,rho_comp)*phi(i,j,k)
        enddo
     enddo
  enddo

end subroutine compute_rhoh_3d

end module thermal_conduct_module
