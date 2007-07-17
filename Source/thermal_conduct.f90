module thermal_conduct_module

  use bl_types
  use bc_module
  use multifab_module
  use boxarray_module
  use stencil_module
  use macproject_module

  implicit none

contains 

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
  type(multifab), allocatable :: kthn(:),kth2(:),cpn(:),cp2(:)
  integer                     :: n,nlevs

  if (parallel_IOProcessor()) print *,'... Entering thermal_conduct ...'

  nlevs = mla%nlevel

  allocate(rh(nlevs),phi(nlevs),alpha(nlevs),beta(nlevs))
  allocate(kthn(nlevs),kth2(nlevs),cpn(nlevs),cp2(nlevs))

  do n = 1,nlevs
     call multifab_build(   rh(n), mla%la(n), 1, 0)
     call multifab_build(  phi(n), mla%la(n), 1, 1)
     call multifab_build(alpha(n), mla%la(n), 1, 1)
     call multifab_build( beta(n), mla%la(n), 1, 1)

     call multifab_build( kthn(n), mla%la(n), 1, 1)
     call multifab_build( kth2(n), mla%la(n), 1, 1)
     call multifab_build(  cpn(n), mla%la(n), 1, 1)
     call multifab_build(  cp2(n), mla%la(n), 1, 1)
  end do

  if (parallel_IOProcessor()) print *,'... Setting alpha = rho ...'
  ! Copy rho directly into alpha
  call multifab_copy_c(alpha(n),1,s2(n),rho_comp,1)

  if (parallel_IOProcessor()) print *,'... Setting beta ...'
  ! Compute k_th^(2) - (temporarily set to 1 until I hook it into the eos)
  do n = 1,nlevs
     call setval(kth2(n),ONE,all=.true.)
  enddo

  ! Compute c_p^(2) - (temporarily set to 1 until I hook it into the eos)
  do n = 1,nlevs
     call setval(cp2(n),ONE,all=.true.)
  enddo

  ! Create beta = \frac{\Delta t k_th^(2)}{2 c_p^(2)}


  if (parallel_IOProcessor()) print *,'... Making RHS ...'
  ! Compute k_th^n
  ! Temporarily set to 1 until I hook into eos
  do n = 1,nlevs
     call setval(kthn(n),ONE,all=.true.)
  enddo

  ! Compute c_p^n - (temporarily set to 1 until I hook it into the eos)
  do n = 1,nlevs
     call setval(cpn(n),ONE,all=.true.)
  enddo

  ! RHS = (\rho h)^(2) + \nabla\cdot(\frac{\Delta t k_th^n}{2 c_p^n}\nabla h)


  if (parallel_IOProcessor()) print *,'... Calling solver ...'
  ! Compute solution to (alpha - \nabla\cdot\beta\nabla)\phi = RHS


  ! Compute updated (\rho h) = \rho^(2)h^(2')


  ! Deallocate memory
  do n = 1,nlevs
     call destroy(rh(n))
     call destroy(phi(n))
     call destroy(alpha(n))
     call destroy(beta(n))

     call destroy(kthn(n))
     call destroy(kth2(n))
     call destroy(cpn(n))
     call destroy(cp2(n))
  enddo

  deallocate(rh,phi,alpha,beta)
  deallocate(kthn,kth2,cpn,cp2)

end subroutine thermal_conduct

end module thermal_conduct_module
