module make_explicit_thermal_module

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
! Compute the explicit thermal conduction term needed for "S"
subroutine make_explicit_thermal(mla,dx,dt,thermal,s,p0, &
                                 mg_verbose,cg_verbose,the_bc_tower)

  type(ml_layout), intent(inout) :: mla
  real(dp_t)     , intent(in   ) :: dx(:,:),dt
  type(multifab) , intent(inout) :: thermal(:)
  type(multifab) , intent(in   ) :: s(:)
  real(kind=dp_t), intent(in   ) :: p0(0:)
  integer        , intent(in   ) :: mg_verbose,cg_verbose
  type(bc_tower) , intent(in   ) :: the_bc_tower

  ! Local
  type(multifab), allocatable :: h(:),alpha(:),beta(:)
  type(multifab), allocatable :: scalefactor(:),ccbeta(:)
  integer                     :: i,n,nlevs,dm,stencil_order,ng_0,ng_1,ng_3
  integer                     :: lo(s(1)%dim),hi(s(1)%dim)
  real(kind=dp_t), pointer    :: thermalp(:,:,:,:),sp(:,:,:,:)
  real(kind=dp_t), pointer    :: hp(:,:,:,:),ccbetap(:,:,:,:)
  real(kind=dp_t), pointer    :: betap(:,:,:,:),scalefactorp(:,:,:,:)

  nlevs = mla%nlevel
  dm      = mla%dim
  stencil_order = 2

  allocate(h(nlevs),alpha(nlevs),beta(nlevs),scalefactor(nlevs),ccbeta(nlevs))

  do n = 1,nlevs
     call multifab_build(          h(n), mla%la(n), 1, 1)
     call multifab_build(      alpha(n), mla%la(n), 1, 0)
     call multifab_build(       beta(n), mla%la(n),dm, 1)
     call multifab_build(     ccbeta(n), mla%la(n), 1, 1)
     call multifab_build(scalefactor(n), mla%la(n), 1, 1)  

     call setval(          h(n), ZERO, all=.true.)
     call setval(      alpha(n), ZERO, all=.true.)
     call setval(       beta(n), ZERO, all=.true.)
     call setval(scalefactor(n), ZERO, all=.true.)
     call setval(     ccbeta(n), ZERO, all=.true.)
  end do

  do n=1,nlevs
     ng_0 = thermal(n)%ng
     ng_1 = h(n)%ng
     ng_3 = s(n)%ng

     do i=1,s(n)%nboxes
        if (multifab_remote(s(n),i)) cycle
        sp => dataptr(s(n),i)
        betap => dataptr(beta(n),i)
        hp => dataptr(h(n),i)
        scalefactorp => dataptr(scalefactor(n),i)
        ccbetap => dataptr(ccbeta(n),i)
        lo =  lwb(get_box(s(n), i))
        hi =  upb(get_box(s(n), i))
        select case (dm)
        case (2)
           call make_thermal_coeffs_2d(lo,hi,dt,dx(n,:),ng_0,ng_1,ng_3, &
                                       p0,sp(:,:,1,:),betap(:,:,1,:), &
                                       hp(:,:,1,1),scalefactorp(:,:,1,1), &
                                       ccbetap(:,:,1,1))
        case (3)
           call make_thermal_coeffs_3d(lo,hi,dt,dx(n,:),ng_0,ng_1,ng_3, &
                                       p0,sp(:,:,:,:),betap(:,:,:,:), &
                                       hp(:,:,:,1),scalefactorp(:,:,:,1), &
                                       ccbetap(:,:,:,1))
        end select
     end do
  enddo


  ! applyop
  call mac_applyop(mla,thermal,h,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
                   stencil_order,mla%mba%rr,mg_verbose,cg_verbose)

  ! scale final answer
  do n=1,nlevs
     ng_0 = thermal(n)%ng

     do i=1,s(n)%nboxes
        if (multifab_remote(s(n),i)) cycle
        thermalp => dataptr(thermal(n),i)
        scalefactorp => dataptr(scalefactor(n),i)
        lo =  lwb(get_box(s(n), i))
        hi =  upb(get_box(s(n), i))
        select case (dm)
        case (2)
           call scale_thermal_2d(lo,hi,ng_0, &
                                 thermalp(:,:,1,1),scalefactorp(:,:,1,1))
        case (3)
           call scale_thermal_3d(lo,hi,ng_0, &
                                 thermalp(:,:,:,1),scalefactorp(:,:,:,1))
        end select
     end do
  enddo

  ! Deallocate memory
  do n = 1,nlevs
     call destroy(h(n))
     call destroy(alpha(n))
     call destroy(beta(n))
     call destroy(ccbeta(n))
     call destroy(scalefactor(n))

  enddo

  deallocate(h,alpha,beta,ccbeta,scalefactor)

end subroutine make_explicit_thermal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine make_thermal_coeffs_2d(lo,hi,dt,dx,ng_0,ng_1,ng_3, &
                                  p0,s,beta,h,scalefactor,ccbeta)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt,dx(:)
  integer        , intent(in   ) :: ng_0,ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,:)
  real(kind=dp_t), intent(inout) :: beta(lo(1)-ng_1:,lo(2)-ng_1:,:)
  real(kind=dp_t), intent(inout) :: h(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: scalefactor(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: ccbeta(lo(1)-ng_1:,lo(2)-ng_1:)

  integer :: i,j
  integer :: nx,ny

! dens, pres, and xmass are inputs
  input_flag = 4
  do_diag = .false.

  nx = size(beta,dim=1) - 2*ng_1
  ny = size(beta,dim=2) - 2*ng_1

  ! Compute c_p^(2), k_th^2, and beta
  do j=lo(2)-1,hi(2)+1
     do i=lo(1)-1,hi(1)+1

        den_row(1) = s(i,j,rho_comp)
        if(j .eq. -1) then
           p_row(1) = p0(0)
        else if(j .eq. hi(2)+1) then
           p_row(1) = p0(hi(2))
        else
           p_row(1) = p0(j)
        endif
        xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

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

        ccbeta(i,j) = -conduct_row(1)/cp_row(1)
        scalefactor(i,j) = dpdt_row(1)/(den_row(1)**2*cp_row(1)*dpdr_row(1))

     enddo
  enddo

  do j = 0,ny-1
     do i = 0,nx
        beta(i,j,1) = (ccbeta(i,j) + ccbeta(i-1,j)) / TWO
     end do
  end do
  
  do j = 0,ny
     do i = 0,nx-1
        beta(i,j,2) = (ccbeta(i,j) + ccbeta(i,j-1)) / TWO
     end do
  end do

  ! set phi = h^n for applyop on RHS
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        h(i,j) = s(i,j,rhoh_comp)/s(i,j,rho_comp)
     enddo
  enddo

end subroutine make_thermal_coeffs_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine make_thermal_coeffs_3d(lo,hi,dt,dx,ng_0,ng_1,ng_3, &
                                  p0,s,beta,h,scalefactor,ccbeta)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt,dx(:)
  integer        , intent(in   ) :: ng_0,ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:,:)
  real(kind=dp_t), intent(inout) :: beta(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:,:)
  real(kind=dp_t), intent(inout) :: h(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: scalefactor(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: ccbeta(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)

  integer :: i,j,k
  integer :: nx,ny,nz
  real(kind=dp_t), allocatable :: p0_cart(:,:,:)

! dens, pres, and xmass are inputs
  input_flag = 4
  do_diag = .false.

  nx = size(beta,dim=1) - 2*ng_1
  ny = size(beta,dim=2) - 2*ng_1
  nz = size(beta,dim=3) - 2*ng_1

  if (spherical .eq. 1) then
     allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
     call fill_3d_data(p0_cart,p0,lo,hi,dx,0)
  end if

  ! Compute c_p^(2), k_th^2, and beta
  do k=lo(3)-1,hi(3)+1
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           
           den_row(1) = s(i,j,k,rho_comp)
           xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

           if(spherical .eq. 0) then
              p_row(1) = p0(k)
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
           
           ccbeta(i,j,k) = -ONE/cp_row(1)
           scalefactor(i,j,k) = dpdt_row(1)/(den_row(1)**2*cp_row(1)*dpdr_row(1))
           
        enddo
     enddo
  enddo
  
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

  ! set phi = h^n for applyop on RHS
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           h(i,j,k) = s(i,j,k,rhoh_comp)/s(i,j,k,rho_comp)
        enddo
     enddo
  enddo

end subroutine make_thermal_coeffs_3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine scale_thermal_2d(lo,hi,ng,thermal,scalefactor)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(  out) :: thermal(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(in   ) :: scalefactor(lo(1)-ng:,lo(2)-ng:)

! Local
  integer :: i,j

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        thermal(i,j) = scalefactor(i,j)*thermal(i,j)
     enddo
  enddo
  
end subroutine scale_thermal_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine scale_thermal_3d(lo,hi,ng,thermal,scalefactor)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(  out) :: thermal(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(in   ) :: scalefactor(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

! Local
  integer :: i,j,k

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           thermal(i,j,k) = scalefactor(i,j,k)*thermal(i,j,k)
        enddo
     enddo
  enddo

end subroutine scale_thermal_3d


end module make_explicit_thermal_module
