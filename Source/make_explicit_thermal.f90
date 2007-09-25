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
  use probin_module

  implicit none

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the explicit thermal conduction term needed for "S"
subroutine make_explicit_thermal(mla,dx,thermal,s,p0,mg_verbose, &
     cg_verbose,the_bc_tower,temperature_diffusion)

  type(ml_layout), intent(inout) :: mla
  real(dp_t)     , intent(in   ) :: dx(:,:)
  type(multifab) , intent(inout) :: thermal(:)
  type(multifab) , intent(in   ) :: s(:)
  real(kind=dp_t), intent(in   ) :: p0(0:)
  integer        , intent(in   ) :: mg_verbose,cg_verbose
  type(bc_tower) , intent(in   ) :: the_bc_tower
  logical       ,  intent(in   ) :: temperature_diffusion

  ! Local
  type(multifab), allocatable :: phi(:),alpha(:),beta(:),xik(:)
  type(multifab), allocatable :: kth(:),kthovercp(:),resid(:)
  integer                     :: i,k,n,nlevs,dm,stencil_order,ng_0,ng_1,ng_3
  integer                     :: lo(s(1)%dim),hi(s(1)%dim)
  real(kind=dp_t), pointer    :: thermalp(:,:,:,:),sp(:,:,:,:)
  real(kind=dp_t), pointer    :: phip(:,:,:,:),betap(:,:,:,:)
  real(kind=dp_t), pointer    :: xikp(:,:,:,:)
  real(kind=dp_t), pointer    :: kthp(:,:,:,:),kthovercpp(:,:,:,:)
  real(kind=dp_t), pointer    :: residp(:,:,:,:)

  nlevs = mla%nlevel
  dm      = mla%dim
  stencil_order = 2

  allocate(phi(nlevs),alpha(nlevs),beta(nlevs),xik(nlevs))
  allocate(kth(nlevs),kthovercp(nlevs),resid(nlevs))

  do n = 1,nlevs
     call multifab_build(         phi(n), mla%la(n), 1, 1)
     call multifab_build(       alpha(n), mla%la(n), 1, 1)
     call multifab_build(        beta(n), mla%la(n),dm, 1)
     call multifab_build(         xik(n), mla%la(n),nspec, 1)
     call multifab_build(         kth(n), mla%la(n), 1, 1)
     call multifab_build(   kthovercp(n), mla%la(n), 1, 1)
     call multifab_build(       resid(n), mla%la(n), 1, 0)

     call setval(         phi(n), ZERO, all=.true.)
     call setval(       alpha(n), ZERO, all=.true.)
     call setval(        beta(n), ZERO, all=.true.)
     call setval(         xik(n), ZERO, all=.true.)
     call setval(         kth(n), ZERO, all=.true.)
     call setval(   kthovercp(n), ZERO, all=.true.)
     call setval(       resid(n), ZERO, all=.true.)
     call setval(     thermal(n), ZERO, all=.true.)
  end do

  ! create kth
  ! create kthovercp
  do n=1,nlevs
     ng_0 = resid(n)%ng
     ng_1 = phi(n)%ng
     ng_3 = s(n)%ng

     do i=1,s(n)%nboxes
        if (multifab_remote(s(n),i)) cycle
        sp            => dataptr(s(n),i)
        kthp          => dataptr(kth(n),i)
        kthovercpp    => dataptr(kthovercp(n),i)
        xikp          => dataptr(xik(n),i)
        lo =  lwb(get_box(s(n), i))
        hi =  upb(get_box(s(n), i))
        select case (dm)
        case (2)
           call make_coeffs_2d(lo,hi,dx(n,:),ng_1,ng_3,p0,sp(:,:,1,:), &
                               kthp(:,:,1,1),kthovercpp(:,:,1,1), &
                               xikp(:,:,1,:))
        case (3)
           call make_coeffs_3d(lo,hi,dx(n,:),ng_1,ng_3,p0,sp(:,:,:,:), &
                               kthp(:,:,:,1),kthovercpp(:,:,:,1), &
                               xikp(:,:,:,:))
        end select
     end do
  enddo


  if(temperature_diffusion) then
     ! load T into phi
     ! setup beta = -k_{th} for T term
     do n=1,nlevs
        ng_1 = phi(n)%ng
        ng_3 = s(n)%ng
        
        do i=1,s(n)%nboxes
           if (multifab_remote(s(n),i)) cycle
           sp => dataptr(s(n),i)
           kthp => dataptr(kth(n),i)
           phip => dataptr(phi(n),i)
           betap => dataptr(beta(n),i)
           lo =  lwb(get_box(s(n), i))
           hi =  upb(get_box(s(n), i))
           select case (dm)
           case (2)
              call setup_T_op_2d(lo,hi,ng_1,ng_3,sp(:,:,1,:), &
                   kthp(:,:,1,1),phip(:,:,1,1),betap(:,:,1,:))
           case (3)
              call setup_T_op_3d(lo,hi,ng_1,ng_3,sp(:,:,:,:), &
                   kthp(:,:,:,1),phip(:,:,:,1),betap(:,:,:,:))
           end select
        end do
     enddo
     
     ! applyop
     call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
          stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
     
     ! scale residual by sigma/rho and add to thermal
     do n=1,nlevs
        call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
     enddo

  else
     ! load h into phi
     ! setup beta = -k_{th}/c_p for h term
     do n=1,nlevs
        ng_1 = phi(n)%ng
        ng_3 = s(n)%ng
        
        do i=1,s(n)%nboxes
           if (multifab_remote(s(n),i)) cycle
           sp => dataptr(s(n),i)
           kthovercpp => dataptr(kthovercp(n),i)
           phip => dataptr(phi(n),i)
           betap => dataptr(beta(n),i)
           lo =  lwb(get_box(s(n), i))
           hi =  upb(get_box(s(n), i))
           select case (dm)
           case (2)
              call setup_h_op_2d(lo,hi,ng_1,ng_3,sp(:,:,1,:), &
                   kthovercpp(:,:,1,1),phip(:,:,1,1),betap(:,:,1,:))
           case (3)
              call setup_h_op_3d(lo,hi,ng_1,ng_3,sp(:,:,:,:), &
                   kthovercpp(:,:,:,1),phip(:,:,:,1),betap(:,:,:,:))
           end select
        end do
     enddo
     
     ! applyop
     call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
          stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
     
     ! scale residual by sigma/rho and add to thermal
     do n=1,nlevs
        call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
     enddo
     
     ! load X_k into phi
     ! compute \xi_k
     ! setup beta = \xi_k * (k_{th}/c_p)
     do k=1,nspec
        do n=1,nlevs
           ng_1 = phi(n)%ng
           ng_3 = s(n)%ng
           
           do i=1,s(n)%nboxes
              if (multifab_remote(s(n),i)) cycle
              sp         => dataptr(s(n),i)
              kthovercpp => dataptr(kthovercp(n),i)
              phip       => dataptr(phi(n),i)
              betap      => dataptr(beta(n),i)
              xikp       => dataptr(xik(n),i)
              lo =  lwb(get_box(s(n), i))
              hi =  upb(get_box(s(n), i))
              select case (dm)
              case (2)
                 call setup_Xk_op_2d(k,lo,hi,dx(n,:),ng_1,ng_3,sp(:,:,1,:), &
                      kthovercpp(:,:,1,1),phip(:,:,1,1),betap(:,:,1,:), &
                      xikp(:,:,1,:))
              case (3)
                 call setup_Xk_op_3d(k,lo,hi,dx(n,:),ng_1,ng_3,sp(:,:,:,:), &
                      kthovercpp(:,:,:,1),phip(:,:,:,1),betap(:,:,:,:), &
                      xikp(:,:,:,:))
              end select
           end do
        enddo
        
        ! applyop
        call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower, &
             dm+spec_comp+k-1,stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
        
        ! scale residual by sigma/rho and add to thermal
        do n=1,nlevs
           call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
           call multifab_fill_boundary(thermal(n))
        enddo
     enddo
  endif

  !  call fabio_ml_multifab_write_d(thermal,mla%mba%rr(:,1),"a_thermal")

  ! Deallocate memory
  do n = 1,nlevs
     call destroy(phi(n))
     call destroy(alpha(n))
     call destroy(beta(n))
     call destroy(xik(n))
     call destroy(kth(n))
     call destroy(kthovercp(n))
     call destroy(resid(n))
  enddo

  deallocate(phi,alpha,beta,xik,kth,kthovercp,resid)

end subroutine make_explicit_thermal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create kth
! create kthovercp
subroutine make_coeffs_2d(lo,hi,dx,ng_1,ng_3,p0,s,kth,kthovercp,xik)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dx(:)
  integer        , intent(in   ) :: ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,:)
  real(kind=dp_t), intent(inout) :: kth(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: kthovercp(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: xik(lo(1)-ng_1:,lo(2)-ng_1:,:)

  ! local
  integer :: i,j,n


  do j=lo(2)-1,hi(2)+1
     do i=lo(1)-1,hi(1)+1

        den_row(1) = s(i,j,rho_comp)
        temp_row(1) = s(i,j,temp_comp)
        xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

        ! dens, temp, and xmass are inputs
        input_flag = 1
        do_diag = .false.

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

        kth(i,j) = conduct_row(1)
        kthovercp(i,j) = conduct_row(1)/cp_row(1)

        do n=1,nspec
           xik(i,j,n) = dhdX_row(1,n)
           if(use_big_h) then
              xik(i,j,n) = xik(i,j,n) + ebin(n)
           endif
        enddo

     enddo
  enddo

end subroutine make_coeffs_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create kth
! create kthovercp
subroutine make_coeffs_3d(lo,hi,dx,ng_1,ng_3,p0,s,kth,kthovercp,xik)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dx(:)
  integer        , intent(in   ) :: ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:,:)
  real(kind=dp_t), intent(inout) :: kth(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: kthovercp(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: xik(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:,:)

  ! local
  integer :: i,j,k,n
  real(kind=dp_t), allocatable :: p0_cart(:,:,:)

  if (spherical .eq. 1) then
     allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
     call fill_3d_data(p0_cart,p0,lo,hi,dx,0)
  end if

  do k=lo(3)-1,hi(3)+1
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           
           den_row(1) = s(i,j,k,rho_comp)
           temp_row(1) = s(i,j,k,temp_comp)
           xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

           ! dens, temp, and xmass are inputs
           input_flag = 1
           do_diag = .false.
        
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

           kth(i,j,k) = conduct_row(1)
           kthovercp(i,j,k) = conduct_row(1)/cp_row(1)

           do n=1,nspec
              xik(i,j,k,n) = dhdX_row(1,n)
              if(use_big_h) then
                 xik(i,j,k,n) = xik(i,j,k,n) + ebin(n)
              endif
           enddo

        enddo
     enddo
  enddo

end subroutine make_coeffs_3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load T into phi
! setup beta for T term
subroutine setup_T_op_2d(lo,hi,ng_1,ng_3,s,kth,phi,beta)

  integer        , intent(in   ) :: lo(:),hi(:)
  integer        , intent(in   ) :: ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,:)
  real(kind=dp_t), intent(in   ) :: kth(-ng_1:,-ng_1:)
  real(kind=dp_t), intent(inout) :: phi(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: beta(-ng_1:,-ng_1:,:)

  integer :: i,j
  integer :: nx,ny

  nx = size(beta,dim=1) - 2*ng_1
  ny = size(beta,dim=2) - 2*ng_1

  do j=lo(2)-1,hi(2)+1
     do i=lo(1)-1,hi(1)+1
        phi(i,j) = s(i,j,temp_comp)
     enddo
  enddo

  ! set beta
  do j = 0,ny-1
     do i = 0,nx
        beta(i,j,1) = -TWO*(kth(i,j)*kth(i-1,j))/(kth(i,j) + kth(i-1,j))
     end do
  end do
  
  do j = 0,ny
     do i = 0,nx-1
        beta(i,j,2) = -TWO*(kth(i,j)*kth(i,j-1))/(kth(i,j) + kth(i,j-1))
     end do
  end do

end subroutine setup_T_op_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load T into phi
! setup beta for T term
subroutine setup_T_op_3d(lo,hi,ng_1,ng_3,s,kth,phi,beta)

  integer        , intent(in   ) :: lo(:),hi(:)
  integer        , intent(in   ) :: ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:,:)
  real(kind=dp_t), intent(in   ) :: kth(-ng_1:,-ng_1:,-ng_1:)
  real(kind=dp_t), intent(inout) :: phi(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: beta(-ng_1:,-ng_1:,-ng_1:,:)

  integer :: i,j,k
  integer :: nx,ny,nz

  nx = size(beta,dim=1) - 2*ng_1
  ny = size(beta,dim=2) - 2*ng_1
  nz = size(beta,dim=3) - 2*ng_1

  do k=lo(3)-1,hi(3)+1
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           phi(i,j,k) = s(i,j,k,temp_comp)
        enddo
     enddo
  enddo
  
  ! set beta
  do k = 0,nz-1
     do j = 0,ny-1
        do i = 0,nx
           beta(i,j,k,1) = -TWO*(kth(i,j,k)*kth(i-1,j,k))/(kth(i,j,k) &
                + kth(i-1,j,k))
        end do
     end do
  end do
  
  do k = 0,nz-1
     do j = 0,ny
        do i = 0,nx-1
           beta(i,j,k,2) = -TWO*(kth(i,j,k)*kth(i,j-1,k))/(kth(i,j,k) &
                + kth(i,j-1,k))
        end do
     end do
  end do
  
  do k = 0,nz
     do j = 0,ny-1
        do i = 0,nx-1
           beta(i,j,k,3) = -TWO*(kth(i,j,k)*kth(i,j,k-1))/(kth(i,j,k) &
                + kth(i,j,k-1))
        end do
     end do
  end do

end subroutine setup_T_op_3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load h into phi
! setup beta for h term
subroutine setup_h_op_2d(lo,hi,ng_1,ng_3,s,kthovercp,phi,beta)

  integer        , intent(in   ) :: lo(:),hi(:)
  integer        , intent(in   ) :: ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,:)
  real(kind=dp_t), intent(in   ) :: kthovercp(-ng_1:,-ng_1:)
  real(kind=dp_t), intent(inout) :: phi(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: beta(-ng_1:,-ng_1:,:)

  integer :: i,j
  integer :: nx,ny

  nx = size(beta,dim=1) - 2*ng_1
  ny = size(beta,dim=2) - 2*ng_1

  do j=lo(2)-1,hi(2)+1
     do i=lo(1)-1,hi(1)+1
        phi(i,j) = s(i,j,rhoh_comp)/s(i,j,rho_comp)
     enddo
  enddo

  ! set beta
  do j = 0,ny-1
     do i = 0,nx
        beta(i,j,1) = -TWO*(kthovercp(i,j)*kthovercp(i-1,j))/(kthovercp(i,j) &
             + kthovercp(i-1,j))
     end do
  end do
  
  do j = 0,ny
     do i = 0,nx-1
        beta(i,j,2) = -TWO*(kthovercp(i,j)*kthovercp(i,j-1))/(kthovercp(i,j) &
             + kthovercp(i,j-1))
     end do
  end do

end subroutine setup_h_op_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load h into phi
! setup beta for h term
subroutine setup_h_op_3d(lo,hi,ng_1,ng_3,s,kthovercp,phi,beta)

  integer        , intent(in   ) :: lo(:),hi(:)
  integer        , intent(in   ) :: ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:,:)
  real(kind=dp_t), intent(in   ) :: kthovercp(-ng_1:,-ng_1:,-ng_1:)
  real(kind=dp_t), intent(inout) :: phi(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: beta(-ng_1:,-ng_1:,-ng_1:,:)

  integer :: i,j,k
  integer :: nx,ny,nz

  nx = size(beta,dim=1) - 2*ng_1
  ny = size(beta,dim=2) - 2*ng_1
  nz = size(beta,dim=3) - 2*ng_1

  do k=lo(3)-1,hi(3)+1
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           phi(i,j,k) = s(i,j,k,rhoh_comp)/s(i,j,k,rho_comp)
        enddo
     enddo
  enddo
  
  ! set beta
  do k = 0,nz-1
     do j = 0,ny-1
        do i = 0,nx
           beta(i,j,k,1) = -TWO*(kthovercp(i,j,k)*kthovercp(i-1,j,k)) &
                /(kthovercp(i,j,k) + kthovercp(i-1,j,k))
        end do
     end do
  end do
  
  do k = 0,nz-1
     do j = 0,ny
        do i = 0,nx-1
           beta(i,j,k,2) = -TWO*(kthovercp(i,j,k)*kthovercp(i,j-1,k)) &
                /(kthovercp(i,j,k) + kthovercp(i,j-1,k))
        end do
     end do
  end do
  
  do k = 0,nz
     do j = 0,ny-1
        do i = 0,nx-1
           beta(i,j,k,3) = -TWO*(kthovercp(i,j,k) + kthovercp(i,j,k-1)) &
                /(kthovercp(i,j,k) + kthovercp(i,j,k-1))
        end do
     end do
  end do

end subroutine setup_h_op_3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load X_k into phi
! compute \xi_k
! setup beta = \xi_k * (k_{th}/c_p)
subroutine setup_Xk_op_2d(spec,lo,hi,dx,ng_1,ng_3,s,kthovercp,phi,beta,xik)

  integer        , intent(in   ) :: spec,lo(:),hi(:)
  real(dp_t)     , intent(in   ) :: dx(:)
  integer        , intent(in   ) :: ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,:)
  real(kind=dp_t), intent(in   ) :: kthovercp(-ng_1:,-ng_1:)
  real(kind=dp_t), intent(inout) :: phi(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: beta(-ng_1:,-ng_1:,:)
  real(kind=dp_t), intent(in   ) :: xik(-ng_1:,-ng_1:,:)

  integer :: i,j
  integer :: nx,ny

  nx = size(beta,dim=1) - 2*ng_1
  ny = size(beta,dim=2) - 2*ng_1

  ! Load X_k into phi
  do j=lo(2)-1,hi(2)+1
     do i=lo(1)-1,hi(1)+1
        phi(i,j) = s(i,j,spec_comp+spec-1)/s(i,j,rho_comp)
     enddo
  enddo

  ! set beta
  do j = 0,ny-1
     do i = 0,nx
        beta(i,j,1) = TWO*(xik(i,j,spec)*kthovercp(i,j)*xik(i-1,j,spec)*kthovercp(i-1,j)) &
             /(xik(i,j,spec)*kthovercp(i,j) + xik(i-1,j,spec)*kthovercp(i-1,j))
     end do
  end do
  
  do j = 0,ny
     do i = 0,nx-1
        beta(i,j,2) = TWO*(xik(i,j,spec)*kthovercp(i,j)*xik(i,j-1,spec)*kthovercp(i,j-1)) &
             /(xik(i,j,spec)*kthovercp(i,j) + xik(i,j-1,spec)*kthovercp(i,j-1))
     end do
  end do

end subroutine setup_Xk_op_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load X_k into phi
! compute \xi_k
! setup beta = \xi_k * (k_{th}/c_p)
subroutine setup_Xk_op_3d(spec,lo,hi,dx,ng_1,ng_3,s,kthovercp,phi,beta,xik)

  integer        , intent(in   ) :: spec,lo(:),hi(:)
  real(dp_t)     , intent(in   ) :: dx(:)
  integer        , intent(in   ) :: ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:,:)
  real(kind=dp_t), intent(in   ) :: kthovercp(-ng_1:,-ng_1:,-ng_1:)
  real(kind=dp_t), intent(inout) :: phi(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: beta(-ng_1:,-ng_1:,-ng_1:,:)
  real(kind=dp_t), intent(in   ) :: xik(-ng_1:,-ng_1:,-ng_1:,:)

  integer :: i,j,k
  integer :: nx,ny,nz

  nx = size(beta,dim=1) - 2*ng_1
  ny = size(beta,dim=2) - 2*ng_1
  nz = size(beta,dim=3) - 2*ng_1

  ! Load X_k into phi
  do k=lo(3)-1,hi(3)+1
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
            phi(i,j,k) = s(i,j,k,spec_comp+spec-1)/s(i,j,k,rho_comp)
        enddo
     enddo
  enddo
  
  ! set beta
  do k = 0,nz-1
     do j = 0,ny-1
        do i = 0,nx
           beta(i,j,k,1) = TWO*(xik(i,j,k,spec)*kthovercp(i,j,k)*xik(i-1,j,k,spec)*kthovercp(i-1,j,k)) &
                /(xik(i,j,k,spec)*kthovercp(i,j,k) + xik(i-1,j,k,spec)*kthovercp(i-1,j,k))
        end do
     end do
  end do
  
  do k = 0,nz-1
     do j = 0,ny
        do i = 0,nx-1
           beta(i,j,k,2) = TWO*(xik(i,j,k,spec)*kthovercp(i,j,k)*xik(i,j-1,k,spec)*kthovercp(i,j-1,k)) &
                /(xik(i,j,k,spec)*kthovercp(i,j,k) + xik(i,j-1,k,spec)*kthovercp(i,j-1,k))
        end do
     end do
  end do
  
  do k = 0,nz
     do j = 0,ny-1
        do i = 0,nx-1
           beta(i,j,k,3) = TWO*(xik(i,j,k,spec)*kthovercp(i,j,k)*xik(i,j,k-1,spec)*kthovercp(i,j,k-1)) &
                /(xik(i,j,k,spec)*kthovercp(i,j,k) + xik(i,j,k-1,spec)*kthovercp(i,j,k-1))
        end do
     end do
  end do

end subroutine setup_Xk_op_3d


end module make_explicit_thermal_module
