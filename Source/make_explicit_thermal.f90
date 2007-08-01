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
  type(multifab), allocatable :: phi(:),alpha(:),beta(:)
  type(multifab), allocatable :: sigmaoverrho(:),ccbeta(:)
  type(multifab), allocatable :: kthovercp(:),temp(:)
  integer                     :: i,k,n,nlevs,dm,stencil_order,ng_0,ng_1,ng_3
  integer                     :: lo(s(1)%dim),hi(s(1)%dim)
  real(kind=dp_t), pointer    :: thermalp(:,:,:,:),sp(:,:,:,:)
  real(kind=dp_t), pointer    :: phip(:,:,:,:),ccbetap(:,:,:,:)
  real(kind=dp_t), pointer    :: betap(:,:,:,:),sigmaoverrhop(:,:,:,:)
  real(kind=dp_t), pointer    :: kthovercpp(:,:,:,:),tempp(:,:,:,:)

  nlevs = mla%nlevel
  dm      = mla%dim
  stencil_order = 2

  allocate(phi(nlevs),alpha(nlevs),beta(nlevs),sigmaoverrho(nlevs),ccbeta(nlevs))
  allocate(kthovercp(nlevs),temp(nlevs))

  do n = 1,nlevs
     call multifab_build(         phi(n), mla%la(n), 1, 1)
     call multifab_build(       alpha(n), mla%la(n), 1, 0)
     call multifab_build(        beta(n), mla%la(n),dm, 1)
     call multifab_build(      ccbeta(n), mla%la(n), 1, 1)
     call multifab_build(sigmaoverrho(n), mla%la(n), 1, 1)
     call multifab_build(   kthovercp(n), mla%la(n), 1, 1)
     call multifab_build(        temp(n), mla%la(n), 1, 0)

     call setval(         phi(n), ZERO, all=.true.)
     call setval(       alpha(n), ZERO, all=.true.)
     call setval(        beta(n), ZERO, all=.true.)
     call setval(sigmaoverrho(n), ZERO, all=.true.)
     call setval(      ccbeta(n), ZERO, all=.true.)
     call setval(   kthovercp(n), ZERO, all=.true.)
     call setval(        temp(n), ZERO, all=.true.)
  end do

  ! create sigmaoverrho = p_T(\rho^2 c_p p_\rho)
  ! create kthovercp
  ! create beta for h diffusion
  ! load h
  do n=1,nlevs
     ng_0 = temp(n)%ng
     ng_1 = phi(n)%ng
     ng_3 = s(n)%ng

     do i=1,s(n)%nboxes
        if (multifab_remote(s(n),i)) cycle
        sp => dataptr(s(n),i)
        betap => dataptr(beta(n),i)
        phip => dataptr(phi(n),i)
        sigmaoverrhop => dataptr(sigmaoverrho(n),i)
        ccbetap => dataptr(ccbeta(n),i)
        kthovercpp => dataptr(kthovercp(n),i)
        lo =  lwb(get_box(s(n), i))
        hi =  upb(get_box(s(n), i))
        select case (dm)
        case (2)
           call make_thermal_coeffs_2d(lo,hi,dt,dx(n,:),ng_0,ng_1,ng_3, &
                                       p0,sp(:,:,1,:),betap(:,:,1,:), &
                                       phip(:,:,1,1),sigmaoverrhop(:,:,1,1), &
                                       ccbetap(:,:,1,1),kthovercpp(:,:,1,1))
        case (3)
           call make_thermal_coeffs_3d(lo,hi,dt,dx(n,:),ng_0,ng_1,ng_3, &
                                       p0,sp(:,:,:,:),betap(:,:,:,:), &
                                       phip(:,:,:,1),sigmaoverrhop(:,:,:,1), &
                                       ccbetap(:,:,:,1),kthovercpp(:,:,:,1))
        end select
     end do
  enddo

  call fabio_ml_multifab_write_d(alpha,mla%mba%rr(:,1),"a_alpha")
  call fabio_ml_multifab_write_d(beta,mla%mba%rr(:,1),"a_beta")
  call fabio_ml_multifab_write_d(phi,mla%mba%rr(:,1),"a_h")

  ! applyop
  call mac_applyop(mla,temp,phi,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
                   stencil_order,mla%mba%rr,mg_verbose,cg_verbose)

  call fabio_ml_multifab_write_d(temp,mla%mba%rr(:,1),"a_resid")

  ! scale residual by sigma/rho
  do n=1,nlevs
     ng_0 = temp(n)%ng

     do i=1,s(n)%nboxes
        if (multifab_remote(s(n),i)) cycle
        tempp => dataptr(temp(n),i)
        sigmaoverrhop => dataptr(sigmaoverrho(n),i)
        lo =  lwb(get_box(s(n), i))
        hi =  upb(get_box(s(n), i))
        select case (dm)
        case (2)
           call scale_resid_2d(lo,hi,ng_0, &
                               tempp(:,:,1,1),sigmaoverrhop(:,:,1,1))
        case (3)
           call scale_resid_3d(lo,hi,ng_0, &
                               tempp(:,:,:,1),sigmaoverrhop(:,:,:,1))
        end select
     end do
  enddo

  call fabio_ml_multifab_write_d(temp,mla%mba%rr(:,1),"a_scaled_resid")
  call fabio_ml_multifab_write_d(thermal,mla%mba%rr(:,1),"a_thermal0")

  ! add resid to final answer
  do n=1,nlevs
     ng_0 = temp(n)%ng

     do i=1,s(n)%nboxes
        if (multifab_remote(s(n),i)) cycle
        tempp => dataptr(temp(n),i)
        thermalp => dataptr(thermal(n),i)
        lo =  lwb(get_box(s(n), i))
        hi =  upb(get_box(s(n), i))
        select case (dm)
        case (2)
           call add_resid_to_thermal_2d(lo,hi,ng_0, &
                                        tempp(:,:,1,1),thermalp(:,:,1,1))
        case (3)
           call add_resid_to_thermal_3d(lo,hi,ng_0, &
                                        tempp(:,:,:,1),thermalp(:,:,:,1))
        end select
     end do
  enddo

  call fabio_ml_multifab_write_d(thermal,mla%mba%rr(:,1),"a_thermal1")
  stop

  ! now for the species terms
  do k=1,nspec
     
     ! compute \xi_k and ccbeta = -k_{th}\xi_k/c_p (note the minus sign!)
     ! compute beta and load X_k
     do n=1,nlevs
        ng_0 = temp(n)%ng
        ng_1 = phi(n)%ng
        ng_3 = s(n)%ng
        
        do i=1,s(n)%nboxes
           if (multifab_remote(s(n),i)) cycle
           sp => dataptr(s(n),i)
           betap => dataptr(beta(n),i)
           phip => dataptr(phi(n),i)
           kthovercpp => dataptr(kthovercp(n),i)
           ccbetap => dataptr(ccbeta(n),i)
           lo =  lwb(get_box(s(n), i))
           hi =  upb(get_box(s(n), i))
           select case (dm)
           case (2)
              call make_species_coeffs_2d(k,lo,hi,dt,dx(n,:),ng_0,ng_1,ng_3, &
                   p0,sp(:,:,1,:),betap(:,:,1,:), &
                   phip(:,:,1,1),kthovercpp(:,:,1,1),ccbetap(:,:,1,1))
           case (3)
              call make_species_coeffs_3d(k,lo,hi,dt,dx(n,:),ng_0,ng_1,ng_3, &
                   p0,sp(:,:,:,:),betap(:,:,:,:), &
                   phip(:,:,:,1),kthovercpp(:,:,:,1),ccbetap(:,:,:,1))
           end select
        end do
     enddo

     ! applyop
     call mac_applyop(mla,temp,phi,alpha,beta,dx,the_bc_tower,dm+spec_comp+k-1, &
          stencil_order,mla%mba%rr,mg_verbose,cg_verbose)

     ! scale residual by sigma/rho
     do n=1,nlevs
        ng_0 = temp(n)%ng
        
        do i=1,s(n)%nboxes
           if (multifab_remote(s(n),i)) cycle
           tempp => dataptr(temp(n),i)
           sigmaoverrhop => dataptr(sigmaoverrho(n),i)
           lo =  lwb(get_box(s(n), i))
           hi =  upb(get_box(s(n), i))
           select case (dm)
           case (2)
              call scale_resid_2d(lo,hi,ng_0, &
                   tempp(:,:,1,1),sigmaoverrhop(:,:,1,1))
           case (3)
              call scale_resid_3d(lo,hi,ng_0, &
                   tempp(:,:,:,1),sigmaoverrhop(:,:,:,1))
           end select
        end do
     enddo

     ! add resid to final answer
     do n=1,nlevs
        ng_0 = temp(n)%ng
        
        do i=1,s(n)%nboxes
           if (multifab_remote(s(n),i)) cycle
           tempp => dataptr(temp(n),i)
           thermalp => dataptr(thermal(n),i)
           lo =  lwb(get_box(s(n), i))
           hi =  upb(get_box(s(n), i))
           select case (dm)
           case (2)
              call add_resid_to_thermal_2d(lo,hi,ng_0, &
                   tempp(:,:,1,1),thermalp(:,:,1,1))
           case (3)
              call add_resid_to_thermal_3d(lo,hi,ng_0, &
                   tempp(:,:,:,1),thermalp(:,:,:,1))
           end select
        end do
     enddo

  enddo

  ! Deallocate memory
  do n = 1,nlevs
     call destroy(phi(n))
     call destroy(alpha(n))
     call destroy(beta(n))
     call destroy(ccbeta(n))
     call destroy(sigmaoverrho(n))
     call destroy(kthovercp(n))
     call destroy(temp(n))
  enddo

  deallocate(phi,alpha,beta,ccbeta,sigmaoverrho,kthovercp,temp)

end subroutine make_explicit_thermal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine make_thermal_coeffs_2d(lo,hi,dt,dx,ng_0,ng_1,ng_3, &
                                  p0,s,beta,h,sigmaoverrho,ccbeta,kthovercp)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt,dx(:)
  integer        , intent(in   ) :: ng_0,ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,:)
  real(kind=dp_t), intent(inout) :: beta(lo(1)-ng_1:,lo(2)-ng_1:,:)
  real(kind=dp_t), intent(inout) :: h(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: sigmaoverrho(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: ccbeta(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: kthovercp(lo(1)-ng_1:,lo(2)-ng_1:)

  integer :: i,j
  integer :: nx,ny

! dens, pres, and xmass are inputs
  input_flag = 4
  do_diag = .false.

  nx = size(beta,dim=1) - 2*ng_1
  ny = size(beta,dim=2) - 2*ng_1

  ! Compute c_p^(2), k_th^2, betacc, and h
  do j=lo(2)-1,hi(2)+1
     do i=lo(1)-1,hi(1)+1

        den_row(1) = s(i,j,rho_comp)
        xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

        if(j .eq. lo(2)-1) then
           p_row(1) = p0(lo(2))
        else if(j .eq. hi(2)+1) then
           p_row(1) = p0(hi(2))
        else
           p_row(1) = p0(j)
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

        sigmaoverrho(i,j) = dpdt_row(1)/(den_row(1)**2*cp_row(1)*dpdr_row(1))
        kthovercp(i,j) = conduct_row(1)/cp_row(1)
        h(i,j) = s(i,j,rhoh_comp)/s(i,j,rho_comp)
        ccbeta(i,j) = -conduct_row(1)/cp_row(1)

     enddo
  enddo

  ! set beta
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

end subroutine make_thermal_coeffs_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine make_thermal_coeffs_3d(lo,hi,dt,dx,ng_0,ng_1,ng_3, &
                                  p0,s,beta,h,sigmaoverrho,ccbeta,kthovercp)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dt,dx(:)
  integer        , intent(in   ) :: ng_0,ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:,:)
  real(kind=dp_t), intent(inout) :: beta(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:,:)
  real(kind=dp_t), intent(inout) :: h(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: sigmaoverrho(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: ccbeta(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: kthovercp(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)

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

  ! Compute c_p^(2), k_th^2, betacc, and h
  do k=lo(3)-1,hi(3)+1
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           
           den_row(1) = s(i,j,k,rho_comp)
           xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

           if(spherical .eq. 0) then
              if(k .eq. lo(3)-1) then
                 p_row(1) = p0(lo(3))
              else if(k .eq. hi(3)+1) then
                 p_row(1) = p0(hi(3))
              else
                 p_row(1) = p0(k)
              endif
           else
              ! This still needs to be rewritten to handle out-of-domain cases!
              print*, "Computation of beta for spherical case not written!"
              stop
              
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
                do_diag, conduct_row)

           sigmaoverrho(i,j,k) = dpdt_row(1)/(den_row(1)**2*cp_row(1)*dpdr_row(1))
           kthovercp(i,j,k) = conduct_row(1)/cp_row(1)
           h(i,j,k) = s(i,j,k,rhoh_comp)/s(i,j,k,rho_comp)           
           ccbeta(i,j,k) = -ONE/cp_row(1)

        enddo
     enddo
  enddo
  
  ! set beta
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

end subroutine make_thermal_coeffs_3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine make_species_coeffs_2d(spec,lo,hi,dt,dx,ng_0,ng_1,ng_3, &
                                  p0,s,beta,Xk,kthovercp,ccbeta)

  integer        , intent(in   ) :: spec,lo(:),hi(:)
  real(dp_t)     , intent(in   ) :: dt,dx(:)
  integer        , intent(in   ) :: ng_0,ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,:)
  real(kind=dp_t), intent(inout) :: beta(lo(1)-ng_1:,lo(2)-ng_1:,:)
  real(kind=dp_t), intent(inout) :: Xk(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(in   ) :: kthovercp(lo(1)-ng_1:,lo(2)-ng_1:)
  real(kind=dp_t), intent(inout) :: ccbeta(lo(1)-ng_1:,lo(2)-ng_1:)

  integer :: i,j
  integer :: nx,ny

! dens, pres, and xmass are inputs
  input_flag = 4
  do_diag = .false.

  nx = size(beta,dim=1) - 2*ng_1
  ny = size(beta,dim=2) - 2*ng_1

  ! Compute xi_k
  do j=lo(2)-1,hi(2)+1
     do i=lo(1)-1,hi(1)+1

        den_row(1) = s(i,j,rho_comp)
        xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

        if(j .eq. lo(2)-1) then
           p_row(1) = p0(lo(2))
        else if(j .eq. hi(2)+1) then
           p_row(1) = p0(hi(2))
        else
           p_row(1) = p0(j)
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

        ccbeta(i,j) = kthovercp(i,j)*dhdX_row(1,spec)
        Xk(i,j) = s(i,j,spec_comp+spec-1)/s(i,j,rho_comp)

     enddo
  enddo

  ! set beta
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

end subroutine make_species_coeffs_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine make_species_coeffs_3d(spec,lo,hi,dt,dx,ng_0,ng_1,ng_3, &
                                  p0,s,beta,Xk,kthovercp,ccbeta)

  integer        , intent(in   ) :: spec,lo(:),hi(:)
  real(dp_t)     , intent(in   ) :: dt,dx(:)
  integer        , intent(in   ) :: ng_0,ng_1,ng_3
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:,:)
  real(kind=dp_t), intent(inout) :: beta(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:,:)
  real(kind=dp_t), intent(inout) :: Xk(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
  real(kind=dp_t), intent(inout) :: kthovercp(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
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

  ! Compute c_p^(2), k_th^2, betacc, and h
  do k=lo(3)-1,hi(3)+1
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           
           den_row(1) = s(i,j,k,rho_comp)
           xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

           if(spherical .eq. 0) then
              if(k .eq. lo(3)-1) then
                 p_row(1) = p0(lo(3))
              else if(k .eq. hi(3)+1) then
                 p_row(1) = p0(hi(3))
              else
                 p_row(1) = p0(k)
              endif
           else
              ! This still needs to be rewritten to handle out-of-domain cases!
              print*, "Computation of beta for spherical case not written!"
              stop
              
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
           
           ccbeta(i,j,k) = kthovercp(i,j,k)*dhdX_row(1,spec)
           Xk(i,j,k) = s(i,j,k,spec_comp+spec-1)/s(i,j,k,rho_comp)

        enddo
     enddo
  enddo
  
  ! set beta
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

end subroutine make_species_coeffs_3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine scale_resid_2d(lo,hi,ng,resid,sigmaoverrho)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(  out) :: resid(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(in   ) :: sigmaoverrho(lo(1)-ng:,lo(2)-ng:)

! Local
  integer :: i,j

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        resid(i,j) = sigmaoverrho(i,j)*resid(i,j)
     enddo
  enddo
  
end subroutine scale_resid_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine scale_resid_3d(lo,hi,ng,resid,sigmaoverrho)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(  out) :: resid(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(in   ) :: sigmaoverrho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

! Local
  integer :: i,j,k

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           resid(i,j,k) = sigmaoverrho(i,j,k)*resid(i,j,k)
        enddo
     enddo
  enddo

end subroutine scale_resid_3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine add_resid_to_thermal_2d(lo,hi,ng,resid,thermal)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(in   ) :: resid(lo(1)-ng:,lo(2)-ng:)
  real(kind=dp_t), intent(  out) :: thermal(lo(1)-ng:,lo(2)-ng:)

! Local
  integer :: i,j

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        thermal(i,j) = thermal(i,j) + resid(i,j)
     enddo
  enddo
  
end subroutine add_resid_to_thermal_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine add_resid_to_thermal_3d(lo,hi,ng,resid,thermal)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(in   ) :: resid(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  real(kind=dp_t), intent(  out) :: thermal(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

! Local
  integer :: i,j,k

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           thermal(i,j,k) = thermal(i,j,k) + resid(i,j,k)
        enddo
     enddo
  enddo

end subroutine add_resid_to_thermal_3d


end module make_explicit_thermal_module
