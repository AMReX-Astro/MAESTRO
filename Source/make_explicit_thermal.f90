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
  use thermal_conduct_module

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
  type(multifab), allocatable :: phi(:),alpha(:),beta(:),Xkcoeff(:)
  type(multifab), allocatable :: Tcoeff(:),hcoeff(:),pcoeff(:),resid(:)
  integer                     :: i,k,n,nlevs,dm,stencil_order
  integer                     :: lo(s(1)%dim),hi(s(1)%dim)
  real(kind=dp_t), pointer    :: thermalp(:,:,:,:),sp(:,:,:,:)
  real(kind=dp_t), pointer    :: phip(:,:,:,:),betap(:,:,:,:)
  real(kind=dp_t), pointer    :: Xkcoeffp(:,:,:,:)
  real(kind=dp_t), pointer    :: Tcoeffp(:,:,:,:),hcoeffp(:,:,:,:)
  real(kind=dp_t), pointer    :: pcoeffp(:,:,:,:),residp(:,:,:,:)

  type(bc_level) ::  bc

  nlevs = mla%nlevel
  dm      = mla%dim
  stencil_order = 2

  allocate(phi(nlevs),alpha(nlevs),beta(nlevs),Xkcoeff(nlevs))
  allocate(Tcoeff(nlevs),hcoeff(nlevs),pcoeff(nlevs),resid(nlevs))

  do n = 1,nlevs
     call multifab_build(         phi(n), mla%la(n), 1, 1)
     call multifab_build(       alpha(n), mla%la(n), 1, 1)
     call multifab_build(        beta(n), mla%la(n),dm, 1)
     call multifab_build(     Xkcoeff(n), mla%la(n),nspec, 1)
     call multifab_build(      Tcoeff(n), mla%la(n), 1, 1)
     call multifab_build(      hcoeff(n), mla%la(n), 1, 1)
     call multifab_build(      pcoeff(n), mla%la(n), 1, 1)
     call multifab_build(       resid(n), mla%la(n), 1, 0)

     call setval(         phi(n), ZERO, all=.true.)
     call setval(       alpha(n), ZERO, all=.true.)
     call setval(        beta(n), ZERO, all=.true.)
     call setval(     Xkcoeff(n), ZERO, all=.true.)
     call setval(      Tcoeff(n), ZERO, all=.true.)
     call setval(      hcoeff(n), ZERO, all=.true.)
     call setval(      pcoeff(n), ZERO, all=.true.)
     call setval(       resid(n), ZERO, all=.true.)
     call setval(     thermal(n), ZERO, all=.true.)
  end do

  ! create Tcoeff = -kth, hcoeff = -kth/cp, Xkcoeff = xik*kth/cp, pcoeff = hp*kth/cp
  do n=1,nlevs
     do i=1,s(n)%nboxes
        if (multifab_remote(s(n),i)) cycle
        sp            => dataptr(s(n),i)
        Tcoeffp       => dataptr(Tcoeff(n),i)
        hcoeffp       => dataptr(hcoeff(n),i)
        Xkcoeffp      => dataptr(Xkcoeff(n),i)
        pcoeffp       => dataptr(pcoeff(n),i)
        lo =  lwb(get_box(s(n), i))
        hi =  upb(get_box(s(n), i))
        select case (dm)
        case (2)
           call make_coeffs_2d(lo,hi,dx(n,:),p0,sp(:,:,1,:), &
                               Tcoeffp(:,:,1,1),hcoeffp(:,:,1,1), &
                               Xkcoeffp(:,:,1,:),pcoeffp(:,:,1,1))
        case (3)
           call make_coeffs_3d(lo,hi,dx(n,:),p0,sp(:,:,:,:), &
                               Tcoeffp(:,:,:,1),hcoeffp(:,:,:,1), &
                               Xkcoeffp(:,:,:,:),pcoeffp(:,:,:,1))
        end select
     end do
  enddo

  if(temperature_diffusion) then
     ! load T into phi
     do n=1,nlevs
        call multifab_copy_c(phi(n),1,s(n),temp_comp,1,1)
     enddo

     ! setup beta = Tcoeff on faces
     do n=1,nlevs
        do i=1,s(n)%nboxes
           if (multifab_remote(s(n),i)) cycle
           Tcoeffp => dataptr(Tcoeff(n),i)
           betap   => dataptr(beta(n),i)
           lo =  lwb(get_box(s(n), i))
           hi =  upb(get_box(s(n), i))
           select case (dm)
           case (2)
              call put_beta_on_faces_2d(lo,hi,Tcoeffp(:,:,1,1),betap(:,:,1,:))
           case (3)
              call put_beta_on_faces_3d(lo,hi,Tcoeffp(:,:,:,1),betap(:,:,:,:))
           end select
        end do
     enddo
     
     ! applyop
     call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
          stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
     
     ! scale residual by sigma/rho and add to thermal
     do n=1,nlevs
        call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
        call multifab_fill_boundary(thermal(n))
     enddo

  else
     ! load h into phi
     do n=1,nlevs
        call multifab_copy_c(phi(n),1,s(n),rhoh_comp,1,1)
        call multifab_div_div_c(phi(n),1,s(n),rho_comp,1,1)
     enddo

     ! setup beta = hcoeff on faces
     do n=1,nlevs
        do i=1,s(n)%nboxes
           if (multifab_remote(s(n),i)) cycle
           hcoeffp    => dataptr(hcoeff(n),i)
           betap      => dataptr(beta(n),i)
           lo =  lwb(get_box(s(n), i))
           hi =  upb(get_box(s(n), i))
           select case (dm)
           case (2)
              call put_beta_on_faces_2d(lo,hi,hcoeffp(:,:,1,1),betap(:,:,1,:))
           case (3)
              call put_beta_on_faces_3d(lo,hi,hcoeffp(:,:,:,1),betap(:,:,:,:))
           end select
        end do
     enddo

     ! applyop
     call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
          stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
     
     ! scale residual by sigma/rho and add to thermal
     do n=1,nlevs
        call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
        call multifab_fill_boundary(thermal(n))
     enddo
     
     ! loop over species
     do k=1,nspec
        ! load X_k into phi
        do n=1,nlevs
           call multifab_copy_c(phi(n),1,s(n),spec_comp+k-1,1,1)
           call multifab_div_div_c(phi(n),1,s(n),rho_comp,1,1)
        enddo

        ! setup beta = +Xkcoeff on faces
        do n=1,nlevs
           do i=1,s(n)%nboxes
              if (multifab_remote(s(n),i)) cycle
              Xkcoeffp   => dataptr(Xkcoeff(n),i)
              betap      => dataptr(beta(n),i)
              lo =  lwb(get_box(s(n), i))
              hi =  upb(get_box(s(n), i))
              select case (dm)
              case (2)
                 call put_beta_on_faces_2d(lo,hi,Xkcoeffp(:,:,1,k),betap(:,:,1,:))
              case (3)
                 call put_beta_on_faces_3d(lo,hi,Xkcoeffp(:,:,:,k),betap(:,:,:,:))
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

     ! load p0 into phi
     do n=1,nlevs
        do i=1,s(n)%nboxes
           if (multifab_remote(phi(n),i)) cycle
           phip => dataptr(phi(n),i)
           lo =  lwb(get_box(phi(n), i))
           hi =  upb(get_box(phi(n), i))
           select case (dm)
           case (2)
              call put_base_state_on_multifab_2d(lo,hi,p0,phip(:,:,1,1))
           case (3)
              call put_base_state_on_multifab_3d(lo,hi,p0,phip(:,:,:,1))
           end select
        end do
     enddo

     ! set the boundary conditions for pressure
     do n=1,nlevs
        call multifab_fill_boundary(phi(n))
        bc = the_bc_tower%bc_tower_array(n)
        do i = 1, phi(n)%nboxes
           if ( multifab_remote(phi(n), i) ) cycle
           phip => dataptr(phi(n), i)
           lo =  lwb(get_box(phi(n), i))
           hi =  upb(get_box(phi(n), i))
           select case (dm)
           case (2)
              call setbc_2d(phip(:,:,1,1), lo, 1, &
                   bc%adv_bc_level_array(i,:,:,neumann_comp), &
                   dx(n,:),neumann_comp)
           case (3)
              call setbc_3d(phip(:,:,:,1), lo, 1, &
                   bc%adv_bc_level_array(i,:,:,neumann_comp), &
                   dx(n,:),neumann_comp)
           end select
        enddo
     enddo

     ! setup beta = pcoeff on faces
!     do n=1,nlevs
!        do i=1,beta(n)%nboxes
!           if (multifab_remote(beta(n),i)) cycle
!           pcoeffp => dataptr(pcoeff(n),i)
!           betap   => dataptr(beta(n),i)
!           lo = lwb(get_box(beta(n), i))
!           hi = upb(get_box(beta(n), i))
!           select case (dm)
!           case (2)
!              call put_beta_on_faces_2d(lo,hi,pcoeffp(:,:,1,1),betap(:,:,1,:))
!           case (3)
!              call put_beta_on_faces_3d(lo,hi,pcoeffp(:,:,:,1),betap(:,:,:,:))
!           end select
!        end do
!     enddo
!
!     ! applyop
!     call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,neumann_comp, &
!          stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
!     
!     ! scale residual by sigma/rho and add to thermal
!     do n=1,nlevs
!        call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
!        call multifab_fill_boundary(thermal(n))
!     enddo
  endif

  ! Deallocate memory
  do n = 1,nlevs
     call destroy(phi(n))
     call destroy(alpha(n))
     call destroy(beta(n))
     call destroy(Xkcoeff(n))
     call destroy(Tcoeff(n))
     call destroy(hcoeff(n))
     call destroy(pcoeff(n))
     call destroy(resid(n))
  enddo

  deallocate(phi,alpha,beta,Xkcoeff,Tcoeff,hcoeff,pcoeff,resid)

end subroutine make_explicit_thermal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create Tcoeff = -kth, hcoeff = -kth/cp, Xkcoeff = xik*kth/cp, pcoeff = hp*kth/cp
subroutine make_coeffs_2d(lo,hi,dx,p0,s,Tcoeff,hcoeff,Xkcoeff,pcoeff)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dx(:)
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-3:,lo(2)-3:,:)
  real(kind=dp_t), intent(inout) :: Tcoeff(lo(1)-1:,lo(2)-1:)
  real(kind=dp_t), intent(inout) :: hcoeff(lo(1)-1:,lo(2)-1:)
  real(kind=dp_t), intent(inout) :: Xkcoeff(lo(1)-1:,lo(2)-1:,:)
  real(kind=dp_t), intent(inout) :: pcoeff(lo(1)-1:,lo(2)-1:)

  ! local
  integer :: i,j,n


  do j=lo(2)-1,hi(2)+1
     do i=lo(1)-1,hi(1)+1

        den_row(1) = s(i,j,rho_comp)
        temp_row(1) = s(i,j,temp_comp)
        xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

        ! dens, temp, and xmass are inputs
        do_diag = .false.

        call conducteos(eos_input_rt, den_row, temp_row, &
             npts, nspec, &
             xn_zone, aion, zion, &
             p_row, h_row, e_row, & 
             cv_row, cp_row, xne_row, eta_row, pele_row, &
             dpdt_row, dpdr_row, dedt_row, dedr_row, &
             dpdX_row, dhdX_row, &
             gam1_row, cs_row, s_row, &
             dsdt_row, dsdr_row, &
             do_diag, conduct_row)

        Tcoeff(i,j) = -conduct_row(1)
        hcoeff(i,j) = -conduct_row(1)/cp_row(1)
        pcoeff(i,j) = (conduct_row(1)/cp_row(1))* &
             ((1.0d0/den_row(1))* &
              (1.0d0-p_row(1)/(den_row(1)*dpdr_row(1)))+dedr_row(1)/dpdr_row(1))

        if(use_big_h) then
           do n=1,nspec
              Xkcoeff(i,j,n) = (conduct_row(1)/cp_row(1))*(dhdX_row(1,n) + ebin(n))
           enddo
        else
           do n=1,nspec
              Xkcoeff(i,j,n) = (conduct_row(1)/cp_row(1))*dhdX_row(1,n)
           enddo
        endif
     enddo
  enddo

end subroutine make_coeffs_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create Tcoeff = -kth, hcoeff = -kth/cp, Xkcoeff = xik*kth/cp, pcoeff = hp*kth/cp
subroutine make_coeffs_3d(lo,hi,dx,p0,s,Tcoeff,hcoeff,Xkcoeff,pcoeff)

  integer        , intent(in   ) :: lo(:),hi(:)
  real(dp_t)    ,  intent(in   ) :: dx(:)
  real(kind=dp_t), intent(in   ) :: p0(0:)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-3:,lo(2)-3:,lo(3)-3:,:)
  real(kind=dp_t), intent(inout) :: Tcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:)
  real(kind=dp_t), intent(inout) :: hcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:)
  real(kind=dp_t), intent(inout) :: Xkcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
  real(kind=dp_t), intent(inout) :: pcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:)

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
           do_diag = .false.
        
           call conducteos(eos_input_rt, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, & 
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                dsdt_row, dsdr_row, &
                do_diag, conduct_row)

           Tcoeff(i,j,k) = -conduct_row(1)
           hcoeff(i,j,k) = -conduct_row(1)/cp_row(1)
           pcoeff(i,j,k) = (conduct_row(1)/cp_row(1))* &
                ((1.0d0/den_row(1))* &
                 (1.0d0-p_row(1)/(den_row(1)*dpdr_row(1)))+dedr_row(1)/dpdr_row(1))

           if(use_big_h) then
              do n=1,nspec
                 Xkcoeff(i,j,k,n) = (conduct_row(1)/cp_row(1))*(dhdX_row(1,n) + ebin(n))
              enddo
           else
              do n=1,nspec
                 Xkcoeff(i,j,k,n) = (conduct_row(1)/cp_row(1))*dhdX_row(1,n)
              enddo
           endif
        enddo
     enddo
  enddo

end subroutine make_coeffs_3d





end module make_explicit_thermal_module
