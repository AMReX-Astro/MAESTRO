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
  use thermal_conduct_module
  use ml_restriction_module
  use multifab_fill_ghost_module
  use ml_layout_module
  use bl_constants_module
  use variables
  use probin_module, ONLY: use_big_h
  use geometry
  use multifab_physbc_module
  
  implicit none

  private
  public :: make_explicit_thermal

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute thermal = del dot kappa grad T (if temperature_diffusion)
! Otherwise, compute thermal with grad h + grad X_k + grad p_0 formulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_explicit_thermal(mla,dx,thermal,s,p0,mg_verbose,cg_verbose, &
                                   the_bc_tower,temperature_diffusion)

    type(ml_layout), intent(inout) :: mla
    real(dp_t)     , intent(in   ) :: dx(:,:)
    type(multifab) , intent(inout) :: thermal(:)
    type(multifab) , intent(in   ) :: s(:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    integer        , intent(in   ) :: mg_verbose,cg_verbose
    type(bc_tower) , intent(in   ) :: the_bc_tower
    logical        , intent(in   ) :: temperature_diffusion

    ! Local
    type(multifab), allocatable :: phi(:),alpha(:),beta(:),Xkcoeff(:)
    type(multifab), allocatable :: Tcoeff(:),hcoeff(:),pcoeff(:),resid(:)
    integer                     :: i,comp,n,nlevs,dm,stencil_order
    integer                     :: lo(s(1)%dim),hi(s(1)%dim)
    real(kind=dp_t), pointer    :: sp(:,:,:,:),phip(:,:,:,:)
    real(kind=dp_t), pointer    :: betap(:,:,:,:),Xkcoeffp(:,:,:,:)
    real(kind=dp_t), pointer    :: Tcoeffp(:,:,:,:),hcoeffp(:,:,:,:)
    real(kind=dp_t), pointer    :: pcoeffp(:,:,:,:),residp(:,:,:,:)
    type(bc_level)              :: bc

    nlevs = mla%nlevel
    dm = mla%dim
    stencil_order = 2

    allocate(phi(nlevs),alpha(nlevs),beta(nlevs),Xkcoeff(nlevs))
    allocate(Tcoeff(nlevs),hcoeff(nlevs),pcoeff(nlevs),resid(nlevs))
    
    do n=1,nlevs
       call multifab_build( phi(n),     mla%la(n), 1,     1)
       call multifab_build( alpha(n),   mla%la(n), 1,     1)
       call multifab_build( beta(n),    mla%la(n), dm,    1)
       call multifab_build( Xkcoeff(n), mla%la(n), nspec, 1)
       call multifab_build( Tcoeff(n),  mla%la(n), 1,     1)
       call multifab_build( hcoeff(n),  mla%la(n), 1,     1)
       call multifab_build( pcoeff(n),  mla%la(n), 1,     1)
       call multifab_build( resid(n),   mla%la(n), 1,     0)
       
       call setval( phi(n),     ZERO, all=.true.)
       call setval( alpha(n),   ZERO, all=.true.)
       call setval( beta(n),    ZERO, all=.true.)
       call setval( Xkcoeff(n), ZERO, all=.true.)
       call setval( Tcoeff(n),  ZERO, all=.true.)
       call setval( hcoeff(n),  ZERO, all=.true.)
       call setval( pcoeff(n),  ZERO, all=.true.)
       call setval( resid(n),   ZERO, all=.true.)
       call setval( thermal(n), ZERO, all=.true.)
    end do
    
    ! create Tcoeff = -kth, hcoeff = -kth/cp, Xkcoeff = xik*kth/cp, pcoeff = hp*kth/cp
    do n=1,nlevs
       do i=1,s(n)%nboxes
          if (multifab_remote(s(n),i)) cycle
          sp       => dataptr(s(n),i)
          Tcoeffp  => dataptr(Tcoeff(n),i)
          hcoeffp  => dataptr(hcoeff(n),i)
          Xkcoeffp => dataptr(Xkcoeff(n),i)
          pcoeffp  => dataptr(pcoeff(n),i)
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call make_coeffs_2d(lo,hi,dx(n,:),p0(n,:),sp(:,:,1,:), &
                                 Tcoeffp(:,:,1,1),hcoeffp(:,:,1,1), &
                                 Xkcoeffp(:,:,1,:),pcoeffp(:,:,1,1))
          case (3)
             call make_coeffs_3d(n,lo,hi,dx(n,:),p0(n,:),sp(:,:,:,:), &
                                 Tcoeffp(:,:,:,1),hcoeffp(:,:,:,1), &
                                 Xkcoeffp(:,:,:,:),pcoeffp(:,:,:,1))
          end select
       end do
    enddo

    if(temperature_diffusion) then

       do n=1,nlevs
          ! load T into phi
          call multifab_copy_c(phi(n),1,s(n),temp_comp,1,1)
       
          ! setup beta = Tcoeff on faces
          do i=1,s(n)%nboxes
             if (multifab_remote(s(n),i)) cycle
             Tcoeffp => dataptr(Tcoeff(n),i)
             betap   => dataptr(beta(n),i)
             lo = lwb(get_box(s(n),i))
             hi = upb(get_box(s(n),i))
             select case (dm)
             case (2)
                call put_beta_on_faces_2d(lo,hi,Tcoeffp(:,:,1,1),betap(:,:,1,:))
             case (3)
                call put_beta_on_faces_3d(lo,hi,Tcoeffp(:,:,:,1),betap(:,:,:,:))
             end select
          end do
       enddo ! end loop over levels
       
       ! applyop to compute resid = del dot Tcoeff grad T
       call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
                        stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
     
       ! add residual to thermal
       do n=1,nlevs
          call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
       enddo

    else ! the if(.not. temperature_diffusion) case
       
       do n=1,nlevs
          ! load h into phi
          call multifab_copy_c(phi(n),1,s(n),rhoh_comp,1,1)
          call multifab_div_div_c(phi(n),1,s(n),rho_comp,1,1)
                 
          ! setup beta = hcoeff on faces
          do i=1,s(n)%nboxes
             if (multifab_remote(s(n),i)) cycle
             hcoeffp => dataptr(hcoeff(n),i)
             betap   => dataptr(beta(n),i)
             lo =  lwb(get_box(s(n),i))
             hi =  upb(get_box(s(n),i))
             select case (dm)
             case (2)
                call put_beta_on_faces_2d(lo,hi,hcoeffp(:,:,1,1),betap(:,:,1,:))
             case (3)
                call put_beta_on_faces_3d(lo,hi,hcoeffp(:,:,:,1),betap(:,:,:,:))
             end select
          end do
       enddo ! end loop over levels
       
       ! applyop to compute resid = del dot hcoeff grad h
       call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
                        stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
       
       ! add residual to thermal
       do n=1,nlevs
          call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
       enddo
     
       ! loop over species
       do comp=1,nspec
          do n=1,nlevs
             ! load X_k into phi
             call multifab_copy_c(phi(n),1,s(n),spec_comp+comp-1,1,1)
             call multifab_div_div_c(phi(n),1,s(n),rho_comp,1,1)
             
             ! setup beta = Xkcoeff on faces
             do i=1,s(n)%nboxes
                if (multifab_remote(s(n),i)) cycle
                Xkcoeffp => dataptr(Xkcoeff(n),i)
                betap    => dataptr(beta(n),i)
                lo = lwb(get_box(s(n),i))
                hi = upb(get_box(s(n),i))
                select case (dm)
                case (2)
                   call put_beta_on_faces_2d(lo,hi,Xkcoeffp(:,:,1,comp),betap(:,:,1,:))
                case (3)
                   call put_beta_on_faces_3d(lo,hi,Xkcoeffp(:,:,:,comp),betap(:,:,:,:))
                end select
             end do
          enddo ! end loop over levels
          
          ! applyop to compute resid = del dot Xkcoeff grad X_k
          call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+spec_comp+comp-1, &
                           stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
          
          ! add residual to thermal
          do n=1,nlevs
             call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
          enddo
       enddo ! end loop over species
       
       ! load p0 into phi
       do n=1,nlevs
          do i=1,s(n)%nboxes
             if (multifab_remote(phi(n),i)) cycle
             phip => dataptr(phi(n),i)
             lo = lwb(get_box(phi(n),i))
             hi = upb(get_box(phi(n),i))
             select case (dm)
             case (2)
                call put_base_state_on_multifab_2d(lo,hi,p0(n,:),phip(:,:,1,1))
             case (3)
                call put_base_state_on_multifab_3d(lo,hi,p0(n,:),phip(:,:,:,1))
             end select
          end do
       enddo
       
       ! set the boundary conditions for pressure
       do n=1,nlevs
          call multifab_fill_boundary(phi(n))
          call multifab_physbc(phi(n),1,foextrap_comp,1,dx(n,:), &
                               the_bc_tower%bc_tower_array(n))
       enddo

       ! setup beta = pcoeff on faces
       do n=1,nlevs
          do i=1,beta(n)%nboxes
             if (multifab_remote(beta(n),i)) cycle
             pcoeffp => dataptr(pcoeff(n),i)
             betap   => dataptr(beta(n),i)
             lo = lwb(get_box(beta(n),i))
             hi = upb(get_box(beta(n),i))
             select case (dm)
             case (2)
                call put_beta_on_faces_2d(lo,hi,pcoeffp(:,:,1,1),betap(:,:,1,:))
             case (3)
                call put_beta_on_faces_3d(lo,hi,pcoeffp(:,:,:,1),betap(:,:,:,:))
             end select
          end do
       enddo
       
       ! applyop to compute resid = del dot pcoeff grad p0
       call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,foextrap_comp, &
                        stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
       
       ! add residual to thermal
       do n=1,nlevs
          call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
       enddo

    endif ! end if(temperature_diffusion) logic
    
    do n=1,nlevs
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary conditions
       call multifab_fill_boundary(thermal(n))

       ! A call to multifab_physbc for thermal (which may be used as a force in the 
       ! Godunov step) is not required.  Even though the Godunov step 
       ! references values of force outside of the domain, the boundary conditions
       ! ensure that the values of force outside of the domain do not influce the result.
    enddo

    do n=nlevs,2,-1
       ! make sure that coarse cells are the average of the fine cells covering it.
       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       call ml_cc_restriction(thermal(n-1),thermal(n),mla%mba%rr(n-1,:))

       ! fill fine ghost cells using interpolation from the underlying coarse data
       bc = the_bc_tower%bc_tower_array(n-1)
       call multifab_fill_ghost_cells(thermal(n),thermal(n-1), &
                                      1,mla%mba%rr(n-1,:), &
                                      the_bc_tower%bc_tower_array(n-1), &
                                      the_bc_tower%bc_tower_array(n  ), &
                                      1,foextrap_comp,1)
    enddo
    
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create Tcoeff = -kth, hcoeff = -kth/cp, Xkcoeff = xik*kth/cp, pcoeff = hp*kth/cp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    integer :: i,j,comp    
    
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          
          den_eos(1) = s(i,j,rho_comp)
          temp_eos(1) = s(i,j,temp_comp)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
          
          ! dens, temp, and xmass are inputs
          do_diag = .false.
          
          call conducteos(eos_input_rt, den_eos, temp_eos, &
                          npts, nspec, &
                          xn_eos, &
                          p_eos, h_eos, e_eos, & 
                          cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                          dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                          dpdX_eos, dhdX_eos, &
                          gam1_eos, cs_eos, s_eos, &
                          dsdt_eos, dsdr_eos, &
                          do_diag, conduct_eos)
          
          Tcoeff(i,j) = -conduct_eos(1)
          hcoeff(i,j) = -conduct_eos(1)/cp_eos(1)
          pcoeff(i,j) = (conduct_eos(1)/cp_eos(1))* &
               ((1.0d0/den_eos(1))* &
               (1.0d0-p_eos(1)/(den_eos(1)*dpdr_eos(1)))+dedr_eos(1)/dpdr_eos(1))
          
          if(use_big_h) then
             do comp=1,nspec
                Xkcoeff(i,j,comp) = (conduct_eos(1)/cp_eos(1))*(dhdX_eos(1,comp) &
                     + ebin(comp))
             enddo
          else
             do comp=1,nspec
                Xkcoeff(i,j,comp) = (conduct_eos(1)/cp_eos(1))*dhdX_eos(1,comp)
             enddo
          endif
       enddo
    enddo
    
  end subroutine make_coeffs_2d
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create Tcoeff = -kth, hcoeff = -kth/cp, Xkcoeff = xik*kth/cp, pcoeff = hp*kth/cp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_coeffs_3d(n,lo,hi,dx,p0,s,Tcoeff,hcoeff,Xkcoeff,pcoeff)
    
    integer        , intent(in   ) :: n,lo(:),hi(:)
    real(dp_t)    ,  intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: p0(0:)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-3:,lo(2)-3:,lo(3)-3:,:)
    real(kind=dp_t), intent(inout) :: Tcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(inout) :: hcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(inout) :: Xkcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
    real(kind=dp_t), intent(inout) :: pcoeff(lo(1)-1:,lo(2)-1:,lo(3)-1:)

    ! local
    integer :: i,j,k,comp
    real(kind=dp_t), allocatable :: p0_cart(:,:,:)
    
    if (spherical .eq. 1) then
       allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
       call fill_3d_data(n,p0_cart,p0,lo,hi,dx,0)
    end if
    
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             
             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = s(i,j,k,temp_comp)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
             
             ! dens, temp, and xmass are inputs
             do_diag = .false.
             
             call conducteos(eos_input_rt, den_eos, temp_eos, &
                             npts, nspec, &
                             xn_eos, &
                             p_eos, h_eos, e_eos, & 
                             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                             dpdX_eos, dhdX_eos, &
                             gam1_eos, cs_eos, s_eos, &
                             dsdt_eos, dsdr_eos, &
                             do_diag, conduct_eos)
             
             Tcoeff(i,j,k) = -conduct_eos(1)
             hcoeff(i,j,k) = -conduct_eos(1)/cp_eos(1)
             pcoeff(i,j,k) = (conduct_eos(1)/cp_eos(1))* &
                  ((1.0d0/den_eos(1))* &
                  (1.0d0-p_eos(1)/(den_eos(1)*dpdr_eos(1)))+dedr_eos(1)/dpdr_eos(1))
             
             if(use_big_h) then
                do comp=1,nspec
                   Xkcoeff(i,j,k,comp) = (conduct_eos(1)/cp_eos(1))*(dhdX_eos(1,comp) &
                        + ebin(comp))
                enddo
             else
                do comp=1,nspec
                   Xkcoeff(i,j,k,comp) = (conduct_eos(1)/cp_eos(1))*dhdX_eos(1,comp)
                enddo
             endif
          enddo
       enddo
    enddo
    
  end subroutine make_coeffs_3d
  
end module make_explicit_thermal_module
