module make_explicit_thermal_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bl_constants_module

  implicit none

  private

  public :: make_explicit_thermal, make_thermal_coeffs

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the quantity: thermal = del dot kappa grad T
!
!  if temp_diffusion_formulation = 1, then we compute this directly.
!  if temp_diffusion_formulation = 2, then we compute the algebraically
!     equivalent form with grad h + grad X_k + grad p_0 formulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_explicit_thermal(mla,dx,thermal,s,Tcoeff,hcoeff,Xkcoeff,pcoeff, &
                                   p0,the_bc_tower)

    use bc_module
    use bl_prof_module
    use cc_stencil_module
    use mac_applyop_module
    use network, only: nspec
    use ml_restriction_module, only : ml_cc_restriction
    use multifab_fill_ghost_module
    use bl_constants_module
    use variables, only: temp_comp, rho_comp, rhoh_comp, spec_comp, foextrap_comp
    use multifab_physbc_module
    use fill_3d_module
    use probin_module, only: temp_diffusion_formulation


    type(ml_layout), intent(inout) :: mla
    real(dp_t)     , intent(in   ) :: dx(:,:)
    type(multifab) , intent(inout) :: thermal(:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: Tcoeff(:)
    type(multifab) , intent(in   ) :: hcoeff(:)
    type(multifab) , intent(in   ) :: Xkcoeff(:)
    type(multifab) , intent(in   ) :: pcoeff(:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! Local
    type(multifab) :: phi(mla%nlevel),alpha(mla%nlevel),beta(mla%nlevel,mla%dim)
    type(multifab) :: resid(mla%nlevel)

    integer                     :: comp,i,n,stencil_order,dm,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_explicit_thermal")

    dm = mla%dim
    nlevs = mla%nlevel

    stencil_order = 2

    do n=1,nlevs
       call setval(thermal(n), ZERO, all=.true.)
    end do

    ! for either diffusion formulation, we will use mac_applyop to
    ! construct the form of the conduction term.  mac_applyop forms
    ! the generic quantity:
    !
    !   (alpha MINUS del dot beta grad) phi = RHS

    if (temp_diffusion_formulation .eq. 1) then

       ! compute del dot (Tcoeff grad T) directly

       do n=1,nlevs
          call multifab_build(phi(n), mla%la(n), 1,  1)
          do i = 1,dm
             call multifab_build_edge(beta(n,i),mla%la(n),1,1,i)
          end do
       end do

       do n=1,nlevs
          ! load T into phi
          call multifab_copy_c(phi(n),1,s(n),temp_comp,1,1)
       end do

       call put_data_on_faces(mla,Tcoeff,1,beta,.true.)

       do n=1,nlevs
          call multifab_build(alpha(n), mla%la(n), 1, 1)
          call multifab_build(resid(n), mla%la(n), 1, 0)
          call setval(alpha(n), ZERO, all=.true.)
       end do

       ! applyop to compute resid = del dot Tcoeff grad T
       call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
                        stencil_order,mla%mba%rr)
     
       do n=1,nlevs
          call destroy(phi(n))
          call destroy(alpha(n))
          do i = 1,dm
             call destroy(beta(n,i))
          end do
       end do

       ! add residual to thermal
       do n=1,nlevs
          do i = 1,dm
             call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
          end do
       enddo

       do n = 1,nlevs
          call destroy(resid(n))
       enddo

    else ! if (temp_diffusion_formulation .eq. 2) case

       ! compute thermal = del dot ( hcoeff grad h) +
       !             sum_k del dot (Xkcoeff grad X_k) +
       !                   del dot ( pcoeff grad p_0)
       
       do n=1,nlevs
          call multifab_build(phi(n),  mla%la(n), 1,  1)
          do i = 1,dm
             call multifab_build_edge(beta(n,i), mla%la(n), 1, 1, i)
          end do
       end do

       do n=1,nlevs
          ! load h into phi
          call multifab_copy_c(phi(n),1,s(n),rhoh_comp,1,1)
          call multifab_div_div_c(phi(n),1,s(n),rho_comp,1,1)
       end do

       call put_data_on_faces(mla,hcoeff,1,beta,.true.)

       do n=1,nlevs
          call multifab_build(alpha(n), mla%la(n), 1, 1)
          call multifab_build(resid(n), mla%la(n), 1, 0)
          call setval(alpha(n), ZERO, all=.true.)
       end do
       
       ! applyop to compute resid = del dot hcoeff grad h
       call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+rhoh_comp, &
                        stencil_order,mla%mba%rr)

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
          end do

          call put_data_on_faces(mla,Xkcoeff,comp,beta,.true.)
          
          ! applyop to compute resid = del dot Xkcoeff grad X_k
          call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,dm+spec_comp+comp-1, &
                           stencil_order,mla%mba%rr)
          
          ! add residual to thermal
          do n=1,nlevs
             call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
          enddo
       enddo ! end loop over species

       call put_1d_array_on_cart(p0,phi,foextrap_comp,.false.,.false., &
                                 dx,the_bc_tower%bc_tower_array,mla)       

       if (nlevs .eq. 1) then

          ! fill ghost cells for two adjacent grids at the same level
          ! this includes periodic domain boundary ghost cells
          call multifab_fill_boundary(phi(nlevs))

          ! fill non-periodic domain boundary ghost cells
          call multifab_physbc(phi(nlevs),1,foextrap_comp,1, &
                               the_bc_tower%bc_tower_array(nlevs))

       else

          do n=nlevs,2,-1

             ! we shouldn't need a call to ml_cc_restriction here as
             ! long as the coarse phi under fine cells is reasonably
             ! valued, the results of mac_applyop are identical

             ! fill level n ghost cells using interpolation from level
             ! n-1 data note that multifab_fill_boundary and
             ! multifab_physbc are called for both levels n-1 and n
             call multifab_fill_ghost_cells(phi(n),phi(n-1),1,mla%mba%rr(n-1,:), &
                                            the_bc_tower%bc_tower_array(n-1), &
                                            the_bc_tower%bc_tower_array(n), &
                                            1,foextrap_comp,1,fill_crse_input=.false.)
          end do

       end if

       call put_data_on_faces(mla,pcoeff,1,beta,.true.)

       ! applyop to compute resid = del dot pcoeff grad p0
       call mac_applyop(mla,resid,phi,alpha,beta,dx,the_bc_tower,foextrap_comp, &
                        stencil_order,mla%mba%rr)
       
       do n=1,nlevs
          call destroy(phi(n))
          call destroy(alpha(n))
          do i = 1,dm 
             call destroy(beta(n,i))
          end do
       end do

       ! add residual to thermal
       do n=1,nlevs
          call multifab_plus_plus_c(thermal(n),1,resid(n),1,1,0)
       enddo
       
       do n = 1,nlevs
          call destroy(resid(n))
       enddo

    endif ! end temp_diffusion_formulation logic
    
    call destroy(bpt)
    
  end subroutine make_explicit_thermal


  subroutine make_thermal_coeffs(s,Tcoeff,hcoeff,Xkcoeff,pcoeff)

    ! create the coefficients for grad{T}, grad{h}, grad{X_k}, and grad{p_0}
    ! for the thermal diffusion term in the enthalpy equation.  
    !
    ! note: we explicitly fill the ghostcells by looping over them directly
    ! in the _2d and _3d routines below.

    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: Tcoeff(:)
    type(multifab) , intent(inout) :: hcoeff(:)
    type(multifab) , intent(inout) :: Xkcoeff(:)
    type(multifab) , intent(inout) :: pcoeff(:)

    ! local
    integer :: n,i
    integer :: ng_s,ng_T,ng_h,ng_X,ng_p
    integer :: lo(get_dim(s(1))),hi(get_dim(s(1))),dm,nlevs

    real(kind=dp_t), pointer    :: sp(:,:,:,:)
    real(kind=dp_t), pointer    :: Tcoeffp(:,:,:,:),hcoeffp(:,:,:,:)
    real(kind=dp_t), pointer    :: Xkcoeffp(:,:,:,:),pcoeffp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    dm = get_dim(s(1))
    nlevs = size(s)

999 format('... Level ', i1, ' create thermal coeffs:')

    call build(bpt, "make_thermal_coeffs")

    ng_s = nghost(s(1))
    ng_T = nghost(Tcoeff(1))
    ng_h = nghost(hcoeff(1))
    ng_X = nghost(Xkcoeff(1))
    ng_p = nghost(pcoeff(1))

    ! create Tcoeff = -kth, 
    !        hcoeff = -kth/cp, 
    !        Xkcoeff = xik*kth/cp, 
    !        pcoeff = hp*kth/cp
    do n=1,nlevs
       if (parallel_IOProcessor()) write (6, 999) n

       do i=1,nboxes(s(n))
          if (multifab_remote(s(n),i)) cycle
          sp       => dataptr(s(n),i)
          Tcoeffp  => dataptr(Tcoeff(n),i)
          hcoeffp  => dataptr(hcoeff(n),i)
          Xkcoeffp => dataptr(Xkcoeff(n),i)
          pcoeffp  => dataptr(pcoeff(n),i)
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (1)
             call make_thermal_coeffs_1d(lo,hi,sp(:,1,1,:),ng_s,Tcoeffp(:,1,1,1),ng_T, &
                                         hcoeffp(:,1,1,1),ng_h,Xkcoeffp(:,1,1,:),ng_X, &
                                         pcoeffp(:,1,1,1),ng_p)
          case (2)
             call make_thermal_coeffs_2d(lo,hi,sp(:,:,1,:),ng_s,Tcoeffp(:,:,1,1),ng_T, &
                                         hcoeffp(:,:,1,1),ng_h,Xkcoeffp(:,:,1,:),ng_X, &
                                         pcoeffp(:,:,1,1),ng_p)
          case (3)
             call make_thermal_coeffs_3d(lo,hi,sp(:,:,:,:),ng_s,Tcoeffp(:,:,:,1),ng_T, &
                                         hcoeffp(:,:,:,1),ng_h,Xkcoeffp(:,:,:,:),ng_X, &
                                         pcoeffp(:,:,:,1),ng_p)
          end select
       end do
    enddo

    call destroy(bpt)

  end subroutine make_thermal_coeffs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create Tcoeff = -kth, 
!        hcoeff = -kth/cp, 
!       Xkcoeff = xik*kth/cp, 
!        pcoeff = hp*kth/cp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_thermal_coeffs_1d(lo,hi,s,ng_s,Tcoeff,ng_T,hcoeff,ng_h, &
                                    Xkcoeff,ng_X,pcoeff,ng_p)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use conductivity_module
    use network, only: nspec
    use probin_module, only: buoyancy_cutoff_factor, base_cutoff_density, limit_conductivity

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_T,ng_h,ng_X,ng_p
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::  Tcoeff(lo(1)-ng_T:)
    real(kind=dp_t), intent(inout) ::  hcoeff(lo(1)-ng_h:)
    real(kind=dp_t), intent(inout) :: Xkcoeff(lo(1)-ng_X:,:)
    real(kind=dp_t), intent(inout) ::  pcoeff(lo(1)-ng_p:)
    
    ! local
    integer :: i,comp    
    
    do i=lo(1)-1,hi(1)+1

       ! if we are outside the star, turn off the conductivity
       if (limit_conductivity .and. &
             s(i,rho_comp) < buoyancy_cutoff_factor*base_cutoff_density) then

           Tcoeff(i) = ZERO
           hcoeff(i) = ZERO
           pcoeff(i) = ZERO
           Xkcoeff(i,:) = ZERO

        else
          
          den_eos = s(i,rho_comp)
          temp_eos = s(i,temp_comp)
          xn_eos(:) = s(i,spec_comp:spec_comp+nspec-1)/den_eos

          ! dens, temp, and xmass are inputs
          call conducteos(eos_input_rt, den_eos, temp_eos, &
                          nspec, &
                          xn_eos, &
                          p_eos, h_eos, e_eos, & 
                          cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                          dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                          dpdX_eos, dhdX_eos, &
                          gam1_eos, cs_eos, s_eos, &
                          dsdt_eos, dsdr_eos, &
                          .false., conduct_eos)

          Tcoeff(i) = -conduct_eos
          hcoeff(i) = -conduct_eos/cp_eos
          pcoeff(i) = (conduct_eos/cp_eos)* &
               ((1.0d0/den_eos)* &
               (1.0d0-p_eos/(den_eos*dpdr_eos))+dedr_eos/dpdr_eos)

          do comp=1,nspec
             Xkcoeff(i,comp) = (conduct_eos/cp_eos)*dhdX_eos(comp)
          enddo

       endif

    enddo
    
  end subroutine make_thermal_coeffs_1d
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create Tcoeff = -kth, 
!        hcoeff = -kth/cp, 
!       Xkcoeff = xik*kth/cp, 
!        pcoeff = hp*kth/cp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_thermal_coeffs_2d(lo,hi,s,ng_s,Tcoeff,ng_T,hcoeff,ng_h, &
                                    Xkcoeff,ng_X,pcoeff,ng_p)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use conductivity_module
    use network, only: nspec
    use probin_module, only: buoyancy_cutoff_factor, base_cutoff_density, limit_conductivity

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_T,ng_h,ng_X,ng_p
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::  Tcoeff(lo(1)-ng_T:,lo(2)-ng_T:)
    real(kind=dp_t), intent(inout) ::  hcoeff(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(inout) :: Xkcoeff(lo(1)-ng_X:,lo(2)-ng_X:,:)
    real(kind=dp_t), intent(inout) ::  pcoeff(lo(1)-ng_p:,lo(2)-ng_p:)
    
    ! local
    integer :: i,j,comp    
    
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1

          ! if we are outside the star, turn off the conductivity
           if (limit_conductivity .and. &
                s(i,j,rho_comp) < buoyancy_cutoff_factor*base_cutoff_density) &
                then

              Tcoeff(i,j) = ZERO
              hcoeff(i,j) = ZERO
              pcoeff(i,j) = ZERO
              Xkcoeff(i,j,:) = ZERO

           else
          
             den_eos = s(i,j,rho_comp)
             temp_eos = s(i,j,temp_comp)
             xn_eos(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos

             ! dens, temp, and xmass are inputs
             call conducteos(eos_input_rt, den_eos, temp_eos, &
                             nspec, &
                             xn_eos, &
                             p_eos, h_eos, e_eos, & 
                             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                             dpdX_eos, dhdX_eos, &
                             gam1_eos, cs_eos, s_eos, &
                             dsdt_eos, dsdr_eos, &
                             .false., conduct_eos)

             Tcoeff(i,j) = -conduct_eos
             hcoeff(i,j) = -conduct_eos/cp_eos
             pcoeff(i,j) = (conduct_eos/cp_eos)* &
                  ((1.0d0/den_eos)* &
                  (1.0d0-p_eos/(den_eos*dpdr_eos))+dedr_eos/dpdr_eos)

             do comp=1,nspec
                Xkcoeff(i,j,comp) = (conduct_eos/cp_eos)*dhdX_eos(comp)
             enddo

          endif

       enddo
 enddo
    
  end subroutine make_thermal_coeffs_2d
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create Tcoeff = -kth, 
!        hcoeff = -kth/cp, 
!       Xkcoeff = xik*kth/cp, 
!        pcoeff = hp*kth/cp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_thermal_coeffs_3d(lo,hi,s,ng_s,Tcoeff,ng_T,hcoeff,ng_h, &
                                    Xkcoeff,ng_X,pcoeff,ng_p)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use conductivity_module
    use network, only: nspec
    use geometry, only: spherical
    use fill_3d_module
    use probin_module, only: buoyancy_cutoff_factor, base_cutoff_density, limit_conductivity
    
    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_T,ng_h,ng_X,ng_p
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::  Tcoeff(lo(1)-ng_T:,lo(2)-ng_T:,lo(3)-ng_T:)
    real(kind=dp_t), intent(inout) ::  hcoeff(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(inout) :: Xkcoeff(lo(1)-ng_X:,lo(2)-ng_X:,lo(3)-ng_X:,:)
    real(kind=dp_t), intent(inout) ::  pcoeff(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)

    ! local
    integer :: i,j,k,comp
    
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

              if (limit_conductivity .and. &
                   s(i,j,k,rho_comp) < &
                   buoyancy_cutoff_factor*base_cutoff_density) then
                
                 Tcoeff(i,j,k) = ZERO
                 hcoeff(i,j,k) = ZERO
                 pcoeff(i,j,k) = ZERO
                 Xkcoeff(i,j,k,:) = ZERO

              else
             
                den_eos = s(i,j,k,rho_comp)
                temp_eos = s(i,j,k,temp_comp)
                xn_eos(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos

                ! dens, temp, and xmass are inputs
                call conducteos(eos_input_rt, den_eos, temp_eos, &
                                nspec, &
                                xn_eos, &
                                p_eos, h_eos, e_eos, & 
                                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                                dpdX_eos, dhdX_eos, &
                                gam1_eos, cs_eos, s_eos, &
                                dsdt_eos, dsdr_eos, &
                                .false., conduct_eos)

                Tcoeff(i,j,k) = -conduct_eos
                hcoeff(i,j,k) = -conduct_eos/cp_eos
                pcoeff(i,j,k) = (conduct_eos/cp_eos)* &
                     ((1.0d0/den_eos)* &
                     (1.0d0-p_eos/(den_eos*dpdr_eos))+dedr_eos/dpdr_eos)

                do comp=1,nspec
                   Xkcoeff(i,j,k,comp) = (conduct_eos/cp_eos)*dhdX_eos(comp)
                enddo

             endif

          enddo
       enddo
    enddo
    
  end subroutine make_thermal_coeffs_3d

end module make_explicit_thermal_module
