module mkflux_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: mkflux
  
contains

  subroutine mkflux(nlevs,sflux,sold,sedge,umac,w0,w0mac, &
                    rho0_old,rho0_edge_old,rho0_old_cart, &
                    rho0_new,rho0_edge_new,rho0_new_cart, &
                    rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                    rhoh0_new,rhoh0_edge_new,rhoh0_new_cart, &
                    startcomp,endcomp,mla)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical, dm
    use ml_restriction_module, only: ml_edge_restriction_c
    use variables, only: nscal

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(in   ) :: sold(:),sedge(:,:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:),rho0_edge_old(:,0:)
    type(multifab) , intent(in   ) :: rho0_old_cart(:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:),rho0_edge_new(:,0:)
    type(multifab) , intent(in   ) :: rho0_new_cart(:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(:,0:),rhoh0_edge_old(:,0:)
    type(multifab) , intent(in   ) :: rhoh0_old_cart(:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new(:,0:),rhoh0_edge_new(:,0:)
    type(multifab) , intent(in   ) :: rhoh0_new_cart(:)
    integer        , intent(in   ) :: startcomp,endcomp
    type(ml_layout), intent(inout) :: mla

    ! local    
    type(box) :: domain

    integer :: domlo(sold(1)%dim),domhi(sold(1)%dim)
    integer :: i,n
    integer :: lo(sold(1)%dim),hi(sold(1)%dim)
    integer :: ng_sf,ng_se,ng_um,ng_ro,ng_rn,ng_ho,ng_hn,ng_w0

    real(kind=dp_t), pointer :: sfxp(:,:,:,:)
    real(kind=dp_t), pointer :: sfyp(:,:,:,:)
    real(kind=dp_t), pointer :: sfzp(:,:,:,:)
    real(kind=dp_t), pointer :: sexp(:,:,:,:)
    real(kind=dp_t), pointer :: seyp(:,:,:,:)
    real(kind=dp_t), pointer :: sezp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    real(kind=dp_t), pointer :: rho0op(:,:,:,:)
    real(kind=dp_t), pointer :: rho0np(:,:,:,:)
    real(kind=dp_t), pointer :: rhoh0op(:,:,:,:)
    real(kind=dp_t), pointer :: rhoh0np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mkflux")

    ng_sf = sflux(1,1)%ng
    ng_se = sedge(1,1)%ng
    ng_um = umac(1,1)%ng
    ng_ro = rho0_old_cart(1)%ng
    ng_rn = rho0_new_cart(1)%ng
    ng_ho = rhoh0_old_cart(1)%ng
    ng_hn = rhoh0_new_cart(1)%ng
    ng_w0 = w0mac(1,1)%ng
    
    do n=1,nlevs

       domain = layout_get_pd(sold(n)%la)
       domlo = lwb(domain)
       domhi = upb(domain)

       do i=1, sold(n)%nboxes
          if ( multifab_remote(sold(n),i) ) cycle
          sfxp => dataptr(sflux(n,1),i)
          sfyp => dataptr(sflux(n,2),i)
          sexp => dataptr(sedge(n,1),i)
          seyp => dataptr(sedge(n,2),i)
          ump  => dataptr(umac(n,1),i)
          vmp  => dataptr(umac(n,2),i)
          lo = lwb(get_box(sold(n),i))
          hi = upb(get_box(sold(n),i))
          select case (dm)
          case (2)
             call mkflux_2d(sfxp(:,:,1,:), sfyp(:,:,1,:), ng_sf, &
                            sexp(:,:,1,:), seyp(:,:,1,:), ng_se, &
                            ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                            rho0_old(n,:), rho0_edge_old(n,:), &
                            rho0_new(n,:), rho0_edge_new(n,:), &
                            rhoh0_old(n,:), rhoh0_edge_old(n,:), &
                            rhoh0_new(n,:), rhoh0_edge_new(n,:), &
                            w0(n,:),startcomp,endcomp,lo,hi)
          case (3)
             sfzp => dataptr(sflux(n,3),i)
             sezp => dataptr(sedge(n,3),i)
             wmp  => dataptr(umac(n,3),i)
             if(spherical .eq. 0) then
                call mkflux_3d_cart(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), ng_sf, &
                                    sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), ng_se, &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                    rho0_old(n,:), rho0_edge_old(n,:), &
                                    rho0_new(n,:), rho0_edge_new(n,:), &
                                    rhoh0_old(n,:), rhoh0_edge_old(n,:), &
                                    rhoh0_new(n,:), rhoh0_edge_new(n,:), &
                                    w0(n,:),startcomp,endcomp,lo,hi)

             else
                rho0op => dataptr(rho0_old_cart(n), i)
                rho0np => dataptr(rho0_new_cart(n), i)
                rhoh0op => dataptr(rhoh0_old_cart(n), i)
                rhoh0np => dataptr(rhoh0_new_cart(n), i)
                w0xp => dataptr(w0mac(n,1),i)
                w0yp => dataptr(w0mac(n,2),i)
                w0zp => dataptr(w0mac(n,3),i)
                call mkflux_3d_sphr(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), ng_sf, &
                                    sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), ng_se, &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                    rho0op(:,:,:,1), ng_ro, rho0np(:,:,:,1), ng_rn, &
                                    rhoh0op(:,:,:,1), ng_ho, rhoh0np(:,:,:,1), ng_hn, &
                                    w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1), &
                                    ng_w0,startcomp,endcomp,lo,hi,domlo,domhi)
             endif
          end select
       end do

    end do ! end loop over levels

    ! synchronize fluxes at coarse-fine interface
    do n = nlevs,2,-1
       do i = 1, dm
          call ml_edge_restriction_c(sflux(n-1,i),1,sflux(n,i),1,mla%mba%rr(n-1,:),i,nscal)
       enddo
    enddo

    call destroy(bpt)
    
  end subroutine mkflux
  
  subroutine mkflux_2d(sfluxx,sfluxy,ng_sf,sedgex,sedgey,ng_se, &
                       umac,vmac,ng_um,rho0_old,rho0_edge_old,rho0_new,rho0_edge_new, &
                       rhoh0_old,rhoh0_edge_old,rhoh0_new,rhoh0_edge_new, &
                       w0,startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : rho_comp, rhoh_comp, trac_comp, ntrac
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:),ng_sf,ng_se,ng_um
    real(kind=dp_t), intent(inout) ::     sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sedgex(lo(1)-ng_se:,lo(2)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::     sedgey(lo(1)-ng_se:,lo(2)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:), rho0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:), rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(0:), rhoh0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new(0:), rhoh0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local
    integer         :: comp
    integer         :: i,j
    real(kind=dp_t) :: rho0_edge, rhoh0_edge
    logical :: test
    
    ! loop over components
    do comp = startcomp, endcomp

       ! test = T means the edge states are NOT in perturbational form
       test = ( (comp.eq.rhoh_comp).and. &
                     ( enthalpy_pred_type.eq.predict_h .or. &
                       enthalpy_pred_type.eq.predict_T_then_h ) ) &
         .or. ( (comp.ge.trac_comp).and.(comp.le.trac_comp+ntrac-1) )

       ! create x-fluxes
       if (test) then

          do j=lo(2),hi(2)
             
             rho0_edge = HALF*(rho0_old(j)+rho0_new(j))
             
             do i=lo(1),hi(1)+1
                
                sfluxx(i,j,comp) = umac(i,j)* &
                     (rho0_edge+sedgex(i,j,rho_comp))*sedgex(i,j,comp)
                
             end do
             
          end do
                
       else
              
          do j=lo(2),hi(2)
             
             rhoh0_edge = HALF*(rhoh0_old(j)+rhoh0_new(j))
             
             do i=lo(1),hi(1)+1
                
                sfluxx(i,j,comp) = umac(i,j)*(rhoh0_edge+sedgex(i,j,comp))
                
             end do
             
          end do
          
       end if
        
       ! create y-fluxes
       if (test) then
       
          do j=lo(2),hi(2)+1
             
             rho0_edge = HALF*(rho0_edge_old(j)+rho0_edge_new(j))
             
             do i=lo(1),hi(1)
                
                sfluxy(i,j,comp) = &
                     (vmac(i,j)+w0(j))*(rho0_edge+sedgey(i,j,rho_comp))*sedgey(i,j,comp)
                
             end do

          end do
          
       else
          
          do j=lo(2),hi(2)+1
             
             rhoh0_edge = HALF*(rhoh0_edge_old(j)+rhoh0_edge_new(j))
             
             do i=lo(1),hi(1)
                
                sfluxy(i,j,comp) = (vmac(i,j)+w0(j))*(sedgey(i,j,comp)+rhoh0_edge)
                
             end do
             
          end do

       end if
             
    end do

  end subroutine mkflux_2d

  subroutine mkflux_3d_cart(sfluxx,sfluxy,sfluxz,ng_sf,&
                            sedgex,sedgey,sedgez,ng_se,umac,vmac,wmac,ng_um, &
                            rho0_old,rho0_edge_old,rho0_new,rho0_edge_new, &
                            rhoh0_old,rhoh0_edge_old,rhoh0_new,rhoh0_edge_new, &
                            w0,startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : rho_comp, rhoh_comp, trac_comp, ntrac
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:),ng_sf,ng_se,ng_um
    real(kind=dp_t), intent(inout) ::     sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sfluxz(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sedgex(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::     sedgey(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::     sedgez(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:), rho0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:), rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(0:), rhoh0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new(0:), rhoh0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp

   ! local
    integer         :: comp
    integer         :: i,j,k
    real(kind=dp_t) :: rho0_edge, rhoh0_edge
    logical         :: test
    
    ! loop over components
    do comp = startcomp, endcomp

       ! test = T means the edge states are NOT in perturbational form
       test = ( (comp.eq.rhoh_comp).and. &
                     ( enthalpy_pred_type.eq.predict_h .or. &
                       enthalpy_pred_type.eq.predict_T_then_h ) ) &
         .or. ( (comp.ge.trac_comp).and.(comp.le.trac_comp+ntrac-1) )
       
       ! create x-fluxes and y-fluxes
       if (test) then

          do k=lo(3),hi(3)
             
             rho0_edge = HALF*(rho0_old(k)+rho0_new(k))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                   
                   ! sedgex is either h or X at edges
                   sfluxx(i,j,k,comp) = &
                        umac(i,j,k)*(rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp)
                   
                end do
             end do
             
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                   
                   ! sedgey is either h or X at edges
                   sfluxy(i,j,k,comp) = &
                        vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp)

                end do
             end do
             
          end do
                
       else
              
          do k=lo(3),hi(3)
             
             rhoh0_edge = HALF*(rhoh0_old(k)+rhoh0_new(k))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                
                   sfluxx(i,j,k,comp) = umac(i,j,k)*(rhoh0_edge+sedgex(i,j,k,comp))
                   
                end do
             end do
             
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                
                   sfluxy(i,j,k,comp) = vmac(i,j,k)*(rhoh0_edge+sedgey(i,j,k,comp))
                   
                end do
             end do
             
          end do
          
       end if
        
       ! create z-fluxes
       if (test) then
       
          do k=lo(3),hi(3)+1
             
             rho0_edge = HALF*(rho0_edge_old(k)+rho0_edge_new(k))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                
                   ! sedgez is either h or X at edges
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))* &
                        (rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp)
                
                end do
             end do

          end do
          
       else
          
          do k=lo(3),hi(3)+1
             
             rhoh0_edge = HALF*(rhoh0_edge_old(k)+rhoh0_edge_new(k))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))*(sedgez(i,j,k,comp)+rhoh0_edge)
                                
                end do
             end do

          end do

       end if
    end do

     
  end subroutine mkflux_3d_cart

  subroutine mkflux_3d_sphr(sfluxx,sfluxy,sfluxz,ng_sf,&
                            sedgex,sedgey,sedgez,ng_se, &
                            umac,vmac,wmac,ng_um, &
                            rho0_old_cart,ng_ro,rho0_new_cart,ng_rn, &
                            rhoh0_old_cart,ng_ho,rhoh0_new_cart,ng_hn, &
                            w0macx,w0macy,w0macz,ng_w0,startcomp,endcomp,lo,hi,domlo,domhi)

    use bl_constants_module
    use network, only: nspec
    use variables, only: rho_comp, rhoh_comp, trac_comp, ntrac
    use pred_parameters
    use probin_module, only: enthalpy_pred_type

    integer        , intent(in   ) :: lo(:),hi(:),domlo(:),domhi(:)
    integer        , intent(in   ) :: ng_sf,ng_se,ng_um,ng_ro,ng_rn,ng_ho,ng_hn,ng_w0
    real(kind=dp_t), intent(inout) ::        sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::        sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::        sfluxz(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::        sedgex(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::        sedgey(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::        sedgez(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::          umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::          vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::          wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) :: rho0_old_cart(lo(1)-ng_ro:,lo(2)-ng_ro:,lo(3)-ng_ro:)
    real(kind=dp_t), intent(in   ) :: rho0_new_cart(lo(1)-ng_rn:,lo(2)-ng_rn:,lo(3)-ng_rn:)
    real(kind=dp_t), intent(in   ) ::rhoh0_old_cart(lo(1)-ng_ho:,lo(2)-ng_ho:,lo(3)-ng_ho:)
    real(kind=dp_t), intent(in   ) ::rhoh0_new_cart(lo(1)-ng_hn:,lo(2)-ng_hn:,lo(3)-ng_hn:)
    real(kind=dp_t), intent(in   ) ::        w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::        w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::        w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local
    integer         :: comp
    integer         :: i,j,k
    real(kind=dp_t) :: bc_lox,bc_loy,bc_loz
    real(kind=dp_t) :: rho0_edge
    logical         :: test
    
    do comp = startcomp, endcomp

       ! test = T means the edge states are NOT in perturbational form
       test = ( (comp.eq.rhoh_comp).and. &
                     ( enthalpy_pred_type.eq.predict_h .or. &
                       enthalpy_pred_type.eq.predict_T_then_h ) ) &
         .or. ( (comp.ge.trac_comp).and.(comp.le.trac_comp+ntrac-1) )

       ! loop for x-fluxes
       if (test) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1
                   
                   if (i.eq.domlo(1)) then
                      rho0_edge = HALF * (rho0_old_cart(i,j,k)+rho0_new_cart(i,j,k))
                   else if (i.eq.domhi(1)+1) then
                      rho0_edge = HALF * (rho0_old_cart(i-1,j,k)+rho0_new_cart(i-1,j,k))
                   else
                      rho0_edge = ( rho0_old_cart(i,j,k)+rho0_old_cart(i-1,j,k) &
                                   +rho0_new_cart(i,j,k)+rho0_new_cart(i-1,j,k) ) * FOURTH
                   end if
                   
                   sfluxx(i,j,k,comp) = (umac(i,j,k) + w0macx(i,j,k)) * &
                        (rho0_edge + sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp)
                
                end do
             end do
          end do

       else

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1


                   if (i.eq.domlo(1)) then
                      bc_lox = HALF * (rhoh0_old_cart(i,j,k)+rhoh0_new_cart(i,j,k))
                   else if (i.eq.domhi(1)+1) then
                      bc_lox = HALF * (rhoh0_old_cart(i-1,j,k)+rhoh0_new_cart(i-1,j,k))
                   else
                      bc_lox = ( rhoh0_old_cart(i,j,k)+rhoh0_old_cart(i-1,j,k) &
                                +rhoh0_new_cart(i,j,k)+rhoh0_new_cart(i-1,j,k) ) * FOURTH
                   end if

                   sfluxx(i,j,k,comp) = &
                        (umac(i,j,k)+w0macx(i,j,k))*(bc_lox+sedgex(i,j,k,comp))
                
                end do
             end do
          end do

       endif


       ! loop for y-fluxes
       if (test) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)
                   
                   if (j.eq.domlo(2)) then
                      rho0_edge = HALF * (rho0_old_cart(i,j,k)+rho0_new_cart(i,j,k))
                   else if (j.eq.domhi(2)+1) then
                      rho0_edge = HALF * (rho0_old_cart(i,j-1,k)+rho0_new_cart(i,j-1,k))
                   else
                      rho0_edge = ( rho0_old_cart(i,j,k)+rho0_old_cart(i,j-1,k) &
                                   +rho0_new_cart(i,j,k)+rho0_new_cart(i,j-1,k) ) * FOURTH
                   end if
                   
                   sfluxy(i,j,k,comp) = (vmac(i,j,k) + w0macy(i,j,k)) * &
                        (rho0_edge + sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp)
                      
                end do
             end do
          end do

       else

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)
                
                   if (j.eq.domlo(2)) then
                      bc_loy = HALF * (rhoh0_old_cart(i,j,k)+rhoh0_new_cart(i,j,k))
                   else if (j.eq.domhi(2)+1) then
                      bc_loy = HALF * (rhoh0_old_cart(i,j-1,k)+rhoh0_new_cart(i,j-1,k))
                   else
                      bc_loy = ( rhoh0_old_cart(i,j,k)+rhoh0_old_cart(i,j-1,k) &
                                +rhoh0_new_cart(i,j,k)+rhoh0_new_cart(i,j-1,k) ) * FOURTH
                   end if
                   
                   sfluxy(i,j,k,comp) = &
                        (vmac(i,j,k)+w0macy(i,j,k))*(sedgey(i,j,k,comp)+bc_loy)

                end do
             end do
          end do

       endif


       ! loop for z-fluxes
       if (test) then

          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   
                   if (k.eq.domlo(3)) then
                      rho0_edge = HALF * (rho0_old_cart(i,j,k)+rho0_new_cart(i,j,k))
                   else if (k.eq.domhi(3)+1) then
                      rho0_edge = HALF * (rho0_old_cart(i,j,k-1)+rho0_new_cart(i,j,k-1))
                   else
                      rho0_edge = ( rho0_old_cart(i,j,k)+rho0_old_cart(i,j,k-1) &
                                   +rho0_new_cart(i,j,k)+rho0_new_cart(i,j,k-1) ) * FOURTH
                   end if
                   
                   sfluxz(i,j,k,comp) = (wmac(i,j,k) + w0macz(i,j,k)) * &
                        (rho0_edge + sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp)

                end do
             end do
          end do

       else

          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)


                   if (k.eq.domlo(3)) then
                      bc_loz = HALF * (rhoh0_old_cart(i,j,k)+rhoh0_new_cart(i,j,k))
                   else if (k.eq.domhi(3)+1) then
                      bc_loz = HALF * (rhoh0_old_cart(i,j,k-1)+rhoh0_new_cart(i,j,k-1))
                   else
                     bc_loz = ( rhoh0_old_cart(i,j,k)+rhoh0_old_cart(i,j,k-1) &
                               +rhoh0_new_cart(i,j,k)+rhoh0_new_cart(i,j,k-1) ) * FOURTH
                   end if
                   
                   sfluxz(i,j,k,comp) = &
                        (wmac(i,j,k)+w0macz(i,j,k))*(sedgez(i,j,k,comp)+bc_loz)

                end do
             end do
          end do

       endif

    end do ! end loop over components
     
  end subroutine mkflux_3d_sphr
   
end module mkflux_module
