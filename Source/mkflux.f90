module mkflux_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: mkflux
  
contains

  subroutine mkflux(nlevs,sflux,etarhoflux,sold,sedge,umac,w0,w0_cart_vec, &
                    rho0_old,rho0_edge_old,rho0_old_cart, &
                    rho0_new,rho0_edge_new,rho0_new_cart, &
                    rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                    rhoh0_new,rhoh0_edge_new,rhoh0_new_cart, &
                    rho0_predicted_edge,startcomp,endcomp,mla,dx)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical
    use ml_restriction_module, only: ml_edge_restriction_c
    use variables, only: nscal

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(inout) :: etarhoflux(:)
    type(multifab) , intent(in   ) :: sold(:),sedge(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:),rho0_edge_old(:,0:)
    type(multifab) , intent(in   ) :: rho0_old_cart(:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:),rho0_edge_new(:,0:)
    type(multifab) , intent(in   ) :: rho0_new_cart(:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(:,0:),rhoh0_edge_old(:,0:)
    type(multifab) , intent(in   ) :: rhoh0_old_cart(:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new(:,0:),rhoh0_edge_new(:,0:)
    type(multifab) , intent(in   ) :: rhoh0_new_cart(:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    integer        , intent(in   ) :: startcomp,endcomp
    type(ml_layout), intent(inout) :: mla
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local    
    type(box) :: domain

    integer :: domlo(sold(1)%dim),domhi(sold(1)%dim)
    integer :: i,dm,n
    integer :: lo(sold(1)%dim),hi(sold(1)%dim)

    real(kind=dp_t), pointer :: sfxp(:,:,:,:)
    real(kind=dp_t), pointer :: sfyp(:,:,:,:)
    real(kind=dp_t), pointer :: sfzp(:,:,:,:)
    real(kind=dp_t), pointer :: efp(:,:,:,:)
    real(kind=dp_t), pointer :: sexp(:,:,:,:)
    real(kind=dp_t), pointer :: seyp(:,:,:,:)
    real(kind=dp_t), pointer :: sezp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)
    real(kind=dp_t), pointer :: rho0op(:,:,:,:)
    real(kind=dp_t), pointer :: rho0np(:,:,:,:)
    real(kind=dp_t), pointer :: rhoh0op(:,:,:,:)
    real(kind=dp_t), pointer :: rhoh0np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mkflux")

    dm = sold(1)%dim
    
    do n=1,nlevs

       domain = layout_get_pd(sold(n)%la)
       domlo = lwb(domain)
       domhi = upb(domain)

       do i=1, sold(n)%nboxes
          if ( multifab_remote(sold(n),i) ) cycle
          sfxp => dataptr(sflux(n,1),i)
          sfyp => dataptr(sflux(n,2),i)
          efp  => dataptr(etarhoflux(n),i)
          sexp => dataptr(sedge(n,1),i)
          seyp => dataptr(sedge(n,2),i)
          ump  => dataptr(umac(n,1),i)
          vmp  => dataptr(umac(n,2),i)
          lo = lwb(get_box(sold(n),i))
          hi = upb(get_box(sold(n),i))
          select case (dm)
          case (2)
             call mkflux_2d(sfxp(:,:,1,:), sfyp(:,:,1,:), &
                            efp(:,:,1,1), &
                            sexp(:,:,1,:), seyp(:,:,1,:), &
                            ump(:,:,1,1), vmp(:,:,1,1), &
                            rho0_old(n,:), rho0_edge_old(n,:), &
                            rho0_new(n,:), rho0_edge_new(n,:), &
                            rhoh0_old(n,:), rhoh0_edge_old(n,:), &
                            rhoh0_new(n,:), rhoh0_edge_new(n,:), &
                            rho0_predicted_edge(n,:), &
                            w0(n,:),startcomp,endcomp,lo,hi)
          case (3)
             sfzp => dataptr(sflux(n,3),i)
             sezp => dataptr(sedge(n,3),i)
             wmp  => dataptr(umac(n,3),i)
             if(spherical .eq. 0) then
                call mkflux_3d_cart(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), &
                                    efp(:,:,:,1), &
                                    sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                    rho0_old(n,:), rho0_edge_old(n,:), &
                                    rho0_new(n,:), rho0_edge_new(n,:), &
                                    rhoh0_old(n,:), rhoh0_edge_old(n,:), &
                                    rhoh0_new(n,:), rhoh0_edge_new(n,:), &
                                    rho0_predicted_edge(n,:), &
                                    w0(n,:),startcomp,endcomp,lo,hi)

             else
                rho0op => dataptr(rho0_old_cart(n), i)
                rho0np => dataptr(rho0_new_cart(n), i)
                rhoh0op => dataptr(rhoh0_old_cart(n), i)
                rhoh0np => dataptr(rhoh0_new_cart(n), i)
                w0p => dataptr(w0_cart_vec(n),i)
                call mkflux_3d_sphr(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), &
                                    efp(:,:,:,1), &
                                    sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                    rho0op(:,:,:,1), rho0np(:,:,:,1), &
                                    rhoh0op(:,:,:,1), rhoh0np(:,:,:,1), &
                                    w0p(:,:,:,:),startcomp,endcomp,lo,hi,domlo,domhi)
             endif
          end select
       end do

    end do ! end loop over levels

    ! synchronize fluxes at coarse-fine interface
    do n = nlevs,2,-1
       do i = 1, dm
          call ml_edge_restriction_c(sflux(n-1,i),1,sflux(n,i),1,mla%mba%rr(n-1,:),i,nscal)
       enddo

       call ml_edge_restriction_c(etarhoflux(n-1),1,etarhoflux(n),1,mla%mba%rr(n-1,:),dm,1)

    enddo

    call destroy(bpt)
    
  end subroutine mkflux
  
  subroutine mkflux_2d(sfluxx,sfluxy,etarhoflux,sedgex,sedgey,umac,vmac, &
                       rho0_old,rho0_edge_old,rho0_new,rho0_edge_new, &
                       rhoh0_old,rhoh0_edge_old,rhoh0_new,rhoh0_edge_new, &
                       rho0_predicted_edge,w0,startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac
    use probin_module, only: enthalpy_pred_type, predict_rho
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::  sfluxx(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) ::  sfluxy(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) :: etarhoflux(lo(1)  :,lo(2)  :)
    real(kind=dp_t), intent(inout) ::  sedgex(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) ::  sedgey(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) ::    umac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) ::    vmac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:), rho0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:), rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(0:), rhoh0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new(0:), rhoh0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local
    integer         :: comp
    integer         :: i,j
    real(kind=dp_t) :: rho_prime, rho0_edge, rhoh0_edge
    logical :: test
    
    ! loop over components
    do comp = startcomp, endcomp

       ! test = T means the edge states are NOT in perturbational form
       test = ( (comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1) ) &
         .or. ( (comp.eq.rhoh_comp).and. &
                     ( enthalpy_pred_type.eq.predict_h .or. &
                       enthalpy_pred_type.eq.predict_T_then_h ) ) &
         .or. ( (comp.ge.trac_comp).and.(comp.le.trac_comp+ntrac-1) )

       ! create x-fluxes
       if (test) then

          do j=lo(2),hi(2)
             
             rho0_edge = HALF*(rho0_old(j)+rho0_new(j))
             
             do i=lo(1),hi(1)+1
                
                rho_prime = sedgex(i,j,rho_comp)

                if (predict_rho) then
                   sfluxx(i,j,comp) = umac(i,j)*rho_prime*sedgex(i,j,comp)
                else
                   sfluxx(i,j,comp) = umac(i,j)*(rho0_edge+rho_prime)*sedgex(i,j,comp)
                end if
                
             end do
             
          end do
                
       else
              
          do j=lo(2),hi(2)
             
             rhoh0_edge = HALF*(rhoh0_old(j)+rhoh0_new(j))
             
             do i=lo(1),hi(1)+1
                
                sfluxx(i,j,comp) = umac(i,j)*(rhoh0_edge + sedgex(i,j,comp))
                
             end do
             
          end do
          
       end if
        
       ! create y-fluxes
       if (test) then
       
          do j=lo(2),hi(2)+1
             
             rho0_edge = HALF*(rho0_edge_old(j)+rho0_edge_new(j))
             
             do i=lo(1),hi(1)
                
                rho_prime = sedgey(i,j,rho_comp)
                
                if (predict_rho) then
                   sfluxy(i,j,comp) = (vmac(i,j)+w0(j))*rho_prime*sedgey(i,j,comp)
                else
                   sfluxy(i,j,comp) = &
                        (vmac(i,j)+w0(j))*(rho0_edge+rho_prime)*sedgey(i,j,comp)
                end if
                
                if ( (comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1) ) then

                   etarhoflux(i,j) = etarhoflux(i,j) + sfluxy(i,j,comp)

                   if ( comp.eq.spec_comp+nspec-1) then
                      etarhoflux(i,j) = etarhoflux(i,j) - w0(j)*rho0_predicted_edge(j)
                   end if

                end if
                
             end do

          end do
          
       else
          
          do j=lo(2),hi(2)+1
             
             rhoh0_edge = HALF*(rhoh0_edge_old(j)+rhoh0_edge_new(j))
             
             do i=lo(1),hi(1)
                
                sfluxy(i,j,comp) = (vmac(i,j)+w0(j))*sedgey(i,j,comp) + vmac(i,j)*rhoh0_edge
                
             end do
             
          end do

       end if
             
    end do

  end subroutine mkflux_2d

  subroutine mkflux_3d_cart(sfluxx,sfluxy,sfluxz,etarhoflux,sedgex,sedgey,sedgez, &
                            umac,vmac,wmac, &
                            rho0_old,rho0_edge_old,rho0_new,rho0_edge_new, &
                            rhoh0_old,rhoh0_edge_old,rhoh0_new,rhoh0_edge_new, &
                            rho0_predicted_edge,w0,startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac
    use probin_module, only: enthalpy_pred_type, predict_rho
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::  sfluxx(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) ::  sfluxy(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) ::  sfluxz(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: etarhoflux(lo(1)  :,lo(2)  :,lo(3)  :)
    real(kind=dp_t), intent(inout) ::  sedgex(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) ::  sedgey(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) ::  sedgez(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) ::    umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::    vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::    wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:), rho0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:), rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(0:), rhoh0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new(0:), rhoh0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp

   ! local
    integer         :: comp
    integer         :: i,j,k
    real(kind=dp_t) :: rho_prime, rho0_edge, rhoh0_edge
    logical         :: test
    
    ! loop over components
    do comp = startcomp, endcomp

       ! test = T means the edge states are NOT in perturbational form
       test = ( (comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1) ) &
         .or. ( (comp.eq.rhoh_comp).and. &
                     ( enthalpy_pred_type.eq.predict_h .or. &
                       enthalpy_pred_type.eq.predict_T_then_h ) ) &
         .or. ( (comp.ge.trac_comp).and.(comp.le.trac_comp+ntrac-1) )
       
       ! create x-fluxes and y-fluxes
       if (test) then

          do k=lo(3),hi(3)
             
             rho0_edge = HALF*(rho0_old(k)+rho0_new(k))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                
                   rho_prime = sedgex(i,j,k,rho_comp)
                   
                   ! sedgex is either h or X at edges
                   if (predict_rho) then
                      sfluxx(i,j,k,comp) = umac(i,j,k)*rho_prime*sedgex(i,j,k,comp)
                   else
                      sfluxx(i,j,k,comp) = &
                           umac(i,j,k)*(rho0_edge+rho_prime)*sedgex(i,j,k,comp)
                   end if
                   
                end do
             end do
             
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                
                   rho_prime = sedgey(i,j,k,rho_comp)
                   
                   ! sedgey is either h or X at edges
                   if (predict_rho) then
                      sfluxy(i,j,k,comp) = vmac(i,j,k)*rho_prime*sedgey(i,j,k,comp)
                   else
                      sfluxy(i,j,k,comp) =&
                           vmac(i,j,k)*(rho0_edge+rho_prime)*sedgey(i,j,k,comp)
                   end if

                end do
             end do
             
          end do
                
       else
              
          do k=lo(3),hi(3)
             
             rhoh0_edge = HALF*(rhoh0_old(k)+rhoh0_new(k))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                
                   sfluxx(i,j,k,comp) = umac(i,j,k)*(rhoh0_edge + sedgex(i,j,k,comp))
                   
                end do
             end do
             
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                
                   sfluxy(i,j,k,comp) = vmac(i,j,k)*(rhoh0_edge + sedgey(i,j,k,comp))
                   
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
                
                   rho_prime = sedgez(i,j,k,rho_comp)
                
                   ! sedgez is either h or X at edges
                   if (predict_rho) then
                      sfluxz(i,j,k,comp) = &
                           (wmac(i,j,k)+w0(k))*rho_prime*sedgez(i,j,k,comp)
                   else
                      sfluxz(i,j,k,comp) = &
                           (wmac(i,j,k)+w0(k))*(rho0_edge+rho_prime)*sedgez(i,j,k,comp)
                   end if

                   if ( (comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1) ) then
                      
                      etarhoflux(i,j,k) = etarhoflux(i,j,k) + sfluxz(i,j,k,comp)
                      
                      if ( comp.eq.spec_comp+nspec-1) then
                         etarhoflux(i,j,k) = &
                              etarhoflux(i,j,k) - w0(k)*rho0_predicted_edge(k)
                      end if
                      
                   end if
                
                end do
             end do

          end do
          
       else
          
          do k=lo(3),hi(3)+1
             
             rhoh0_edge = HALF*(rhoh0_edge_old(k)+rhoh0_edge_new(k))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                
                sfluxz(i,j,k,comp) = &
                     (wmac(i,j,k)+w0(k))*sedgez(i,j,k,comp) + wmac(i,j,k)*rhoh0_edge
                                
                end do
             end do

          end do

       end if
    end do

     
  end subroutine mkflux_3d_cart

  subroutine mkflux_3d_sphr(sfluxx,sfluxy,sfluxz,etarhoflux,sedgex,sedgey,sedgez, &
                            umac,vmac,wmac, &
                            rho0_old_cart,rho0_new_cart, &
                            rhoh0_old_cart,rhoh0_new_cart, &
                            w0_cart,startcomp,endcomp,lo,hi,domlo,domhi)

    use bl_constants_module
    use network, only: nspec
    use variables, only: spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac
    use pred_parameters
    use probin_module, only: enthalpy_pred_type, predict_rho

    integer        , intent(in   ) :: lo(:),hi(:),domlo(:),domhi(:)
    real(kind=dp_t), intent(inout) :: sfluxx(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: sfluxy(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: sfluxz(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: etarhoflux(lo(1)  :,lo(2)  :,lo(3)  :)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: sedgez(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: rho0_old_cart(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: rho0_new_cart(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old_cart(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new_cart(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::     w0_cart(lo(1)-2:,lo(2)-2:,lo(3)-2:,:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local
    integer         :: comp
    integer         :: i,j,k
    real(kind=dp_t) :: bc_lox,bc_loy,bc_loz
    real(kind=dp_t) :: w0_edgex, w0_edgey, w0_edgez
    real(kind=dp_t) :: rho0_edge
    logical         :: test
    
    do comp = startcomp, endcomp

       ! test = T means the edge states are NOT in perturbational form
       test = ( (comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1) ) &
         .or. ( (comp.eq.rhoh_comp).and. &
                     ( enthalpy_pred_type.eq.predict_h .or. &
                       enthalpy_pred_type.eq.predict_T_then_h ) ) &
         .or. ( (comp.ge.trac_comp).and.(comp.le.trac_comp+ntrac-1) )

       ! loop for x-fluxes
       if (test) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1
                   
                   w0_edgex = HALF * ( w0_cart(i  ,j,k,1) +w0_cart(i-1,j,k,1) )
                   
                   if (predict_rho) then
                   
                      sfluxx(i,j,k,comp) = (umac(i,j,k) + w0_edgex) * &
                           (sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp)

                   else

                      if (i.eq.domlo(1)) then
                         rho0_edge = HALF * (rho0_old_cart(i,j,k)+rho0_new_cart(i,j,k))
                      else if (i.eq.domhi(1)+1) then
                         rho0_edge = HALF * (rho0_old_cart(i-1,j,k)+rho0_new_cart(i-1,j,k))
                      else
                         rho0_edge = ( rho0_old_cart(i,j,k)+rho0_old_cart(i-1,j,k) &
                                      +rho0_new_cart(i,j,k)+rho0_new_cart(i-1,j,k) ) * FOURTH
                      end if
                      
                      sfluxx(i,j,k,comp) = (umac(i,j,k) + w0_edgex) * &
                           (rho0_edge + sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp)
                      
                   end if
                
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
                      bc_lox = (rhoh0_old_cart(i,j,k)+rhoh0_old_cart(i-1,j,k) + &
                                rhoh0_new_cart(i,j,k)+rhoh0_new_cart(i-1,j,k) ) * FOURTH
                   end if

                   w0_edgex = HALF * ( w0_cart(i  ,j,k,1) +w0_cart(i-1,j,k,1) )

                   sfluxx(i,j,k,comp) = bc_lox*umac(i,j,k) + &
                        (umac(i,j,k) + w0_edgex)*sedgex(i,j,k,comp)
                
                end do
             end do
          end do

       endif


       ! loop for y-fluxes
       if (test) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)
                   
                   w0_edgey = HALF * ( w0_cart(i,j  ,k,2) + w0_cart(i,j-1,k,2) )
                   
                   if (predict_rho) then

                      sfluxy(i,j,k,comp) = (vmac(i,j,k) + w0_edgey) * &
                           (sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp)

                   else

                      if (j.eq.domlo(2)) then
                         rho0_edge = HALF * (rho0_old_cart(i,j,k)+rho0_new_cart(i,j,k))
                      else if (j.eq.domhi(2)+1) then
                         rho0_edge = HALF * (rho0_old_cart(i,j-1,k)+rho0_new_cart(i,j-1,k))
                      else
                         rho0_edge = (rho0_old_cart(i,j,k)+rho0_old_cart(i,j-1,k) + &
                                      rho0_new_cart(i,j,k)+rho0_new_cart(i,j-1,k) ) * FOURTH
                      end if

                      sfluxy(i,j,k,comp) = (vmac(i,j,k) + w0_edgey) * &
                           (rho0_edge + sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp)
                      
                   end if

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
                      bc_loy = (rhoh0_old_cart(i,j,k)+rhoh0_old_cart(i,j-1,k) + &
                                rhoh0_new_cart(i,j,k)+rhoh0_new_cart(i,j-1,k) ) * FOURTH
                   end if
                   
                   w0_edgey = HALF * ( w0_cart(i,j  ,k,2) + w0_cart(i,j-1,k,2) )
                   
                   sfluxy(i,j,k,comp) = bc_loy*vmac(i,j,k) + &
                        (vmac(i,j,k) + w0_edgey)*sedgey(i,j,k,comp)

                end do
             end do
          end do

       endif


       ! loop for z-fluxes
       if (test) then

          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   w0_edgez = HALF * ( w0_cart(i,j,k  ,3) + w0_cart(i,j,k-1,3) )
                   
                   if (predict_rho) then

                      sfluxz(i,j,k,comp) = (wmac(i,j,k) + w0_edgez) * &
                           (sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp)

                   else

                      
                      if (k.eq.domlo(3)) then
                         rho0_edge = HALF * (rho0_old_cart(i,j,k)+rho0_new_cart(i,j,k))
                      else if (k.eq.domhi(3)+1) then
                         rho0_edge = HALF * (rho0_old_cart(i,j,k-1)+rho0_new_cart(i,j,k-1))
                      else
                         rho0_edge = ( rho0_old_cart(i,j,k)+rho0_old_cart(i,j,k-1) &
                                      +rho0_new_cart(i,j,k)+rho0_new_cart(i,j,k-1) ) * FOURTH
                      end if
                      
                      sfluxz(i,j,k,comp) = (wmac(i,j,k) + w0_edgez) * &
                           (rho0_edge + sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp)
                      
                   end if

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
                     bc_loz = (rhoh0_old_cart(i,j,k)+rhoh0_old_cart(i,j,k-1) + &
                               rhoh0_new_cart(i,j,k)+rhoh0_new_cart(i,j,k-1) ) * FOURTH
                   end if
                   
                   w0_edgez = HALF * ( w0_cart(i,j,k  ,3) + w0_cart(i,j,k-1,3) )
                   
                   sfluxz(i,j,k,comp) = bc_loz*wmac(i,j,k) + &
                        (wmac(i,j,k) + w0_edgez)*sedgez(i,j,k,comp)

                end do
             end do
          end do

       endif

    end do ! end loop over components
     
  end subroutine mkflux_3d_sphr
   
end module mkflux_module
