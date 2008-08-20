module mk_rhoX_flux_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: mk_rhoX_flux
  
contains

  subroutine mk_rhoX_flux(mla,sflux,etarhoflux,sold,sedge,umac,w0,w0mac, &
                          rho0_old,rho0_edge_old,rho0_old_cart, &
                          rho0_new,rho0_edge_new,rho0_new_cart, &
                          rho0_predicted_edge,startcomp,endcomp)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical
    use ml_restriction_module, only: ml_edge_restriction_c
    use variables, only: nscal

    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(inout) :: etarhoflux(:)
    type(multifab) , intent(in   ) :: sold(:),sedge(:,:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:),rho0_edge_old(:,0:)
    type(multifab) , intent(in   ) :: rho0_old_cart(:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:),rho0_edge_new(:,0:)
    type(multifab) , intent(in   ) :: rho0_new_cart(:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local    
    type(box) :: domain

    integer :: domlo(sold(1)%dim),domhi(sold(1)%dim)
    integer :: i,dm,n,nlevs
    integer :: lo(sold(1)%dim),hi(sold(1)%dim)
    integer :: ng_sf,ng_ef,ng_se,ng_um,ng_ro,ng_rn,ng_w0


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
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    real(kind=dp_t), pointer :: rho0op(:,:,:,:)
    real(kind=dp_t), pointer :: rho0np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mk_rhoX_flux")

    dm    = mla%dim
    nlevs = mla%nlevel
    
    ng_sf = sflux(1,1)%ng
    ng_ef = etarhoflux(1)%ng
    ng_se = sedge(1,1)%ng
    ng_um = umac(1,1)%ng
    ng_ro = rho0_old_cart(1)%ng
    ng_rn = rho0_new_cart(1)%ng
    ng_w0 = w0mac(1,1)%ng
    
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
             call mk_rhoX_flux_2d(sfxp(:,:,1,:), sfyp(:,:,1,:), ng_sf, &
                            efp(:,:,1,1), ng_ef, &
                            sexp(:,:,1,:), seyp(:,:,1,:), ng_se, &
                            ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                            rho0_old(n,:), rho0_edge_old(n,:), &
                            rho0_new(n,:), rho0_edge_new(n,:), &
                            rho0_predicted_edge(n,:), &
                            w0(n,:),startcomp,endcomp,lo,hi)
          case (3)
             sfzp => dataptr(sflux(n,3),i)
             sezp => dataptr(sedge(n,3),i)
             wmp  => dataptr(umac(n,3),i)
             if(spherical .eq. 0) then
                call mk_rhoX_flux_3d_cart(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), ng_sf, &
                                    efp(:,:,:,1), ng_ef, &
                                    sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), ng_se, &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                    rho0_old(n,:), rho0_edge_old(n,:), &
                                    rho0_new(n,:), rho0_edge_new(n,:), &
                                    rho0_predicted_edge(n,:), &
                                    w0(n,:),startcomp,endcomp,lo,hi)

             else
                rho0op => dataptr(rho0_old_cart(n), i)
                rho0np => dataptr(rho0_new_cart(n), i)
                w0xp => dataptr(w0mac(n,1),i)
                w0yp => dataptr(w0mac(n,2),i)
                w0zp => dataptr(w0mac(n,3),i)
                call mk_rhoX_flux_3d_sphr(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), ng_sf, &
                                    sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), ng_se, &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                    rho0op(:,:,:,1), ng_ro, rho0np(:,:,:,1), ng_rn, &
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

       if (spherical .eq. 0) &
         call ml_edge_restriction_c(etarhoflux(n-1),1,etarhoflux(n),1,mla%mba%rr(n-1,:),dm,1)

    enddo

    call destroy(bpt)
    
  end subroutine mk_rhoX_flux
  
  subroutine mk_rhoX_flux_2d(sfluxx,sfluxy,ng_sf,etarhoflux,ng_ef,sedgex,sedgey,ng_se, &
                       umac,vmac,ng_um,rho0_old,rho0_edge_old,rho0_new,rho0_edge_new, &
                       rho0_predicted_edge,w0,startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : spec_comp, rho_comp
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:),ng_sf,ng_ef,ng_se,ng_um
    real(kind=dp_t), intent(inout) ::     sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real(kind=dp_t), intent(inout) :: etarhoflux(lo(1)-ng_ef:,lo(2)-ng_ef:)
    real(kind=dp_t), intent(inout) ::     sedgex(lo(1)-ng_se:,lo(2)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::     sedgey(lo(1)-ng_se:,lo(2)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:), rho0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:), rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local
    integer         :: comp
    integer         :: i,j
    real(kind=dp_t) :: rho0_edge
    
    ! loop over components -- note the edges states are X
    do comp = startcomp, endcomp

       ! create x-fluxes
       do j=lo(2),hi(2)
          rho0_edge = HALF*(rho0_old(j)+rho0_new(j))
          do i=lo(1),hi(1)+1
             sfluxx(i,j,comp) = umac(i,j)* &
                  (rho0_edge+sedgex(i,j,rho_comp))*sedgex(i,j,comp)
          end do
       end do
             
       ! create y-fluxes
       do j = lo(2),hi(2)+1
          rho0_edge = HALF*(rho0_edge_old(j)+rho0_edge_new(j))
          do i = lo(1),hi(1)

             sfluxy(i,j,comp) = &
                  (vmac(i,j)+w0(j))*(rho0_edge+sedgey(i,j,rho_comp))*sedgey(i,j,comp)

             if (comp .ge. spec_comp .and. comp .le. spec_comp+nspec-1) &
                etarhoflux(i,j) = etarhoflux(i,j) + sfluxy(i,j,comp)

             if ( comp.eq.spec_comp+nspec-1) &
                etarhoflux(i,j) = etarhoflux(i,j) - w0(j)*rho0_predicted_edge(j)

          end do
       end do
    end do

  end subroutine mk_rhoX_flux_2d

  subroutine mk_rhoX_flux_3d_cart(sfluxx,sfluxy,sfluxz,ng_sf,etarhoflux,ng_ef, &
                            sedgex,sedgey,sedgez,ng_se,umac,vmac,wmac,ng_um, &
                            rho0_old,rho0_edge_old,rho0_new,rho0_edge_new, &
                            rho0_predicted_edge,w0,startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : spec_comp, rho_comp
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:),ng_sf,ng_ef,ng_se,ng_um
    real(kind=dp_t), intent(inout) ::     sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sfluxz(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) :: etarhoflux(lo(1)-ng_ef:,lo(2)-ng_ef:,lo(3)-ng_ef:)
    real(kind=dp_t), intent(inout) ::     sedgex(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::     sedgey(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::     sedgez(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:), rho0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:), rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp

   ! local
    integer         :: comp
    integer         :: i,j,k
    real(kind=dp_t) :: rho0_edge
    
    ! loop over components -- note the edges states are X
    do comp = startcomp, endcomp

       ! create x-fluxes and y-fluxes
       do k=lo(3),hi(3)
             
          rho0_edge = HALF*(rho0_old(k)+rho0_new(k))
             
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                
                sfluxx(i,j,k,comp) = &
                     umac(i,j,k)*(rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp)
                
             end do
          end do
          
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                
                sfluxy(i,j,k,comp) = &
                     vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp)

             end do
          end do
          
       end do
        
       ! create z-fluxes
       do k=lo(3),hi(3)+1
             
          rho0_edge = HALF*(rho0_edge_old(k)+rho0_edge_new(k))
          
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
             
                sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))* &
                     (rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp)

                if (comp .ge. spec_comp .and. comp .le. spec_comp+nspec-1) &
                   etarhoflux(i,j,k) = etarhoflux(i,j,k) + sfluxz(i,j,k,comp)
                   
                if ( comp.eq.spec_comp+nspec-1) &
                   etarhoflux(i,j,k) = &
                        etarhoflux(i,j,k) - w0(k)*rho0_predicted_edge(k)
             end do
          end do

       end do
    end do
     
  end subroutine mk_rhoX_flux_3d_cart

  subroutine mk_rhoX_flux_3d_sphr(sfluxx,sfluxy,sfluxz,ng_sf,&
                            sedgex,sedgey,sedgez,ng_se, &
                            umac,vmac,wmac,ng_um, &
                            rho0_old_cart,ng_ro,rho0_new_cart,ng_rn, &
                            w0macx,w0macy,w0macz,ng_w0,startcomp,endcomp,lo,hi,domlo,domhi)

    use bl_constants_module
    use network, only: nspec
    use variables, only: spec_comp, rho_comp
    use pred_parameters
    use probin_module, only: enthalpy_pred_type

    integer        , intent(in   ) :: lo(:),hi(:),domlo(:),domhi(:)
    integer        , intent(in   ) :: ng_sf,ng_se,ng_um,ng_ro,ng_rn,ng_w0
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
    real(kind=dp_t), intent(in   ) ::        w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::        w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::        w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local
    integer         :: comp
    integer         :: i,j,k
    real(kind=dp_t) :: rho0_edge
    
    ! loop over components -- note the edges states are X
    do comp = startcomp, endcomp

       ! loop for x-fluxes
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

       ! loop for y-fluxes
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

       ! loop for z-fluxes
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

    end do ! end loop over components
     
  end subroutine mk_rhoX_flux_3d_sphr
   
end module mk_rhoX_flux_module
