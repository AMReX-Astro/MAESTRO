! The mkflux routines take the predicted edges states of the scalars
! and the MAC velocities and compute the fluxes through the
! interfaces.

! For the species fluxes, the construction of the fluxes depends on
! what form the incoming edge states take.  This depends on
! species_pred_type:
!
! predict_rhoprime_and_X: 
!    We have rho' and X, and need a edge-centered base state to
!    make the final fluxes
!
! predict_rhoprime_and_rhoX:
!    We use the (rho X) edge state directly to compute the fluxes.
!    No base state input needed.
!
! predict_rho_and_X:
!   The fluxes are computed from the product of the rho and X 
!   edge states, again, no base state input needed.
!
!
! For enthalpy, there are a wide range of quantities that we predict,
! but they fall into 2 categories.  The enthalpy edge states either
! contain predictions of h or (rho h)'.  (There is limited support for
! h' prediction, but it is not well tested).  If we have h, then we
! construct a rho depending on the species states (i.e. species_pred_type).
! If we have (rho h)', then we use the base state to make (rho h)_0 on
! edges.

module mkflux_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: mk_rhoX_flux, mk_rhoh_flux
  
contains

  !**************************************************************************
  ! density and tracer fluxes
  !**************************************************************************

  subroutine mk_rhoX_flux(mla,sflux,etarhoflux,sold,sedge,umac,w0,w0mac, &
                          rho0_old,rho0_edge_old,rho0mac_old, &
                          rho0_new,rho0_edge_new,rho0mac_new, &
                          rho0_predicted_edge,startcomp,endcomp)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical
    use ml_restriction_module, only: ml_edge_restriction_c
    use variables, only: spec_comp
    use network, only: nspec

    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(inout) :: etarhoflux(:)
    type(multifab) , intent(in   ) :: sold(:),sedge(:,:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:),rho0_edge_old(:,0:)
    type(multifab) , intent(in   ) :: rho0mac_old(:,:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:),rho0_edge_new(:,0:)
    type(multifab) , intent(in   ) :: rho0mac_new(:,:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local    
    integer :: i,n,lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: ng_sf,ng_ef,ng_se,ng_um,ng_w0,ng_ro,ng_rn

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
    real(kind=dp_t), pointer :: r0xo(:,:,:,:)
    real(kind=dp_t), pointer :: r0yo(:,:,:,:)
    real(kind=dp_t), pointer :: r0zo(:,:,:,:)
    real(kind=dp_t), pointer :: r0xn(:,:,:,:)
    real(kind=dp_t), pointer :: r0yn(:,:,:,:)
    real(kind=dp_t), pointer :: r0zn(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mk_rhoX_flux")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_sf = nghost(sflux(1,1))
    ng_ef = nghost(etarhoflux(1))
    ng_se = nghost(sedge(1,1))
    ng_um = nghost(umac(1,1))
    ng_w0 = nghost(w0mac(1,1))
    ng_ro = nghost(rho0mac_old(1,1))
    ng_rn = nghost(rho0mac_new(1,1))
    
    do n=1,nlevs
       do i=1, nboxes(sold(n))
          if ( multifab_remote(sold(n),i) ) cycle
          sfxp => dataptr(sflux(n,1),i)
          efp  => dataptr(etarhoflux(n),i)
          sexp => dataptr(sedge(n,1),i)
          ump  => dataptr(umac(n,1),i)
          lo = lwb(get_box(sold(n),i))
          hi = upb(get_box(sold(n),i))
          select case (dm)
          case (1)
             call mk_rhoX_flux_1d(sfxp(:,1,1,:), ng_sf, &
                                   efp(:,1,1,1), ng_ef, &
                                  sexp(:,1,1,:), ng_se, &
                                   ump(:,1,1,1), ng_um, &
                                  rho0_edge_old(n,:), &
                                  rho0_edge_new(n,:), &
                                  rho0_predicted_edge(n,:), &
                                  w0(n,:),startcomp,endcomp,lo,hi)
          case (2)
             sfyp => dataptr(sflux(n,2),i)
             seyp => dataptr(sedge(n,2),i)
             vmp  => dataptr(umac(n,2),i)
             call mk_rhoX_flux_2d(sfxp(:,:,1,:), sfyp(:,:,1,:), ng_sf, &
                                  efp(:,:,1,1), ng_ef, &
                                  sexp(:,:,1,:), seyp(:,:,1,:), ng_se, &
                                  ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                  rho0_old(n,:), rho0_edge_old(n,:), &
                                  rho0_new(n,:), rho0_edge_new(n,:), &
                                  rho0_predicted_edge(n,:), &
                                  w0(n,:),startcomp,endcomp,lo,hi)
          case (3)
             sfyp => dataptr(sflux(n,2),i)
             sfzp => dataptr(sflux(n,3),i)
             seyp => dataptr(sedge(n,2),i)
             sezp => dataptr(sedge(n,3),i)
             vmp  => dataptr(umac(n,2),i)
             wmp  => dataptr(umac(n,3),i)
             if(spherical .eq. 0) then
                call mk_rhoX_flux_3d_cart(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), &
                                          ng_sf, efp(:,:,:,1), ng_ef, &
                                          sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), & 
                                          ng_se, ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                          ng_um, rho0_old(n,:), rho0_edge_old(n,:), &
                                          rho0_new(n,:), rho0_edge_new(n,:), &
                                          rho0_predicted_edge(n,:), &
                                          w0(n,:),startcomp,endcomp,lo,hi)

             else
                w0xp => dataptr(w0mac(n,1),i)
                w0yp => dataptr(w0mac(n,2),i)
                w0zp => dataptr(w0mac(n,3),i)
                r0xo => dataptr(rho0mac_old(n,1),i)
                r0yo => dataptr(rho0mac_old(n,2),i)
                r0zo => dataptr(rho0mac_old(n,3),i)
                r0xn => dataptr(rho0mac_new(n,1),i)
                r0yn => dataptr(rho0mac_new(n,2),i)
                r0zn => dataptr(rho0mac_new(n,3),i)
                call mk_rhoX_flux_3d_sphr(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), &
                                          ng_sf, sexp(:,:,:,:), seyp(:,:,:,:), &
                                          sezp(:,:,:,:), ng_se, ump(:,:,:,1), vmp(:,:,:,1), &
                                          wmp(:,:,:,1), ng_um, &
                                          w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_w0, &
                                          r0xo(:,:,:,1),r0yo(:,:,:,1),r0zo(:,:,:,1),ng_ro, &
                                          r0xn(:,:,:,1),r0yn(:,:,:,1),r0zn(:,:,:,1),ng_rn, &
                                          startcomp,endcomp,lo,hi)
             endif
          end select
       end do
    end do ! end loop over levels

    ! synchronize fluxes at coarse-fine interface
    do n = nlevs,2,-1
       do i = 1, dm
          call ml_edge_restriction_c(sflux(n-1,i),spec_comp,sflux(n,i),spec_comp, &
                                     mla%mba%rr(n-1,:),i,nspec)
       enddo

       if (spherical .eq. 0) then
          call ml_edge_restriction_c(etarhoflux(n-1),1,etarhoflux(n),1,mla%mba%rr(n-1,:), &
                                     dm,1)
       end if

    enddo

    call destroy(bpt)
    
  end subroutine mk_rhoX_flux


  !----------------------------------------------------------------------------
  ! mk_rhoX_flux_1d
  !----------------------------------------------------------------------------
  subroutine mk_rhoX_flux_1d(sfluxx,ng_sf,etarhoflux,ng_ef,sedgex,ng_se, &
                             umac,ng_um,rho0_edge_old,rho0_edge_new, &
                             rho0_predicted_edge,w0,startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : spec_comp, rho_comp
    use pred_parameters
    use probin_module, only: species_pred_type

    integer        , intent(in   ) :: lo(:),hi(:),ng_sf,ng_ef,ng_se,ng_um
    real(kind=dp_t), intent(inout) ::     sfluxx(lo(1)-ng_sf:,:)
    real(kind=dp_t), intent(inout) :: etarhoflux(lo(1)-ng_ef:)
    real(kind=dp_t), intent(inout) ::     sedgex(lo(1)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(0:), rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local
    integer         :: comp
    integer         :: i
    real(kind=dp_t) :: rho0_edge
    

    do comp = startcomp, endcomp

       ! create x-fluxes
       do i = lo(1),hi(1)+1

          if (species_pred_type == predict_rhoprime_and_X) then
             ! edge states are rho' and X.  To make the (rho X) flux,
             ! we need the edge state of rho0
             rho0_edge = HALF*(rho0_edge_old(i)+rho0_edge_new(i))
             sfluxx(i,comp) = &
                  (umac(i)+w0(i))*(rho0_edge+sedgex(i,rho_comp))*sedgex(i,comp)

          else if (species_pred_type == predict_rhoprime_and_rhoX) then
             ! edge states are (rho X)
             sfluxx(i,comp) = &
                  (umac(i)+w0(i))*sedgex(i,comp)             

          else if (species_pred_type == predict_rho_and_X) then
             ! edge states are rho and X
             sfluxx(i,comp) = &
                  (umac(i)+w0(i))*sedgex(i,rho_comp)*sedgex(i,comp)
          endif

          if (comp .ge. spec_comp .and. comp .le. spec_comp+nspec-1) then
             etarhoflux(i) = etarhoflux(i) + sfluxx(i,comp)
          end if

          if ( comp.eq.spec_comp+nspec-1) then
             etarhoflux(i) = etarhoflux(i) - w0(i)*rho0_predicted_edge(i)
          end if
       end do
    end do

  end subroutine mk_rhoX_flux_1d
  
  !----------------------------------------------------------------------------
  ! mk_rhoX_flux_2d
  !----------------------------------------------------------------------------
  subroutine mk_rhoX_flux_2d(sfluxx,sfluxy,ng_sf,etarhoflux,ng_ef,sedgex,sedgey,ng_se, &
                             umac,vmac,ng_um,rho0_old,rho0_edge_old,rho0_new,rho0_edge_new, &
                             rho0_predicted_edge,w0,startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : spec_comp, rho_comp
    use pred_parameters
    use probin_module, only: species_pred_type

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
    
    do comp = startcomp, endcomp

       ! create x-fluxes
       do j=lo(2),hi(2)
          rho0_edge = HALF*(rho0_old(j)+rho0_new(j))
          do i=lo(1),hi(1)+1

             if (species_pred_type == predict_rhoprime_and_X) then
                ! edge states are rho' and X.  To make the (rho X) flux,
                ! we need the edge state of rho0
                sfluxx(i,j,comp) = umac(i,j)* &
                     (rho0_edge+sedgex(i,j,rho_comp))*sedgex(i,j,comp)

             else if (species_pred_type == predict_rhoprime_and_rhoX) then
                ! edge states are (rho X)
                sfluxx(i,j,comp) = umac(i,j)*sedgex(i,j,comp)

             else if (species_pred_type == predict_rho_and_X) then
                ! edge states are rho and X
                sfluxx(i,j,comp) = umac(i,j)* &
                     sedgex(i,j,rho_comp)*sedgex(i,j,comp)
                
             endif

          end do
       end do
             
       ! create y-fluxes
       do j = lo(2),hi(2)+1
          rho0_edge = HALF*(rho0_edge_old(j)+rho0_edge_new(j))
          do i = lo(1),hi(1)

             if (species_pred_type == predict_rhoprime_and_X) then
                ! edge states are rho' and X.  To make the (rho X) flux,
                ! we need the edge state of rho0
                sfluxy(i,j,comp) = &
                     (vmac(i,j)+w0(j))*(rho0_edge+sedgey(i,j,rho_comp))*sedgey(i,j,comp)

             else if (species_pred_type == predict_rhoprime_and_rhoX) then
                ! edge states are (rho X)
                sfluxy(i,j,comp) = &
                     (vmac(i,j)+w0(j))*sedgey(i,j,comp)

             else if (species_pred_type == predict_rho_and_X) then
                ! edge state are rho and X
                sfluxy(i,j,comp) = &
                     (vmac(i,j)+w0(j))*sedgey(i,j,rho_comp)*sedgey(i,j,comp)

             endif

             if (comp .ge. spec_comp .and. comp .le. spec_comp+nspec-1) then
                etarhoflux(i,j) = etarhoflux(i,j) + sfluxy(i,j,comp)
             end if

             if ( comp.eq.spec_comp+nspec-1) then
                etarhoflux(i,j) = etarhoflux(i,j) - w0(j)*rho0_predicted_edge(j)
             end if
          end do
       end do
    end do

  end subroutine mk_rhoX_flux_2d

  !----------------------------------------------------------------------------
  ! mk_rhoX_flux_3d_cart
  !----------------------------------------------------------------------------
  subroutine mk_rhoX_flux_3d_cart(sfluxx,sfluxy,sfluxz,ng_sf,etarhoflux,ng_ef, &
                                  sedgex,sedgey,sedgez,ng_se,umac,vmac,wmac,ng_um, &
                                  rho0_old,rho0_edge_old,rho0_new,rho0_edge_new, &
                                  rho0_predicted_edge,w0,startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : spec_comp, rho_comp
    use pred_parameters
    use probin_module, only: species_pred_type

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
    
    do comp = startcomp, endcomp

       ! create x-fluxes and y-fluxes
       do k=lo(3),hi(3)
          rho0_edge = HALF*(rho0_old(k)+rho0_new(k))

          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   sfluxx(i,j,k,comp) = &
                        umac(i,j,k)*(rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp)

                else if (species_pred_type == predict_rhoprime_and_rhoX) then
                   ! edge states are (rho X)
                   sfluxx(i,j,k,comp) = &
                        umac(i,j,k)*sedgex(i,j,k,comp)          

                else if (species_pred_type == predict_rho_and_X) then         
                   ! edge states are rho and X
                   sfluxx(i,j,k,comp) = &
                        umac(i,j,k)*sedgex(i,j,k,rho_comp)*sedgex(i,j,k,comp)

                endif

             end do
          end do
          
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   sfluxy(i,j,k,comp) = &
                        vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp)

                else if (species_pred_type == predict_rhoprime_and_rhoX) then
                   ! edge states are (rho X)
                   sfluxy(i,j,k,comp) = &
                        vmac(i,j,k)*sedgey(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxy(i,j,k,comp) = &
                        vmac(i,j,k)*sedgey(i,j,k,rho_comp)*sedgey(i,j,k,comp)

                endif

             end do
          end do
       end do
        
       ! create z-fluxes
       do k=lo(3),hi(3)+1
          rho0_edge = HALF*(rho0_edge_old(k)+rho0_edge_new(k))
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))* &
                     (rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp)

                else if (species_pred_type == predict_rhoprime_and_rhoX) then
                   ! edge states are (rho X)
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))*sedgez(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))* &
                        sedgez(i,j,k,rho_comp)*sedgez(i,j,k,comp)

                endif

                if (comp .ge. spec_comp .and. comp .le. spec_comp+nspec-1) then
                   etarhoflux(i,j,k) = etarhoflux(i,j,k) + sfluxz(i,j,k,comp)
                end if
                   
                if ( comp.eq.spec_comp+nspec-1) then
                   etarhoflux(i,j,k) = etarhoflux(i,j,k) - w0(k)*rho0_predicted_edge(k)
                end if
             end do
          end do
       end do

    end do
     
  end subroutine mk_rhoX_flux_3d_cart

  !----------------------------------------------------------------------------
  ! mk_rhoX_flux_3d_sphr
  !----------------------------------------------------------------------------
  subroutine mk_rhoX_flux_3d_sphr(sfluxx,sfluxy,sfluxz,ng_sf,&
                                  sedgex,sedgey,sedgez,ng_se, &
                                  umac,vmac,wmac,ng_um, &
                                  w0macx,w0macy,w0macz,ng_w0, &
                                  rho0macx_old,rho0macy_old,rho0macz_old,ng_ro, &
                                  rho0macx_new,rho0macy_new,rho0macz_new,ng_rn, &
                                  startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only: nspec
    use variables, only: spec_comp, rho_comp
    use pred_parameters
    use probin_module, only: species_pred_type

    integer        , intent(in   ) :: lo(:),hi(:)
    integer        , intent(in   ) :: ng_sf,ng_se,ng_um,ng_w0,ng_ro,ng_rn
    real(kind=dp_t), intent(inout) ::        sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::        sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::        sfluxz(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::        sedgex(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::        sedgey(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::        sedgez(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::          umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::          vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::          wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::        w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::        w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::        w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::  rho0macx_old(lo(1)-ng_ro:,lo(2)-ng_ro:,lo(3)-ng_ro:)
    real(kind=dp_t), intent(in   ) ::  rho0macy_old(lo(1)-ng_ro:,lo(2)-ng_ro:,lo(3)-ng_ro:)
    real(kind=dp_t), intent(in   ) ::  rho0macz_old(lo(1)-ng_ro:,lo(2)-ng_ro:,lo(3)-ng_ro:)
    real(kind=dp_t), intent(in   ) ::  rho0macx_new(lo(1)-ng_rn:,lo(2)-ng_rn:,lo(3)-ng_rn:)
    real(kind=dp_t), intent(in   ) ::  rho0macy_new(lo(1)-ng_rn:,lo(2)-ng_rn:,lo(3)-ng_rn:)
    real(kind=dp_t), intent(in   ) ::  rho0macz_new(lo(1)-ng_rn:,lo(2)-ng_rn:,lo(3)-ng_rn:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local
    integer         :: comp
    integer         :: i,j,k
    real(kind=dp_t) :: rho0_edge
    

    do comp = startcomp, endcomp

       ! loop for x-fluxes
       !$OMP PARALLEL PRIVATE(i,j,k,rho0_edge)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   rho0_edge = HALF*(rho0macx_old(i,j,k)+rho0macx_new(i,j,k))
                   sfluxx(i,j,k,comp) = (umac(i,j,k) + w0macx(i,j,k)) * &
                        (rho0_edge + sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp)

                else if (species_pred_type == predict_rhoprime_and_rhoX) then
                   ! edge states are (rho X)
                   sfluxx(i,j,k,comp) = (umac(i,j,k) + w0macx(i,j,k)) * sedgex(i,j,k,comp)                   

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxx(i,j,k,comp) = (umac(i,j,k) + w0macx(i,j,k)) * &
                        sedgex(i,j,k,rho_comp)*sedgex(i,j,k,comp)

                endif

             end do
          end do
       end do
       !$OMP END DO NOWAIT

       ! loop for y-fluxes
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   rho0_edge = HALF*(rho0macy_old(i,j,k)+rho0macy_new(i,j,k))
                   sfluxy(i,j,k,comp) = (vmac(i,j,k) + w0macy(i,j,k)) * &
                        (rho0_edge + sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp)

                else if (species_pred_type == predict_rhoprime_and_rhoX) then
                   ! edge states are (rho X)
                   sfluxy(i,j,k,comp) = (vmac(i,j,k) + w0macy(i,j,k)) * sedgey(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxy(i,j,k,comp) = (vmac(i,j,k) + w0macy(i,j,k)) * &
                        sedgey(i,j,k,rho_comp)*sedgey(i,j,k,comp)

                endif

             end do
          end do
       end do
       !$OMP END DO NOWAIT

       ! loop for z-fluxes
       !$OMP DO
       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   rho0_edge = HALF*(rho0macz_old(i,j,k)+rho0macz_new(i,j,k))
                   sfluxz(i,j,k,comp) = (wmac(i,j,k) + w0macz(i,j,k)) * &
                     (rho0_edge + sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp)

                else if (species_pred_type == predict_rhoprime_and_rhoX) then
                   ! edge states are (rho X)
                   sfluxz(i,j,k,comp) = (wmac(i,j,k) + w0macz(i,j,k)) * sedgez(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxz(i,j,k,comp) = (wmac(i,j,k) + w0macz(i,j,k)) * &
                     sedgez(i,j,k,rho_comp)*sedgez(i,j,k,comp)

                endif

             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL

    end do ! end loop over components
     
  end subroutine mk_rhoX_flux_3d_sphr



  !**************************************************************************
  ! enthalpy fluxes
  !**************************************************************************

  subroutine mk_rhoh_flux(mla,sflux,sold,sedge,umac,w0,w0mac, &
                          rho0_old,rho0_edge_old,rho0mac_old, &
                          rho0_new,rho0_edge_new,rho0mac_new, &
                          rhoh0_old,rhoh0_edge_old,rhoh0mac_old, &
                          rhoh0_new,rhoh0_edge_new,rhoh0mac_new, &
                          h0mac_old,h0mac_new)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical
    use ml_restriction_module, only: ml_edge_restriction_c
    use variables, only: rhoh_comp
    use probin_module, only: enthalpy_pred_type, species_pred_type
    use pred_parameters

    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(in   ) :: sold(:),sedge(:,:),umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:),rho0_edge_old(:,0:)
    type(multifab) , intent(in   ) :: rho0mac_old(:,:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:),rho0_edge_new(:,0:)
    type(multifab) , intent(in   ) :: rho0mac_new(:,:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(:,0:),rhoh0_edge_old(:,0:)
    type(multifab) , intent(in   ) :: rhoh0mac_old(:,:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new(:,0:),rhoh0_edge_new(:,0:)
    type(multifab) , intent(in   ) :: rhoh0mac_new(:,:)
    type(multifab) , intent(in   ) :: h0mac_old(:,:),h0mac_new(:,:)

    ! local    
    integer :: i,n,lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: ng_sf,ng_se,ng_um,ng_w0,ng_0m

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
    real(kind=dp_t), pointer :: r0mxop(:,:,:,:)
    real(kind=dp_t), pointer :: rh0mxop(:,:,:,:)
    real(kind=dp_t), pointer :: h0mxop(:,:,:,:)
    real(kind=dp_t), pointer :: r0mxnp(:,:,:,:)
    real(kind=dp_t), pointer :: rh0mxnp(:,:,:,:)
    real(kind=dp_t), pointer :: h0mxnp(:,:,:,:)
    real(kind=dp_t), pointer :: r0myop(:,:,:,:)
    real(kind=dp_t), pointer :: rh0myop(:,:,:,:)
    real(kind=dp_t), pointer :: h0myop(:,:,:,:)
    real(kind=dp_t), pointer :: r0mynp(:,:,:,:)
    real(kind=dp_t), pointer :: rh0mynp(:,:,:,:)
    real(kind=dp_t), pointer :: h0mynp(:,:,:,:)
    real(kind=dp_t), pointer :: r0mzop(:,:,:,:)
    real(kind=dp_t), pointer :: rh0mzop(:,:,:,:)
    real(kind=dp_t), pointer :: h0mzop(:,:,:,:)
    real(kind=dp_t), pointer :: r0mznp(:,:,:,:)
    real(kind=dp_t), pointer :: rh0mznp(:,:,:,:)
    real(kind=dp_t), pointer :: h0mznp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mk_rhoh_flux")

    dm = mla%dim
    nlevs = mla%nlevel

    if (enthalpy_pred_type == predict_hprime) then
       if (dm /= 3) then
          call bl_error("ERROR: predict_hprime not supported in mk_rhoh_flux for 1- or 2-d")
       else if (spherical == 0 .or. species_pred_type == predict_rho_and_X) then
          call bl_error("ERROR: predict_hprime not supported for this geometry/species_pred_type")
       endif
    endif

    ng_sf = nghost(sflux(1,1))
    ng_se = nghost(sedge(1,1))
    ng_um = nghost(umac(1,1))
    ng_w0 = nghost(w0mac(1,1))
    ng_0m = nghost(rho0mac_old(1,1))
    
    do n=1,nlevs
       do i=1, nboxes(sold(n))
          if ( multifab_remote(sold(n),i) ) cycle
          sfxp => dataptr(sflux(n,1),i)
          sexp => dataptr(sedge(n,1),i)
          ump  => dataptr(umac(n,1),i)
          lo = lwb(get_box(sold(n),i))
          hi = upb(get_box(sold(n),i))
          select case (dm)
          case (1)
             call mk_rhoh_flux_1d(sfxp(:,1,1,:), ng_sf, &
                                  sexp(:,1,1,:), ng_se, &
                                   ump(:,1,1,1), ng_um, &
                                  rho0_edge_old(n,:), &
                                  rho0_edge_new(n,:), &
                                  rhoh0_edge_old(n,:), &
                                  rhoh0_edge_new(n,:), &
                                  w0(n,:),lo,hi)
          case (2)
             sfyp => dataptr(sflux(n,2),i)
             seyp => dataptr(sedge(n,2),i)
             vmp  => dataptr(umac(n,2),i)
             call mk_rhoh_flux_2d(sfxp(:,:,1,:), sfyp(:,:,1,:), ng_sf, &
                                  sexp(:,:,1,:), seyp(:,:,1,:), ng_se, &
                                  ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                  rho0_old(n,:), rho0_edge_old(n,:), &
                                  rho0_new(n,:), rho0_edge_new(n,:), &
                                  rhoh0_old(n,:), rhoh0_edge_old(n,:), &
                                  rhoh0_new(n,:), rhoh0_edge_new(n,:), &
                                  w0(n,:),lo,hi)
          case (3)
             sfyp => dataptr(sflux(n,2),i)
             sfzp => dataptr(sflux(n,3),i)
             seyp => dataptr(sedge(n,2),i)
             sezp => dataptr(sedge(n,3),i)
             vmp  => dataptr(umac(n,2),i)
             wmp  => dataptr(umac(n,3),i)
             if(spherical .eq. 0) then
                call mk_rhoh_flux_3d_cart(sfxp(:,:,:,:),sfyp(:,:,:,:),sfzp(:,:,:,:),ng_sf, &
                                          sexp(:,:,:,:),seyp(:,:,:,:),sezp(:,:,:,:),ng_se, &
                                          ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                          rho0_old(n,:), rho0_edge_old(n,:), &
                                          rho0_new(n,:), rho0_edge_new(n,:), &
                                          rhoh0_old(n,:), rhoh0_edge_old(n,:), &
                                          rhoh0_new(n,:), rhoh0_edge_new(n,:), &
                                          w0(n,:),lo,hi)

             else
                w0xp => dataptr(w0mac(n,1),i)
                w0yp => dataptr(w0mac(n,2),i)
                w0zp => dataptr(w0mac(n,3),i)
                r0mxop => dataptr(rho0mac_old(n,1),i)
                r0myop => dataptr(rho0mac_old(n,2),i)
                r0mzop => dataptr(rho0mac_old(n,3),i)
                rh0mxop => dataptr(rhoh0mac_old(n,1),i)
                rh0myop => dataptr(rhoh0mac_old(n,2),i)
                rh0mzop => dataptr(rhoh0mac_old(n,3),i)
                h0mxop => dataptr(h0mac_old(n,1),i)
                h0myop => dataptr(h0mac_old(n,2),i)
                h0mzop => dataptr(h0mac_old(n,3),i)
                r0mxnp => dataptr(rho0mac_new(n,1),i)
                r0mynp => dataptr(rho0mac_new(n,2),i)
                r0mznp => dataptr(rho0mac_new(n,3),i)
                rh0mxnp => dataptr(rhoh0mac_new(n,1),i)
                rh0mynp => dataptr(rhoh0mac_new(n,2),i)
                rh0mznp => dataptr(rhoh0mac_new(n,3),i)
                h0mxnp => dataptr(h0mac_new(n,1),i)
                h0mynp => dataptr(h0mac_new(n,2),i)
                h0mznp => dataptr(h0mac_new(n,3),i)
                call mk_rhoh_flux_3d_sphr(sfxp(:,:,:,:),sfyp(:,:,:,:),sfzp(:,:,:,:), ng_sf, &
                                          sexp(:,:,:,:),seyp(:,:,:,:),sezp(:,:,:,:), ng_se, &
                                          ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                          w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_w0, &
                                          r0mxop(:,:,:,1),r0myop(:,:,:,1),r0mzop(:,:,:,1), &
                                          h0mxop(:,:,:,1),h0myop(:,:,:,1),h0mzop(:,:,:,1), &
                                          r0mxnp(:,:,:,1),r0mynp(:,:,:,1),r0mznp(:,:,:,1), &
                                          h0mxnp(:,:,:,1),h0mynp(:,:,:,1),h0mznp(:,:,:,1), &
                                          ng_0m,lo,hi)
             endif
          end select
       end do
    end do ! end loop over levels

    ! synchronize fluxes at coarse-fine interface
    do n = nlevs,2,-1
       do i = 1, dm
          call ml_edge_restriction_c(sflux(n-1,i),rhoh_comp,sflux(n,i),rhoh_comp, &
                                     mla%mba%rr(n-1,:),i,1)
       enddo
    enddo

    call destroy(bpt)
    
  end subroutine mk_rhoh_flux

  !----------------------------------------------------------------------------
  ! mk_rhoh_flux_1d
  !----------------------------------------------------------------------------
  subroutine mk_rhoh_flux_1d(sfluxx,ng_sf,sedgex,ng_se, &
                             umac,ng_um,rho0_edge_old,rho0_edge_new, &
                             rhoh0_edge_old,rhoh0_edge_new, &
                             w0,lo,hi)
    
    use bl_constants_module
    use network, only : nspec
    use variables, only : rho_comp, rhoh_comp
    use probin_module, only: enthalpy_pred_type, species_pred_type
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:),ng_sf,ng_se,ng_um
    real(kind=dp_t), intent(inout) ::     sfluxx(lo(1)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sedgex(lo(1)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)

    ! local
    integer         :: i
    real(kind=dp_t) :: rho0_edge,rhoh0_edge
    logical         :: have_h, have_hprime

    have_h = enthalpy_pred_type.eq.predict_h .or. &
             enthalpy_pred_type.eq.predict_T_then_h .or. &
             enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime
    
    ! create x-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho'
          do i=lo(1),hi(1)+1
             rho0_edge = HALF*(rho0_edge_old(i)+rho0_edge_new(i))
             sfluxx(i,rhoh_comp) = &
                  (umac(i)+w0(i))*(rho0_edge+sedgex(i,rho_comp))*sedgex(i,rhoh_comp)
          end do
          
       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho
          do i=lo(1),hi(1)+1
             sfluxx(i,rhoh_comp) = &
                  (umac(i)+w0(i))*sedgex(i,rho_comp)*sedgex(i,rhoh_comp)
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call bl_error("mk_rhoh_flux_1d : predict_hprime not coded yet")

    else

       ! enthalpy edge state is (rho h)'
       do i=lo(1),hi(1)+1
          rhoh0_edge = HALF*(rhoh0_edge_old(i)+rhoh0_edge_new(i))
          sfluxx(i,rhoh_comp) = (umac(i)+w0(i))*(sedgex(i,rhoh_comp)+rhoh0_edge)
       end do

    end if

  end subroutine mk_rhoh_flux_1d

  !----------------------------------------------------------------------------
  ! mk_rhoh_flux_2d
  !----------------------------------------------------------------------------
  subroutine mk_rhoh_flux_2d(sfluxx,sfluxy,ng_sf,sedgex,sedgey,ng_se, &
                             umac,vmac,ng_um,rho0_old,rho0_edge_old,rho0_new,rho0_edge_new, &
                             rhoh0_old,rhoh0_edge_old,rhoh0_new,rhoh0_edge_new, &
                             w0,lo,hi)
    
    use bl_constants_module
    use network, only : nspec
    use variables, only : rho_comp, rhoh_comp
    use probin_module, only: enthalpy_pred_type, species_pred_type
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:),ng_sf,ng_se,ng_um
    real(kind=dp_t), intent(inout) ::     sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::     sedgex(lo(1)-ng_se:,lo(2)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::     sedgey(lo(1)-ng_se:,lo(2)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:),  rho0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:),  rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(0:), rhoh0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new(0:), rhoh0_edge_new(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)

    ! local
    integer         :: i,j
    real(kind=dp_t) :: rho0_edge,rhoh0_edge
    logical         :: have_h, have_hprime

    have_h = enthalpy_pred_type.eq.predict_h .or. &
             enthalpy_pred_type.eq.predict_T_then_h .or. &
             enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime
    
    ! create x-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho'          
          do j=lo(2),hi(2)
             rho0_edge = HALF*(rho0_old(j)+rho0_new(j))
             do i=lo(1),hi(1)+1
                sfluxx(i,j,rhoh_comp) = &
                     umac(i,j)*(rho0_edge+sedgex(i,j,rho_comp))*sedgex(i,j,rhoh_comp)
             end do
          end do

       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                sfluxx(i,j,rhoh_comp) = &
                     umac(i,j)*sedgex(i,j,rho_comp)*sedgex(i,j,rhoh_comp)
             end do
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call bl_error("mk_rhoh_flux_2d : predict_hprime not coded yet")

    else

       ! enthalpy edge state is (rho h)'
       do j=lo(2),hi(2)
          rhoh0_edge = HALF*(rhoh0_old(j)+rhoh0_new(j))
          do i=lo(1),hi(1)+1
             sfluxx(i,j,rhoh_comp) = umac(i,j)*(rhoh0_edge+sedgex(i,j,rhoh_comp))
          end do
       end do

    end if

    ! create y-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho'
          do j=lo(2),hi(2)+1
             rho0_edge = HALF*(rho0_edge_old(j)+rho0_edge_new(j))
             do i=lo(1),hi(1)
                sfluxy(i,j,rhoh_comp) = &
                     (vmac(i,j)+w0(j))*(rho0_edge+sedgey(i,j,rho_comp))*sedgey(i,j,rhoh_comp)
             end do
          end do
          
       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                sfluxy(i,j,rhoh_comp) = &
                     (vmac(i,j)+w0(j))*sedgey(i,j,rho_comp)*sedgey(i,j,rhoh_comp)
             end do
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call bl_error("mk_rhoh_flux_2d : predict_hprime not coded yet")

    else

       ! enthalpy edge state is (rho h)'
       do j=lo(2),hi(2)+1
          rhoh0_edge = HALF*(rhoh0_edge_old(j)+rhoh0_edge_new(j))
          do i=lo(1),hi(1)
             sfluxy(i,j,rhoh_comp) = (vmac(i,j)+w0(j))*(sedgey(i,j,rhoh_comp)+rhoh0_edge)
          end do
       end do

    end if

  end subroutine mk_rhoh_flux_2d

  !----------------------------------------------------------------------------
  ! mk_rhoh_flux_3d_cart
  !----------------------------------------------------------------------------
  subroutine mk_rhoh_flux_3d_cart(sfluxx,sfluxy,sfluxz,ng_sf,&
                                  sedgex,sedgey,sedgez,ng_se,umac,vmac,wmac,ng_um, &
                                  rho0_old,rho0_edge_old,rho0_new,rho0_edge_new, &
                                  rhoh0_old,rhoh0_edge_old,rhoh0_new,rhoh0_edge_new, &
                                  w0,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : rho_comp, rhoh_comp
    use probin_module, only: species_pred_type, enthalpy_pred_type
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

   ! local
    integer         :: i,j,k
    real(kind=dp_t) :: rho0_edge,rhoh0_edge
    logical         :: have_h, have_hprime
    
    have_h = enthalpy_pred_type.eq.predict_h .or. &
             enthalpy_pred_type.eq.predict_T_then_h .or. &
             enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime

    ! create x-fluxes and y-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho'
          do k=lo(3),hi(3)
             rho0_edge = HALF*(rho0_old(k)+rho0_new(k))
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                   sfluxx(i,j,k,rhoh_comp) = &
                        umac(i,j,k)*(rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,rhoh_comp)
                end do
             end do

             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                   sfluxy(i,j,k,rhoh_comp) = &
                        vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,rhoh_comp)
                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho
          do k=lo(3),hi(3)

             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                   sfluxx(i,j,k,rhoh_comp) = &
                        umac(i,j,k)*sedgex(i,j,k,rho_comp)*sedgex(i,j,k,rhoh_comp)
                end do
             end do

             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                   sfluxy(i,j,k,rhoh_comp) = &
                        vmac(i,j,k)*sedgey(i,j,k,rho_comp)*sedgey(i,j,k,rhoh_comp)
                end do
             end do
          end do

       endif

    else if (have_hprime) then
       
       ! enthalpy edge state is h'
       call bl_error("mk_rhoh_flux_3d_cart : predict_hprime not coded yet")

    else

       ! enthalpy edge state is (rho h)'
       do k=lo(3),hi(3)
          rhoh0_edge = HALF*(rhoh0_old(k)+rhoh0_new(k))
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                sfluxx(i,j,k,rhoh_comp) = umac(i,j,k)*(rhoh0_edge+sedgex(i,j,k,rhoh_comp))
             end do
          end do

          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*(rhoh0_edge+sedgey(i,j,k,rhoh_comp))
             end do
          end do
       end do

    end if

    ! create z-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho'
          do k=lo(3),hi(3)+1
             rho0_edge = HALF*(rho0_edge_old(k)+rho0_edge_new(k))
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k)+w0(k))* &
                        (rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,rhoh_comp)
                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho
          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k)+w0(k))* &
                        sedgez(i,j,k,rho_comp)*sedgez(i,j,k,rhoh_comp)
                end do
             end do
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call bl_error("mk_rhoh_flux_3d_cart : predict_hprime not coded yet")

    else

       ! enthalpy edge state is (rho h)'
       do k=lo(3),hi(3)+1
          rhoh0_edge = HALF*(rhoh0_edge_old(k)+rhoh0_edge_new(k))
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sfluxz(i,j,k,rhoh_comp) = &
                     (wmac(i,j,k)+w0(k))*(sedgez(i,j,k,rhoh_comp)+rhoh0_edge)
             end do
          end do
       end do

    end if

  end subroutine mk_rhoh_flux_3d_cart

  !----------------------------------------------------------------------------
  ! mk_rhoh_flux_3d_sphr
  !----------------------------------------------------------------------------
  subroutine mk_rhoh_flux_3d_sphr(sfluxx,sfluxy,sfluxz,ng_sf,&
                                  sedgex,sedgey,sedgez,ng_se, &
                                  umac,vmac,wmac,ng_um, &
                                  w0macx,w0macy,w0macz,ng_w0, &
                                  rho0macx_old,rho0macy_old,rho0macz_old, &
                                  h0macx_old,h0macy_old,h0macz_old, &
                                  rho0macx_new,rho0macy_new,rho0macz_new, &
                                  h0macx_new,h0macy_new,h0macz_new, &
                                  ng_0m,lo,hi)

    use bl_constants_module
    use network, only: nspec
    use variables, only: rho_comp, rhoh_comp
    use pred_parameters
    use probin_module, only: species_pred_type, enthalpy_pred_type

    integer        , intent(in   ) :: lo(:),hi(:)
    integer        , intent(in   ) :: ng_sf,ng_se,ng_um,ng_w0,ng_0m
    real(kind=dp_t), intent(inout) ::        sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::        sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::        sfluxz(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real(kind=dp_t), intent(inout) ::        sedgex(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::        sedgey(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) ::        sedgez(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::          umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::          vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::          wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::        w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::        w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::        w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::  rho0macx_old(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::  rho0macy_old(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::  rho0macz_old(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::    h0macx_old(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::    h0macy_old(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::    h0macz_old(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::  rho0macx_new(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::  rho0macy_new(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::  rho0macz_new(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::    h0macx_new(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::    h0macy_new(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)
    real(kind=dp_t), intent(in   ) ::    h0macz_new(lo(1)-ng_0m:,lo(2)-ng_0m:,lo(3)-ng_0m:)

    ! local
    integer         :: i,j,k
    real(kind=dp_t) :: rho0_edge,h0_edge
!   real(kind=dp_t) :: rhoh0_edge
    logical         :: have_h, have_hprime
    
    have_h = enthalpy_pred_type.eq.predict_h .or. &
             enthalpy_pred_type.eq.predict_T_then_h .or. &
             enthalpy_pred_type.eq.predict_Tprime_then_h
    
    have_hprime = enthalpy_pred_type.eq.predict_hprime
    
    ! create x-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho' 

          !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1

                   rho0_edge = HALF*(rho0macx_old(i,j,k)+rho0macx_new(i,j,k))

                   sfluxx(i,j,k,rhoh_comp) = (umac(i,j,k) + w0macx(i,j,k)) * &
                        (rho0_edge + sedgex(i,j,k,rho_comp))*sedgex(i,j,k,rhoh_comp)

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho 

          !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1

                   sfluxx(i,j,k,rhoh_comp) = (umac(i,j,k) + w0macx(i,j,k)) * &
                        sedgex(i,j,k,rho_comp)*sedgex(i,j,k,rhoh_comp)

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       endif


    else if (have_hprime) then

       ! enthalpy edge state is h'

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho'

          ! (rho h)_edge = (h' + h_0) * (rho' + rho_0) where h0 is
          ! computed from (rho h)_0 / rho_0 
          ! sfluxx = (umac(i,j,k)+w0macx(i,j,k)) * (rho h)_edge

          !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge,h0_edge)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1
                   
                   rho0_edge = HALF*(rho0macx_old(i,j,k)+rho0macx_new(i,j,k))
                   
                   h0_edge = HALF*(h0macx_old(i,j,k)+h0macx_new(i,j,k))
                   
                   sfluxx(i,j,k,rhoh_comp) = (umac(i,j,k)+w0macx(i,j,k)) * &
                        (sedgex(i,j,k,rho_comp)+rho0_edge) * (sedgex(i,j,k,rhoh_comp)+h0_edge)
                   
                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho
          call bl_error("ERROR: predict_rho_and_X and predict_hprime not supported together")

       endif

    else

       ! enthalpy edge state is (rho h)'

       !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge,h0_edge)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.  
                !  (rho h)_edge = (rho h)' + (rho_0 * h_0) 
                ! where h_0 is computed from (rho h)_0 / rho_0
                rho0_edge = HALF*(rho0macx_old(i,j,k)+rho0macx_new(i,j,k))

                h0_edge = HALF*(h0macx_old(i,j,k)+h0macx_new(i,j,k))

                sfluxx(i,j,k,rhoh_comp) = &
                     (umac(i,j,k)+w0macx(i,j,k))*(rho0_edge*h0_edge+sedgex(i,j,k,rhoh_comp))

                ! alternate options that needs further testing
                ! rhoh0_edge = HALF*(rhoh0macx_old(i,j,k)+rhoh0macx_new(i,j,k))
                ! sfluxx(i,j,k,rhoh_comp) = &
                !    (umac(i,j,k)+w0macx(i,j,k))*(rhoh0_edge+sedgex(i,j,k,rhoh_comp))

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    endif

    ! create y-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho'

          !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)

                   rho0_edge = HALF*(rho0macy_old(i,j,k)+rho0macy_new(i,j,k))

                   sfluxy(i,j,k,rhoh_comp) = (vmac(i,j,k) + w0macy(i,j,k)) * &
                        (rho0_edge + sedgey(i,j,k,rho_comp))*sedgey(i,j,k,rhoh_comp)

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho

          !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)

                   sfluxy(i,j,k,rhoh_comp) = (vmac(i,j,k) + w0macy(i,j,k)) * &
                        sedgey(i,j,k,rho_comp)*sedgey(i,j,k,rhoh_comp)

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho'

          ! (rho h)_edge = (h' + h_0) * (rho' + rho_0) where h0 is
          ! computed from (rho h)_0 / rho_0 
          ! sfluxy = (vmac(i,j,k)+w0macy(i,j,k)) * (rho h)_edge

          !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge,h0_edge)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)
                   
                   rho0_edge = HALF*(rho0macy_old(i,j,k)+rho0macy_new(i,j,k))
                   
                   h0_edge = HALF*(h0macy_old(i,j,k)+h0macy_new(i,j,k))
                   
                   sfluxy(i,j,k,rhoh_comp) = (vmac(i,j,k)+w0macy(i,j,k)) * &
                        (sedgey(i,j,k,rho_comp)+rho0_edge) * (sedgey(i,j,k,rhoh_comp)+h0_edge)
                   
                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho
          call bl_error("ERROR: predict_rho_and_X and predict_hprime not supported together")
          
       endif

    else

       ! enthalpy edge state is (rho h)'

       !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge,h0_edge)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.  
                !  (rho h)_edge = (rho h)' + (rho_0 * h_0) 
                ! where h_0 is computed from (rho h)_0 / rho_0
                rho0_edge = HALF*(rho0macy_old(i,j,k)+rho0macy_new(i,j,k))

                h0_edge = HALF*(h0macy_old(i,j,k)+h0macy_new(i,j,k))

                sfluxy(i,j,k,rhoh_comp) = &
                     (vmac(i,j,k)+w0macy(i,j,k))*(rho0_edge*h0_edge+sedgey(i,j,k,rhoh_comp))

                ! alternate options that needs further testing
                ! rhoh0_edge = HALF*(rhoh0macy_old(i,j,k)+rhoh0macy_new(i,j,k))
                ! sfluxy(i,j,k,rhoh_comp) = &
                !   (vmac(i,j,k)+w0macy(i,j,k))*(rhoh0_edge+sedgey(i,j,k,rhoh_comp))

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    endif


    ! create z-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho'

          !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge)
          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   rho0_edge = HALF*(rho0macz_old(i,j,k)+rho0macz_new(i,j,k))

                   sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k) + w0macz(i,j,k)) * &
                        (rho0_edge + sedgez(i,j,k,rho_comp))*sedgez(i,j,k,rhoh_comp)

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho

          !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge)
          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k) + w0macz(i,j,k)) * &
                        sedgez(i,j,k,rho_comp)*sedgez(i,j,k,rhoh_comp)

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'

       if (species_pred_type == predict_rhoprime_and_X .or. &
           species_pred_type == predict_rhoprime_and_rhoX) then

          ! density edge state is rho'

          ! (rho h)_edge = (h' + h_0) * (rho' + rho_0)
          ! where h0 is computed from (rho h)_0 / rho_0
          ! sfluxz = (wmac(i,j,k)+w0macz(i,j,k)) * (rho h)_edge
          
          !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge,h0_edge)
          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   
                   rho0_edge = HALF*(rho0macz_old(i,j,k)+rho0macz_new(i,j,k))
                   
                   h0_edge = HALF*(h0macz_old(i,j,k)+h0macz_new(i,j,k))
                   
                   sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k)+w0macz(i,j,k)) * &
                        (sedgez(i,j,k,rho_comp)+rho0_edge) * (sedgez(i,j,k,rhoh_comp)+h0_edge)
                   
                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else if (species_pred_type == predict_rho_and_X) then

          ! density edge state is rho
          call bl_error("ERROR: predict_rho_and_X and predict_hprime not supported together")

       endif

    else

       ! enthalpy edge state is (rho h)'

       !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge,h0_edge)
       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.  
                !  (rho h)_edge = (rho h)' + (rho_0 * h_0) 
                ! where h_0 is computed from (rho h)_0 / rho_0
                rho0_edge = HALF*(rho0macz_old(i,j,k)+rho0macz_new(i,j,k))

                h0_edge = HALF*(h0macz_old(i,j,k)+h0macz_new(i,j,k))

                sfluxz(i,j,k,rhoh_comp) = &
                     (wmac(i,j,k)+w0macz(i,j,k))*(rho0_edge*h0_edge+sedgez(i,j,k,rhoh_comp))

                ! alternate options that needs further testing
                ! rhoh0_edge = HALF*(rhoh0macz_old(i,j,k)+rhoh0macz_new(i,j,k))
                ! sfluxz(i,j,k,rhoh_comp) = &
                !   (wmac(i,j,k)+w0macz(i,j,k))*(rhoh0_edge+sedgez(i,j,k,rhoh_comp))

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    endif

  end subroutine mk_rhoh_flux_3d_sphr

end module mkflux_module
