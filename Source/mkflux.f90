module mkflux_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: mkflux
  
contains

  subroutine mkflux(nlevs,sflux,etaflux,sold,sedge,umac,w0,w0_cart_vec,s0_old,s0_edge_old, &
                    s0_old_cart,s0_new,s0_edge_new,s0_new_cart, &
                    s0_predicted_edge,startcomp,endcomp,mla,dx,dt)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical
    use ml_restriction_module, only: ml_edge_restriction_c
    use variables, only: nscal

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(inout) :: etaflux(:)
    type(multifab) , intent(in   ) :: sold(:),sedge(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,0:,:),s0_edge_old(:,0:,:)
    type(multifab) , intent(in   ) :: s0_old_cart(:)
    real(kind=dp_t), intent(in   ) :: s0_new(:,0:,:),s0_edge_new(:,0:,:)
    type(multifab) , intent(in   ) :: s0_new_cart(:)
    real(kind=dp_t), intent(in   ) :: s0_predicted_edge(:,0:,:)
    integer        , intent(in   ) :: startcomp,endcomp
    type(ml_layout), intent(inout) :: mla
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt

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
    real(kind=dp_t), pointer :: s0op(:,:,:,:)
    real(kind=dp_t), pointer :: s0np(:,:,:,:)

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
          efp  => dataptr(etaflux(n),i)
          sexp => dataptr(sedge(n,1),i)
          seyp => dataptr(sedge(n,2),i)
          ump  => dataptr(umac(n,1),i)
          vmp  => dataptr(umac(n,2),i)
          lo = lwb(get_box(sold(n),i))
          hi = upb(get_box(sold(n),i))
          select case (dm)
          case (2)
             call mkflux_2d(sfxp(:,:,1,:), sfyp(:,:,1,:), &
                            efp(:,:,1,:), &
                            sexp(:,:,1,:), seyp(:,:,1,:), &
                            ump(:,:,1,1), vmp(:,:,1,1), &
                            s0_old(n,:,:), s0_edge_old(n,:,:), &
                            s0_new(n,:,:), s0_edge_new(n,:,:), &
                            s0_predicted_edge(n,:,:), &
                            w0(n,:),startcomp,endcomp,lo,hi,dx(n,:),dt)
          case (3)
             sfzp => dataptr(sflux(n,3),i)
             sezp => dataptr(sedge(n,3),i)
             wmp  => dataptr(umac(n,3),i)
             if(spherical .eq. 0) then
                call mkflux_3d_cart(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), &
                                    efp(:,:,:,:), &
                                    sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                    s0_old(n,:,:), s0_edge_old(n,:,:), &
                                    s0_new(n,:,:), s0_edge_new(n,:,:), &
                                    s0_predicted_edge(n,:,:), &
                                    w0(n,:),startcomp,endcomp,lo,hi)

             else
                s0op => dataptr(s0_old_cart(n), i)
                s0np => dataptr(s0_new_cart(n), i)
                w0p => dataptr(w0_cart_vec(n),i)
                call mkflux_3d_sphr(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), &
                                    efp(:,:,:,:), &
                                    sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                    s0_old(n,:,:), s0_edge_old(n,:,:), s0op(:,:,:,:), &
                                    s0_new(n,:,:), s0_edge_new(n,:,:), s0np(:,:,:,:), &
                                    w0(n,:),w0p(:,:,:,:),startcomp,endcomp,lo,hi,domlo,domhi)
             endif
          end select
       end do

    end do ! end loop over levels

    ! synchronize fluxes at coarse-fine interface
    do n = nlevs,2,-1
       do i = 1, dm
          call ml_edge_restriction_c(sflux(n-1,i),1,sflux(n,i),1,mla%mba%rr(n-1,:),i,nscal)
       enddo

       call ml_edge_restriction_c(etaflux(n-1),1,etaflux(n),1,mla%mba%rr(n-1,:),dm,nscal)

    enddo

    call destroy(bpt)
    
  end subroutine mkflux
  
  subroutine mkflux_2d(sfluxx,sfluxy,etaflux,sedgex,sedgey,umac,vmac,s0_old,s0_edge_old, &
                       s0_new,s0_edge_new,s0_pred_edge,w0,startcomp,endcomp, &
                       lo,hi,dx,dt)

    use bl_constants_module
    use network, only : nspec
    use variables, only : spec_comp, rho_comp, rhoh_comp
    use probin_module, only: predict_X_at_edges, enthalpy_pred_type

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::  sfluxx(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) ::  sfluxy(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) :: etaflux(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) ::  sedgex(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) ::  sedgey(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) ::    umac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) ::    vmac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_pred_edge(0:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: comp,spec
    integer :: i,j
    real(kind=dp_t) :: s0_edge
    real(kind=dp_t) :: rho_prime, rho0_edge
    logical :: test,needrhoprime
    
    ! loop over components
    do comp = startcomp, endcomp

       test = ((comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1).and.predict_X_at_edges) &
         .or. ((comp.eq.rhoh_comp).and. &
                     (enthalpy_pred_type.eq.2 .or. &
                     (enthalpy_pred_type.eq.3.and.predict_X_at_edges)))
       
       needrhoprime = ((comp.eq.rhoh_comp).and. &
            enthalpy_pred_type.eq.2.and.(.not.predict_X_at_edges))

       if (needrhoprime) then

          ! compute rho' on x-faces
          sedgex(lo(1):hi(1)+1,lo(2):hi(2),rho_comp) = ZERO
          do spec = 1,nspec    
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1
                   sedgex(i,j,rho_comp) = sedgex(i,j,rho_comp)+sedgex(i,j,spec_comp+spec-1)
                end do
             end do
          end do
          
          ! compute rho' on y-faces
          sedgey(lo(1):hi(1),lo(2):hi(2)+1,rho_comp) = ZERO
          do spec = 1,nspec    
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)
                   sedgey(i,j,rho_comp) = sedgey(i,j,rho_comp) + sedgey(i,j,spec_comp+spec-1)
                end do
             end do
          end do
          
       end if

       ! create x-fluxes
       if (test) then

          do j=lo(2),hi(2)
             
             rho0_edge = HALF*(s0_old(j,rho_comp)+s0_new(j,rho_comp))
             
             do i=lo(1),hi(1)+1
                
                rho_prime = sedgex(i,j,rho_comp)
                
                ! sedgex is either h or X at edges
                sfluxx(i,j,comp) = umac(i,j)*(rho0_edge+rho_prime)*sedgex(i,j,comp)
                
             end do
             
          end do
                
       else
              
          do j=lo(2),hi(2)
             
             s0_edge = HALF*(s0_old(j,comp)+s0_new(j,comp))
             
             do i=lo(1),hi(1)+1
                
                ! s0_edge is either (rho h)_0, (rho X)_0, or (rho trac)_0 at edges
                ! sedgex is either (rho h)', (rho X)', or (rho trac)' at edges
                sfluxx(i,j,comp) = umac(i,j)*(s0_edge + sedgex(i,j,comp))
                
             end do
             
          end do
          
       end if
        
       ! create y-fluxes
       if (test) then
       
          do j=lo(2),hi(2)+1
             
             rho0_edge = HALF*(s0_edge_old(j,rho_comp)+s0_edge_new(j,rho_comp))
             
             do i=lo(1),hi(1)
                
                rho_prime = sedgey(i,j,rho_comp)
                
                ! sedgey is either h or X at edges
                sfluxy(i,j,comp) = (vmac(i,j)+w0(j))*(rho0_edge+rho_prime)*sedgey(i,j,comp)
                
                etaflux(i,j,comp) = &
                     sfluxy(i,j,comp) - w0(j)*s0_pred_edge(j,comp)*s0_pred_edge(j,rho_comp)
                
             end do

          end do
          
       else
          
          do j=lo(2),hi(2)+1
             
             s0_edge = HALF*(s0_edge_old(j,comp)+s0_edge_new(j,comp))
             
             do i=lo(1),hi(1)
                
                ! s0_edge is either (rho h)_0, (rho X)_0, or (rho trac)_0 at edges
                ! sedgey is either (rho h)', (rho X)', or (rho trac)' at edges
                sfluxy(i,j,comp) = (vmac(i,j)+w0(j))*sedgey(i,j,comp) + vmac(i,j)*s0_edge
                
                etaflux(i,j,comp) = sfluxy(i,j,comp)
                
             end do
             
          end do

       end if
             
    end do

  end subroutine mkflux_2d

  subroutine mkflux_3d_cart(sfluxx,sfluxy,sfluxz,etaflux,sedgex,sedgey,sedgez, &
                            umac,vmac,wmac,s0_old,s0_edge_old,s0_new,s0_edge_new, &
                            s0_pred_edge,w0,startcomp,endcomp,lo,hi)

    use bl_constants_module
    use network, only : nspec
    use variables, only : spec_comp, rho_comp, rhoh_comp
    use probin_module, only: predict_X_at_edges, enthalpy_pred_type

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::  sfluxx(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) ::  sfluxy(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) ::  sfluxz(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: etaflux(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) ::  sedgex(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) ::  sedgey(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) ::  sedgez(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) ::    umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::    vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::    wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_pred_edge(0:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp

   ! local
    real(kind=dp_t) :: s0_edge
    real(kind=dp_t) :: rho_prime, rho0_edge
    integer         :: comp,i,j,k
    logical         :: test
    
    ! loop over components
    do comp = startcomp, endcomp

       test = ((comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1).and.predict_X_at_edges) &
         .or. ((comp.eq.rhoh_comp).and.enthalpy_pred_type.eq.2)
       
       ! create x-fluxes and y-fluxes
       if (test) then

          do k=lo(3),hi(3)
             
             rho0_edge = HALF*(s0_old(k,rho_comp)+s0_new(k,rho_comp))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                
                   rho_prime = sedgex(i,j,k,rho_comp)
                   
                   ! sedgex is either h or X at edges
                   sfluxx(i,j,k,comp) = umac(i,j,k)*(rho0_edge+rho_prime)*sedgex(i,j,k,comp)
                   
                end do
             end do
             
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                
                   rho_prime = sedgey(i,j,k,rho_comp)
                   
                   ! sedgey is either h or X at edges
                   sfluxy(i,j,k,comp) = vmac(i,j,k)*(rho0_edge+rho_prime)*sedgey(i,j,k,comp)
                   
                end do
             end do
             
          end do
                
       else
              
          do k=lo(3),hi(3)
             
             s0_edge = HALF*(s0_old(k,comp)+s0_new(k,comp))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                
                   ! s0_edge is either (rho h)_0, (rho X)_0, or (rho trac)_0 at edges
                   ! sedgex is either (rho h)', (rho X)', or (rho trac)' at edges
                   sfluxx(i,j,k,comp) = umac(i,j,k)*(s0_edge + sedgex(i,j,k,comp))
                   
                end do
             end do
             
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                
                   ! s0_edge is either (rho h)_0, (rho X)_0, or (rho trac)_0 at edges
                   ! sedgex is either (rho h)', (rho X)', or (rho trac)' at edges
                   sfluxy(i,j,k,comp) = vmac(i,j,k)*(s0_edge + sedgey(i,j,k,comp))
                   
                end do
             end do
             
          end do
          
       end if
        
       ! create z-fluxes
       if (test) then
       
          do k=lo(3),hi(3)+1
             
             rho0_edge = HALF*(s0_edge_old(k,rho_comp)+s0_edge_new(k,rho_comp))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                
                   rho_prime = sedgez(i,j,k,rho_comp)
                
                   ! sedgez is either h or X at edges
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))*(rho0_edge+rho_prime)*sedgez(i,j,k,comp)
                   
                   etaflux(i,j,k,comp) = &
                        sfluxz(i,j,k,comp) - w0(k)*s0_pred_edge(k,comp)*s0_pred_edge(k,rho_comp)
                
                end do
             end do

          end do
          
       else
          
          do k=lo(3),hi(3)+1
             
             s0_edge = HALF*(s0_edge_old(k,comp)+s0_edge_new(k,comp))
             
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                
                ! s0_edge is either (rho h)_0, (rho X)_0, or (rho trac)_0 at edges
                ! sedgey is either (rho h)', (rho X)', or (rho trac)' at edges
                sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))*sedgez(i,j,k,comp) + wmac(i,j,k)*s0_edge
                
                etaflux(i,j,k,comp) = sfluxz(i,j,k,comp)
                
                end do
             end do

          end do

       end if
    end do

     
  end subroutine mkflux_3d_cart

  subroutine mkflux_3d_sphr(sfluxx,sfluxy,sfluxz,etaflux,sedgex,sedgey,sedgez, &
                            umac,vmac,wmac, &
                            s0_old,s0_edge_old,s0_old_cart,s0_new,s0_edge_new,s0_new_cart, &
                            w0,w0_cart,startcomp,endcomp,lo,hi,domlo,domhi)

    use bl_constants_module
    use addw0_module

    integer        , intent(in   ) :: lo(:),hi(:),domlo(:),domhi(:)
    real(kind=dp_t), intent(inout) :: sfluxx(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: sfluxy(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: sfluxz(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: etaflux(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) :: sedgex(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) :: sedgey(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) :: sedgez(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(inout) ::   wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::      s0_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_old_cart(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
    real(kind=dp_t), intent(in   ) ::      s0_new(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_edge_new(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new_cart(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
    real(kind=dp_t), intent(in   ) ::          w0(0:)
    real(kind=dp_t), intent(in   ) ::     w0_cart(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
    integer        , intent(in   ) :: startcomp,endcomp

    ! local
    integer         :: i,j,k,comp
    real(kind=dp_t) :: mult
    real(kind=dp_t) :: bc_lox,bc_loy,bc_loz
    real(kind=dp_t) :: w0_edgex, w0_edgey, w0_edgez

    do comp = startcomp, endcomp

       ! loop for x-fluxes
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                bc_lox = (s0_old_cart(i,j,k,comp)+s0_old_cart(i-1,j,k,comp) &
                     +s0_new_cart(i,j,k,comp)+s0_new_cart(i-1,j,k,comp) ) * FOURTH

                if (i.eq.domlo(1)) then
                   bc_lox = HALF * (s0_old_cart(i,j,k,comp)+s0_new_cart(i,j,k,comp))
                end if
                if (i.eq.domhi(1)+1) then
                   bc_lox = HALF * (s0_old_cart(i-1,j,k,comp)+s0_new_cart(i-1,j,k,comp))
                end if

                w0_edgex = HALF * ( w0_cart(i  ,j,k,1) +w0_cart(i-1,j,k,1) )

                sfluxx(i,j,k,comp) = bc_lox*umac(i,j,k) + &
                     (umac(i,j,k) + w0_edgex)*sedgex(i,j,k,comp)
                
             end do
          end do
       end do


       ! loop for y-fluxes
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                bc_loy = (s0_old_cart(i,j,k,comp)+s0_old_cart(i,j-1,k,comp) &
                     +s0_new_cart(i,j,k,comp)+s0_new_cart(i,j-1,k,comp) ) * FOURTH
                
                if (j.eq.domlo(2)) then
                   bc_loy = HALF * (s0_old_cart(i,j,k,comp)+s0_new_cart(i,j,k,comp))
                end if
                if (j.eq.domhi(2)+1) then
                   bc_loy = HALF * (s0_old_cart(i,j-1,k,comp)+s0_new_cart(i,j-1,k,comp))
                end if

                w0_edgey = HALF * ( w0_cart(i,j  ,k,2) + w0_cart(i,j-1,k,2) )

                sfluxy(i,j,k,comp) = bc_loy*vmac(i,j,k) + &
                     (vmac(i,j,k) + w0_edgey)*sedgey(i,j,k,comp)

             end do
          end do
       end do


       ! loop for z-fluxes
       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                bc_loz = (s0_old_cart(i,j,k,comp)+s0_old_cart(i,j,k-1,comp) &
                     +s0_new_cart(i,j,k,comp)+s0_new_cart(i,j,k-1,comp) ) * FOURTH

                if (k.eq.domlo(3)) then
                   bc_loz = HALF * (s0_old_cart(i,j,k,comp)+s0_new_cart(i,j,k,comp))
                end if
                if (k.eq.domhi(3)+1) then
                   bc_loz = HALF * (s0_old_cart(i,j,k-1,comp)+s0_new_cart(i,j,k-1,comp))
                end if

                w0_edgez = HALF * ( w0_cart(i,j,k  ,3) + w0_cart(i,j,k-1,3) )

                sfluxz(i,j,k,comp) = bc_loz*wmac(i,j,k) + &
                     (wmac(i,j,k) + w0_edgez)*sedgez(i,j,k,comp)

             end do
          end do
       end do

    end do ! end loop over components
     
  end subroutine mkflux_3d_sphr
   
end module mkflux_module
