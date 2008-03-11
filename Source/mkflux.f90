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
    use variables, only : spec_comp, rho_comp, rhoh_comp
    use network, only : nspec
    use probin_module, only: predict_X_at_edges, predict_h_at_edges

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::  sfluxx(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) ::  sfluxy(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) :: etaflux(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) ::  sedgex(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) ::  sedgey(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) ::    umac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) ::    vmac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_pred_edge(0:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: comp
    integer :: i,j
    real(kind=dp_t) :: vel_avg,w0_avg,s0_edge
    real(kind=dp_t) :: rho_prime, rho0_edge
    logical :: test
    
    ! loop over components
    do comp = startcomp, endcomp

       test = ((comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1).and.predict_X_at_edges) &
         .or. ((comp.eq.rhoh_comp).and.predict_h_at_edges)
       
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
                
                etaflux(i,j,comp) = &
                     sfluxy(i,j,comp) - w0(j)*s0_pred_edge(j,comp)*s0_pred_edge(j,rho_comp)
                
             end do
             
          end do

       end if
             
    end do

  end subroutine mkflux_2d

  subroutine make_edge_spec_1d(lo,s,sedgex,umac,dx,dt)

     use probin_module, only: slope_order
     use bl_constants_module
     
     integer        , intent(in   ) :: lo
     real(kind=dp_t), intent(in   ) ::      s( 0  :)
     real(kind=dp_t), intent(inout) :: sedgex(lo  :)
     real(kind=dp_t), intent(in   ) ::   umac(lo-1:)
     real(kind=dp_t), intent(in   ) :: dx,dt
     
     real(kind=dp_t), allocatable::  slopex(:)
     real(kind=dp_t), allocatable::  s_l(:),s_r(:)
     real(kind=dp_t), allocatable:: dxscr(:,:)
     real(kind=dp_t) :: dmin,dpls,ds,del,slim,sflag
     real(kind=dp_t) :: ubardth, dth, savg
     real(kind=dp_t) :: abs_eps, eps, umax, u
     
     integer :: i,is,ie,hi,nx,nr
     integer :: istart,iend
     integer        , parameter :: cen = 1, lim = 2, flag = 3, fromm = 4
     real(kind=dp_t), parameter :: fourthirds = 4.0_dp_t / 3.0_dp_t
     
     nr = size(s,dim=1)
     nx = size(sedgex,dim=1)-1
     hi = lo + (nx-1)
     
     allocate(s_l(lo-1:hi+2),s_r(lo-1:hi+2))
     allocate(slopex(lo-1:hi+1))
     allocate(dxscr(lo-2:hi+2,4))

     ! Default to zero for physical boundaries and slope_order = 0.
     slopex = ZERO
     sedgex = ZERO

     ! Default to zero for physical boundaries
     dxscr(:,:) = ZERO

     abs_eps = 1.0d-8
     
     dth = HALF*dt
     
     is = lo
     ie = hi
     
     umax = ZERO
     do i = is,ie+1
        umax = max(umax,abs(umac(i)))
     end do
     
     eps = abs_eps * umax

     if (slope_order .eq. 2) then

        if (is .eq. 0) then
          istart = is+1
        else 
          istart = is-1
        end if
        if (hi .eq. nr-1) then
          iend = ie-1
        else 
          iend = ie+1
        end if

        do i = istart,iend
           del = half*(s(i+1) - s(i-1))
           dpls = two*(s(i+1) - s(i  ))
           dmin = two*(s(i  ) - s(i-1))
           slim = min(abs(dpls), abs(dmin))
           slim = merge(slim, zero, dpls*dmin.gt.ZERO)
           sflag = sign(one,del)
           slopex(i)= sflag*min(slim,abs(del))
        enddo

     else if (slope_order .eq. 4) then

        if (is .eq. 0) then
          istart = is+1
        else 
          istart = is-2
        end if
        if (hi .eq. nr-1) then
          iend = ie-1
        else 
          iend = ie+2
        end if
     
        do i = istart,iend
           dxscr(i,cen) = half*(s(i+1)-s(i-1))
           dpls = two*(s(i+1)-s(i  ))
           dmin = two*(s(i  )-s(i-1))
           dxscr(i,lim)= min(abs(dmin),abs(dpls))
           dxscr(i,lim) = merge(dxscr(i,lim),zero,dpls*dmin.gt.ZERO)
           dxscr(i,flag) = sign(one,dxscr(i,cen))
           dxscr(i,fromm)= dxscr(i,flag)*min(dxscr(i,lim),abs(dxscr(i,cen)))
        enddo
     
        istart = min(istart+1,is+1)
        iend   = max(iend  -1,ie-1)
        do i = istart,iend
           ds = fourthirds * dxscr(i,cen) - sixth * (dxscr(i+1,fromm) + dxscr(i-1,fromm))
           slopex(i) = dxscr(i,flag)*min(abs(ds),dxscr(i,lim))
        enddo

     end if
        
     ! Compute edge values using slopes and forcing terms.
     if (is .eq. 0) then
       istart = is
     else 
       istart = is-1
     end if
     if (hi .eq. nr-1) then
       iend = ie
     else 
       iend = ie+1
     end if

     do i = istart,iend
        
        u = HALF * (umac(i) + umac(i+1))
        ubardth = dth*u/dx
        
        s_l(i+1)= s(i) + (HALF-ubardth)*slopex(i)
        s_r(i  )= s(i) - (HALF+ubardth)*slopex(i)
        
     enddo
     
     if (is .eq.    0) sedgex(is  ) = s_r(is  )
     if (hi .eq. nr-1) sedgex(ie+1) = s_l(ie+1)
     
     do i = istart+1, iend
        sedgex(i)=merge(s_l(i),s_r(i),umac(i).gt.ZERO)
        savg = HALF*(s_r(i) + s_l(i))
        sedgex(i)=merge(savg,sedgex(i),abs(umac(i)) .lt. eps)
     enddo
     
  end subroutine make_edge_spec_1d
  
  subroutine mkflux_3d_cart(sfluxx,sfluxy,sfluxz,etaflux,sedgex,sedgey,sedgez, &
                            umac,vmac,wmac,s0_old,s0_edge_old,s0_new,s0_edge_new, &
                            w0,startcomp,endcomp,lo,hi)

    use bl_constants_module

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
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp

   ! local
    integer :: comp
    integer :: i,j,k

    ! loop over components
    do comp = startcomp, endcomp

       ! create x-fluxes
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                
                sfluxx(i,j,k,comp) = &
                     umac(i,j,k)*(sedgex(i,j,k,comp) + HALF*(s0_old(k,comp)+s0_new(k,comp)))

             end do
          end do
       end do

       ! create y-fluxes
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                
                sfluxy(i,j,k,comp) = &
                     vmac(i,j,k)*(sedgey(i,j,k,comp) + HALF*(s0_old(k,comp)+s0_new(k,comp)))

             end do
          end do
       end do

       ! create z-fluxes
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                
                sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))*sedgez(i,j,k,comp) &
                     + wmac(i,j,k)*HALF*(s0_edge_old(k,comp)+s0_edge_new(k,comp))

                etaflux(i,j,k,comp) = wmac(i,j,k)*sedgez(i,j,k,comp)
                
             end do
          end do
       end do

    end do ! end loop over components
     
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

    ! Note the umac here does NOT have w0 in it

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

                sfluxx(i,j,k,comp) = bc_lox*umac(i,j,k)
                
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

                sfluxy(i,j,k,comp) = bc_loy*vmac(i,j,k)

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

                sfluxz(i,j,k,comp) = bc_loz*wmac(i,j,k)

             end do
          end do
       end do

    end do ! end loop over components

    mult = ONE
    call addw0_3d_sphr(umac,vmac,wmac,w0_cart,lo,hi,mult)

    ! Note the umac here DOES have w0 in it

    do comp = startcomp, endcomp

       ! loop for x-fluxes
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1
                sfluxx(i,j,k,comp) = sfluxx(i,j,k,comp) + umac(i,j,k)*sedgex(i,j,k,comp)
             end do
          end do
       end do

       ! loop for y-fluxes
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)
                sfluxy(i,j,k,comp) = sfluxy(i,j,k,comp) + vmac(i,j,k)*sedgey(i,j,k,comp)
             end do
          end do
       end do

       ! loop for z-fluxes
       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                sfluxz(i,j,k,comp) = sfluxz(i,j,k,comp) + wmac(i,j,k)*sedgez(i,j,k,comp)
             end do
          end do
       end do

    end do ! end loop over components

    mult = -ONE
    call addw0_3d_sphr(umac,vmac,wmac,w0_cart,lo,hi,mult)
     
  end subroutine mkflux_3d_sphr
   
end module mkflux_module
