module mkflux_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use slope_module
  use fill_3d_module
  use geometry
  use define_bc_module
  use ml_layout_module
  use ml_restriction_module
  use variables
  use network, ONLY: nspec
  use addw0_module

  implicit none

  private
  public :: mkflux
  
contains

  subroutine mkflux(nlevs,sflux,sold,sedge,umac,w0,w0_cart_vec,s0_old,s0_edge_old, &
                    s0_old_cart,s0_new,s0_edge_new,s0_new_cart, &
                    startcomp,endcomp,which_step,dx,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(in   ) :: sold(:),sedge(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,0:,:),s0_edge_old(:,0:,:)
    type(multifab) , intent(in   ) :: s0_old_cart(:)
    real(kind=dp_t), intent(in   ) :: s0_new(:,0:,:),s0_edge_new(:,0:,:)
    type(multifab) , intent(in   ) :: s0_new_cart(:)
    integer        , intent(in   ) :: startcomp,endcomp,which_step
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(inout) :: mla

    ! local    
    type(box) :: domain

    integer :: domlo(sold(1)%dim),domhi(sold(1)%dim)
    integer :: i,dm,n,comp
    integer :: lo(sold(1)%dim),hi(sold(1)%dim)

    real(kind=dp_t), pointer :: sfxp(:,:,:,:)
    real(kind=dp_t), pointer :: sfyp(:,:,:,:)
    real(kind=dp_t), pointer :: sfzp(:,:,:,:)
    real(kind=dp_t), pointer :: sexp(:,:,:,:)
    real(kind=dp_t), pointer :: seyp(:,:,:,:)
    real(kind=dp_t), pointer :: sezp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)
    real(kind=dp_t), pointer :: s0op(:,:,:,:)
    real(kind=dp_t), pointer :: s0np(:,:,:,:)

    dm = sold(1)%dim
    
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
             call mkflux_2d(sfxp(:,:,1,:), sfyp(:,:,1,:), &
                            sexp(:,:,1,:), seyp(:,:,1,:), &
                            ump(:,:,1,1), vmp(:,:,1,1), &
                            s0_old(n,:,:), s0_edge_old(n,:,:), &
                            s0_new(n,:,:), s0_edge_new(n,:,:), &
                            w0(n,:), &
                            startcomp,endcomp,which_step,lo,hi)
          case (3)
             sfzp => dataptr(sflux(n,3),i)
             sezp => dataptr(sedge(n,3),i)
             wmp  => dataptr(umac(n,3),i)
             if(spherical .eq. 0) then
                call mkflux_3d_cart(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), &
                                    sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                    s0_old(n,:,:), s0_edge_old(n,:,:), &
                                    s0_new(n,:,:), s0_edge_new(n,:,:), &
                                    w0(n,:), &
                                    startcomp,endcomp,which_step,lo,hi)

             else
                s0op => dataptr(s0_old_cart(n), i)
                s0np => dataptr(s0_new_cart(n), i)
                w0p => dataptr(w0_cart_vec(n),i)
                call mkflux_3d_sphr(sfxp(:,:,:,:), sfyp(:,:,:,:), sfzp(:,:,:,:), &
                                    sexp(:,:,:,:), seyp(:,:,:,:), sezp(:,:,:,:), &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                    s0_old(n,:,:), s0_edge_old(n,:,:), s0op(:,:,:,:), &
                                    s0_new(n,:,:), s0_edge_new(n,:,:), s0np(:,:,:,:), &
                                    w0(n,:), w0p(:,:,:,:), &
                                    startcomp,endcomp,which_step,lo,hi,domlo,domhi)
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
    
  end subroutine mkflux
  
  subroutine mkflux_2d(sfluxx,sfluxy,sedgex,sedgey,umac,vmac,s0_old,s0_edge_old, &
                       s0_new,s0_edge_new,w0,startcomp,endcomp,which_step,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sfluxx(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) :: sfluxy(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) :: sedgex(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) :: sedgey(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp,which_step

    ! local
    integer :: comp
    integer :: i,j

    ! loop over components
    do comp = startcomp, endcomp

       ! create x-fluxes
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             if(which_step .eq. 1) then
                sfluxx(i,j,comp) = umac(i,j)* &
                     (sedgex(i,j,comp) + s0_old(j,comp))
             else
                sfluxx(i,j,comp) = umac(i,j)* &
                     (sedgex(i,j,comp) + HALF*(s0_old(j,comp)+s0_new(j,comp)))
             endif

          end do
       end do

       ! create y-fluxes
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             if(which_step .eq. 1) then
                sfluxy(i,j,comp) = (vmac(i,j)+w0(j))*sedgey(i,j,comp) &
                     + vmac(i,j)*(s0_edge_old(j,comp))
             else
                sfluxy(i,j,comp) = (vmac(i,j)+w0(j))*sedgey(i,j,comp) &
                     + vmac(i,j)*(HALF*(s0_edge_old(j,comp)+s0_edge_new(j,comp)))
             end if

          end do
       end do

    end do ! end loop over components

  end subroutine mkflux_2d
  
  subroutine mkflux_3d_cart(sfluxx,sfluxy,sfluxz,sedgex,sedgey,sedgez,umac,vmac,wmac, &
                            s0_old,s0_edge_old,s0_new,s0_edge_new,w0,startcomp,endcomp, &
                            which_step,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sfluxx(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: sfluxy(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: sfluxz(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) :: sedgex(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) :: sedgey(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) :: sedgez(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: startcomp,endcomp,which_step

   ! local
    integer :: comp
    integer :: i,j,k

    ! loop over components
    do comp = startcomp, endcomp

       ! create x-fluxes
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                
                if(which_step .eq. 1) then
                   sfluxx(i,j,k,comp) = umac(i,j,k)* &
                        (sedgex(i,j,k,comp) + s0_old(k,comp))
                else
                   sfluxx(i,j,k,comp) = umac(i,j,k)* &
                        (sedgex(i,j,k,comp) + HALF*(s0_old(k,comp)+s0_new(k,comp)))
                endif
                
             end do
          end do
       end do

       ! create y-fluxes
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                
                if(which_step .eq. 1) then
                   sfluxy(i,j,k,comp) = vmac(i,j,k)* &
                        (sedgey(i,j,k,comp) + s0_old(k,comp))
                else
                   sfluxy(i,j,k,comp) = vmac(i,j,k)* &
                        (sedgey(i,j,k,comp) + HALF*(s0_old(k,comp)+s0_new(k,comp)))
                endif

             end do
          end do
       end do

       ! create z-fluxes
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                
                if(which_step .eq. 1) then
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))*sedgez(i,j,k,comp) &
                        + wmac(i,j,k)*s0_edge_old(k,comp)
                else
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(k))*sedgez(i,j,k,comp) &
                        + wmac(i,j,k)*HALF*(s0_edge_old(k,comp)+s0_edge_new(k,comp))
                end if
                
             end do
          end do
       end do

    end do ! end loop over components
     
  end subroutine mkflux_3d_cart

  subroutine mkflux_3d_sphr(sfluxx,sfluxy,sfluxz,sedgex,sedgey,sedgez,umac,vmac,wmac, &
                            s0_old,s0_edge_old,s0_old_cart,s0_new,s0_edge_new,s0_new_cart, &
                            w0,w0_cart,startcomp,endcomp,which_step,lo,hi,domlo,domhi)

    integer        , intent(in   ) :: lo(:),hi(:),domlo(:),domhi(:)
    real(kind=dp_t), intent(inout) :: sfluxx(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: sfluxy(lo(1)  :,lo(2)  :,lo(3)  :,:)
    real(kind=dp_t), intent(inout) :: sfluxz(lo(1)  :,lo(2)  :,lo(3)  :,:)
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
    integer        , intent(in   ) :: startcomp,endcomp,which_step

    ! local
    integer         :: i,j,k,comp
    real(kind=dp_t) :: mult
    real(kind=dp_t) :: bc_lox,bc_loy,bc_loz
    real(kind=dp_t) :: bc_hix,bc_hiy,bc_hiz

    ! Note the umac here does NOT have w0 in it

    do comp = startcomp, endcomp

       ! loop for x-fluxes
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                if (which_step .eq. 1) then

                   bc_lox = (s0_old_cart(i,j,k,comp)+s0_old_cart(i-1,j,k,comp)) * HALF
                   bc_hix = (s0_old_cart(i,j,k,comp)+s0_old_cart(i+1,j,k,comp)) * HALF
                   
                   if (i.eq.domlo(1)) bc_lox = s0_old_cart(i,j,k,comp)
                   if (i.eq.domhi(1)+1) bc_hix = s0_old_cart(i,j,k,comp)

                else if (which_step .eq. 2) then

                   bc_lox = (s0_old_cart(i,j,k,comp)+s0_old_cart(i-1,j,k,comp) &
                        +s0_new_cart(i,j,k,comp)+s0_new_cart(i-1,j,k,comp) ) * FOURTH
                   bc_hix = (s0_old_cart(i,j,k,comp)+s0_old_cart(i+1,j,k,comp) &
                        +s0_new_cart(i,j,k,comp)+s0_new_cart(i+1,j,k,comp) ) * FOURTH

                   if (i.eq.domlo(1)) bc_lox = &
                        HALF * (s0_old_cart(i,j,k,comp)+s0_new_cart(i,j,k,comp))
                   if (i.eq.domhi(1)+1) bc_hix = &
                        HALF * (s0_old_cart(i,j,k,comp)+s0_new_cart(i,j,k,comp))

                end if

                sfluxx(i,j,k,comp) = bc_hix*umac(i+1,j,k) - bc_lox*umac(i,j,k)
                
             end do
          end do
       end do

       ! loop for y-fluxes
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                if (which_step .eq. 1) then

                   bc_loy = (s0_old_cart(i,j,k,comp)+s0_old_cart(i,j-1,k,comp)) * HALF
                   bc_hiy = (s0_old_cart(i,j,k,comp)+s0_old_cart(i,j+1,k,comp)) * HALF
                   
                   if (j.eq.domlo(2)) bc_loy = s0_old_cart(i,j,k,comp)
                   if (j.eq.domhi(2)+1) bc_hiy = s0_old_cart(i,j,k,comp)

                else if (which_step .eq. 2) then

                   bc_loy = (s0_old_cart(i,j,k,comp)+s0_old_cart(i,j-1,k,comp) &
                        +s0_new_cart(i,j,k,comp)+s0_new_cart(i,j-1,k,comp) ) * FOURTH
                   bc_hiy = (s0_old_cart(i,j,k,comp)+s0_old_cart(i,j+1,k,comp) &
                        +s0_new_cart(i,j,k,comp)+s0_new_cart(i,j+1,k,comp) ) * FOURTH

                   if (j.eq.domlo(2)) bc_loy = &
                        HALF * (s0_old_cart(i,j,k,comp)+s0_new_cart(i,j,k,comp))
                   if (j.eq.domhi(2)+1) bc_hiy = &
                        HALF * (s0_old_cart(i,j,k,comp)+s0_new_cart(i,j,k,comp))

                end if

                sfluxy(i,j,k,comp) = bc_hiy*vmac(i,j+1,k) - bc_loy*vmac(i,j,k)

             end do
          end do
       end do

       ! loop for z-fluxes
       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (which_step .eq. 1) then

                   bc_loz = (s0_old_cart(i,j,k,comp)+s0_old_cart(i,j,k-1,comp)) * HALF
                   bc_hiz = (s0_old_cart(i,j,k,comp)+s0_old_cart(i,j,k+1,comp)) * HALF
                   
                   if (k.eq.domlo(3)) bc_loz = s0_old_cart(i,j,k,comp)
                   if (k.eq.domhi(3)+1) bc_hiz = s0_old_cart(i,j,k,comp)

                else if (which_step .eq. 2) then

                   bc_loz = (s0_old_cart(i,j,k,comp)+s0_old_cart(i,j,k-1,comp) &
                        +s0_new_cart(i,j,k,comp)+s0_new_cart(i,j,k-1,comp) ) * FOURTH
                   bc_hiz = (s0_old_cart(i,j,k,comp)+s0_old_cart(i,j,k+1,comp) &
                        +s0_new_cart(i,j,k,comp)+s0_new_cart(i,j,k+1,comp) ) * FOURTH

                   if (k.eq.domlo(3)) bc_loz = &
                        HALF * (s0_old_cart(i,j,k,comp)+s0_new_cart(i,j,k,comp))
                   if (k.eq.domhi(3)+1) bc_hiz = &
                        HALF * (s0_old_cart(i,j,k,comp)+s0_new_cart(i,j,k,comp))

                end if

                sfluxz(i,j,k,comp) = bc_hiz*wmac(i,j,k+1) - bc_loz*wmac(i,j,k)

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
