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
  use cell_to_edge_module
  use network, ONLY: nspec

  implicit none

  private
  public :: mkflux
  
contains

  subroutine mkflux(nlevs,sflux,sold,sedge,umac,w0,w0_cart_vec,s0_old,s0_new,dx,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(in   ) :: sold(:),sedge(:,:),umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,0:,:),s0_new(:,0:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(inout) :: mla

    ! local
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
    real(kind=dp_t), pointer :: w0op(:,:,:,:)
    real(kind=dp_t), pointer :: w0np(:,:,:,:)

    real(kind=dp_t), allocatable :: s0_edge_old(:,:,:)
    real(kind=dp_t), allocatable :: s0_edge_new(:,:,:)

    ! for spherical case
    type(multifab), allocatable :: s0_old_cart(:)
    type(multifab), allocatable :: s0_new_cart(:)

    allocate(s0_edge_old(nlevs,0:nr(nlevs),nscal))
    allocate(s0_edge_new(nlevs,0:nr(nlevs),nscal))

    allocate(s0_old_cart(nlevs))
    allocate(s0_new_cart(nlevs))

    dm = sold(1)%dim

    do n = 1, nlevs
       call cell_to_edge_allcomps(n,s0_old(n,:,:),s0_edge_old(n,:,:))
       call cell_to_edge_allcomps(n,s0_new(n,:,:),s0_edge_new(n,:,:))
    enddo

    do n = 1, nlevs
       call build(s0_old_cart(n), sold(n)%la, nscal, 1)
       call build(s0_new_cart(n), sold(n)%la, nscal, 1)
       call setval(s0_old_cart(n),ZERO,all=.true.)
       call setval(s0_new_cart(n),ZERO,all=.true.)
    enddo

    ! Define s0_old_cart and s0_new_cart
    if (spherical .eq. 1) then
       call fill_3d_data_wrapper(nlevs,s0_old_cart,s0_old(:,:,rhoh_comp),dx,rhoh_comp)
       call fill_3d_data_wrapper(nlevs,s0_new_cart,s0_new(:,:,rhoh_comp),dx,rhoh_comp)
       do comp = spec_comp, spec_comp+nspec-1
          call fill_3d_data_wrapper(nlevs,s0_old_cart,s0_old(:,:,comp),dx,comp)
          call fill_3d_data_wrapper(nlevs,s0_new_cart,s0_new(:,:,comp),dx,comp)
       enddo
       do comp = trac_comp, trac_comp+ntrac-1
          call fill_3d_data_wrapper(nlevs,s0_old_cart,s0_old(:,:,comp),dx,comp)
          call fill_3d_data_wrapper(nlevs,s0_new_cart,s0_new(:,:,comp),dx,comp)
       enddo
    end if
    
    do n=1,nlevs

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
                            lo,hi)
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
                                    lo,hi)

             else
                ! call mkflux_3d_sphr not written yet
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

    deallocate(s0_edge_old,s0_edge_new)
    
  end subroutine mkflux
  
  subroutine mkflux_2d(sfluxx,sfluxy,sedgex,sedgey,umac,vmac, &
                       s0_old,s0_edge_old,s0_new,s0_edge_new,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sfluxx(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(inout) :: sfluxy(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) :: sedgex(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) :: sedgey(lo(1)  :,lo(2)  :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)





  end subroutine mkflux_2d
  
  subroutine mkflux_3d_cart(sfluxx,sfluxy,sfluxz,sedgex,sedgey,sedgez,umac,vmac,wmac, &
                            s0_old,s0_edge_old,s0_new,s0_edge_new,lo,hi)

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




     
  end subroutine mkflux_3d_cart
   
end module mkflux_module
