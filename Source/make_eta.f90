! Compute eta_{X_k}, eta_{rho h}, and eta_{rho}

module make_eta_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: make_eta

contains

  subroutine make_eta(nlevs,eta,sold,sflux,dx,mla)

    use bl_constants_module
    use geometry, only: spherical, nr
    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer           , intent(in   ) :: nlevs
    real(kind=dp_t)   , intent(inout) :: eta(:,0:,:)
    type(multifab)    , intent(in   ) :: sold(:)
    type(multifab)    , intent(in   ) :: sflux(:,:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:)
    type(ml_layout)   , intent(inout) :: mla

    ! local
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    
    real(kind=dp_t), allocatable :: ncell_proc(:,:)
    real(kind=dp_t), allocatable :: ncell(:,:)
    real(kind=dp_t), allocatable :: etasum_proc(:,:,:)
    real(kind=dp_t), allocatable :: etasum(:,:,:)
    real(kind=dp_t), allocatable :: source_buffer(:)
    real(kind=dp_t), allocatable :: target_buffer(:)

    type(box) :: domain

    integer :: domlo(sold(1)%dim),domhi(sold(1)%dim)
    integer :: lo(sold(1)%dim),hi(sold(1)%dim)
    integer :: i,k,n,dm,rr,comp

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_eta")

    dm = sold(1)%dim

    allocate( ncell_proc(nlevs,0:nr(nlevs)))
    allocate( ncell     (nlevs,0:nr(nlevs)))
    allocate(etasum_proc(nlevs,0:nr(nlevs),nscal))
    allocate(etasum     (nlevs,0:nr(nlevs),nscal))

    allocate(source_buffer(0:nr(nlevs)))
    allocate(target_buffer(0:nr(nlevs)))

    ncell_proc  = ZERO
    ncell       = ZERO
    etasum_proc = ZERO
    etasum      = ZERO

    if (spherical .eq. 0) then
    
       domain = layout_get_pd(sold(1)%la)
       domlo = lwb(domain)
       domhi = upb(domain)

       if (dm .eq. 2) then
          ncell(1,:) = domhi(1)-domlo(1)+1
       else if(dm .eq. 3) then
          ncell(1,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
       end if
       
       ! the first step is to compute eta assuming the coarsest level 
       ! is the only level in existence
       do i=1,sold(1)%nboxes
          if ( multifab_remote(sold(1), i) ) cycle
          fluxrp => dataptr(sflux(1,dm), i)
          lo =  lwb(get_box(sold(1), i))
          hi =  upb(get_box(sold(1), i))
          select case (dm)
          case (2)
            call sum_eta_coarsest_2d(lo,hi,fluxrp(:,:,1,:),etasum_proc(1,:,:))
          case (3)
            call sum_eta_coarsest_3d(lo,hi,fluxrp(:,:,:,:),etasum_proc(1,:,:))
          end select
       end do

       source_buffer = etasum_proc(1,:,rhoh_comp)
       call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
       etasum(1,:,rhoh_comp) = target_buffer

       do comp=spec_comp,spec_comp+nspec-1
          source_buffer = etasum_proc(1,:,comp)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          etasum(1,:,comp) = target_buffer
       end do

       do comp=spec_comp,spec_comp+nspec-1
          etasum(1,:,rho_comp) = etasum(1,:,rho_comp) + etasum(1,:,comp)
       end do

       do k=0,nr(1)
          eta(1,k,rhoh_comp) = etasum(1,k,rhoh_comp) / dble(ncell(1,k))
          do comp=spec_comp,spec_comp+nspec-1
             eta(1,k,comp) = etasum(1,k,comp) / dble(ncell(1,k))
          end do
          eta(1,k,rho_comp) = etasum(1,k,rho_comp) / dble(ncell(1,k))
       end do

       ! now we compute eta at the finer levels
       do n=2,nlevs
          
       end do

    else

       ! haven't written the spherical case yet

    end if


    deallocate(ncell_proc,ncell)

    call destroy(bpt)

  end subroutine make_eta

  subroutine sum_eta_coarsest_2d(lo,hi,fluxy,etasum)

    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer         , intent(in   ) :: lo(:), hi(:)
    real (kind=dp_t), intent(in   ) :: fluxy(lo(1):,lo(2):,:)
    real (kind=dp_t), intent(inout) :: etasum(0:,:)

    ! local
    integer :: i,j,comp

    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          etasum(j,rhoh_comp) = etasum(j,rhoh_comp) + fluxy(i,j,rhoh_comp)
          do comp=spec_comp,spec_comp+nspec-1
             etasum(j,comp) = etasum(j,comp) + fluxy(i,j,comp)
          end do
       end do
    end do

  end subroutine sum_eta_coarsest_2d

  subroutine sum_eta_coarsest_3d(lo,hi,fluxy,etasum)

    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer         , intent(in   ) :: lo(:), hi(:)
    real (kind=dp_t), intent(in   ) :: fluxz(lo(1):,lo(2):,lo(3):,:)
    real (kind=dp_t), intent(inout) :: etasum(0:,:)

    ! local
    integer :: i,j,k,comp

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etasum(k,rhoh_comp) = etasum(k,rhoh_comp) + fluxz(i,j,k,rhoh_comp)
             do comp=spec_comp,spec_comp+nspec-1
                etasum(k,comp) = etasum(k,comp) + fluxz(i,j,k,comp)
             end do
          end do
       end do
    end do

  end subroutine sum_eta_coarsest_3d

end module make_eta_module
