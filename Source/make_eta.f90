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
    type(multifab)    , intent(inout) :: sflux(:,:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:)
    type(ml_layout)   , intent(inout) :: mla

    ! local
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: fpc(:,:,:,:)
    
    real(kind=dp_t), allocatable :: ncell_proc(:,:)
    real(kind=dp_t), allocatable :: ncell(:,:)
    real(kind=dp_t), allocatable :: etasum_proc(:,:,:)
    real(kind=dp_t), allocatable :: etasum(:,:,:)
    real(kind=dp_t), allocatable :: etapert_proc(:,:,:)
    real(kind=dp_t), allocatable :: etapert(:,:,:)
    real(kind=dp_t), allocatable :: source_buffer(:)
    real(kind=dp_t), allocatable :: target_buffer(:)

    type(box) :: domain

    integer :: domlo(sold(1)%dim),domhi(sold(1)%dim)
    integer :: lo(sold(1)%dim),hi(sold(1)%dim)
    integer :: i,k,n,dm,rr,comp

    type(layout) :: lasfluxcoarse
    type(multifab) :: sfluxcoarse

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_eta")

    dm = sold(1)%dim

    allocate(ncell_proc(nlevs,0:nr(nlevs)))
    allocate(ncell     (nlevs,0:nr(nlevs)))

    allocate(etasum_proc(nlevs,0:nr(nlevs),nscal))
    allocate(etasum     (nlevs,0:nr(nlevs),nscal))

    allocate(etapert_proc(nlevs,0:nr(nlevs),nscal))
    allocate(etapert     (nlevs,0:nr(nlevs),nscal))

    allocate(source_buffer(0:nr(nlevs)))
    allocate(target_buffer(0:nr(nlevs)))

    ncell_proc   = ZERO
    ncell        = ZERO
    etasum_proc  = ZERO
    etasum       = ZERO
    etapert_proc = ZERO
    etapert      = ZERO

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
          fp => dataptr(sflux(1,dm), i)
          lo =  lwb(get_box(sold(1), i))
          hi =  upb(get_box(sold(1), i))
          select case (dm)
          case (2)
            call sum_eta_coarsest_2d(lo,hi,domhi,fp(:,:,1,:),etasum_proc(1,:,:))
          case (3)
            call sum_eta_coarsest_3d(lo,hi,domhi,fp(:,:,:,:),etasum_proc(1,:,:))
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

       do k=0,nr(1)
          eta(1,k,rhoh_comp) = etasum(1,k,rhoh_comp) / dble(ncell(1,k))
          do comp=spec_comp,spec_comp+nspec-1
             eta(1,k,comp) = etasum(1,k,comp) / dble(ncell(1,k))
             eta(1,k,rho_comp) = eta(1,k,rho_comp) + eta(1,k,comp)
          end do
       end do

       ! now we compute eta at the finer levels
       do n=2,nlevs
          
          rr = mla%mba%rr(n-1,dm)
          if(rr .ne. 2) then
             print*,"Error: in make_eta, refinement ratio must be 2"
             stop
          end if

          domain = layout_get_pd(sold(n)%la)
          domlo  = lwb(domain)
          domhi  = upb(domain)

          if (dm .eq. 2) then
             ncell(n,:) = domhi(1)-domlo(1)+1
          else if (dm .eq. 3) then
             ncell(n,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
          end if
          
          ! on faces that exist at the next coarser level, eta is the same since
          ! the fluxes have been restricted.  We copy these values of eta directly and
          ! copy scaled values of etasum directly.
          do k=0,nr(n-1)
             eta(n,k*rr,rhoh_comp) = eta(n-1,k,rhoh_comp)
             etasum(n,k*rr,rhoh_comp) = etasum(n-1,k,rhoh_comp)*rr**(dm-1)
             do comp=spec_comp,spec_comp+nspec-1
                eta(n,k*rr,comp) = eta(n-1,k,comp)
                etasum(n,k*rr,comp) = etasum(n-1,k,comp)*rr**(dm-1)
             end do
          end do

          ! on faces that do not exist at the next coarser level, we use linear
          ! interpolation to get etasum at these faces.
          do k=1,nr(n)-1,2
             etasum(n,k,rhoh_comp) = HALF*(etasum(n,k-1,rhoh_comp)+etasum(n,k+1,rhoh_comp))
             do comp=spec_comp,spec_comp+nspec-1
                etasum(n,k,comp) = HALF*(etasum(n,k-1,comp)+etasum(n,k+1,comp))
             end do
          end do

          ! create a temporary coarse multifab that corresponds directly to regions 
          ! where the finer flux, sflux(n,dm), exists
          call layout_build_coarse(lasfluxcoarse, sflux(n,dm)%la, mla%mba%rr(n-1,:))
          call multifab_build(sfluxcoarse, lasfluxcoarse, nc=sflux(n,dm)%nc, ng=0, &
                              nodal=sflux(n,dm)%nodal)

          ! copy data from the coarse multifab into my temporary coarse multifab
          call copy(sfluxcoarse, 1, sflux(n-1,dm), 1, nscal)

          ! compute etapert_proc on faces that do not exist at the coarser level
          do i=1,sold(n)%nboxes
             if ( multifab_remote(sold(n), i) ) cycle
             fp  => dataptr(sflux(n,dm), i)
             fpc  => dataptr(sfluxcoarse, i)
             lo =  lwb(get_box(sold(n), i))
             hi =  upb(get_box(sold(n), i))
             select case (dm)
             case (2)
                call compute_etapert_2d(lo,hi,fp(:,:,1,:),fpc(:,:,1,:), &
                                        etasum_proc(1,:,:),rr)
             case (3)
                call compute_etapert_3d(lo,hi,fp(:,:,:,:),fpc(:,:,:,:), &
                                        etasum_proc(1,:,:),rr)
             end select
          end do

          ! gather etapert for rhoh
          source_buffer = etapert_proc(n,:,rhoh_comp)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          etapert(n,:,rhoh_comp) = target_buffer

          ! gather etapert for rhoX
          do comp=spec_comp,spec_comp+nspec-1
             source_buffer = etapert_proc(n,:,comp)
             call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
             etapert(n,:,comp) = target_buffer
          end do

          ! update etasum on faces that do not exist at the coarser level
          ! then recompute eta on these faces
          do k=1,nr(n)-1,2
             etasum(n,k,rhoh_comp) = etasum(n,k,rhoh_comp) + etapert(n,k,rhoh_comp)
             eta(n,k,rhoh_comp) = etasum(n,k,rhoh_comp) / dble(ncell(n,k))
             do comp=spec_comp,spec_comp+nspec-1
                etasum(n,k,comp) = etasum(n,k,comp) + etapert(n,k,comp)
                eta(n,k,comp) = etasum(n,k,comp) / dble(ncell(n,k))
                eta(n,k,rho_comp) = eta(n,k,rho_comp) + eta(n,k,comp)
             end do
          end do

       end do ! end loop over levels

    else

       ! haven't written the spherical case yet

    end if

    deallocate(ncell_proc,ncell)
    deallocate(etasum_proc,etasum)
    deallocate(etapert_proc,etapert)
    deallocate(source_buffer,target_buffer)

    call destroy(bpt)

  end subroutine make_eta

  subroutine sum_eta_coarsest_2d(lo,hi,domhi,fluxr,etasum)

    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer         , intent(in   ) :: lo(:), hi(:), domhi(:)
    real (kind=dp_t), intent(in   ) :: fluxr(lo(1):,lo(2):,:)
    real (kind=dp_t), intent(inout) :: etasum(0:,:)

    ! local
    integer :: i,j,comp

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          etasum(j,rhoh_comp) = etasum(j,rhoh_comp) + fluxr(i,j,rhoh_comp)
          do comp=spec_comp,spec_comp+nspec-1
             etasum(j,comp) = etasum(j,comp) + fluxr(i,j,comp)
          end do
       end do
    end do

    ! we only add the contribution at the top edge if we are at the top of the domain
    ! this prevents double counting
    if(hi(2) .eq. domhi(2)) then
       j=hi(2)+1
       do i=lo(1),hi(1)
          etasum(j,rhoh_comp) = etasum(j,rhoh_comp) + fluxr(i,j,rhoh_comp)
          do comp=spec_comp,spec_comp+nspec-1
             etasum(j,comp) = etasum(j,comp) + fluxr(i,j,comp)
          end do
       end do
    end if

  end subroutine sum_eta_coarsest_2d

  subroutine sum_eta_coarsest_3d(lo,hi,domhi,fluxr,etasum)

    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer         , intent(in   ) :: lo(:), hi(:), domhi(:)
    real (kind=dp_t), intent(in   ) :: fluxr(lo(1):,lo(2):,lo(3):,:)
    real (kind=dp_t), intent(inout) :: etasum(0:,:)

    ! local
    integer :: i,j,k,comp

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etasum(k,rhoh_comp) = etasum(k,rhoh_comp) + fluxr(i,j,k,rhoh_comp)
             do comp=spec_comp,spec_comp+nspec-1
                etasum(k,comp) = etasum(k,comp) + fluxr(i,j,k,comp)
             end do
          end do
       end do
    end do

    ! we only add the contribution at the top edge if we are at the top of the domain
    ! this prevents double counting
    if(hi(3) .eq. domhi(3)) then
       k=hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etasum(k,rhoh_comp) = etasum(k,rhoh_comp) + fluxr(i,j,k,rhoh_comp)
             do comp=spec_comp,spec_comp+nspec-1
                etasum(k,comp) = etasum(k,comp) + fluxr(i,j,k,comp)
             end do
          end do
       end do
    end if

  end subroutine sum_eta_coarsest_3d

  subroutine compute_etapert_2d(lo,hi,fluxr,fluxrc,etapert,rr)

    use variables, only: rhoh_comp, spec_comp
    use network, only: nspec
    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), rr
    real (kind=dp_t), intent(in   ) ::  fluxr(lo(1):,lo(2):,:)
    real (kind=dp_t), intent(in   ) :: fluxrc(lo(1):,lo(2):,:)
    real (kind=dp_t), intent(inout) :: etapert(0:,:)
    
    ! local
    integer         :: i,j,comp
    real(kind=dp_t) :: crseval

    do j=lo(2)+1,hi(2)-1,2
       do i=lo(1),hi(1)

          crseval = (fluxrc(i/rr,(j-1)/rr,rhoh_comp) + fluxrc(i/rr,(j-1)/rr+1,rhoh_comp))/TWO
          etapert(j,rhoh_comp) = fluxr(i,j,rhoh_comp) - crseval

          do comp=spec_comp,spec_comp+nspec-1
             crseval = (fluxrc(i/rr,(j-1)/rr,comp) + fluxrc(i/rr,(j-1)/rr+1,comp))/TWO
             etapert(j,comp) = fluxr(i,j,comp) - crseval
          end do

       end do
    end do


  end subroutine compute_etapert_2d

  subroutine compute_etapert_3d(lo,hi,fluxr,fluxrc,etapert,rr)

    use variables, only: rhoh_comp, spec_comp
    use network, only: nspec
    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), rr
    real (kind=dp_t), intent(in   ) ::  fluxr(lo(1):,lo(2):,lo(3):,:)
    real (kind=dp_t), intent(in   ) :: fluxrc(lo(1):,lo(2):,lo(3):,:)
    real (kind=dp_t), intent(inout) :: etapert(0:,:)
    
    ! local
    integer         :: i,j,k,comp
    real(kind=dp_t) :: crseval

    do k=lo(3)+1,hi(3)-1,2
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             crseval = (fluxrc(i/rr,j/rr,(k-1)/rr,rhoh_comp) + &
                  fluxrc(i/rr,j/rr,(k-1)/rr+1,rhoh_comp))/TWO
             etapert(k,rhoh_comp) = fluxr(i,j,k,rhoh_comp) - crseval
             
             do comp=spec_comp,spec_comp+nspec-1
                crseval = (fluxrc(i/rr,j/rr,(k-1)/rr,comp) + &
                     fluxrc(i/rr,j/rr,(k-1)/rr+1,comp))/TWO
                etapert(k,comp) = fluxr(i,j,k,comp) - crseval
             end do
             
          end do
       end do
    end do

  end subroutine compute_etapert_3d

end module make_eta_module
