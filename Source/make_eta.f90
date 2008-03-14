! Compute eta_{X_k}, eta_{rho h}, and eta_{rho}

module make_eta_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: make_eta

contains

  subroutine make_eta(nlevs,eta,etaflux,mla)

    use bl_constants_module
    use geometry, only: spherical, nr
    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer           , intent(in   ) :: nlevs
    real(kind=dp_t)   , intent(inout) :: eta(:,0:,:)
    type(multifab)    , intent(inout) :: etaflux(:)
    type(ml_layout)   , intent(inout) :: mla

    ! local
    real(kind=dp_t), pointer :: efp(:,:,:,:)
    
    real(kind=dp_t), allocatable :: ncell(:,:)
    real(kind=dp_t), allocatable :: etasum_proc(:,:,:)
    real(kind=dp_t), allocatable :: etasum(:,:,:)
    real(kind=dp_t), allocatable :: etapert_proc(:,:,:)
    real(kind=dp_t), allocatable :: etapert(:,:,:)
    real(kind=dp_t), allocatable :: source_buffer(:)
    real(kind=dp_t), allocatable :: target_buffer(:)

    type(box) :: domain

    integer :: domlo(mla%dim),domhi(mla%dim)
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: i,r,rpert,n,dm,rr,comp

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_eta")

    dm = mla%dim

    allocate(ncell       (nlevs,0:nr(nlevs))) ! ncell is a function of r only for spherical
    allocate(etasum_proc (nlevs,0:nr(nlevs),nscal))
    allocate(etasum      (nlevs,0:nr(nlevs),nscal))
    allocate(etapert_proc(nlevs,0:nr(nlevs),nscal))
    allocate(etapert     (nlevs,0:nr(nlevs),nscal))

    allocate(source_buffer(0:nr(nlevs)))
    allocate(target_buffer(0:nr(nlevs)))

    ncell        = ZERO
    etasum_proc  = ZERO
    etasum       = ZERO
    etapert_proc = ZERO
    etapert      = ZERO
    eta          = ZERO

    if (spherical .eq. 0) then
    
       domain = layout_get_pd(mla%la(1))
       domlo = lwb(domain)
       domhi = upb(domain)

       if (dm .eq. 2) then
          ncell(1,:) = domhi(1)-domlo(1)+1
       else if(dm .eq. 3) then
          ncell(1,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
       end if
       
       ! the first step is to compute eta assuming the coarsest level 
       ! is the only level in existence
       do i=1,layout_nboxes(mla%la(1))
          if ( multifab_remote(etaflux(1), i) ) cycle
          efp => dataptr(etaflux(1), i)
          lo =  lwb(get_box(mla%la(1), i))
          hi =  upb(get_box(mla%la(1), i))
          select case (dm)
          case (2)
            call sum_eta_coarsest_2d(lo,hi,domhi,efp(:,:,1,:),etasum_proc(1,:,:))
          case (3)
            call sum_eta_coarsest_3d(lo,hi,domhi,efp(:,:,:,:),etasum_proc(1,:,:))
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

       do r=0,nr(1)
          eta(1,r,rhoh_comp) = etasum(1,r,rhoh_comp) / dble(ncell(1,r))
          do comp=spec_comp,spec_comp+nspec-1
             eta(1,r,comp) = etasum(1,r,comp) / dble(ncell(1,r))
             eta(1,r,rho_comp) = eta(1,r,rho_comp) + eta(1,r,comp)
          end do
       end do

       ! now we compute eta at the finer levels
       do n=2,nlevs
          
          rr = mla%mba%rr(n-1,dm)

          if (mla%mba%rr(n-1,1) .ne. mla%mba%rr(n-1,dm) .or. &
               mla%mba%rr(n-1,2) .ne. mla%mba%rr(n-1,dm)) then
             print*,"ERROR: In make_eta, refinement ratio in each direction must match"
             stop
          endif

          domain = layout_get_pd(mla%la(n))
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
          do r=0,nr(n-1)
             eta   (n,r*rr,rhoh_comp) = eta   (n-1,r,rhoh_comp)
             etasum(n,r*rr,rhoh_comp) = etasum(n-1,r,rhoh_comp)*rr**(dm-1)
             do comp=spec_comp,spec_comp+nspec-1
                eta   (n,r*rr,comp)     = eta   (n-1,r,comp)
                etasum(n,r*rr,comp)     = etasum(n-1,r,comp)*rr**(dm-1)
                eta   (n,r*rr,rho_comp) = eta   (n,r*rr,rho_comp) + eta(n,r*rr,comp)
             end do
          end do

          ! on faces that do not exist at the next coarser level, we use linear
          ! interpolation to get etasum at these faces.
          do r=0,nr(n-1)-1
             do rpert=1,rr-1
                etasum(n,r*rr+rpert,rhoh_comp) = &
                     dble(rpert)/dble(rr)*(etasum(n,r*rr,rhoh_comp)) + &
                     dble(rr-rpert)/dble(rr)*(etasum(n,(r+1)*rr,rhoh_comp))
                do comp=spec_comp,spec_comp+nspec-1
                   etasum(n,r*rr+rpert,comp) = &
                        dble(rpert)/dble(rr)*(etasum(n,r*rr,comp)) + &
                        dble(rr-rpert)/dble(rr)*(etasum(n,(r+1)*rr,comp))
                end do
             end do
          end do

          ! compute etapert_proc on faces that do not exist at the coarser level
          do i=1,layout_nboxes(mla%la(n))
             if ( multifab_remote(etaflux(n), i) ) cycle
             efp  => dataptr(etaflux(n), i)
             lo =  lwb(get_box(mla%la(n), i))
             hi =  upb(get_box(mla%la(n), i))
             select case (dm)
             case (2)
                call compute_etapert_2d(lo,hi,efp(:,:,1,:),etasum_proc(1,:,:),rr)
             case (3)
                call compute_etapert_3d(lo,hi,efp(:,:,:,:),etasum_proc(1,:,:),rr)
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
          do r=0,nr(n-1)-1
             do rpert=1,rr-1
                etasum(n,r*rr+rpert,rhoh_comp) = &
                     etasum(n,r*rr+rpert,rhoh_comp) + etapert(n,r*rr+rpert,rhoh_comp)
                eta(n,r*rr+rpert,rhoh_comp) = &
                     etasum(n,r*rr+rpert,rhoh_comp)/dble(ncell(n,r*rr+rpert))
                do comp=spec_comp,spec_comp+nspec-1
                   etasum(n,r*rr+rpert,comp) = &
                        etasum(n,r*rr+rpert,comp) + etapert(n,r*rr+rpert,comp)
                   eta(n,r*rr+rpert,comp) = &
                        etasum(n,r*rr+rpert,comp)/dble(ncell(n,r*rr+rpert))
                   eta(n,r*rr+rpert,rho_comp) = &
                        eta(n,r*rr+rpert,rho_comp) + eta(n,r*rr+rpert,comp)
                end do
             end do
          end do

       end do ! end loop over levels

    else

       ! haven't written the spherical case yet

    end if

    deallocate(ncell)
    deallocate(etasum_proc,etasum)
    deallocate(etapert_proc,etapert)
    deallocate(source_buffer,target_buffer)

    call destroy(bpt)

  end subroutine make_eta

  subroutine sum_eta_coarsest_2d(lo,hi,domhi,etaflux,etasum)

    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer         , intent(in   ) :: lo(:), hi(:), domhi(:)
    real (kind=dp_t), intent(in   ) :: etaflux(lo(1):,lo(2):,:)
    real (kind=dp_t), intent(inout) :: etasum(0:,:)

    ! local
    integer :: i,j,comp

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          etasum(j,rhoh_comp) = etasum(j,rhoh_comp) + etaflux(i,j,rhoh_comp)
          do comp=spec_comp,spec_comp+nspec-1
             etasum(j,comp) = etasum(j,comp) + etaflux(i,j,comp)
          end do
       end do
    end do

    ! we only add the contribution at the top edge if we are at the top of the domain
    ! this prevents double counting
    if(hi(2) .eq. domhi(2)) then
       j=hi(2)+1
       do i=lo(1),hi(1)
          etasum(j,rhoh_comp) = etasum(j,rhoh_comp) + etaflux(i,j,rhoh_comp)
          do comp=spec_comp,spec_comp+nspec-1
             etasum(j,comp) = etasum(j,comp) + etaflux(i,j,comp)
          end do
       end do
    end if

  end subroutine sum_eta_coarsest_2d

  subroutine sum_eta_coarsest_3d(lo,hi,domhi,etaflux,etasum)

    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer         , intent(in   ) :: lo(:), hi(:), domhi(:)
    real (kind=dp_t), intent(in   ) :: etaflux(lo(1):,lo(2):,lo(3):,:)
    real (kind=dp_t), intent(inout) :: etasum(0:,:)

    ! local
    integer :: i,j,k,comp

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etasum(k,rhoh_comp) = etasum(k,rhoh_comp) + etaflux(i,j,k,rhoh_comp)
             do comp=spec_comp,spec_comp+nspec-1
                etasum(k,comp) = etasum(k,comp) + etaflux(i,j,k,comp)
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
             etasum(k,rhoh_comp) = etasum(k,rhoh_comp) + etaflux(i,j,k,rhoh_comp)
             do comp=spec_comp,spec_comp+nspec-1
                etasum(k,comp) = etasum(k,comp) + etaflux(i,j,k,comp)
             end do
          end do
       end do
    end if

  end subroutine sum_eta_coarsest_3d

  subroutine compute_etapert_2d(lo,hi,etaflux,etapert,rr)

    use variables, only: rhoh_comp, spec_comp
    use network, only: nspec
    use bl_constants_module
    use geometry, only: nr

    integer         , intent(in   ) :: lo(:), hi(:), rr
    real (kind=dp_t), intent(in   ) :: etaflux(lo(1):,lo(2):,:)
    real (kind=dp_t), intent(inout) :: etapert(0:,:)
    
    ! local
    integer         :: i,j,ipert,jpert,comp
    real(kind=dp_t) :: loavg,hiavg,crseavg

    do j=lo(2),hi(2)-rr,rr
       do jpert=1,rr-1

          do i=lo(1),hi(1)-rr,rr

             loavg = ZERO
             hiavg = ZERO
             do ipert=0,rr-1
                loavg = loavg + etaflux(i+ipert,j,   rhoh_comp)
                hiavg = hiavg + etaflux(i+ipert,j+rr,rhoh_comp)
             end do
             loavg = loavg / dble(rr)
             hiavg = hiavg / dble(rr)
             crseavg = dble(jpert)/dble(rr)*loavg + dble(rr-jpert)/dble(rr)*hiavg
             do ipert=0,rr-1
                etapert(j+jpert,rhoh_comp) = etaflux(i+ipert,j+jpert,rhoh_comp) - crseavg
             end do

             do comp=spec_comp,spec_comp+nspec-1
                loavg = ZERO
                hiavg = ZERO
                do ipert=0,rr-1
                   loavg = loavg + etaflux(i+ipert,j,   comp)
                   hiavg = hiavg + etaflux(i+ipert,j+rr,comp)
                end do
                loavg = loavg / dble(rr)
                hiavg = hiavg / dble(rr)
                crseavg = dble(jpert)/dble(rr)*loavg + dble(rr-jpert)/dble(rr)*hiavg
                do ipert=0,rr-1
                   etapert(j+jpert,comp) = etaflux(i+ipert,j+jpert,comp) - crseavg
                end do
             end do
             
          end do

       end do
    end do

  end subroutine compute_etapert_2d

  subroutine compute_etapert_3d(lo,hi,etaflux,etapert,rr)

    use variables, only: rhoh_comp, spec_comp
    use network, only: nspec
    use bl_constants_module
    use geometry, only: nr

    integer         , intent(in   ) :: lo(:), hi(:), rr
    real (kind=dp_t), intent(in   ) :: etaflux(lo(1):,lo(2):,lo(3):,:)
    real (kind=dp_t), intent(inout) :: etapert(0:,:)
    
    ! local
    integer         :: i,j,k,ipert,jpert,kpert,comp
    real(kind=dp_t) :: loavg,hiavg,crseavg

    do k=lo(3),hi(3)-rr,rr
       do kpert=1,rr-1

          do j=lo(2),hi(2)-rr,rr
             do i=lo(1),hi(1)-rr,rr

                loavg = ZERO
                hiavg = ZERO
                do ipert=0,rr-1
                   do jpert=0,rr-1
                      loavg = loavg + etaflux(i+ipert,j+jpert,k,   rhoh_comp)
                      hiavg = hiavg + etaflux(i+ipert,j+jpert,k+rr,rhoh_comp)
                   end do
                end do
                loavg = loavg / dble(rr**2)
                hiavg = hiavg / dble(rr**2)
                crseavg = dble(kpert)/dble(rr)*loavg + dble(rr-kpert)/dble(rr)*hiavg
                do ipert=0,rr-1
                   do jpert=0,rr-1
                      etapert(k+kpert,rhoh_comp) = &
                           etaflux(i+ipert,j+jpert,k+kpert,rhoh_comp) - crseavg
                   end do
                end do

                do comp=spec_comp,spec_comp+nspec-1
                   loavg = ZERO
                   hiavg = ZERO
                   do ipert=0,rr-1
                      do jpert=0,rr-1
                         loavg = loavg + etaflux(i+ipert,j+jpert,k,   comp)
                         hiavg = hiavg + etaflux(i+ipert,j+jpert,k+rr,comp)
                      end do
                   end do
                   loavg = loavg / dble(rr**2)
                   hiavg = hiavg / dble(rr**2)
                   crseavg = dble(kpert)/dble(rr)*loavg + dble(rr-kpert)/dble(rr)*hiavg
                   do ipert=0,rr-1
                      do jpert=0,rr-1
                         etapert(k+kpert,comp) = &
                              etaflux(i+ipert,j+jpert,k+kpert,comp) - crseavg
                      end do
                   end do
                end do

             end do
          end do

       end do
    end do

  end subroutine compute_etapert_3d

end module make_eta_module
