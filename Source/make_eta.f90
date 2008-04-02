! Compute etarho

module make_eta_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: make_etarho

contains

  subroutine make_etarho(nlevs,etarho,etarhoflux,mla)

    use bl_constants_module
    use geometry, only: spherical, nr
    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer           , intent(in   ) :: nlevs
    real(kind=dp_t)   , intent(inout) :: etarho(:,0:)
    type(multifab)    , intent(inout) :: etarhoflux(:)
    type(ml_layout)   , intent(inout) :: mla

    ! local
    real(kind=dp_t), pointer :: efp(:,:,:,:)
    
    real(kind=dp_t), allocatable :: ncell(:,:)
    real(kind=dp_t), allocatable :: etarhosum_proc(:,:)
    real(kind=dp_t), allocatable :: etarhosum(:,:)
    real(kind=dp_t), allocatable :: etarhopert_proc(:,:)
    real(kind=dp_t), allocatable :: etarhopert(:,:)
    real(kind=dp_t), allocatable :: source_buffer(:)
    real(kind=dp_t), allocatable :: target_buffer(:)

    type(box) :: domain

    integer :: domlo(mla%dim),domhi(mla%dim)
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: i,r,rpert,n,dm,rr,comp

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_etarho")

    dm = mla%dim

    allocate(ncell       (nlevs,0:nr(nlevs))) ! ncell is a function of r only for spherical
    allocate(etarhosum_proc (nlevs,0:nr(nlevs)))
    allocate(etarhosum      (nlevs,0:nr(nlevs)))
    allocate(etarhopert_proc(nlevs,0:nr(nlevs)))
    allocate(etarhopert     (nlevs,0:nr(nlevs)))

    allocate(source_buffer(0:nr(nlevs)))
    allocate(target_buffer(0:nr(nlevs)))

    ncell        = ZERO
    etarhosum_proc  = ZERO
    etarhosum       = ZERO
    etarhopert_proc = ZERO
    etarhopert      = ZERO
    etarho          = ZERO

    if (spherical .eq. 0) then
    
       domain = layout_get_pd(mla%la(1))
       domlo = lwb(domain)
       domhi = upb(domain)

       if (dm .eq. 2) then
          ncell(1,:) = domhi(1)-domlo(1)+1
       else if(dm .eq. 3) then
          ncell(1,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
       end if
       
       ! the first step is to compute etarho assuming the coarsest level 
       ! is the only level in existence
       do i=1,layout_nboxes(mla%la(1))
          if ( multifab_remote(etarhoflux(1), i) ) cycle
          efp => dataptr(etarhoflux(1), i)
          lo =  lwb(get_box(mla%la(1), i))
          hi =  upb(get_box(mla%la(1), i))
          select case (dm)
          case (2)
            call sum_etarho_coarsest_2d(lo,hi,domhi,efp(:,:,1,1),etarhosum_proc(1,:))
          case (3)
            call sum_etarho_coarsest_3d(lo,hi,domhi,efp(:,:,:,1),etarhosum_proc(1,:))
          end select
       end do

       source_buffer = etarhosum_proc(1,:)
       call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
       etarhosum(1,:) = target_buffer

       do r=0,nr(1)
          etarho(1,r) = etarhosum(1,r) / dble(ncell(1,r))
       end do

       ! now we compute etarho at the finer levels
       do n=2,nlevs
          
          rr = mla%mba%rr(n-1,dm)

          if (mla%mba%rr(n-1,1) .ne. mla%mba%rr(n-1,dm) .or. &
               mla%mba%rr(n-1,2) .ne. mla%mba%rr(n-1,dm)) then
             print*,"ERROR: In make_etarho, refinement ratio in each direction must match"
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
          
          ! on faces that exist at the next coarser level, etarho is the same since
          ! the fluxes have been restricted.  We copy these values of etarho directly and
          ! copy scaled values of etarhosum directly.
          do r=0,nr(n-1)
             etarho   (n,r*rr) = etarho   (n-1,r)
             etarhosum(n,r*rr) = etarhosum(n-1,r)*rr**(dm-1)
          end do

          ! on faces that do not exist at the next coarser level, we use linear
          ! interpolation to get etarhosum at these faces.
          do r=0,nr(n-1)-1
             do rpert=1,rr-1
                etarhosum(n,r*rr+rpert) = &
                     dble(rpert)/dble(rr)*(etarhosum(n,r*rr)) + &
                     dble(rr-rpert)/dble(rr)*(etarhosum(n,(r+1)*rr))
             end do
          end do

          ! compute etarhopert_proc on faces that do not exist at the coarser level
          do i=1,layout_nboxes(mla%la(n))
             if ( multifab_remote(etarhoflux(n), i) ) cycle
             efp  => dataptr(etarhoflux(n), i)
             lo =  lwb(get_box(mla%la(n), i))
             hi =  upb(get_box(mla%la(n), i))
             select case (dm)
             case (2)
                call compute_etarhopert_2d(lo,hi,efp(:,:,1,1),etarhosum_proc(1,:),rr)
             case (3)
                call compute_etarhopert_3d(lo,hi,efp(:,:,:,1),etarhosum_proc(1,:),rr)
             end select
          end do

          ! gather etarhopert for rhoh
          source_buffer = etarhopert_proc(n,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          etarhopert(n,:) = target_buffer

          ! update etasum on faces that do not exist at the coarser level
          ! then recompute eta on these faces
          do r=0,nr(n-1)-1
             do rpert=1,rr-1
                etarhosum(n,r*rr+rpert) = &
                     etarhosum(n,r*rr+rpert) + etarhopert(n,r*rr+rpert)
                etarho(n,r*rr+rpert) = &
                     etarhosum(n,r*rr+rpert)/dble(ncell(n,r*rr+rpert))
             end do
          end do

       end do ! end loop over levels

    else

       ! haven't written the spherical case yet

    end if

    deallocate(ncell)
    deallocate(etarhosum_proc,etarhosum)
    deallocate(etarhopert_proc,etarhopert)
    deallocate(source_buffer,target_buffer)

    call destroy(bpt)

  end subroutine make_etarho

  subroutine sum_etarho_coarsest_2d(lo,hi,domhi,etarhoflux,etarhosum)

    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer         , intent(in   ) :: lo(:), hi(:), domhi(:)
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1):,lo(2):)
    real (kind=dp_t), intent(inout) :: etarhosum(0:)

    ! local
    integer :: i,j,comp

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          etarhosum(j) = etarhosum(j) + etarhoflux(i,j)
       end do
    end do

    ! we only add the contribution at the top edge if we are at the top of the domain
    ! this prevents double counting
    if(hi(2) .eq. domhi(2)) then
       j=hi(2)+1
       do i=lo(1),hi(1)
          etarhosum(j) = etarhosum(j) + etarhoflux(i,j)
       end do
    end if

  end subroutine sum_etarho_coarsest_2d

  subroutine sum_etarho_coarsest_3d(lo,hi,domhi,etarhoflux,etarhosum)

    use variables, only: nscal, rho_comp, rhoh_comp, spec_comp
    use network, only: nspec

    integer         , intent(in   ) :: lo(:), hi(:), domhi(:)
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1):,lo(2):,lo(3):)
    real (kind=dp_t), intent(inout) :: etarhosum(0:)

    ! local
    integer :: i,j,k,comp

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etarhosum(k) = etarhosum(k) + etarhoflux(i,j,k)
          end do
       end do
    end do

    ! we only add the contribution at the top edge if we are at the top of the domain
    ! this prevents double counting
    if(hi(3) .eq. domhi(3)) then
       k=hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etarhosum(k) = etarhosum(k) + etarhoflux(i,j,k)
          end do
       end do
    end if

  end subroutine sum_etarho_coarsest_3d

  subroutine compute_etarhopert_2d(lo,hi,etarhoflux,etarhopert,rr)

    use variables, only: rhoh_comp, spec_comp
    use network, only: nspec
    use bl_constants_module
    use geometry, only: nr

    integer         , intent(in   ) :: lo(:), hi(:), rr
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1):,lo(2):)
    real (kind=dp_t), intent(inout) :: etarhopert(0:)
    
    ! local
    integer         :: i,j,ipert,jpert,comp
    real(kind=dp_t) :: loavg,hiavg,crseavg

    do j=lo(2),hi(2)-rr,rr
       do jpert=1,rr-1

          do i=lo(1),hi(1)-rr,rr

             loavg = ZERO
             hiavg = ZERO
             do ipert=0,rr-1
                loavg = loavg + etarhoflux(i+ipert,j)
                hiavg = hiavg + etarhoflux(i+ipert,j+rr)
             end do
             loavg = loavg / dble(rr)
             hiavg = hiavg / dble(rr)
             crseavg = dble(jpert)/dble(rr)*loavg + dble(rr-jpert)/dble(rr)*hiavg
             do ipert=0,rr-1
                etarhopert(j+jpert) = etarhoflux(i+ipert,j+jpert) - crseavg
             end do
             
          end do

       end do
    end do

  end subroutine compute_etarhopert_2d

  subroutine compute_etarhopert_3d(lo,hi,etarhoflux,etarhopert,rr)

    use variables, only: rhoh_comp, spec_comp
    use network, only: nspec
    use bl_constants_module
    use geometry, only: nr

    integer         , intent(in   ) :: lo(:), hi(:), rr
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1):,lo(2):,lo(3):)
    real (kind=dp_t), intent(inout) :: etarhopert(0:)
    
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
                      loavg = loavg + etarhoflux(i+ipert,j+jpert,k)
                      hiavg = hiavg + etarhoflux(i+ipert,j+jpert,k+rr)
                   end do
                end do
                loavg = loavg / dble(rr**2)
                hiavg = hiavg / dble(rr**2)
                crseavg = dble(kpert)/dble(rr)*loavg + dble(rr-kpert)/dble(rr)*hiavg
                do ipert=0,rr-1
                   do jpert=0,rr-1
                      etarhopert(k+kpert) = &
                           etarhoflux(i+ipert,j+jpert,k+kpert) - crseavg
                   end do
                end do

             end do
          end do

       end do
    end do

  end subroutine compute_etarhopert_3d

end module make_eta_module
