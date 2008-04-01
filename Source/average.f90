! given a multifab of data (phi), average down to a base state quantity,
! phibar.  If we are in plane-parallel, the averaging is at constant
! height.  If we are spherical, then the averaging is done at constant
! radius.  

module average_module
  
  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private
  public :: average

contains

  subroutine average(mla,phi,phibar,dx,incomp,outcomp,numcomp)

    use geometry, only: nr, spherical, center, dr
    use bl_prof_module
    use bl_constants_module

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: incomp,outcomp,numcomp
    type(multifab) , intent(inout) :: phi(:)   ! Need the out so layout_aveassoc() can 
                                               ! modify the layout.
    real(kind=dp_t), intent(inout) :: phibar(:,0:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    real(kind=dp_t), pointer     :: pp(:,:,:,:)
    logical, pointer             :: mp(:,:,:,:)
    type(box)                    :: domain
    integer                      :: domlo(phi(1)%dim),domhi(phi(1)%dim)
    integer                      :: lo(phi(1)%dim),hi(phi(1)%dim)
    integer                      :: i,r,n,nlevs,ng,dm,rr,nsub,comp
    real(kind=dp_t), allocatable :: ncell_grid(:,:)
    real(kind=dp_t), allocatable :: ncell_proc(:,:)
    real(kind=dp_t), allocatable :: ncell(:,:)
    real(kind=dp_t), allocatable :: phisum_proc(:,:,:)
    real(kind=dp_t), allocatable :: phisum(:,:,:)
    real(kind=dp_t), allocatable :: phipert_proc(:,:,:)
    real(kind=dp_t), allocatable :: phipert(:,:,:)
    real(kind=dp_t), allocatable :: source_buffer(:)
    real(kind=dp_t), allocatable :: target_buffer(:)

    type(aveassoc) :: avasc

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average")

    dm = phi(1)%dim
    ng = phi(1)%ng
    nlevs = size(dx,dim=1)

    phibar = ZERO    

    if (spherical .eq. 1) allocate(ncell_grid(nlevs,0:nr(nlevs)-1))

    allocate(ncell_proc(nlevs,0:nr(nlevs)-1))
    allocate(     ncell(nlevs,0:nr(nlevs)-1))

    allocate(phisum_proc(nlevs,0:nr(nlevs)-1,numcomp))
    allocate(     phisum(nlevs,0:nr(nlevs)-1,numcomp))

    allocate(phipert_proc(nlevs,0:nr(nlevs)-1,numcomp))
    allocate(     phipert(nlevs,0:nr(nlevs)-1,numcomp))

    allocate(source_buffer(nr(nlevs)))
    allocate(target_buffer(nr(nlevs)))

    ncell        = ZERO
    ncell_proc   = ZERO
    phisum       = ZERO       
    phisum_proc  = ZERO
    phipert      = ZERO
    phipert_proc = ZERO

    if (spherical .eq. 0) then

       domain = layout_get_pd(phi(1)%la)
       domlo  = lwb(domain)
       domhi  = upb(domain)

       if (dm .eq. 2) then
          ncell(1,:) = domhi(1)-domlo(1)+1
       else if (dm .eq. 3) then
          ncell(1,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
       end if

       ! the first step is to compute phibar assuming the coarsest level 
       ! is the only level in existence
       do i=1,phi(1)%nboxes
          if ( multifab_remote(phi(1), i) ) cycle
          pp => dataptr(phi(1), i)
          lo =  lwb(get_box(phi(1), i))
          hi =  upb(get_box(phi(1), i))
          select case (dm)
          case (2)
             call sum_phi_coarsest_2d(pp(:,:,1,:),phisum_proc(1,:,:),lo,hi,ng, &
                                      incomp,numcomp)
          case (3)
             call sum_phi_coarsest_3d(pp(:,:,:,:),phisum_proc(1,:,:),lo,hi,ng, &
                                      incomp,numcomp)
          end select
       end do

       do comp=1,numcomp
          ! gather phisum
          source_buffer = phisum_proc(1,:,comp)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          phisum(1,:,comp) = target_buffer
       end do
       do comp=1,numcomp
          do r=0,nr(1)-1
             phibar(1,r,outcomp+comp-1) = phisum(1,r,comp) / dble(ncell(1,r))
          end do
       end do

       ! now we compute the phibar at the finer levels
       do n=2,nlevs

          rr = mla%mba%rr(n-1,dm)

          domain = layout_get_pd(phi(n)%la)
          domlo  = lwb(domain)
          domhi  = upb(domain)

          if (dm .eq. 2) then
             ncell(n,:) = domhi(1)-domlo(1)+1
          else if (dm .eq. 3) then
             ncell(n,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
          end if

          ! compute phisum at next finer level
          ! begin by assuming piecewise constant interpolation
          do comp=1,numcomp
             do r=0,nr(n)-1
                phisum(n,r,comp) = phisum(n-1,r/rr,comp)*rr**(dm-1)
             end do
          end do

          ! compute phipert_proc
          do i=1,phi(n)%nboxes
             if ( multifab_remote(phi(n), i) ) cycle
             pp => dataptr(phi(n), i)
             lo =  lwb(get_box(phi(n), i))
             hi =  upb(get_box(phi(n), i))
             select case (dm)
             case (2)
                call compute_phipert_2d(pp(:,:,1,:),phipert_proc(n,:,:),lo,hi,ng,incomp, &
                                        numcomp,rr)
             case (3)
                call compute_phipert_3d(pp(:,:,:,:),phipert_proc(n,:,:),lo,hi,ng,incomp, &
                                        numcomp,rr)
             end select
          end do

          do comp=1,numcomp
             ! gather phipert
             source_buffer = phipert_proc(n,:,comp)
             call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
             phipert(n,:,comp) = target_buffer
          end do
          ! update phisum and compute phibar
          do comp=1,numcomp
             do r=0,nr(n)-1
                phisum(n,r,comp) = phisum(n,r,comp) + phipert(n,r,comp)
                phibar(n,r,outcomp+comp-1) = phisum(n,r,comp) / dble(ncell(n,r))
             end do
          end do

       end do

    else if(spherical .eq. 1) then

       ! The spherical case is tricky because the base state only exists at one level
       ! as defined by dr_base in the inputs file.
       ! Therefore, the goal here is to compute phibar(nlevs,:).
       ! phisum(nlevs,:,:) will be the volume weighted sum over all levels.
       ! ncell(nlevs,:) will be the volume weighted number of cells over all levels.
       ! We make sure to use mla%mask to not double count cells, i.e.,
       ! we only sum up cells that are not covered by finer cells.
       ! we use the convention that a cell volume of 1 corresponds to dx(n=1)**3

       ! First we compute ncell(nlevs,:) and phisum(nlevs,:,:) as if the finest level
       ! were the only level in existence.
       ! Then, we add contributions from each coarser cell that is not covered by 
       ! a finer cell.

       do n=nlevs,1,-1

          ! This MUST match the nsub in average_3d_sphr().
          nsub = int(dx(n,1)/dr(nlevs)) + 1  

          avasc = layout_aveassoc(phi(n)%la,nsub,phi(n)%nodal,dx(n,:),center,dr(nlevs))

          do i=1,phi(n)%nboxes
             if ( multifab_remote(phi(n), i) ) cycle
             pp => dataptr(phi(n), i)
             lo =  lwb(get_box(phi(n), i))
             hi =  upb(get_box(phi(n), i))
             ncell_grid(n,:) = ZERO
             if (n .eq. nlevs) then
                call average_3d_sphr(n,nlevs,pp(:,:,:,:),phisum_proc(n,:,:),avasc%fbs(i), &
                                     lo,hi,ng,dx(n,:),ncell_grid(n,:),incomp,numcomp,mla)
             else
                mp => dataptr(mla%mask(n), i)
                call average_3d_sphr(n,nlevs,pp(:,:,:,:),phisum_proc(n,:,:),avasc%fbs(i), &
                                     lo,hi,ng,dx(n,:),ncell_grid(n,:),incomp,numcomp,mla, &
                                     mp(:,:,:,1))
             end if

             ncell_proc(n,:) = ncell_proc(n,:) + ncell_grid(n,:)
          end do

          call parallel_reduce(ncell(n,:), ncell_proc(n,:), MPI_SUM)

          do comp=1,numcomp
             source_buffer = phisum_proc(n,:,comp)
             call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
             phisum(n,:,comp) = target_buffer
          end do
          if (n .ne. nlevs) then
             ncell(nlevs,:) = ncell(nlevs,:) + ncell(n,:)
             do comp=1,numcomp
                do r=0,nr(nlevs)-1
                   phisum(nlevs,r,comp) = phisum(nlevs,r,comp) + phisum(n,r,comp)
                end do
             end do
          end if

       end do

       ! now divide the total phisum by the number of cells to get phibar
       do comp=1,numcomp
          do r=0,nr(nlevs)-1
             if (ncell(nlevs,r) .gt. ZERO) then
                phibar(nlevs,r,outcomp+comp-1) = phisum(nlevs,r,comp) / ncell(nlevs,r)
             else
                phibar(nlevs,r,outcomp+comp-1) = ZERO
             end if
          end do

          ! temporary hack for the case where the outermost radial bin average to zero
          ! because there is no contribution from any Cartesian cell that lies in this bin.
          ! this needs to be addressed - perhaps in the definition of nr_fine in varden.f90
          ! for spherical problems.
          if (ncell(nlevs,nr(nlevs)-1) .eq. ZERO) then
             phibar(nlevs,nr(nlevs)-1,outcomp+comp-1) = &
                  phibar(nlevs,nr(nlevs)-2,outcomp+comp-1)
          end if

       end do

       deallocate(ncell_grid)

       deallocate(ncell_proc,ncell)
       deallocate(phisum_proc,phisum)
       deallocate(phipert_proc,phipert)
       deallocate(source_buffer,target_buffer)

    endif

    call destroy(bpt)

  end subroutine average

  subroutine sum_phi_coarsest_2d(phi,phisum,lo,hi,ng,incomp,numcomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp, numcomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:,:)

    integer :: i,j,comp

    do comp=1,numcomp
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phisum(j,comp) = phisum(j,comp) + phi(i,j,incomp+comp-1)
          end do
       end do
    end do

  end subroutine sum_phi_coarsest_2d

  subroutine sum_phi_coarsest_3d(phi,phisum,lo,hi,ng,incomp,numcomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp, numcomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:,:)

    integer :: i,j,k,comp

    do comp=1,numcomp
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                phisum(k,comp) = phisum(k,comp) + phi(i,j,k,incomp+comp-1)
             end do
          end do
       end do
    end do

  end subroutine sum_phi_coarsest_3d

  subroutine compute_phipert_2d(phi,phipert,lo,hi,ng,incomp,numcomp,rr)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp, numcomp, rr
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: phipert(0:,:)

    ! Local variables
    integer          :: i, j, icrse, jcrse, comp
    real (kind=dp_t) :: crseval

    do comp=1,numcomp

       ! loop over coarse cell index
       do jcrse=lo(2)/rr,hi(2)/rr
          do icrse=lo(1)/rr,hi(1)/rr

             crseval = ZERO

             ! compute coarse cell value by taking average of fine cells
             do j=0,rr-1
                do i=0,rr-1
                   crseval = crseval + phi(icrse*rr+i,jcrse*rr+j,incomp+comp-1)
                end do
             end do
             crseval = crseval / dble(rr**2)

             ! compute phipert
             do j=0,rr-1
                do i=0,rr-1
                   phipert(jcrse*rr+j,comp) = phipert(jcrse*rr+j,comp) &
                        + phi(icrse*rr+i,jcrse*rr+j,incomp+comp-1) - crseval
                end do
             end do

          end do
       end do

    end do

  end subroutine compute_phipert_2d

  subroutine compute_phipert_3d(phi,phipert,lo,hi,ng,incomp,numcomp,rr)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp, numcomp, rr
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phipert(0:,:)

    ! Local variables
    integer          :: i, j, k, icrse, jcrse, kcrse, comp
    real (kind=dp_t) :: crseval

    do comp=1,numcomp

       ! loop over coarse cell index
       do kcrse=lo(3)/rr,hi(3)/rr
          do jcrse=lo(2)/rr,hi(2)/rr
             do icrse=lo(1)/rr,hi(1)/rr

                crseval = ZERO

                ! compute coarse cell value by taking average of fine cells
                do k=0,rr-1
                   do j=0,rr-1
                      do i=0,rr-1
                         crseval = crseval &
                              + phi(icrse*rr+i,jcrse*rr+j,kcrse*rr+k,incomp+comp-1)
                      end do
                   end do
                end do
                crseval = crseval / dble(rr**3)

                ! compute phipert
                do k=0,rr-1
                   do j=0,rr-1
                      do i=0,rr-1
                         phipert(kcrse*rr+k,comp) = phipert(kcrse*rr+k,comp) &
                              + phi(icrse*rr+i,jcrse*rr+j,kcrse*rr+k,incomp+comp-1) - crseval
                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

  end subroutine compute_phipert_3d

  subroutine average_3d_sphr(n,nlevs,phi,phisum,avfab,lo,hi,ng,dx,ncell,incomp,numcomp, &
       mla,mask)

    use geometry, only: spherical, dr, center, nr, base_cc_loc
    use ml_layout_module
    use bl_constants_module

    integer         , intent(in   ) :: n, nlevs
    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp, numcomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    type(avefab)    , intent(in   ) :: avfab
    real (kind=dp_t), intent(inout) :: phisum(0:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)
    real (kind=dp_t), intent(inout) :: ncell(0:)
    type(ml_layout) , intent(in   ) :: mla
    logical         , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    integer          :: i, j, k, l, comp, idx, cnt, nsub
    real (kind=dp_t) :: cell_weight
    logical          :: cell_valid
    !
    ! Compute nsub such that we are always guaranteed to fill each of
    ! the base state radial bins.
    !
    nsub = int(dx(1)/dr(nlevs)) + 1

    cell_weight = 1.d0 / nsub**3
    do i=2,n
       cell_weight = cell_weight / (mla%mba%rr(i-1,1))**3
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             cell_valid = .true.
             if ( present(mask) ) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             end if
             if (cell_valid) then
                do l=1,size(avfab%p(i,j,k)%v,dim=1)
                   idx = avfab%p(i,j,k)%v(l,1)
                   cnt = avfab%p(i,j,k)%v(l,2)
                   do comp=1,numcomp
                      phisum(idx,comp) = phisum(idx,comp) &
                           + cnt*cell_weight*phi(i,j,k,incomp+comp-1)
                   end do
                   ncell(idx) = ncell(idx) + cnt*cell_weight
                end do
             end if
          end do
       end do
    end do

  end subroutine average_3d_sphr

end module average_module
