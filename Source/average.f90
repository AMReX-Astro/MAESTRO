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

  subroutine average(mla,phi,phibar,dx,incomp)

    use geometry, only: nr_fine, r_start_coord, r_end_coord, spherical, numdisjointchunks, &
         dm
    use bl_prof_module
    use bl_constants_module
    use restrict_base_module
    use probin_module, only: nlevs, drdxfac

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: incomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    real(kind=dp_t), pointer     :: pp(:,:,:,:)
    logical, pointer             :: mp(:,:,:,:)
    type(box)                    :: domain
    integer                      :: domlo(dm),domhi(dm)
    integer                      :: lo(dm),hi(dm)
    integer                      :: i,r,n,ng,rr
    real(kind=dp_t), allocatable :: ncell_grid(:,:)
    real(kind=dp_t), allocatable :: ncell_proc(:,:)
    real(kind=dp_t), allocatable :: ncell(:,:)
    real(kind=dp_t), allocatable :: phisum_proc(:,:)
    real(kind=dp_t), allocatable :: phisum(:,:)
    real(kind=dp_t), allocatable :: phipert_proc(:,:)
    real(kind=dp_t), allocatable :: phipert(:,:)
    real(kind=dp_t), allocatable :: source_buffer(:)
    real(kind=dp_t), allocatable :: target_buffer(:)
    logical                      :: fine_grids_span_domain_width

    integer                      :: j
    integer                      :: nr_crse
    real(kind=dp_t), allocatable :: ncell_crse(:)
    real(kind=dp_t), allocatable :: phibar_crse(:)
    real(kind=dp_t)              :: w_lo, w_hi, del_w, wsix, theta

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average")

    ng = phi(1)%ng

    phibar = ZERO

    if (spherical .eq. 1) then
       allocate(ncell_grid(nlevs,0:nr_fine-1))
    end if

    allocate(ncell_proc   (nlevs,0:nr_fine-1))
    allocate(ncell        (nlevs,0:nr_fine-1))
    allocate(phisum_proc  (nlevs,0:nr_fine-1))
    allocate(phisum       (nlevs,0:nr_fine-1))
    allocate(phipert_proc (nlevs,0:nr_fine-1))
    allocate(phipert      (nlevs,0:nr_fine-1))
    allocate(source_buffer      (0:nr_fine-1))
    allocate(target_buffer      (0:nr_fine-1))

    ncell        = ZERO
    ncell_proc   = ZERO
    phisum       = ZERO       
    phisum_proc  = ZERO
    phipert      = ZERO
    phipert_proc = ZERO

    if (spherical .eq. 0) then

       fine_grids_span_domain_width = .true.

       if (fine_grids_span_domain_width) then

          do n=1,nlevs

             domain = layout_get_pd(phi(n)%la)
             domlo  = lwb(domain)
             domhi  = upb(domain)
             
             if (dm .eq. 2) then
                ncell(n,:) = domhi(1)-domlo(1)+1
             else if (dm .eq. 3) then
                ncell(n,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
             end if
             
             do i=1,phi(n)%nboxes
                if ( multifab_remote(phi(n), i) ) cycle
                pp => dataptr(phi(n), i)
                lo =  lwb(get_box(phi(n), i))
                hi =  upb(get_box(phi(n), i))
                select case (dm)
                case (2)
                   call sum_phi_coarsest_2d(pp(:,:,1,:),phisum_proc(n,:),lo,hi,ng,incomp)
                case (3)
                   call sum_phi_coarsest_3d(pp(:,:,:,:),phisum_proc(n,:),lo,hi,ng,incomp)
                end select
             end do
             
             ! gather phisum
             source_buffer = phisum_proc(n,:)
             call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
             phisum(n,:) = target_buffer
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n,i),r_end_coord(n,i)
                   phibar(n,r) = phisum(n,r) / dble(ncell(n,r))
                end do
             end do
             
          end do

          if (nlevs .ge. 2) then
             call fill_ghost_base(phibar,.true.)
             call restrict_base(phibar,.true.)
          end if

       else

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
                call sum_phi_coarsest_2d(pp(:,:,1,:),phisum_proc(1,:),lo,hi,ng,incomp)
             case (3)
                call sum_phi_coarsest_3d(pp(:,:,:,:),phisum_proc(1,:),lo,hi,ng,incomp)
             end select
          end do
          
          ! gather phisum
          source_buffer = phisum_proc(1,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          phisum(1,:) = target_buffer
          do r=0,r_end_coord(1,1)
             phibar(1,r) = phisum(1,r) / dble(ncell(1,r))
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
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n,i),r_end_coord(n,i)
                   phisum(n,r) = phisum(n-1,r/rr)*rr**(dm-1)
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
                   call compute_phipert_2d(pp(:,:,1,:),phipert_proc(n,:),lo,hi,ng,incomp,rr)
                case (3)
                   call compute_phipert_3d(pp(:,:,:,:),phipert_proc(n,:),lo,hi,ng,incomp,rr)
                end select
             end do
             
             ! gather phipert
             source_buffer = phipert_proc(n,:)
             call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
             phipert(n,:) = target_buffer
             
             ! update phisum and compute phibar
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n,i),r_end_coord(n,i)
                   phisum(n,r) = phisum(n,r) + phipert(n,r)
                   phibar(n,r) = phisum(n,r) / dble(ncell(n,r))
                end do
             end do

          end do

       end if

    else if(spherical .eq. 1) then

       ! The spherical case is tricky because the base state only exists at 
       ! one level with dr = dx / drdxfac.

       ! Therefore, the goal here is to compute phibar(nlevs,:).
       ! phisum(nlevs,:,:) will be the volume weighted sum over all levels.
       ! ncell(nlevs,:) will be the volume weighted number of cells over 
       ! all levels.

       ! We make sure to use mla%mask to not double count cells, i.e.,
       ! we only sum up cells that are not covered by finer cells.
       ! we use the convention that a cell volume of 1 corresponds to 
       ! dx(n=1)**3

       ! First we compute ncell(nlevs,:) and phisum(nlevs,:,:) as if the 
       ! finest level were the only level in existence.
       ! Then, we add contributions from each coarser cell that is not 
       ! covered by a finer cell.

       do n=nlevs,1,-1

          do i=1,phi(n)%nboxes
             if ( multifab_remote(phi(n), i) ) cycle
             pp => dataptr(phi(n), i)
             lo =  lwb(get_box(phi(n), i))
             hi =  upb(get_box(phi(n), i))
             ncell_grid(n,:) = ZERO
             if (n .eq. nlevs) then
                call average_3d_sphr(n,nlevs,pp(:,:,:,:),phisum_proc(n,:), &
                                     lo,hi,ng,dx(n,:),ncell_grid(n,:),incomp,mla)
             else
                mp => dataptr(mla%mask(n), i)
                call average_3d_sphr(n,nlevs,pp(:,:,:,:),phisum_proc(n,:), &
                                     lo,hi,ng,dx(n,:),ncell_grid(n,:),incomp,mla, &
                                     mp(:,:,:,1))
             end if

             ncell_proc(n,:) = ncell_proc(n,:) + ncell_grid(n,:)
          end do

          call parallel_reduce(ncell(n,:), ncell_proc(n,:), MPI_SUM)

          source_buffer = phisum_proc(n,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          phisum(n,:) = target_buffer

          if (n .ne. nlevs) then
             ncell(nlevs,:) = ncell(nlevs,:) + ncell(n,:)
             do r=0,nr_fine-1
                phisum(nlevs,r) = phisum(nlevs,r) + phisum(n,r)
             end do
          end if

       end do

       if (drdxfac .ne. 1) then

          nr_crse = nr_fine / drdxfac + 1
          allocate( ncell_crse(0:nr_crse-1))
          allocate(phibar_crse(-2:nr_crse+1))

          do r = 0, nr_crse-1

             phibar_crse(r) = 0.d0
              ncell_crse(r) = 0.d0
   
             ! Sum fine data onto the crse grid
             do j = drdxfac*r, min(drdxfac*r+(drdxfac-1),nr_fine-1)
                phibar_crse(r) = phibar_crse(r) + phisum(nlevs,j)
                 ncell_crse(r) =  ncell_crse(r) +  ncell(nlevs,j)
             end do

             ! Now compute the average
             if (ncell_crse(r) .gt. ZERO) then
                phibar_crse(r) = phibar_crse(r) / ncell_crse(r)
             else if (r .eq. nr_crse-1) then
                phibar_crse(r) = phibar_crse(nr_crse-2)
             else 
                phibar_crse(r) = ZERO
             end if
          end do

          ! Reflect (even) across origin
          phibar_crse(-1) = phibar_crse(0)
          phibar_crse(-2) = phibar_crse(1)

          ! Extend at high r
          phibar_crse(nr_crse  ) = phibar_crse(nr_crse-1)
          phibar_crse(nr_crse+1) = phibar_crse(nr_crse-1)

          ! Put the average back onto the fine grid
          do r = 0, nr_crse-1
   
             w_lo = ( 7.d0 * (phibar_crse(r  ) + phibar_crse(r-1)) &
                     -1.d0 * (phibar_crse(r+1) + phibar_crse(r-2)) ) / 12.d0
             w_hi = ( 7.d0 * (phibar_crse(r  ) + phibar_crse(r+1)) &
                     -1.d0 * (phibar_crse(r-1) + phibar_crse(r+2)) ) / 12.d0

             del_w = w_hi - w_lo

             wsix = 6.d0 * ( phibar_crse(r) - 0.5d0 * (w_lo + w_hi) )

             do j = 0, min(drdxfac-1,nr_fine-drdxfac*r-1)
                ! piecewise constant
!               phibar(nlevs,drdxfac*r+j) = phibar_crse(r)
   
                ! parabolic interpolation
                theta = dble(j) / dble(drdxfac)
                phibar(nlevs,drdxfac*r+j) = w_lo + theta * del_w + &
                                       theta * (1.d0 - theta) * wsix
             end do

          end do

       else 

          ! if drdxfac = 1 then divide the total phisum by the number of cells to get phibar
          do r=0,nr_fine-1
             if (ncell(nlevs,r) .gt. ZERO) then
                phibar(nlevs,r) = phisum(nlevs,r) / ncell(nlevs,r)
             else
                phibar(nlevs,r) = ZERO
             end if
          end do

       end if
       
       ! temporary hack for the case where the outermost radial bin average 
       ! to zero because there is no contribution from any Cartesian cell 
       ! that lies in this bin.  This needs to be addressed - perhaps in the 
       ! definition of nr_fine in varden.f90 for spherical problems.
       if (ncell(nlevs,nr_fine-1) .eq. ZERO) then
          phibar(nlevs,nr_fine-1) = phibar(nlevs,nr_fine-2)
       end if

       deallocate(ncell_grid)

       deallocate(ncell_proc,ncell)
       deallocate(phisum_proc,phisum)
       deallocate(phipert_proc,phipert)
       deallocate(source_buffer,target_buffer)

       if (drdxfac .ne. 1) then
          deallocate( ncell_crse)
          deallocate(phibar_crse)
       end if

    endif

    call destroy(bpt)

  end subroutine average

  subroutine sum_phi_coarsest_2d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i,j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          phisum(j) = phisum(j) + phi(i,j,incomp)
       end do
    end do

  end subroutine sum_phi_coarsest_2d

  subroutine sum_phi_coarsest_3d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i,j,k

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phisum(k) = phisum(k) + phi(i,j,k,incomp)
          end do
       end do
    end do

  end subroutine sum_phi_coarsest_3d

  subroutine compute_phipert_2d(phi,phipert,lo,hi,ng,incomp,rr)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp, rr
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: phipert(0:)

    ! Local variables
    integer          :: i, j, icrse, jcrse
    real (kind=dp_t) :: crseval

    ! loop over coarse cell index
    do jcrse=lo(2)/rr,hi(2)/rr
       do icrse=lo(1)/rr,hi(1)/rr
          
          crseval = ZERO
          
          ! compute coarse cell value by taking average of fine cells
          do j=0,rr-1
             do i=0,rr-1
                crseval = crseval + phi(icrse*rr+i,jcrse*rr+j,incomp)
             end do
          end do
          crseval = crseval / dble(rr**2)
          
          ! compute phipert
          do j=0,rr-1
             do i=0,rr-1
                phipert(jcrse*rr+j) = phipert(jcrse*rr+j) &
                     + phi(icrse*rr+i,jcrse*rr+j,incomp) - crseval
             end do
          end do
          
       end do
    end do

  end subroutine compute_phipert_2d

  subroutine compute_phipert_3d(phi,phipert,lo,hi,ng,incomp,rr)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp, rr
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phipert(0:)

    ! Local variables
    integer          :: i, j, k, icrse, jcrse, kcrse
    real (kind=dp_t) :: crseval

    ! loop over coarse cell index
    do kcrse=lo(3)/rr,hi(3)/rr
       do jcrse=lo(2)/rr,hi(2)/rr
          do icrse=lo(1)/rr,hi(1)/rr
             
             crseval = ZERO
             
             ! compute coarse cell value by taking average of fine cells
             do k=0,rr-1
                do j=0,rr-1
                   do i=0,rr-1
                      crseval = crseval + phi(icrse*rr+i,jcrse*rr+j,kcrse*rr+k,incomp)
                   end do
                end do
             end do
             crseval = crseval / dble(rr**3)
             
             ! compute phipert
             do k=0,rr-1
                do j=0,rr-1
                   do i=0,rr-1
                      phipert(kcrse*rr+k) = phipert(kcrse*rr+k) &
                           + phi(icrse*rr+i,jcrse*rr+j,kcrse*rr+k,incomp) - crseval
                   end do
                end do
             end do
             
          end do
       end do
    end do
    
  end subroutine compute_phipert_3d

  subroutine average_3d_sphr(n,nlevs,phi,phisum,lo,hi,ng,dx,ncell,incomp,mla,mask)

    use geometry, only: spherical, dr, center
    use ml_layout_module
    use bl_constants_module

    integer         , intent(in   ) :: n, nlevs
    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)
    real (kind=dp_t), intent(inout) :: ncell(0:)
    type(ml_layout) , intent(in   ) :: mla
    logical         , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)
    real (kind=dp_t) :: x, y, z, radius
    integer          :: i, j, k, l, idx, cnt
    real (kind=dp_t) :: cell_weight
    logical          :: cell_valid

    cell_weight = ONE
    do i=2,n
       cell_weight = cell_weight / (mla%mba%rr(i-1,1))**3
    end do

    do k=lo(3),hi(3)
       z = (dble(k) + HALF)*dx(3) - center(3)

       do j=lo(2),hi(2)
          y = (dble(j) + HALF)*dx(2) - center(2)

          do i=lo(1),hi(1)
             x = (dble(i) + HALF)*dx(1) - center(1)

             cell_valid = .true.
             if ( present(mask) ) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             end if

             if (cell_valid) then
                radius = sqrt(x**2 + y**2 + z**2)
                idx  = int(radius / dr(n))
                
                phisum(idx) = phisum(idx) + cell_weight*phi(i,j,k,incomp)
                ncell(idx) = ncell(idx) + cell_weight
             end if

          end do
       end do
    end do

  end subroutine average_3d_sphr

end module average_module
