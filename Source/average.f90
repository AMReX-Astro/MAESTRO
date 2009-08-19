! Given a multifab of data (phi), average down to a base state quantity, phibar.
! If we are in plane-parallel, the averaging is at constant height.  
! If we are spherical, then the averaging is done at constant radius.  

module average_module
  
  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private
  public :: average, average_one_level

contains

  subroutine average(mla,phi,phibar,dx,incomp)

    use geometry, only: nr_fine, r_start_coord, r_end_coord, spherical, numdisjointchunks, &
         dm, nlevs, dr
    use bl_prof_module
    use bl_constants_module
    use restrict_base_module
    use probin_module, only: n_cellx

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: incomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    real(kind=dp_t), pointer     :: pp(:,:,:,:)
    logical, pointer             :: mp(:,:,:,:)

    type(box)                    :: domain
    integer                      :: domlo(dm),domhi(dm),lo(dm),hi(dm)
    integer                      :: i,j,r,rcoord,n,ng,nr_crse
    integer                      :: max_radial,r_inner
    real(kind=dp_t)              :: w_lo,w_hi,del_w,wsix,theta,w_min,w_max
    real(kind=dp_t)              :: radius

    real(kind=dp_t), allocatable ::   ncell_proc(:,:)
    real(kind=dp_t), allocatable ::        ncell(:,:)
    real(kind=dp_t), allocatable ::  phisum_proc(:,:)
    real(kind=dp_t), allocatable ::       phisum(:,:)

    real(kind=dp_t), allocatable :: radii(:)

    real(kind=dp_t), allocatable :: source_buffer(:)
    real(kind=dp_t), allocatable :: target_buffer(:)

    real(kind=dp_t), allocatable :: ncell_crse(:)
    real(kind=dp_t), allocatable :: phibar_crse(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average")


    if (spherical .eq. 1) then

       max_radial = (3*(n_cellx/2-0.5d0)**2-0.75d0)/2.d0

       allocate(ncell_proc (nlevs, 0:max_radial))
       allocate(ncell      (nlevs, 0:max_radial))
       allocate(phisum_proc(nlevs, 0:max_radial))
       allocate(phisum     (nlevs,-1:max_radial))

       allocate(radii(-1:max_radial))

       allocate(source_buffer(0:max_radial))
       allocate(target_buffer(0:max_radial))

       ! radii contains every possible distance that a cell-center at the finest
       ! level can map into
       do r=0,max_radial
          radii(r) = sqrt(0.75d0+2.d0*r)*dx(nlevs,1)
       end do

       ! this refers to the center of the star
       radii(-1) = 0.d0

    else

       allocate(ncell_proc (nlevs,0:nr_fine-1))
       allocate(ncell      (nlevs,0:nr_fine-1))
       allocate(phisum_proc(nlevs,0:nr_fine-1))
       allocate(phisum     (nlevs,0:nr_fine-1))

       allocate(source_buffer(0:nr_fine-1))
       allocate(target_buffer(0:nr_fine-1))

    end if

    ng = phi(1)%ng

    phibar       = ZERO
    ncell        = ZERO
    ncell_proc   = ZERO
    phisum       = ZERO       
    phisum_proc  = ZERO

    if (spherical .eq. 0) then
       
       ! The plane-parallel case is straightforward.
       ! Simply average all the cells values at a particular height.

       do n=1,nlevs

          domain = layout_get_pd(phi(n)%la)
          domlo  = lwb(domain)
          domhi  = upb(domain)

          if (dm .eq. 1) then
             ncell(n,:) = 1
          else if (dm .eq. 2) then
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
             case (1)
                call sum_phi_1d(pp(:,1,1,:),phisum_proc(n,:),lo,hi,ng,incomp)
             case (2)
                call sum_phi_2d(pp(:,:,1,:),phisum_proc(n,:),lo,hi,ng,incomp)
             case (3)
                call sum_phi_3d(pp(:,:,:,:),phisum_proc(n,:),lo,hi,ng,incomp)
             end select
          end do

          ! gather phisum
          source_buffer = phisum_proc(n,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          phisum(n,:) = target_buffer

          ! compute phibar by normalizing phisum
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)
                phibar(n,r) = phisum(n,r) / dble(ncell(n,r))
             end do
          end do

       end do

       call restrict_base(phibar,.true.)
       call fill_ghost_base(phibar,.true.)

    else if(spherical .eq. 1) then

       ! For spherical, we construct a 1D array, phisum, that has space
       ! allocated for every possible radius that a cell-center at the finest
       ! level can map into.  The radius locations have been precomputed and stored
       ! in radii(:).

       ! For cells at the non-finest level, map a weighted contribution into the nearest
       ! bin in phisum.

       do n=nlevs,1,-1

          do i=1,phi(n)%nboxes
             if ( multifab_remote(phi(n), i) ) cycle
             pp => dataptr(phi(n), i)
             lo =  lwb(get_box(phi(n), i))
             hi =  upb(get_box(phi(n), i))

             if (n .eq. nlevs) then
                call sum_phi_3d_sphr(n,radii(0:),max_radial,pp(:,:,:,:),phisum_proc(n,:), &
                                     lo,hi,ng,dx(n,:),ncell_proc(n,:),incomp)
             else
                ! we include the mask so we don't double count; i.e., we only consider
                ! cells that we can "see" when constructing the sum
                mp => dataptr(mla%mask(n), i)
                call sum_phi_3d_sphr(n,radii(0:),max_radial,pp(:,:,:,:),phisum_proc(n,:), &
                                     lo,hi,ng,dx(n,:),ncell_proc(n,:),incomp, &
                                     mp(:,:,:,1))
             end if
          end do

          source_buffer = ncell_proc(n,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          ncell(n,:) = target_buffer

          source_buffer = phisum_proc(n,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          phisum(n,0:) = target_buffer

       end do

       ! now gather ncell and phisum from all the levels and store
       ! them in ncell(1,:) and phisum(1,:).  We use 0 based indexing for
       ! phisum so we don't mess with the center point of the star.  We will
       ! compute phisum(1,-1) later.
       do n=2,nlevs
          ncell (1, :) = ncell (1, :) + ncell (n, :)
          phisum(1,0:) = phisum(1,0:) + phisum(n,0:)
       end do

       ! normalize the contributions
       do r=0,max_radial
         if (ncell(1,r) .ne. 0.d0) then
             phisum(1,r) = phisum(1,r) / ncell(1,r)
          end if
       end do

       ! fill the empty bins using piecewise linear functions
       do r=0,max_radial

          if (ncell(1,r) .eq. 0.d0) then
             ! compute the upper bounding point, rcoord.
             rcoord = max_radial+1
             do j=r+1,max_radial
                if (ncell(1,j) .ne. 0) then
                   rcoord = j
                   exit
                end if
             end do

             if (rcoord .ne. max_radial+1) then
                call lin_interp(radii(r), &
                                radii(r-1),radii(rcoord), &
                                phisum(1,r), &
                                phisum(1,r-1),phisum(1,rcoord))
             else
                ! if there is no upper bounding point, linearly extrapolate
                call lin_interp(radii(r), &
                                radii(r-2),radii(r-1), &
                                phisum(1,r), &
                                phisum(1,r-2),phisum(1,r-1))
             end if
          end if

       end do

       ! use quadratic interpolation with homogeneous neumann bc at center of star
       ! to compute value at center of star, indicated with coordinate r=-1
       phisum(1,-1) = (11.d0/8.d0)*phisum(1,0) - (3.d0/8.d0)*phisum(1,1)

       ! compute phibar
       do r=0,nr_fine-1

          ! compute the radius in physical coordinates
          radius = (dble(r)+HALF)*dr(1)

          ! compute the coordinate, rcoord, such that
          ! radius lies in between the radii(r) and radii(r+1)
          rcoord = ((radius/dx(nlevs,1))**2 - 0.75d0)/2.0d0

          ! now overwrite rcoord to correspond to the lo stencil point
          ! for the quadratic interpolation.
          ! compare rcoord-1 and rcoord+2 and see which is closer.
          ! if rcoord-1 is closer, set rcoord=roord-1
          if (rcoord+2 .le. max_radial) then
             if (abs(radii(rcoord-1)-radius) .lt. abs(radii(rcoord+2)-radius)) then
                rcoord = rcoord-1
             end if
          end if

          ! use first order extrapolation if the stencil would have used a point
          ! outside of max_radial
          if (rcoord+2 .gt. max_radial) then
             ! extrapolate the quadratic from the 3 outermost points
             call quad_interp(radius, &
                              radii(max_radial-2),radii(max_radial-1), &
                              radii(max_radial), &
                              phibar(1,r), &
                              phisum(1,max_radial-2),phisum(1,max_radial-1), &
                              phisum(1,max_radial))
          else
             ! interpolate from the three nearest points
             call quad_interp(radius, &
                              radii(rcoord),radii(rcoord+1),radii(rcoord+2), &
                              phibar(1,r), &
                              phisum(1,rcoord),phisum(1,rcoord+1),phisum(1,rcoord+2))
          end if

       end do

    endif

    call destroy(bpt)

  end subroutine average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_1d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i

    do i=lo(1),hi(1)
       phisum(i) = phi(i,incomp)
    end do

  end subroutine sum_phi_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_2d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i,j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          phisum(j) = phisum(j) + phi(i,j,incomp)
       end do
    end do

  end subroutine sum_phi_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_3d(phi,phisum,lo,hi,ng,incomp)

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

  end subroutine sum_phi_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_3d_sphr(n,radii,max_radial,phi,phisum,lo,hi,ng,dx,ncell,incomp,mask)

    use geometry, only: dr, center, nlevs, dm
    use ml_layout_module
    use bl_constants_module

    integer         , intent(in   )           :: n, lo(:), hi(:), ng, incomp, max_radial
    real (kind=dp_t), intent(in   )           :: radii(0:)
    real (kind=dp_t), intent(in   )           :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout)           :: phisum(0:)
    real (kind=dp_t), intent(in   )           :: dx(:)
    real (kind=dp_t), intent(inout)           :: ncell(0:)
    logical         , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    ! local
    real (kind=dp_t) :: x, y, z, radius, dx_fine(dm), weight, diff
    integer          :: i, j, k, ii, jj, kk, index, factor
    logical          :: cell_valid

    if (n .eq. nlevs) then
       
       do k=lo(3),hi(3)
          z = (dble(k) + HALF)*dx(3) - center(3)

          do j=lo(2),hi(2)
             y = (dble(j) + HALF)*dx(2) - center(2)

             do i=lo(1),hi(1)
                x = (dble(i) + HALF)*dx(1) - center(1)

                ! compute distance to center
                radius = sqrt(x**2 + y**2 + z**2)

                ! figure out which radii index this point maps into
                index = ((radius / dx(1))**2 - 0.75d0) / 2.d0

                ! due to roundoff error, need to ensure that we are in the proper radial bin
                if (index .lt. max_radial) then
                   if (abs(radius-radii(index)) .gt. abs(radius-radii(index+1))) then
                      index = index+1
                   end if
                end if

                ! update phisum and ncell
                phisum(index) = phisum(index) + phi(i,j,k,incomp)
                ncell(index)  = ncell(index) + 1.d0

             end do
          end do
       end do

    else

       factor = 2**(nlevs-n)
       weight = 8.d0**(nlevs-n)

       do i=1,dm
          dx_fine(i) = dx(i)/dble(factor)
       end do

       do k=lo(3),hi(3)
          z = (dble(k) + HALF)*dx(3) - center(3)

          do j=lo(2),hi(2)
             y = (dble(j) + HALF)*dx(2) - center(2)

             do i=lo(1),hi(1)
                x = (dble(i) + HALF)*dx(1) - center(1)

                ! make sure the cell isn't covered by finer cells
                cell_valid = .true.
                if ( present(mask) ) then
                   if ( (.not. mask(i,j,k)) ) cell_valid = .false.
                end if

                if (cell_valid) then

                   ! compute distance to center
                   radius = sqrt(x**2 + y**2 + z**2)

                   ! figure out which radii index this point maps into
                   index = ((radius / dx_fine(1))**2 - 0.75d0) / 2.d0

                   ! we won't map exactly onto a location in radii, 
                   ! so we use the index of the closest point
                   if (index .lt. max_radial) then
                      if (abs(radius-radii(index)) .gt. abs(radius-radii(index+1))) then
                         index = index+1
                      end if
                   end if

                   ! update phisum and ncell
                   phisum(index) = phisum(index) + weight*phi(i,j,k,incomp)
                   ncell(index)  = ncell(index) + weight

                end if

             end do
          end do
       end do

    end if

  end subroutine sum_phi_3d_sphr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_one_level(n,phi,phibar,incomp)

    use geometry, only: nr_fine, nr, spherical, dm
    use bl_prof_module
    use bl_constants_module

    integer        , intent(in   ) :: n,incomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar(:,0:)

    ! local
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    type(box)                :: domain
    integer                  :: domlo(dm),domhi(dm)
    integer                  :: lo(dm),hi(dm)
    integer                  :: i,r,ng,ncell

    real(kind=dp_t) ::   phisum_proc(0:nr_fine-1)
    real(kind=dp_t) ::        phisum(0:nr_fine-1)
    real(kind=dp_t) :: source_buffer(0:nr_fine-1)
    real(kind=dp_t) :: target_buffer(0:nr_fine-1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average_one_level")

    ng = phi(1)%ng

    phibar = ZERO

    ncell        = ZERO
    phisum       = ZERO       
    phisum_proc  = ZERO

    if (spherical .eq. 0) then
       
       domain = layout_get_pd(phi(n)%la)
       domlo  = lwb(domain)
       domhi  = upb(domain)

       if (dm .eq. 1) then
          ncell = 1
       else if (dm .eq. 2) then
          ncell = domhi(1)-domlo(1)+1
       else if (dm .eq. 3) then
          ncell = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
       end if

       do i=1,phi(n)%nboxes
          if ( multifab_remote(phi(n), i) ) cycle
          pp => dataptr(phi(n), i)
          lo =  lwb(get_box(phi(n), i))
          hi =  upb(get_box(phi(n), i))
          select case (dm)
          case (1)
             call sum_phi_1d(pp(:,1,1,:),phisum_proc,lo,hi,ng,incomp)
          case (2)
             call sum_phi_2d(pp(:,:,1,:),phisum_proc,lo,hi,ng,incomp)
          case (3)
             call sum_phi_3d(pp(:,:,:,:),phisum_proc,lo,hi,ng,incomp)
          end select
       end do

       ! gather phisum
       source_buffer = phisum_proc
       call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
       phisum = target_buffer
       do r=0,nr(n)-1
          phibar(n,r) = phisum(r) / dble(ncell)
       end do

    else if(spherical .eq. 1) then

       call bl_error("average_one_level not written for multilevel spherical")

    end if

    call destroy(bpt)

  end subroutine average_one_level

  subroutine quad_interp(x,x0,x1,x2,y,y0,y1,y2)

    real(kind=dp_t), intent(in   ) :: x,x0,x1,x2,y0,y1,y2
    real(kind=dp_t), intent(  out) :: y
    
    y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
           + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1)

    if (y .gt. max(y0,y1,y2)) y = max(y0,y1,y2)
    if (y .lt. min(y0,y1,y2)) y = min(y0,y1,y2)

  end subroutine quad_interp

  subroutine lin_interp(x,x0,x1,y,y0,y1)

    real(kind=dp_t), intent(in   ) :: x,x0,x1,y0,y1
    real(kind=dp_t), intent(  out) :: y
    
    y = y0 + (y1-y0)/(x1-x0)*(x-x0)

  end subroutine lin_interp

end module average_module
