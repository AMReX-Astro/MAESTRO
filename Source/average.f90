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

    use geometry, only: nr_fine, nr_irreg, r_start_coord, r_end_coord, spherical, &
         numdisjointchunks, dm, nlevs, dr
    use bl_prof_module
    use bl_constants_module
    use restrict_base_module
    use probin_module, only: prob_hi

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
    integer                      :: max_rcoord(nlevs),rcoord(nlevs),stencil_coord(nlevs)
    integer                      :: i,j,r,n,ng
    real(kind=dp_t)              :: radius

    integer, allocatable ::  ncell_proc(:,:)
    integer, allocatable ::       ncell(:,:)

    integer, allocatable :: which_lev(:)

    real(kind=dp_t), allocatable :: phisum_proc(:,:)
    real(kind=dp_t), allocatable ::      phisum(:,:)

    real(kind=dp_t), allocatable :: radii(:,:)

    real(kind=dp_t), allocatable :: source_buffer(:)
    real(kind=dp_t), allocatable :: target_buffer(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average")

    if (spherical .eq. 1) then

       allocate(ncell_proc (nlevs, 0:nr_irreg))
       allocate(ncell      (nlevs,-1:nr_irreg))

       allocate(which_lev(0:nr_fine-1))

       allocate(phisum_proc(nlevs, 0:nr_irreg))
       allocate(phisum     (nlevs,-1:nr_irreg))
       allocate(radii      (nlevs,-1:nr_irreg+1))

       allocate(source_buffer(0:nr_irreg))
       allocate(target_buffer(0:nr_irreg))

       ! radii contains every possible distance that a cell-center at the finest
       ! level can map into
       do n=1,nlevs
          do r=0,nr_irreg
             radii(n,r) = sqrt(0.75d0+2.d0*r)*dx(n,1)
          end do
       end do
       radii(:,nr_irreg+1) = 1.d99

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
    ncell        = 0
    ncell_proc   = 0
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
                call sum_phi_3d_sphr(radii(n,0:),nr_irreg,pp(:,:,:,:),phisum_proc(n,:), &
                                     lo,hi,ng,dx(n,:),ncell_proc(n,:),incomp)
             else
                ! we include the mask so we don't double count; i.e., we only consider
                ! cells that we can "see" when constructing the sum
                mp => dataptr(mla%mask(n), i)
                call sum_phi_3d_sphr(radii(n,0:),nr_irreg,pp(:,:,:,:),phisum_proc(n,:), &
                                     lo,hi,ng,dx(n,:),ncell_proc(n,:),incomp, &
                                     mp(:,:,:,1))
             end if
          end do

          source_buffer = ncell_proc(n,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          ncell(n,0:) = target_buffer

          source_buffer = phisum_proc(n,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          phisum(n,0:) = target_buffer

       end do

       ! normalize phisum so it actually stores the average at a radius
       do n=1,nlevs
          do r=0,nr_irreg
             if (ncell(n,r) .ne. 0.d0) then
                phisum(n,r) = phisum(n,r) / dble(ncell(n,r))
             end if
          end do
       end do

       ! compute center point for the finest level
       phisum(nlevs,-1) = (11.d0/8.d0)*phisum(nlevs,0) - (3.d0/8.d0)*phisum(nlevs,1)
       radii (nlevs,-1) = 0.d0
       ncell (nlevs,-1) = 1

       ! choose which level to interpolate from
       do n=1,nlevs
          rcoord(n) = 0
       end do

       do r=0,nr_fine-1

         radius = (dble(r)+HALF)*dr(1)

         ! for each level, find the closest coordinate
         do n=1,nlevs
            do j=rcoord(n),nr_irreg
               if (abs(radius-radii(n,j)) .lt. abs(radius-radii(n,j+1))) then
                  rcoord(n) = j
                  exit
               end if
            end do
         end do

         ! make sure closest coordinate is in bounds
         do n=1,nlevs-1
            rcoord(n) = max(rcoord(n),1)
         end do
         do n=1,nlevs
            rcoord(n) = min(rcoord(n),nr_irreg-1)
         end do
         
         ! choose the level with the largest min over the ncell interpolation points
         which_lev(r)=1
         do n=2,nlevs
            if (min(ncell(n,rcoord(n)-1), &
                    ncell(n,rcoord(n)), &
                    ncell(n,rcoord(n)+1)) &
                    .gt. &
                min(ncell(which_lev(r),rcoord(which_lev(r))-1),&
                    ncell(which_lev(r),rcoord(which_lev(r))), &
                    ncell(which_lev(r),rcoord(which_lev(r))+1))) then
               which_lev(r) = n
            end if
         end do

      end do

       ! squish the list at each level down to exclude points with no contribution
       do n=1,nlevs
          j=0
          do r=0,nr_irreg
             do while(ncell(n,j) .eq. 0)
                j = j+1
                if (j .gt. nr_irreg) then
                   exit
                end if
             end do
             if (j .gt. nr_irreg) then
                phisum(n,r:nr_irreg)   = 1.d99
                radii (n,r:nr_irreg+1) = 1.d99
                max_rcoord(n) = r-1
                exit
             end if
             phisum(n,r) = phisum(n,j)
             radii (n,r) = radii (n,j)
             ncell (n,r) = ncell (n,j)
             j = j+1
             if (j .gt. nr_irreg) then
                max_rcoord(n) = r
                exit
             end if
          end do
       end do

       ! compute phibar
       stencil_coord(:) = 0

       do r=0,nr_fine-1

         radius = (dble(r)+HALF)*dr(1)

         ! find the closest coordinate
         do j=stencil_coord(which_lev(r)),max_rcoord(which_lev(r))
            if (abs(radius-radii(which_lev(r),j  )) .lt. &
                abs(radius-radii(which_lev(r),j+1))) then
               stencil_coord(which_lev(r)) = j
               exit
            end if
         end do

         ! make sure the interpolation points will be in bounds
         if (which_lev(r) .ne. nlevs) then
            stencil_coord(which_lev(r)) = max(stencil_coord(which_lev(r)),1)
         end if
         stencil_coord(which_lev(r)) = min(stencil_coord(which_lev(r)),max_rcoord(which_lev(r))-1)

         call quad_interp(radius, &
                          radii(which_lev(r),stencil_coord(which_lev(r))-1), &
                          radii(which_lev(r),stencil_coord(which_lev(r))  ), &
                          radii(which_lev(r),stencil_coord(which_lev(r))+1), &
                          phibar(1,r), &
                          phisum(which_lev(r),stencil_coord(which_lev(r))-1), &
                          phisum(which_lev(r),stencil_coord(which_lev(r))  ), &
                          phisum(which_lev(r),stencil_coord(which_lev(r))+1))

      end do

   end if

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

  subroutine sum_phi_3d_sphr(radii,nr_irreg,phi,phisum,lo,hi,ng,dx,ncell,incomp,mask)

    use geometry, only: dr, center, nlevs, dm
    use ml_layout_module
    use bl_constants_module

    integer         , intent(in   )           :: lo(:), hi(:), ng, incomp, nr_irreg
    real (kind=dp_t), intent(in   )           :: radii(0:)
    real (kind=dp_t), intent(in   )           :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout)           :: phisum(0:)
    real (kind=dp_t), intent(in   )           :: dx(:)
    integer         , intent(inout)           :: ncell(0:)
    logical         , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    ! local
    real (kind=dp_t) :: x, y, z, radius
    integer          :: i, j, k, index
    logical          :: cell_valid

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
                index = ((radius / dx(1))**2 - 0.75d0) / 2.d0
                
                ! due to roundoff error, need to ensure that we are in the proper radial bin
                if (index .lt. nr_irreg) then
                   if (abs(radius-radii(index)) .gt. abs(radius-radii(index+1))) then
                      index = index+1
                   end if
                end if
                
                ! update phisum and ncell
                phisum(index) = phisum(index) + phi(i,j,k,incomp)
                ncell(index)  = ncell(index) + 1

             end if
             
          end do
       end do
    end do

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

    ncell        = 0
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

  subroutine cubic_interp(x,x0,x1,x2,x3,y,y0,y1,y2,y3)

    real(kind=dp_t), intent(in   ) :: x,x0,x1,x2,x3,y0,y1,y2,y3
    real(kind=dp_t), intent(  out) :: y
    
    y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
           + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1) &
           + ( ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0) &
              -((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1))/(x3-x1) ) / (x3-x0) &
            *(x-x0)*(x-x1)*(x-x2)

    if (y .gt. max(y0,y1,y2,y3)) y = max(y0,y1,y2,y3)
    if (y .lt. min(y0,y1,y2,y3)) y = min(y0,y1,y2,y3)

  end subroutine cubic_interp

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
