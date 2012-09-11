! a module for storing the geometric information so we don't have to pass it
!
! This module provides the coordinate value for the left edge of a base-state
! zone (r_edge_loc) and the zone center (r_cc_loc).  As always, it is assumed that 
! the base state arrays begin with index 0, not 1.

module geometry

  use bl_types
  use ml_layout_module

  implicit none

  integer   , save :: nlevs_radial
  integer   , save :: spherical
  real(dp_t), save :: center(3)
  integer   , save :: nr_fine, nr_irreg
  real(dp_t), save :: dr_fine
  real(dp_t), allocatable, save :: dr(:), r_cc_loc(:,:), r_edge_loc(:,:)
  integer   , allocatable, save :: numdisjointchunks(:)
  integer   , allocatable, save :: r_start_coord(:,:), r_end_coord(:,:), nr(:)
  integer   , allocatable, save :: anelastic_cutoff_coord(:)
  integer   , allocatable, save :: base_cutoff_density_coord(:)
  integer   , allocatable, save :: burning_cutoff_density_coord(:)
  real(dp_t), save :: sin_theta, cos_theta, omega

  private

  public :: nlevs_radial, spherical, center, nr_fine, dr_fine, nr_irreg
  public :: dr, r_cc_loc, r_edge_loc
  public :: numdisjointchunks
  public :: r_start_coord, r_end_coord, nr
  public :: anelastic_cutoff_coord
  public :: base_cutoff_density_coord
  public :: burning_cutoff_density_coord
  public :: sin_theta, cos_theta, omega

  public :: init_spherical, init_center, init_radial, init_cutoff, &
       compute_cutoff_coords, init_multilevel, init_rotation, destroy_geometry

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_spherical()

    use probin_module, only: spherical_in

    spherical = spherical_in

  end subroutine init_spherical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_center()

    use probin_module, only: prob_lo, prob_hi, octant, dm_in
    use bl_constants_module

    if (.not. octant) then
       center(1:dm_in) = HALF * (prob_lo(1:dm_in) + prob_hi(1:dm_in))
    else
       if (.not. (spherical == 1 .and. dm_in == 3 .and. &
                  prob_lo(1) == ZERO .and. &
                  prob_lo(2) == ZERO .and. &
                  prob_lo(3) == ZERO)) then
          call bl_error("ERROR: octant requires spherical with prob_lo = 0.0")
       endif
       center(1:dm_in) = ZERO
    endif

  end subroutine init_center

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_radial(num_levs,mba)

    ! computes dr, nr, r_cc_loc, r_edge_loc

    use probin_module, only: prob_lo
    use bl_constants_module

    integer          , intent(in   ) :: num_levs
    type(ml_boxarray), intent(inout) :: mba

    ! local
    integer :: n,i

    if (spherical .eq. 0) then
       
       allocate(dr(num_levs))
       allocate(nr(num_levs))

       allocate(  r_cc_loc(num_levs,0:nr_fine-1))
       allocate(r_edge_loc(num_levs,0:nr_fine))
       
       nr(num_levs) = nr_fine
       dr(num_levs) = dr_fine
       do n=num_levs-1,1,-1
          nr(n) = nr(n+1)/mba%rr(n,mba%dim)
          dr(n) = dr(n+1)*dble(mba%rr(n,mba%dim))
       enddo
       
       do n=1,num_levs
          do i = 0,nr(n)-1
             r_cc_loc(n,i) = prob_lo(mba%dim) + (dble(i)+HALF)*dr(n)
          end do
          do i = 0,nr(n)
             r_edge_loc(n,i) = prob_lo(mba%dim) + (dble(i))*dr(n)
          end do
       enddo

    else

       allocate(dr(1))
       allocate(nr(1))

       allocate(  r_cc_loc(1,0:nr_fine-1))
       allocate(r_edge_loc(1,0:nr_fine))
       
       nr(1) = nr_fine
       dr(1) = dr_fine

       do i=0,nr_fine-1
          r_cc_loc(1,i) = (dble(i)+HALF)*dr(1)
       end do
       
       do i=0,nr_fine
          r_edge_loc(1,i) = (dble(i))*dr(1)
       end do

    end if

  end subroutine init_radial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_cutoff(num_levs)

    ! allocates the cutoff coordinate arrays

    integer, intent(in)    :: num_levs

    if (spherical .eq. 0) then
       allocate(      anelastic_cutoff_coord(num_levs))
       allocate(   base_cutoff_density_coord(num_levs))
       allocate(burning_cutoff_density_coord(num_levs))
    else
       allocate(      anelastic_cutoff_coord(1))
       allocate(   base_cutoff_density_coord(1))
       allocate(burning_cutoff_density_coord(1))
    end if

  end subroutine init_cutoff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_cutoff_coords(rho0)

    use probin_module, only: anelastic_cutoff, base_cutoff_density, burning_cutoff_density

    real(kind=dp_t), intent(in   ) :: rho0(:,0:)

    ! local
    integer :: i,n,r,which_lev
    logical :: found

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute the coordinates of the anelastic cutoff
    found = .false.

    ! find the finest level containing the anelastic cutoff density,
    ! and set the anelastic cutoff coord for this level
    do n=nlevs_radial,1,-1
       do i=1,numdisjointchunks(n)

          if (.not. found) then
             do r=r_start_coord(n,i),r_end_coord(n,i)
                if (rho0(n,r) .le. anelastic_cutoff) then
                   anelastic_cutoff_coord(n) = r
                   which_lev = n
                   found = .true.
                   exit
                end if
             end do
          end if

       end do
    end do

    ! if the anelastic cutoff density was not found anywhere, then set
    ! it to above the top of the domain on the finest level
    if (.not. found) then
       which_lev = nlevs_radial
       anelastic_cutoff_coord(nlevs_radial) = nr(nlevs_radial)
    endif

    ! set the anelastic cutoff coordinate on the finer levels
    do n=which_lev+1,nlevs_radial
       anelastic_cutoff_coord(n) = 2*anelastic_cutoff_coord(n-1)+1
    end do

    ! set the anelastic cutoff coordinate on the coarser levels
    do n=which_lev-1,1,-1
       if (mod(anelastic_cutoff_coord(n+1),2) .eq. 0) then
          anelastic_cutoff_coord(n) = anelastic_cutoff_coord(n+1) / 2
       else
          anelastic_cutoff_coord(n) = anelastic_cutoff_coord(n+1) / 2 + 1
       end if
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute the coordinates of the base cutoff density
    found = .false.

    ! find the finest level containing the base cutoff density,
    ! and set the base cutoff coord for this level
    do n=nlevs_radial,1,-1
       do i=1,numdisjointchunks(n)

          if (.not. found) then
             do r=r_start_coord(n,i),r_end_coord(n,i)
                if (rho0(n,r) .le. base_cutoff_density) then
                   base_cutoff_density_coord(n) = r
                   which_lev = n
                   found = .true.
                   exit
                end if
             end do
          end if

       end do
    end do

    ! if the base cutoff density was not found anywhere, then set
    ! it to above the top of the domain on the finest level
    if (.not. found) then
       which_lev = nlevs_radial
       base_cutoff_density_coord(nlevs_radial) = nr(nlevs_radial)
    endif

    ! set the base cutoff coordinate on the finer levels
    do n=which_lev+1,nlevs_radial
       base_cutoff_density_coord(n) = 2*base_cutoff_density_coord(n-1)+1
    end do

    ! set the base cutoff coordinate on the coarser levels
    do n=which_lev-1,1,-1
       if (mod(base_cutoff_density_coord(n+1),2) .eq. 0) then
          base_cutoff_density_coord(n) = base_cutoff_density_coord(n+1) / 2
       else
          base_cutoff_density_coord(n) = base_cutoff_density_coord(n+1) / 2 + 1
       end if
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute the coordinates of the burning cutoff density
    found = .false.

    ! find the finest level containing the burning cutoff density,
    ! and set the burning cutoff coord for this level
    do n=nlevs_radial,1,-1
       do i=1,numdisjointchunks(n)

          if (.not. found) then
             do r=r_start_coord(n,i),r_end_coord(n,i)
                if (rho0(n,r) .le. burning_cutoff_density) then
                   burning_cutoff_density_coord(n) = r
                   which_lev = n
                   found = .true.
                   exit
                end if
             end do
          end if

       end do
    end do
 
    ! if the burning cutoff density was not found anywhere, then set
    ! it to above the top of the domain on the finest level
    if (.not. found) then
       which_lev = nlevs_radial
       burning_cutoff_density_coord(nlevs_radial) = nr(nlevs_radial)
    endif

    ! set the burning cutoff coordinate on the finer levels 
    do n=which_lev+1,nlevs_radial
       burning_cutoff_density_coord(n) = 2*burning_cutoff_density_coord(n-1)+1
    end do

    ! set the burning cutoff coordinate on the coarser levels
    do n=which_lev-1,1,-1
       if (mod(burning_cutoff_density_coord(n+1),2) .eq. 0) then
          burning_cutoff_density_coord(n) = burning_cutoff_density_coord(n+1) / 2
       else
          burning_cutoff_density_coord(n) = burning_cutoff_density_coord(n+1) / 2 + 1
       end if
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine compute_cutoff_coords

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_multilevel(mf)

    ! computes numdisjointchunks, r_start_coord, r_end_coord

    use bl_constants_module

    type(multifab) , intent(in   ) :: mf(:)

    ! local
    integer :: i,j,n,maxdisjointchunks,temp,dm,nlevs

    integer :: lo(get_dim(mf(1))),hi(get_dim(mf(1)))

    type(boxarray) ::  validboxarr(size(mf))
    type(boxarray) :: diffboxarray(size(mf))
    type(box)      ::  boundingbox(size(mf))
    
    dm = get_dim(mf(1))
    nlevs = size(mf)

    if (spherical .eq. 0) then
    
       ! create a "bounding box" for each level
       ! this the smallest possible box that fits every grid at a particular level
       ! this even includes the empty spaces if there are gaps between grids
       do n=1,nlevs
          boundingbox(n) = get_box(mf(n)%la,1)
          do i=2, nboxes(mf(n)%la)
             boundingbox(n) = box_bbox(boundingbox(n),get_box(mf(n)%la,i))
          end do
       end do
       
       ! compute diffboxarray
       ! each box in diffboxarray corresponds to an "empty space" between valid regions
       ! at each level, excluding the coarsest level.
       ! I am going to use this to compute all of the intermediate r_start_coord and 
       ! r_end_coord
       do n=1,nlevs
          call boxarray_build_copy(validboxarr(n),get_boxarray(mf(n)))
          call boxarray_boxarray_diff(diffboxarray(n),boundingbox(n),validboxarr(n))
          call boxarray_simplify(diffboxarray(n))
       end do
       
       if (allocated(numdisjointchunks)) then
          deallocate(numdisjointchunks)
       end if
       allocate(numdisjointchunks(nlevs))
       
       do n=1,nlevs
          numdisjointchunks(n) = nboxes(diffboxarray(n)) + 1
       end do
       
       maxdisjointchunks = 1
       do n=2,nlevs
          maxdisjointchunks = max(maxdisjointchunks,numdisjointchunks(n))
       end do
       
       if (allocated(r_start_coord)) then
          deallocate(r_start_coord)
       end if
       allocate(r_start_coord(nlevs,maxdisjointchunks))
       if (allocated(r_end_coord)) then
          deallocate(r_end_coord)
       end if
       allocate(r_end_coord(nlevs,maxdisjointchunks))
       
       do n=1,nlevs

          lo = lwb(boundingbox(n))
          hi = upb(boundingbox(n))
          r_start_coord(n,1) = lo(dm)
          r_end_coord(n,1)   = hi(dm)

          if (nboxes(diffboxarray(n)) .gt. 0) then
             do i=1,nboxes(diffboxarray(n))
                lo = lwb(boxarray_get_box(diffboxarray(n),i))
                hi = upb(boxarray_get_box(diffboxarray(n),i))
                r_start_coord(n,i+1) = hi(dm)+1
                r_end_coord(n,i+1)   = lo(dm)-1
             end do

             ! sort start and end coords
             do i=1,nboxes(diffboxarray(n))+1
                do j=1,nboxes(diffboxarray(n))+1-i

                   if (r_start_coord(n,j) .gt. r_start_coord(n,j+1)) then
                      temp = r_start_coord(n,j+1)
                      r_start_coord(n,j+1) = r_start_coord(n,j)
                      r_start_coord(n,j)   = temp
                   end if
                   if (r_end_coord(n,j) .gt. r_end_coord(n,j+1)) then
                      temp = r_end_coord(n,j+1)
                      r_end_coord(n,j+1) = r_end_coord(n,j)
                      r_end_coord(n,j)   = temp
                   end if

                end do
             end do
          end if

       end do

       do n=1,nlevs
          call destroy(validboxarr(n))
          call destroy(diffboxarray(n))
       end do

    else

       if (allocated(numdisjointchunks)) deallocate(numdisjointchunks)
       if (allocated(r_start_coord)) deallocate(r_start_coord)
       if (allocated(r_end_coord)) deallocate(r_end_coord)

       allocate(numdisjointchunks(1))
       allocate(r_start_coord(1,1))
       allocate(r_end_coord(1,1))

       numdisjointchunks(1) = 1
       r_start_coord(1,1) = 0
       r_end_coord(1,1) = nr_fine-1

    end if
  end subroutine init_multilevel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! rotation
! 
! co_latitude, rotation_radius, theta_in_rad, sin_theta and cos_theta are only
! important when wanting to rotate a plane-parallel patch which lives at an 
! angle co_latitude from the rotation axis and at a distance rotation_radius 
! from center().  mk_vel_force should calculate the rotational forcing terms
! for the points within the patch.  
!
! for spherical problems, the only important variable from init_rotation() is 
! omega, the angular frequency - mk_vel_force will calculate the appropriate
! terms for a given coordinate
!
! NOTE: it is currently unclear how to handle BC's with a plane-parallel patch
!       and it is not advisable to utilize rotation for such problems.
!

  subroutine init_rotation()

    use probin_module, only: rotational_frequency, co_latitude
    use bl_constants_module, only: M_PI, ZERO

    real(dp_t) :: theta_in_rad

    theta_in_rad = M_PI * co_latitude / 180_dp_t
    
    sin_theta = sin(theta_in_rad)
    cos_theta = cos(theta_in_rad)

    omega = 2 * M_PI * rotational_frequency

  end subroutine init_rotation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine destroy_geometry()

    deallocate(dr,r_cc_loc,r_edge_loc,r_start_coord,r_end_coord,nr,numdisjointchunks)
    deallocate(anelastic_cutoff_coord,base_cutoff_density_coord,burning_cutoff_density_coord)

  end subroutine destroy_geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module geometry
