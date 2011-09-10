module restrict_base_module

  use bl_types

  implicit none

  private

  public :: restrict_base, fill_ghost_base

contains

  subroutine restrict_base(s0,is_cell_centered)

    use bl_prof_module
    use geometry, only: r_start_coord, r_end_coord, numdisjointchunks, nlevs_radial

    real(kind=dp_t), intent(inout) :: s0(:,0:)
    logical        , intent(in   ) :: is_cell_centered

    ! local
    integer :: n, r, i

    type(bl_prof_timer), save :: bpt

    call build(bpt, "restrict_base")

    if (is_cell_centered) then

       do n=nlevs_radial,2,-1
          ! for level n, make the coarser cells underneath simply the average of the fine
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)-1,2
                s0(n-1,r/2) = 0.5d0 * (s0(n,r) + s0(n,r+1))
             end do
          end do
       end do

    else

       do n=nlevs_radial,2,-1
          ! for level n, make the coarse edge underneath equal to the fine edge value
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,1),r_end_coord(n,1)+1,2
                s0(n-1,r/2) = s0(n,r)
             end do
          end do
       end do

    end if

    call destroy(bpt)

  end subroutine restrict_base

  subroutine fill_ghost_base(s0,is_cell_centered)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: r_start_coord, r_end_coord, nr, numdisjointchunks, nlevs_radial
    use bl_error_module

    real(kind=dp_t), intent(inout) :: s0(:,0:)
    logical        , intent(in   ) :: is_cell_centered

    ! local
    integer         :: n,i,r_crse
    real(kind=dp_t) :: del,dpls,dmin,slim,slope

    type(bl_prof_timer), save :: bpt

    call build(bpt, "fill_ghost_base")

    if (is_cell_centered) then

       ! compute limited slopes at the coarse level and set 4 fine ghost cells 
       ! maintaining conservation
       do n=nlevs_radial,2,-1
          do i=1,numdisjointchunks(n)

             ! lo side
             if (r_start_coord(n,i) .ne. 0) then
                r_crse = r_start_coord(n,i)/2-1
                del = HALF*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1))
                dpls = TWO*(s0(n-1,r_crse+1)-s0(n-1,r_crse  ))
                dmin = TWO*(s0(n-1,r_crse  )-s0(n-1,r_crse-1))
                slim = min(abs(dpls),abs(dmin))
                slim = merge(slim,ZERO,dpls*dmin.gt.ZERO)
                slope=sign(ONE,del)*min(slim,abs(del))
                s0(n,r_start_coord(n,i)-1) = s0(n-1,r_crse) + FOURTH*slope
                s0(n,r_start_coord(n,i)-2) = s0(n-1,r_crse) - FOURTH*slope

                r_crse = r_start_coord(n,i)/2-2
                del = HALF*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1))
                dpls = TWO*(s0(n-1,r_crse+1)-s0(n-1,r_crse  ))
                dmin = TWO*(s0(n-1,r_crse  )-s0(n-1,r_crse-1))
                slim = min(abs(dpls),abs(dmin))
                slim = merge(slim,ZERO,dpls*dmin.gt.ZERO)
                slope=sign(ONE,del)*min(slim,abs(del))
                s0(n,r_start_coord(n,i)-3) = s0(n-1,r_crse) + FOURTH*slope
                s0(n,r_start_coord(n,i)-4) = s0(n-1,r_crse) - FOURTH*slope
             end if
             
             ! hi side
             if (r_end_coord(n,i) .ne. nr(n)-1) then
                r_crse = (r_end_coord(n,i)+1)/2
                del = HALF*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1))
                dpls = TWO*(s0(n-1,r_crse+1)-s0(n-1,r_crse  ))
                dmin = TWO*(s0(n-1,r_crse  )-s0(n-1,r_crse-1))
                slim = min(abs(dpls),abs(dmin))
                slim = merge(slim,ZERO,dpls*dmin.gt.ZERO)
                slope=sign(ONE,del)*min(slim,abs(del))
                s0(n,r_end_coord(n,i)+1) = s0(n-1,r_crse) - FOURTH*slope
                s0(n,r_end_coord(n,i)+2) = s0(n-1,r_crse) + FOURTH*slope

                r_crse = (r_end_coord(n,i)+1)/2+1
                del = HALF*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1))
                dpls = TWO*(s0(n-1,r_crse+1)-s0(n-1,r_crse  ))
                dmin = TWO*(s0(n-1,r_crse  )-s0(n-1,r_crse-1))
                slim = min(abs(dpls),abs(dmin))
                slim = merge(slim,ZERO,dpls*dmin.gt.ZERO)
                slope=sign(ONE,del)*min(slim,abs(del))
                s0(n,r_end_coord(n,i)+3) = s0(n-1,r_crse) - FOURTH*slope
                s0(n,r_end_coord(n,i)+4) = s0(n-1,r_crse) + FOURTH*slope
             end if

          end do
       end do

    else

       do n=nlevs_radial,2,-1
          do i=1,numdisjointchunks(n)

             if (r_start_coord(n,i) .ne. 0) then
                ! quadratic interpolation from the three closest points
                s0(n,r_start_coord(n,i)-1) = -THIRD*s0(n,r_start_coord(n,i)+1) &
                     + s0(n,r_start_coord(n,i)) + THIRD*s0(n-1,r_start_coord(n,i)/2-1)
                ! copy the next ghost cell value directly in from the coarser level
                s0(n,r_start_coord(n,i)-2) = s0(n-1,(r_start_coord(n,i)-2)/2) 
             end if
             
             if (r_end_coord(n,i)+1 .ne. nr(n)) then
                ! quadratic interpolation from the three closest points
                s0(n,r_end_coord(n,i)+2) = -THIRD*s0(n,r_end_coord(n,i)) &
                     + s0(n,r_end_coord(n,i)+1) + THIRD*s0(n-1,(r_end_coord(n,i)+3)/2)
                ! copy the next ghost cell value directly in from the coarser level
                s0(n,r_end_coord(n,i)+3) = s0(n-1,(r_end_coord(n,i)+3)/2)
             end if

          end do
       end do

    end if

    call destroy(bpt)

  end subroutine fill_ghost_base


end module restrict_base_module
