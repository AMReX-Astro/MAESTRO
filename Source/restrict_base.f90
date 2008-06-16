module restrict_base_module

  use bl_types

  implicit none

  private

  public :: restrict_base, fill_ghost_base

contains

  subroutine restrict_base(nlevs,s0,is_cell_centered)

    use bl_prof_module
    use geometry, only: r_start_coord, r_end_coord

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(inout) :: s0(:,0:)
    logical        , intent(in   ) :: is_cell_centered

    ! local
    integer :: n, r

    type(bl_prof_timer), save :: bpt

    call build(bpt, "restrict_base")

    if (is_cell_centered) then

       do n=nlevs,2,-1
          ! for level n, make the coarser cells underneath simply the average of the fine
          do r=r_start_coord(n),r_end_coord(n)-1,2
             s0(n-1,r/2) = 0.5d0 * (s0(n,r) + s0(n,r+1))
          end do
       end do

    else

       do n=nlevs,2,-1
          ! for level n, make the coarse edge underneath equal to the fine edge value
          do r=r_start_coord(n),r_end_coord(n)+1,2
             s0(n-1,r/2) = s0(n,r)
          end do
       end do

    end if

    call destroy(bpt)

  end subroutine restrict_base

  subroutine fill_ghost_base(nlevs,s0,is_cell_centered)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: r_start_coord, r_end_coord, nr

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(inout) :: s0(:,0:)
    logical        , intent(in   ) :: is_cell_centered

    ! local
    integer :: n, r

    type(bl_prof_timer), save :: bpt

    call build(bpt, "fill_ghost_base")

    if (is_cell_centered) then

       do n=nlevs,2,-1

          if (r_start_coord(n) .ne. 0) then
             s0(n,r_start_coord(n)-1) = -(1.d0/5.d0)*s0(n,r_start_coord(n)+1) &
                  + (2.d0/3.d0)*s0(n,r_start_coord(n)) &
                  + (8.d0/15.d0)*s0(n-1,r_start_coord(n)/2-1)
          end if

          if (r_end_coord(n) .ne. nr(n)-1) then
             s0(n,r_end_coord(n)+1) = -(1.d0/5.d0)*s0(n,r_end_coord(n)-1) &
                  + (2.d0/3.d0)*s0(n,r_end_coord(n)) &
                  + (8.d0/15.d0)*s0(n-1,(r_end_coord(n)+1)/2)
          end if

       end do

    else

       do n=nlevs,2,-1

          if (r_start_coord(n) .ne. 0) then
             s0(n,r_start_coord(n)-1) = -THIRD*s0(n,r_start_coord(n)+1) &
                  + s0(n,r_start_coord(n)) + THIRD*s0(n-1,r_start_coord(n)/2-1)
          end if
          
          if (r_end_coord(n)+1 .ne. nr(n)) then
             s0(n,r_end_coord(n)+1) = -THIRD*s0(n,r_end_coord(n)-1) &
                  + s0(n,r_end_coord(n)) + THIRD*s0(n-1,(r_end_coord(n)+1)/2)
          end if

       end do

    end if

  end subroutine fill_ghost_base


end module restrict_base_module
