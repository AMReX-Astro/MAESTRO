module tag_boxes_module

  use BoxLib
  use omp_module
  use f2kcli
  use list_box_module
  use boxarray_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use box_util_module
  use bl_IO_module
  use cluster_module
  use ml_layout_module

  implicit none 

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev)

    type( multifab), intent(in   ) :: mf
    type(lmultifab), intent(inout) :: tagboxes
    real(dp_t)     , intent(in   ) :: dx
    integer        , intent(in   ) :: lev

    real(kind = dp_t), pointer :: sp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(mf%dim), ng

    ng = mf%ng

    do i = 1, mf%nboxes
       if ( multifab_remote(mf, i) ) cycle
       sp => dataptr(mf, i)
       tp => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))

       select case (mf%dim)
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),sp(:,:,1,1),lo,ng,lev)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),sp(:,:,:,1),lo,ng,lev)
       end select
    end do

  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,mf,lo,ng,lev)

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:)
    integer, optional, intent(in   ) :: lev
    integer :: i,j,nx,ny,llev

    llev = 1; if (present(lev)) llev = lev
    nx = size(tagbox,dim=1)
    ny = size(tagbox,dim=2)

    tagbox = .false.

    ! tag all boxes with 3.9d7 <= density <= 2.8d7
    ! if you tag a cell, tag all other cells at that height
    do j = lo(2),lo(2)+ny-1
       do i = lo(1),lo(1)+nx-1
          if (mf(i,j) .ge. 2.8d7 .and. mf(i,j) .le. 3.9d7) then
             tagbox(:,j) = .true.
          end if
       end do
    enddo

  end subroutine tag_boxes_2d

  subroutine tag_boxes_3d(tagbox,mf,lo,ng,lev)

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer, optional, intent(in   ) :: lev

    integer :: i,j,k,nx,ny,nz,llev

    llev = 1; if (present(lev)) llev = lev
    nx = size(tagbox,dim=1)
    ny = size(tagbox,dim=2)
    nz = size(tagbox,dim=3)

    tagbox = .false.

    ! tag all boxes with 3.9d7 <= density <= 2.8d7
    ! if you tag a cell, tag all other cells at that height
    do k = lo(3),lo(3)+nz-1
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             if (mf(i,j,k) .ge. 2.8d7 .and. mf(i,j,k) .le. 3.9e7) then
                tagbox(:,:,k) = .true.
             end if
          end do
       enddo
    end do

  end subroutine tag_boxes_3d

end module tag_boxes_module
