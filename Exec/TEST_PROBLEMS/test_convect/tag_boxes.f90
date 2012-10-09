module tag_boxes_module

  use BoxLib
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

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    type( multifab)          , intent(in   ) :: mf
    type(lmultifab)          , intent(inout) :: tagboxes
    real(dp_t)               , intent(in   ) :: dx
    integer                  , intent(in   ) :: lev
    type( multifab), optional, intent(in   ) :: aux_tag_mf

    real(kind = dp_t), pointer :: sp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(mf%dim), hi(mf%dim), ng

    ng = mf%ng

    do i = 1, nfabs(mf)
       sp => dataptr(mf, i)
       tp => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))
       select case (mf%dim)
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),sp(:,:,1,1),lo,hi,ng,lev)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),sp(:,:,:,1),lo,hi,ng,lev)
       end select
    end do

  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,mf,lo,hi,ng,lev)

    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:)
    integer, optional, intent(in   ) :: lev
    integer :: i,j,llev

    llev = 1; if (present(lev)) llev = lev

    tagbox = .false.

    select case(llev)
    case (1)
       do j = lo(2),hi(2)
          if (j .ge. 128 .and. j .le. 191) then
             tagbox(:,j) = .true.
          end if
       enddo
    case (2)
       do j = lo(2),hi(2)
          if (j .ge. 256 .and. j .le. 383) then
             tagbox(:,j) = .true.
          end if
       end do
    case (3)
       do j = lo(2),hi(2)
          if (j .ge. 512 .and. j .le. 767) then
             tagbox(:,j) = .true.
          end if
       end do
    case default
       call bl_error("tag_boxes.f90: Need to write tagging condition for this level")
    end select

  end subroutine tag_boxes_2d

  subroutine tag_boxes_3d(tagbox,mf,lo,hi,ng,lev)

    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer, optional, intent(in   ) :: lev

    integer :: i,j,k,llev

    llev = 1; if (present(lev)) llev = lev

    tagbox = .false.

    select case(llev)
    case (1)
       do k = lo(3),hi(3)
          if (k .ge. 128 .and. k .le. 191) then
             tagbox(:,:,k) = .true.
          end if
       end do
    case (2)
       do k = lo(3),hi(3)
          if (k .ge. 256 .and. k .le. 383) then
             tagbox(:,:,k) = .true.
          end if
       end do
    case (3)
       do k = lo(3),hi(3)
          if (k .ge. 512 .and. k .le. 767) then
             tagbox(:,:,k) = .true.
          end if
       end do
    case default
       call bl_error("tag_boxes.f90: Need to write tagging condition for this level")
    end select

  end subroutine tag_boxes_3d

end module tag_boxes_module
