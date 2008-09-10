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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes(mf,tagboxes,lev)

    use variables, only: rho_comp, spec_comp
    use geometry, only: dm, nr_fine

    type( multifab), intent(in   ) :: mf
    type(lmultifab), intent(inout) :: tagboxes
    integer        , intent(in   ) :: lev

    real(kind = dp_t), pointer :: sp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(dm), ng_s
    logical           ::      radialtag(0:nr_fine-1)
    logical           :: radialtag_proc(0:nr_fine-1)

    radialtag = .false.
    radialtag_proc = .false.

    ng_s = mf%ng

    do i = 1, mf%nboxes
       if ( multifab_remote(mf, i) ) cycle
       sp => dataptr(mf, i)
       lo =  lwb(get_box(tagboxes, i))
       select case (dm)
       case (2)
          call radialtag_2d(radialtag_proc,sp(:,:,1,spec_comp),sp(:,:,1,rho_comp),lo,ng_s,lev)
       case  (3)
          call radialtag_3d(radialtag_proc,sp(:,:,:,spec_comp),sp(:,:,:,rho_comp),lo,ng_s,lev)
       end select
    end do

    ! gather radialtag
    call parallel_reduce(radialtag, radialtag_proc, MPI_LOR)

    do i = 1, mf%nboxes
       if ( multifab_remote(mf, i) ) cycle
       tp => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       select case (dm)
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),radialtag,lo)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),radialtag,lo)
       end select
    end do

  end subroutine tag_boxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine radialtag_2d(radialtag,he,rho,lo,ng,lev)

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) ::  he(lo(1)-ng:,lo(2)-ng:)
    real(kind = dp_t), intent(in   ) :: rho(lo(1)-ng:,lo(2)-ng:)
    integer, optional, intent(in   ) :: lev

    ! local
    integer :: i,j,nx,ny,llev
    real(kind=dp_t) :: Xhe

    llev = 1; if (present(lev)) llev = lev
    nx = size(he,dim=1) - 2*ng
    ny = size(he,dim=2) - 2*ng

    select case(llev)
    case (1)
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             Xhe = he(i,j)/rho(i,j)
             if (Xhe .gt. 1.d-4 .and. Xhe .lt. 0.99d0) then
                radialtag(j) = .true.
             end if
          end do
       enddo
    case (2)
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             Xhe = he(i,j)/rho(i,j)
             if (Xhe .gt. 1.d-1 .and. Xhe .lt. 0.99d0) then
                radialtag(j) = .true.
             end if
          end do
       end do
    case default
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             Xhe = he(i,j)/rho(i,j)
             if (Xhe .gt. 1.d-1 .and. Xhe .lt. 0.99d0) then
                radialtag(j) = .true.
             end if
          end do
       end do
    end select

  end subroutine radialtag_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine radialtag_3d(radialtag,he,rho,lo,ng,lev)

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) ::  he(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind = dp_t), intent(in   ) :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer, optional, intent(in   ) :: lev

    ! local
    integer :: i,j,k,nx,ny,nz,llev
    real(kind=dp_t) :: Xhe

    llev = 1; if (present(lev)) llev = lev
    nx = size(he,dim=1) - 2*ng
    ny = size(he,dim=2) - 2*ng
    nz = size(he,dim=3) - 2*ng

    select case(llev)
    case (1)
       do k = lo(3),lo(3)+nz-1
          do j = lo(2),lo(2)+ny-1
             do i = lo(1),lo(1)+nx-1
                Xhe = he(i,j,k)/rho(i,j,k)
                if (Xhe .gt. 1.d-4 .and. Xhe .lt. 0.99d0) then
                   radialtag(k) = .true.
                end if
             end do
          enddo
       end do
    case (2)
       do k = lo(3),lo(3)+nz-1
          do j = lo(2),lo(2)+ny-1
             do i = lo(1),lo(1)+nx-1
                Xhe = he(i,j,k)/rho(i,j,k)
                if (Xhe .gt. 1.d-1 .and. Xhe .lt. 0.99d0) then
                   radialtag(k) = .true.
                end if
             end do
          end do
       end do
    case default
       do k = lo(3),lo(3)+nz-1
          do j = lo(2),lo(2)+ny-1
             do i = lo(1),lo(1)+nx-1
                Xhe = he(i,j,k)/rho(i,j,k)
                if (Xhe .gt. 1.d-1 .and. Xhe .lt. 0.99d0) then
                   radialtag(k) = .true.
                end if
             end do
          end do
       end do
    end select

  end subroutine radialtag_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes_2d(tagbox,radialtag,lo)

    integer          , intent(in   ) :: lo(:)
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
    logical          , intent(in   ) :: radialtag(0:)

    integer :: j,ny

    ny = size(tagbox,dim=2)

    tagbox = .false.

    ! tag all boxes with radialtag = .true
    do j = lo(2),lo(2)+ny-1
       tagbox(:,j) = radialtag(j)
    enddo

  end subroutine tag_boxes_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes_3d(tagbox,radialtag,lo)

    integer          , intent(in   ) :: lo(:)
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    logical          , intent(in   ) :: radialtag(0:)

    integer :: k,nz

    nz = size(tagbox,dim=3)

    tagbox = .false.

    ! tag all boxes with radialtag = .true.
    do k = lo(3),lo(3)+nz-1
       tagbox(:,:,k) = radialtag(k)
    end do

  end subroutine tag_boxes_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tag_boxes_module
