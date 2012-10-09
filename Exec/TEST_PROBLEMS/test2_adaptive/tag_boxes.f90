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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes(mf,tagboxes,dx,lev,tempbar,aux_tag_mf)

    use variables, only: temp_comp
    use geometry, only: nr_fine, spherical

    type( multifab)          , intent(in   ) :: mf
    type(lmultifab)          , intent(inout) :: tagboxes
    real(dp_t)               , intent(in   ) :: dx
    integer                  , intent(in   ) :: lev
    real(dp_t)               , intent(in   ) :: tempbar(:,0:)
    type( multifab), optional, intent(in   ) :: aux_tag_mf

    real(kind = dp_t), pointer :: sp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(get_dim(mf)), hi(get_dim(mf)), ng_s
    logical           ::      radialtag(0:nr_fine-1)
    logical           :: radialtag_proc(0:nr_fine-1)
    integer           :: dm

    dm  = get_dim(mf)

    radialtag = .false.
    radialtag_proc = .false.

    ng_s = mf%ng

    if (spherical .eq. 0) then
       do i = 1, nfabs(mf)
          sp => dataptr(mf, i)
          lo =  lwb(get_box(tagboxes, i))
          select case (dm)
          case (2)
             call radialtag_2d(radialtag_proc,sp(:,:,1,temp_comp),lo,ng_s,tempbar,lev)
          case (3)
             call radialtag_3d(radialtag_proc,sp(:,:,:,temp_comp),lo,ng_s,tempbar,lev)
          end select
       end do
       
       ! gather radialtag
       call parallel_reduce(radialtag, radialtag_proc, MPI_LOR)
    end if

    do i = 1, nfabs(mf)
       tp => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))
       select case (dm)
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),radialtag,lo,lev)
       case (3)
          if (spherical .eq. 0) then
             call tag_boxes_3d(tp(:,:,:,1),radialtag,lo,lev)
          else
             sp => dataptr(mf, i)
             call tag_boxes_3d_sphr(tp(:,:,:,1),sp(:,:,:,temp_comp),ng_s,tempbar(lev,:), &
                                    lo,hi,dx,lev)
          end if
       end select
    end do

  end subroutine tag_boxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine radialtag_2d(radialtag,mf,lo,ng,tempbar,lev)

    use geometry, only: dr

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:)
    real(dp_t)       , intent(in   ) :: tempbar(:,0:)
    integer, optional, intent(in   ) :: lev

    ! local
    integer         :: i,j,nx,ny,llev
    real(kind=dp_t) :: height

    llev = 1; if (present(lev)) llev = lev
    nx = size(mf,dim=1) - 2*ng
    ny = size(mf,dim=2) - 2*ng

    do j = lo(2),lo(2)+ny-1
       height = (dble(j)+0.d50)*dr(llev)
       do i = lo(1),lo(1)+nx-1
          if (height .gt. 5.4d7 .and. height .lt. 1.8d8) then
             if (abs(mf(i,j)-tempbar(llev,j)) .gt. 3.d7) then
                radialtag(j) = .true.
             end if
          end if
       end do
    enddo

  end subroutine radialtag_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine radialtag_3d(radialtag,mf,lo,ng,tempbar,lev)

    use geometry, only: dr

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(dp_t)       , intent(in   ) :: tempbar(:,0:)
    integer, optional, intent(in   ) :: lev

    ! local
    integer         :: i,j,k,nx,ny,nz,llev
    real(kind=dp_t) :: height

    llev = 1; if (present(lev)) llev = lev
    nx = size(mf,dim=1) - 2*ng
    ny = size(mf,dim=2) - 2*ng
    nz = size(mf,dim=3) - 2*ng

    do k = lo(3),lo(3)+nz-1
       height = (dble(k)+0.5d0)*dr(llev)
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             if (height .gt. 5.4d7 .and. height .lt. 1.8d8) then
                if (abs(mf(i,j,k)-tempbar(llev,k)) .gt. 3.d7) then
                   radialtag(k) = .true.
                end if
             end if
          end do
       enddo
    end do

  end subroutine radialtag_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes_2d(tagbox,radialtag,lo,lev)

    integer          , intent(in   ) :: lo(:)
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
    logical          , intent(in   ) :: radialtag(0:)
    integer, optional, intent(in   ) :: lev
    integer :: j,ny,llev

    llev = 1; if (present(lev)) llev = lev
    ny = size(tagbox,dim=2)

    tagbox = .false.

    ! tag all boxes with radialtag = .true
    select case(llev)
    case (1)
       do j = lo(2),lo(2)+ny-1
          tagbox(:,j) = radialtag(j)
       enddo
    case (2)
       do j = lo(2),lo(2)+ny-1
          tagbox(:,j) = radialtag(j)
       end do
    case default
       do j = lo(2),lo(2)+ny-1
          tagbox(:,j) = radialtag(j)
       end do
    end select

  end subroutine tag_boxes_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes_3d(tagbox,radialtag,lo,lev)

    integer          , intent(in   ) :: lo(:)
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    logical          , intent(in   ) :: radialtag(0:)
    integer, optional, intent(in   ) :: lev

    integer :: k,nz,llev

    llev = 1; if (present(lev)) llev = lev
    nz = size(tagbox,dim=3)

    tagbox = .false.

    ! tag all boxes with radialtag = .true.
    select case(llev)
    case (1)
       do k = lo(3),lo(3)+nz-1
          tagbox(:,:,k) = radialtag(k)
       end do
    case (2)
       do k = lo(3),lo(3)+nz-1
          tagbox(:,:,k) = radialtag(k)
       end do
    case default
       do k = lo(3),lo(3)+nz-1
          tagbox(:,:,k) = radialtag(k)
       end do
    end select

  end subroutine tag_boxes_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes_3d_sphr(tagbox,mf,ng,tempbar,lo,hi,dx,lev)

    use fill_3d_module, only: put_1d_array_on_cart_3d_sphr

    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1)   :,lo(2)   :,lo(3)   :)
    real(kind = dp_t), intent(in   ) ::     mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(dp_t)       , intent(in   ) :: tempbar(0:)
    real(kind=dp_t)  , intent(in   ) :: dx
    integer, optional, intent(in   ) :: lev

    integer    :: i,j,k,llev
    real(dp_t) :: dx_vec(size(lo)),x,y,z,dist

    real(kind=dp_t), allocatable :: tempbar_cart(:,:,:,:)

    dx_vec(:) = dx

    allocate(tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,tempbar,tempbar_cart,lo,hi,dx_vec,0)

    llev = 1; if (present(lev)) llev = lev

    tagbox = .false.

    ! tag all boxes with radialtag = .true.
    select case(llev)
    case default
       do k = lo(3),hi(3)
          z = (dble(k)+0.5d0) * dx
          do j = lo(2),hi(2)
             y = (dble(j)+0.5d0) * dx
             do i = lo(1),hi(1)
                x = (dble(i)+0.5d0) * dx

                dist = sqrt((x-2.5d8)**2 + (y-2.5d8)**2 + (z-2.5d8)**2)

                if (abs(mf(i,j,k)-tempbar_cart(i,j,k,1)).gt.3.d7 .and. dist.lt.1.5d8) then
                   tagbox(i,j,k) = .true.
                end if

             end do
          end do
       end do
    end select

  end subroutine tag_boxes_3d_sphr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tag_boxes_module
