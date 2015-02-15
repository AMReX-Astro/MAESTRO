! This tagging routine assumes that use_tpert_in_tagging = T, so that the
! quantity that comes through in aux_tag_mf is tpert (not rho_Hnuc).

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MUST set this to .true. if tagging uses ghost cells (e.g., tagging on 
  ! gradients). 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical, save :: tagging_needs_ghost_cells = .false.

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    use variables, only: temp_comp
    use geometry, only: nr_fine, spherical

    type( multifab)          , intent(in   ) :: mf
    type(lmultifab)          , intent(inout) :: tagboxes
    real(dp_t)               , intent(in   ) :: dx
    integer                  , intent(in   ) :: lev
    type( multifab), optional, intent(in   ) :: aux_tag_mf

    real(kind = dp_t), pointer :: sp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(get_dim(mf)), hi(get_dim(mf)), ng_s
    logical           ::      radialtag(0:nr_fine-1)
    logical           :: radialtag_proc(0:nr_fine-1)
    integer           :: dm

    dm  = get_dim(aux_tag_mf)

    radialtag = .false.
    radialtag_proc = .false.

    ng_s = aux_tag_mf%ng

    if (spherical .eq. 0) then
       do i = 1, nfabs(aux_tag_mf)
          sp => dataptr(aux_tag_mf, i)
          lo =  lwb(get_box(tagboxes, i))

          select case (dm)
          case (2)
             call radialtag_2d(radialtag_proc,sp(:,:,1,1),lo,ng_s,lev)

          case (3)
             call radialtag_3d(radialtag_proc,sp(:,:,:,1),lo,ng_s,lev)
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
             call bl_error("ERROR: spherical tagging not implemented")
          end if
       end select
    end do

  end subroutine tag_boxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine radialtag_2d(radialtag,mf,lo,ng,lev)

    use geometry, only: dr

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:)
    integer, optional, intent(in   ) :: lev

    ! local
    integer         :: i,j,nx,ny,llev
    real(kind=dp_t) :: height

    llev = 1; if (present(lev)) llev = lev
    nx = size(mf,dim=1) - 2*ng
    ny = size(mf,dim=2) - 2*ng

    ! mf here is tpert

    do j = lo(2),lo(2)+ny-1
       height = (dble(j)+0.5d0)*dr(llev)

       do i = lo(1),lo(1)+nx-1

          if (height > 5.4d7 .and. height < 1.8d8) then
             if (abs(mf(i,j)) > 3.d7) then
                radialtag(j) = .true.
             endif
          endif

       enddo
    enddo

  end subroutine radialtag_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine radialtag_3d(radialtag,mf,lo,ng,lev)

    use geometry, only: dr

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer, optional, intent(in   ) :: lev

    ! local
    integer         :: i,j,k,nx,ny,nz,llev
    real(kind=dp_t) :: height

    llev = 1; if (present(lev)) llev = lev
    nx = size(mf,dim=1) - 2*ng
    ny = size(mf,dim=2) - 2*ng
    nz = size(mf,dim=3) - 2*ng

    ! mf here is tpert

    do k = lo(3),lo(3)+nz-1
       height = (dble(k)+0.5d0)*dr(llev)

       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             if (height > 5.4d7 .and. height < 1.8d8) then
                if (abs(mf(i,j,k)) > 3.d7) then
                   radialtag(k) = .true.
                endif
             endif

          enddo
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

end module tag_boxes_module
