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
    use geometry, only: dm, nr_fine, spherical

    type( multifab)          , intent(in   ) :: mf
    type(lmultifab)          , intent(inout) :: tagboxes
    real(dp_t)               , intent(in   ) :: dx
    integer                  , intent(in   ) :: lev
    real(dp_t)               , intent(in   ) :: tempbar(:,0:)
    type( multifab), optional, intent(in   ) :: aux_tag_mf

    real(kind = dp_t), pointer :: sp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(dm), hi(dm), ng_s
    logical           ::      radialtag(0:nr_fine-1)
    logical           :: radialtag_proc(0:nr_fine-1)

    radialtag = .false.
    radialtag_proc = .false.

    ng_s = mf%ng

    do i = 1, nfabs(mf)
       tp => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))
       select case (dm)
       case (2)

       case (3)
          if (spherical .eq. 0) then

          else
             sp => dataptr(mf, i)
             call tag_boxes_3d_sphr(tp(:,:,:,1),sp(:,:,:,temp_comp),ng_s,tempbar(1,:), &
                                    lo,hi,dx,lev)
          end if
       end select
    end do

  end subroutine tag_boxes

  subroutine tag_boxes_3d_sphr(tagbox,mf,ng,tempbar,lo,hi,dx,lev)

    use fill_3d_module, only: put_1d_array_on_cart_3d_sphr
    use geometry, only: dm
    use probin_module, only: n_cellx, regrid_int

    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1)   :,lo(2)   :,lo(3)   :)
    real(kind = dp_t), intent(in   ) ::     mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(dp_t)       , intent(in   ) :: tempbar(0:)
    real(kind=dp_t)  , intent(in   ) :: dx
    integer, optional, intent(in   ) :: lev

    integer    :: i,j,k,llev,lotag,hitag,n_cellx_lev
    real(dp_t) :: dx_vec(dm),x,y,z,dist

    real(kind=dp_t), allocatable :: tempbar_cart(:,:,:,:)

    dx_vec(:) = dx

    allocate(tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,tempbar,tempbar_cart,lo,hi,dx_vec,0)

    llev = 1; if (present(lev)) llev = lev

    n_cellx_lev = n_cellx*2**(lev-1)

    ! want a 16^3 refined area
    ! subtract off regrid_int since it will be grown later
    lotag = n_cellx_lev/2 - 4 + min(3,regrid_int)
    hitag = n_cellx_lev/2 + 3 - min(3,regrid_int)

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

                if (i .ge. lotag .and. i .le. hitag .and. &
                    j .ge. lotag .and. j .le. hitag .and. &
                    k .ge. lotag .and. k .le. hitag) then
                   tagbox(i,j,k) = .true.
                end if

             end do
          end do
       end do
    end select

  end subroutine tag_boxes_3d_sphr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tag_boxes_module
