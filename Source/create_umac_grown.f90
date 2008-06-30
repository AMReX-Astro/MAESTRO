module create_umac_grown_module

  use bl_constants_module
  use sparse_solve_module

  private
  public :: create_umac_grown, create_umac_grown_onesided

contains

  subroutine create_umac_grown()




  end subroutine create_umac_grown

  subroutine create_umac_grown_onesided(nlevs,umac,do_normal,do_trans)

    integer       , intent(in   ) :: nlevs
    type(multifab), intent(inout) :: umac(:,:)
    logical       , intent(in   ) :: do_normal, do_trans

    integer :: i,n,dm
    integer :: lo(umac(1,1)%dim),hi(umac(1,1)%dim)

    real(kind=dp_t), pointer :: ump(:,:,:,:) 
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)

    dm = size(umac,dim=2)

    ! we only need to do this for fine levels
    do n=2,nlevs
       do i=1,umac(n,1)%nboxes
          if ( multifab_remote(umac(n,1), i) ) cycle
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)
          lo = lwb(get_box(umac(n,1), i))
          hi = upb(get_box(umac(n,1), i))
          select case (dm)
          case (2)
             call create_umac_grown_onesided_2d(ump(:,:,1,1),vmp(:,:,1,1),do_normal, &
                                                do_trans,lo,hi)
          case (3)
             call create_umac_grown_onesided_3d(ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1), &
                                                do_normal,do_trans,lo,hi)
          end select
       end do
    end do

  end subroutine create_umac_grown_onesided

  subroutine create_umac_grown_onesided_2d(umac,vmac,do_normal,do_trans,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: umac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) :: vmac(lo(1)-1:,lo(2)-1:)
    logical        , intent(in   ) :: do_normal, do_trans

    integer i,j

    if (do_normal) then

       do j=lo(2),hi(2)
          umac(lo(1)-1,j) = TWO*umac(lo(1),j) - umac(lo(1)+1,j)
          umac(hi(1)+2,j) = TWO*umac(hi(1)+1,j) - umac(hi(1),j)
       end do

       do i=lo(1),hi(1)
          vmac(i,lo(2)-1) = TWO*vmac(i,lo(2)) - vmac(i,lo(2)+1)
          vmac(i,hi(2)+2) = TWO*vmac(i,hi(2)+1) - vmac(i,hi(2))
       end do

    end if

    if (do_trans) then

       ! use linear interpolation to fill fine level ghost cells needed for 
       ! transverse derivatives
       do i=lo(1),hi(1)+1
          umac(i,lo(2)-1) = TWO*umac(i,lo(2)) - umac(i,lo(2)+1)
          umac(i,hi(2)+1) = TWO*umac(i,hi(2)) - umac(i,hi(2)-1)
       end do
       
       do j=lo(2),hi(2)+1
          vmac(lo(1)-1,j) = TWO*vmac(lo(1),j) - vmac(lo(1)+1,j)
          vmac(hi(1)+1,j) = TWO*vmac(hi(1),j) - vmac(hi(1)-1,j)
       end do

    end if

  end subroutine create_umac_grown_onesided_2d

  subroutine create_umac_grown_onesided_3d(umac,vmac,wmac,do_normal,do_trans,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(inout) :: vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(inout) :: wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    logical        , intent(in   ) :: do_normal, do_trans

    integer i,j,k

    if (do_normal) then

       do j=lo(2),hi(2)
          do k=lo(3),hi(3)
             umac(lo(1)-1,j,k) = TWO*umac(lo(1),j,k) - umac(lo(1)+1,j,k)
             umac(hi(1)+2,j,k) = TWO*umac(hi(1)+1,j,k) - umac(hi(1),j,k)
          end do
       end do

       do i=lo(1),hi(1)
          do k=lo(3),hi(3)
             vmac(i,lo(2)-1,k) = TWO*vmac(i,lo(2),k) - vmac(i,lo(2)+1,k)
             vmac(i,hi(2)+2,k) = TWO*vmac(i,hi(2)+1,k) - vmac(i,hi(2),k)
          end do
       end do

       do i=lo(1),hi(1)
          do j=lo(2),hi(2)
             wmac(i,j,lo(3)-1) = TWO*wmac(i,j,lo(3)) - wmac(i,j,lo(3)+1)
             wmac(i,j,hi(3)+2) = TWO*wmac(i,j,hi(3)+1) - wmac(i,j,hi(3))
          end do
       end do

    end if

    if (do_trans) then

       ! use linear interpolation to fill fine level ghost cells needed for 
       ! transverse derivatives
       do i=lo(1),hi(1)+1
          do j=lo(2),hi(2)
             umac(i,j,lo(3)-1) = TWO*umac(i,j,lo(3)) - umac(i,j,lo(3)+1)
             umac(i,j,hi(3)+1) = TWO*umac(i,j,hi(3)) - umac(i,j,hi(3)-1)
          end do
       end do
       
       do i=lo(1),hi(1)+1
          do k=lo(3),hi(3)
             umac(i,lo(2)-1,k) = TWO*umac(i,lo(2),k) - umac(i,lo(2)+1,k)
             umac(i,hi(2)+1,k) = TWO*umac(i,hi(2),k) - umac(i,hi(2)-1,k)
          end do
       end do
       
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             vmac(i,j,lo(3)-1) = TWO*vmac(i,j,lo(3)) - vmac(i,j,lo(3)+1)
             vmac(i,j,hi(3)+1) = TWO*vmac(i,j,hi(3)) - vmac(i,j,hi(3)-1)
          end do
       end do
       
       do j=lo(2),hi(2)+1
          do k=lo(3),hi(3)
             vmac(lo(1)-1,j,k) = TWO*vmac(lo(1),j,k) - vmac(lo(1)+1,j,k)
             vmac(hi(1)+1,j,k) = TWO*vmac(hi(1),j,k) - vmac(hi(1)-1,j,k)
          end do
       end do
       
       do k=lo(3),hi(3)+1
          do i=lo(1),hi(1)
             wmac(i,lo(2)-1,k) = TWO*wmac(i,lo(2),k) - wmac(i,lo(2)+1,k)
             wmac(i,hi(2)+1,k) = TWO*wmac(i,hi(2),k) - wmac(i,hi(2)-1,k)
          end do
       end do
       
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             wmac(lo(1)-1,j,k) = TWO*wmac(lo(1),j,k) - wmac(lo(1)+1,j,k)
             wmac(hi(1)+1,j,k) = TWO*wmac(hi(1),j,k) - wmac(hi(1)-1,j,k)
          end do
       end do

    end if
       
  end subroutine create_umac_grown_onesided_3d

end module create_umac_grown_module
