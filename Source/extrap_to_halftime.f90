module extraphalf_module

  use bl_types
  use bl_constants_module
  use variables
  use multifab_module

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine extrap_to_halftime(Source_nph,Source_nm1,Source_old, &
                                dtold,dt)

    type(multifab) , intent(inout) :: Source_nph
    type(multifab) , intent(in   ) :: Source_nm1, Source_old
    real(kind=dp_t), intent(in   ) :: dt, dtold

    real(kind=dp_t), pointer:: Snphp(:,:,:,:)
    real(kind=dp_t), pointer:: Soldp(:,:,:,:)
    real(kind=dp_t), pointer:: Snm1p(:,:,:,:)
    integer :: lo(Source_nph%dim),hi(Source_nph%dim),ng_h,ng_o,dm
    integer :: i,bc_comp
      
    dm = Source_nph%dim
    ng_h = Source_nph%ng
    ng_o = Source_old%ng

    do i = 1, Source_nph%nboxes
       if ( multifab_remote(Source_nph, i) ) cycle
       Snphp => dataptr(Source_nph, i)
       Soldp => dataptr(Source_old, i)
       Snm1p => dataptr(Source_nm1, i)

       lo =  lwb(get_box(Source_nph, i))
       hi =  upb(get_box(Source_nph, i))

       select case (dm)
       case (2)
          call extrap_to_halftime_2d(Snphp(:,:,1,1), &
                                     Snm1p(:,:,1,1), Soldp(:,:,1,1),&
                                     dt, dtold, &
                                     lo,hi,ng_h,ng_o)

       case (3)
          call extrap_to_halftime_3d(Snphp(:,:,:,1), &
                                     Snm1p(:,:,:,1), Soldp(:,:,:,1),&
                                     dt, dtold, &
                                     lo,hi,ng_h,ng_o)

       end select
    end do

  end subroutine extrap_to_halftime


  subroutine extrap_to_halftime_2d(Source_nph, &
                                   Source_nm1, Source_old, &
                                   dt, dtold, &
                                   lo,hi,ng_h,ng_o)

    implicit none

    integer         , intent(in ) :: lo(:), hi(:), ng_h, ng_o
    real (kind=dp_t), intent(out) :: Source_nph(lo(1)-ng_h:,lo(2)-ng_h:)
    real (kind=dp_t), intent(in ) :: Source_nm1(lo(1)-ng_o:,lo(2)-ng_o:)
    real (kind=dp_t), intent(in ) :: Source_old(lo(1)-ng_o:,lo(2)-ng_o:)
    real (kind=dp_t) :: dt, dtold

    ! Local variables
    integer          :: i, j

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          Source_nph(i,j) = Source_old(i,j) + &
               HALF * dt * (Source_old(i,j) - Source_nm1(i,j))/dtold
       end do
    end do
 
  end subroutine extrap_to_halftime_2d


  subroutine extrap_to_halftime_3d(Source_nph, &
                                 Source_nm1, Source_old, &
                                 dt, dtold, &
                                 lo,hi,ng_h,ng_o)

    implicit none

    integer         , intent(in ) :: lo(:), hi(:), ng_h, ng_o
    real (kind=dp_t), intent(out) :: Source_nph(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real (kind=dp_t), intent(in ) :: Source_nm1(lo(1)-ng_o:,lo(2)-ng_o:,lo(3)-ng_o:)
    real (kind=dp_t), intent(in ) :: Source_old(lo(1)-ng_o:,lo(2)-ng_o:,lo(3)-ng_o:)
    real (kind=dp_t) :: dt, dtold

    ! Local variables
    integer          :: i, j, k

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             Source_nph(i,j,k) = Source_old(i,j,k) + &
                  HALF * dt * (Source_old(i,j,k) - Source_nm1(i,j,k))/dtold
          end do
       end do
    end do
 
  end subroutine extrap_to_halftime_3d

end module extraphalf_module
