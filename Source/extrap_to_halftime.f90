! extrap_to_halftime is used to extrapolate S to the half time for the
! Step 1 in the algorithm, when we don't yet have a Source_new

module extraphalf_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: extrap_to_halftime

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine extrap_to_halftime(mla,Source_nph,dSdt,Source_old,dt,the_bc_level)
    
    use variables, only: foextrap_comp
    use ml_restrict_fill_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: Source_nph(:), dSdt(:), Source_old(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    real(kind=dp_t), pointer:: Snphp(:,:,:,:)
    real(kind=dp_t), pointer:: Soldp(:,:,:,:)
    real(kind=dp_t), pointer:: dSdtp(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng_h,ng_o,ng_dS
    integer :: i,n,dm,nlevs
      
    dm = mla%dim
    nlevs = mla%nlevel

    ng_h  = nghost(Source_nph(1))
    ng_o  = nghost(Source_old(1))
    ng_dS = nghost(dSdt(1))

    do n = 1, nlevs

       do i = 1, nfabs(Source_nph(n))
          Snphp => dataptr(Source_nph(n), i)
          Soldp => dataptr(Source_old(n), i)
          dSdtp => dataptr(dSdt(n), i)
          
          lo =  lwb(get_box(Source_nph(n), i))
          hi =  upb(get_box(Source_nph(n), i))
          
          select case (dm)
          case (1)
             call extrap_to_halftime_1d(Snphp(:,1,1,1), &
                                        dSdtp(:,1,1,1), &
                                        Soldp(:,1,1,1), &
                                        dt,lo,hi,ng_h,ng_o,ng_dS)
          case (2)
             call extrap_to_halftime_2d(Snphp(:,:,1,1), &
                                        dSdtp(:,:,1,1), &
                                        Soldp(:,:,1,1), &
                                        dt,lo,hi,ng_h,ng_o,ng_dS)
          case (3)
             call extrap_to_halftime_3d(Snphp(:,:,:,1), &
                                        dSdtp(:,:,:,1), &
                                        Soldp(:,:,:,1), &
                                        dt,lo,hi,ng_h,ng_o,ng_dS)
          end select
       end do
       
    enddo

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,Source_nph,mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=Source_nph(1)%ng)

  end subroutine extrap_to_halftime

  subroutine extrap_to_halftime_1d(Source_nph,dSdt,Source_old, &
                                   dt,lo,hi,ng_h,ng_o,ng_dS)

    use bl_constants_module

    integer         , intent(in ) :: lo(:), hi(:), ng_h, ng_o, ng_dS
    real (kind=dp_t), intent(out) :: Source_nph(lo(1)-ng_h :)
    real (kind=dp_t), intent(in ) ::       dSdt(lo(1)-ng_dS:)
    real (kind=dp_t), intent(in ) :: Source_old(lo(1)-ng_o :)
    real (kind=dp_t) :: dt

    ! Local variables
    integer          :: i

    do i = lo(1),hi(1)
       Source_nph(i) = Source_old(i) + HALF*dt*dSdt(i)
    end do
 
  end subroutine extrap_to_halftime_1d

  subroutine extrap_to_halftime_2d(Source_nph,dSdt,Source_old, &
                                   dt,lo,hi,ng_h,ng_o,ng_dS)

    use bl_constants_module

    integer         , intent(in ) :: lo(:), hi(:), ng_h, ng_o, ng_dS
    real (kind=dp_t), intent(out) :: Source_nph(lo(1)-ng_h :,lo(2)-ng_h :)
    real (kind=dp_t), intent(in ) ::       dSdt(lo(1)-ng_dS:,lo(2)-ng_dS:)
    real (kind=dp_t), intent(in ) :: Source_old(lo(1)-ng_o :,lo(2)-ng_o :)
    real (kind=dp_t) :: dt

    ! Local variables
    integer          :: i, j

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          Source_nph(i,j) = Source_old(i,j) + HALF*dt*dSdt(i,j)
       end do
    end do
 
  end subroutine extrap_to_halftime_2d


  subroutine extrap_to_halftime_3d(Source_nph,dSdt,Source_old, &
                                   dt,lo,hi,ng_h,ng_o,ng_dS)
    use bl_constants_module

    integer         , intent(in ) :: lo(:), hi(:), ng_h, ng_o, ng_dS
    real (kind=dp_t), intent(out) :: Source_nph(lo(1)-ng_h :,lo(2)-ng_h :,lo(3)-ng_h :)
    real (kind=dp_t), intent(in ) ::       dSdt(lo(1)-ng_dS:,lo(2)-ng_dS:,lo(3)-ng_dS:)
    real (kind=dp_t), intent(in ) :: Source_old(lo(1)-ng_o :,lo(2)-ng_o :,lo(3)-ng_o :)
    real (kind=dp_t) :: dt

    ! Local variables
    integer          :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             Source_nph(i,j,k) = Source_old(i,j,k) + HALF*dt*dSdt(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
 
  end subroutine extrap_to_halftime_3d

end module extraphalf_module
