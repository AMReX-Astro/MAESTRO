! extrap_to_halftime is used to extrapolate S to the half time for the
! Step 1 in the algorithm, when we don't yet have a S_cc_new

module extraphalf_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: extrap_to_halftime

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine extrap_to_halftime(mla,S_cc_nph,dSdt,S_cc_old,dt,the_bc_level)
    
    use variables, only: foextrap_comp
    use ml_restrict_fill_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: S_cc_nph(:), dSdt(:), S_cc_old(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    real(kind=dp_t), pointer:: Snphp(:,:,:,:)
    real(kind=dp_t), pointer:: Soldp(:,:,:,:)
    real(kind=dp_t), pointer:: dSdtp(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng_h,ng_o,ng_dS
    integer :: i,n,dm,nlevs
      
    dm = mla%dim
    nlevs = mla%nlevel

    ng_h  = nghost(S_cc_nph(1))
    ng_o  = nghost(S_cc_old(1))
    ng_dS = nghost(dSdt(1))

    do n = 1, nlevs

       do i = 1, nfabs(S_cc_nph(n))
          Snphp => dataptr(S_cc_nph(n), i)
          Soldp => dataptr(S_cc_old(n), i)
          dSdtp => dataptr(dSdt(n), i)
          
          lo =  lwb(get_box(S_cc_nph(n), i))
          hi =  upb(get_box(S_cc_nph(n), i))
          
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

    ! restrict data (S_cc_nph has no ghost cells)
    call ml_restrict_and_fill(nlevs,S_cc_nph,mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=S_cc_nph(1)%ng)

  end subroutine extrap_to_halftime

  subroutine extrap_to_halftime_1d(S_cc_nph,dSdt,S_cc_old, &
                                   dt,lo,hi,ng_h,ng_o,ng_dS)

    use bl_constants_module

    integer         , intent(in ) :: lo(:), hi(:), ng_h, ng_o, ng_dS
    real (kind=dp_t), intent(out) :: S_cc_nph(lo(1)-ng_h :)
    real (kind=dp_t), intent(in ) ::     dSdt(lo(1)-ng_dS:)
    real (kind=dp_t), intent(in ) :: S_cc_old(lo(1)-ng_o :)
    real (kind=dp_t) :: dt

    ! Local variables
    integer          :: i

    do i = lo(1),hi(1)
       S_cc_nph(i) = S_cc_old(i) + HALF*dt*dSdt(i)
    end do
 
  end subroutine extrap_to_halftime_1d

  subroutine extrap_to_halftime_2d(S_cc_nph,dSdt,S_cc_old, &
                                   dt,lo,hi,ng_h,ng_o,ng_dS)

    use bl_constants_module

    integer         , intent(in ) :: lo(:), hi(:), ng_h, ng_o, ng_dS
    real (kind=dp_t), intent(out) :: S_cc_nph(lo(1)-ng_h :,lo(2)-ng_h :)
    real (kind=dp_t), intent(in ) ::     dSdt(lo(1)-ng_dS:,lo(2)-ng_dS:)
    real (kind=dp_t), intent(in ) :: S_cc_old(lo(1)-ng_o :,lo(2)-ng_o :)
    real (kind=dp_t) :: dt

    ! Local variables
    integer          :: i, j

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          S_cc_nph(i,j) = S_cc_old(i,j) + HALF*dt*dSdt(i,j)
       end do
    end do
 
  end subroutine extrap_to_halftime_2d


  subroutine extrap_to_halftime_3d(S_cc_nph,dSdt,S_cc_old, &
                                   dt,lo,hi,ng_h,ng_o,ng_dS)
    use bl_constants_module

    integer         , intent(in ) :: lo(:), hi(:), ng_h, ng_o, ng_dS
    real (kind=dp_t), intent(out) :: S_cc_nph(lo(1)-ng_h :,lo(2)-ng_h :,lo(3)-ng_h :)
    real (kind=dp_t), intent(in ) ::     dSdt(lo(1)-ng_dS:,lo(2)-ng_dS:,lo(3)-ng_dS:)
    real (kind=dp_t), intent(in ) :: S_cc_old(lo(1)-ng_o :,lo(2)-ng_o :,lo(3)-ng_o :)
    real (kind=dp_t) :: dt

    ! Local variables
    integer          :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             S_cc_nph(i,j,k) = S_cc_old(i,j,k) + HALF*dt*dSdt(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
 
  end subroutine extrap_to_halftime_3d

end module extraphalf_module
