! extrap_to_halftime is used to extrapolate S to the half time for the very
! first step, when we don't yet have a Source_new

module extraphalf_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: extrap_to_halftime

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine extrap_to_halftime(nlevs,mla,Source_nph,dSdt,Source_old,dt,the_bc_level)
    
    use variables, only: foextrap_comp
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: Source_nph(:), dSdt(:), Source_old(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    real(kind=dp_t), pointer:: Snphp(:,:,:,:)
    real(kind=dp_t), pointer:: Soldp(:,:,:,:)
    real(kind=dp_t), pointer:: dSdtp(:,:,:,:)
    integer :: lo(Source_nph(1)%dim),hi(Source_nph(1)%dim),ng_h,ng_o,ng_dS,dm
    integer :: i,n
      
    dm = Source_nph(1)%dim
    ng_h = Source_nph(1)%ng
    ng_o = Source_old(1)%ng
    ng_dS = dSdt(1)%ng

    do n = 1, nlevs

       do i = 1, Source_nph(n)%nboxes
          if ( multifab_remote(Source_nph(n), i) ) cycle
          Snphp => dataptr(Source_nph(n), i)
          Soldp => dataptr(Source_old(n), i)
          dSdtp => dataptr(dSdt(n), i)
          
          lo =  lwb(get_box(Source_nph(n), i))
          hi =  upb(get_box(Source_nph(n), i))
          
          select case (dm)
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

    ! fill the ghostcells on the new time-centered S
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(Source_nph(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(Source_nph(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(Source_nph(n-1)    ,Source_nph(n)    ,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(Source_nph(n),Source_nph(n-1), &
                                         ng_h,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         1,foextrap_comp,1)
       enddo

    end if


  end subroutine extrap_to_halftime


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

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             Source_nph(i,j,k) = Source_old(i,j,k) + HALF*dt*dSdt(i,j,k)
          end do
       end do
    end do
 
  end subroutine extrap_to_halftime_3d

end module extraphalf_module
