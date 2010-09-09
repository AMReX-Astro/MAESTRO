module phihalf_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: make_S_at_halftime, make_at_halftime

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_S_at_halftime(mla,shalf,sold,snew,the_bc_level)

    use bl_prof_module
    use multifab_physbc_module
    use ml_restriction_module, only : ml_cc_restriction
    use multifab_fill_ghost_module
    use variables, only: foextrap_comp
    use geometry, only: dm, nlevs

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: shalf(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(in   ) :: snew(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    real(kind=dp_t), pointer:: shp(:,:,:,:)
    real(kind=dp_t), pointer:: sop(:,:,:,:)
    real(kind=dp_t), pointer:: snp(:,:,:,:)
    integer :: lo(dm),hi(dm),ng_h,ng_o,ng_n
    integer :: i,in_comp,out_comp,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_S_at_halftime")

    ng_h = nghost(shalf(1))
    ng_o = nghost(sold(1))
    ng_n = nghost(snew(1))

    in_comp = 1
    out_comp = 1

    do n = 1, nlevs

       do i = 1, nboxes(shalf(n))
          if ( multifab_remote(shalf(n), i) ) cycle
          shp => dataptr(shalf(n), i)
          sop => dataptr(sold(n), i)
          snp => dataptr(snew(n), i)
          lo =  lwb(get_box(shalf(n), i))
          hi =  upb(get_box(shalf(n), i))
          select case (dm)
          case (1)
             call make_at_halftime_1d(shp(:,1,1,out_comp),sop(:,1,1,in_comp), &
                                      snp(:,1,1,in_comp),lo,hi,ng_h,ng_o,ng_n)
          case (2)
             call make_at_halftime_2d(shp(:,:,1,out_comp),sop(:,:,1,in_comp), &
                                      snp(:,:,1,in_comp),lo,hi,ng_h,ng_o,ng_n)
          case (3)
             call make_at_halftime_3d(shp(:,:,:,out_comp),sop(:,:,:,in_comp), &
                                      snp(:,:,:,in_comp),lo,hi,ng_h,ng_o,ng_n)
          end select
       end do

    end do

    ! fill the ghostcells on the new time-centered S
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(shalf(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(shalf(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(shalf(n-1)    ,shalf(n)    ,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(shalf(n),shalf(n-1), &
                                         ng_h,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         1,foextrap_comp,1,fill_crse_input=.false.)
       enddo

    end if


    call destroy(bpt)

  end subroutine make_S_at_halftime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime(phihalf,sold,snew,in_comp,out_comp,the_bc_level,mla)

    use multifab_physbc_module
    use ml_restriction_module, only : ml_cc_restriction
    use multifab_fill_ghost_module
    use geometry, only: dm, nlevs

    type(multifab) , intent(inout) :: phihalf(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(in   ) :: snew(:)
    integer        , intent(in   ) :: in_comp,out_comp
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: rhp(:,:,:,:)
    real(kind=dp_t), pointer:: rop(:,:,:,:)
    real(kind=dp_t), pointer:: rnp(:,:,:,:)
    integer   :: lo(dm),hi(dm)
    integer   :: ng_h,ng_o,ng_n,i,n

    ng_h = nghost(phihalf(1))
    ng_o = nghost(sold(1))
    ng_n = nghost(snew(1))

    do n = 1, nlevs
       do i = 1, nboxes(phihalf(n))
          if ( multifab_remote(phihalf(n), i) ) cycle
          rhp => dataptr(phihalf(n), i)
          rop => dataptr(sold(n), i)
          rnp => dataptr(snew(n), i)
          lo =  lwb(get_box(phihalf(n), i))
          hi =  upb(get_box(phihalf(n), i))
          select case (dm)
          case (1)
             call make_at_halftime_1d(rhp(:,1,1,out_comp),rop(:,1,1,in_comp), &
                                      rnp(:,1,1,in_comp),lo,hi,ng_h,ng_o,ng_n)
          case (2)
             call make_at_halftime_2d(rhp(:,:,1,out_comp),rop(:,:,1,in_comp), &
                                      rnp(:,:,1,in_comp),lo,hi,ng_h,ng_o,ng_n)
          case (3)
             call make_at_halftime_3d(rhp(:,:,:,out_comp),rop(:,:,:,in_comp), &
                                      rnp(:,:,:,in_comp),lo,hi,ng_h,ng_o,ng_n)
          end select
       end do
    end do

    if (nlevs .eq. 1) then
       
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(phihalf(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(phihalf(nlevs),out_comp,dm+in_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(phihalf(n-1),phihalf(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(phihalf(n),phihalf(n-1), &
                                         ng_h,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n  ), &
                                         1,dm+in_comp,1,fill_crse_input=.false.)
       end do

    end if

  end subroutine make_at_halftime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime_1d(phihalf,phiold,phinew,lo,hi,ng_half,ng_old,ng_new)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:),hi(:),ng_half,ng_old,ng_new
    real (kind=dp_t), intent(  out) :: phihalf(lo(1)-ng_half:)
    real (kind=dp_t), intent(in   ) ::  phiold(lo(1)-ng_old :)
    real (kind=dp_t), intent(in   ) ::  phinew(lo(1)-ng_new :)

    !  Local variables
    integer :: i

    do i = lo(1),hi(1)
       phihalf(i) = HALF * (phiold(i) + phinew(i))
    end do

  end subroutine make_at_halftime_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime_2d(phihalf,phiold,phinew,lo,hi,ng_half,ng_old,ng_new)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:),hi(:),ng_half,ng_old,ng_new
    real (kind=dp_t), intent(  out) :: phihalf(lo(1)-ng_half:,lo(2)-ng_half:)
    real (kind=dp_t), intent(in   ) ::  phiold(lo(1)-ng_old :,lo(2)-ng_old :)
    real (kind=dp_t), intent(in   ) ::  phinew(lo(1)-ng_new :,lo(2)-ng_new :)

    !  Local variables
    integer :: i, j

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          phihalf(i,j) = HALF * (phiold(i,j) + phinew(i,j))
       end do
    end do

  end subroutine make_at_halftime_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime_3d(phihalf,phiold,phinew,lo,hi,ng_half,ng_old,ng_new)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:),hi(:),ng_half,ng_old,ng_new
    real (kind=dp_t), intent(  out) :: phihalf(lo(1)-ng_half:,lo(2)-ng_half:,lo(3)-ng_half:)
    real (kind=dp_t), intent(in   ) ::  phiold(lo(1)-ng_old :,lo(2)-ng_old :,lo(3)-ng_old :)
    real (kind=dp_t), intent(in   ) ::  phinew(lo(1)-ng_new :,lo(2)-ng_new :,lo(3)-ng_new :)

    ! Local variables
    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             phihalf(i,j,k) = HALF * (phiold(i,j,k) + phinew(i,j,k))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine make_at_halftime_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module phihalf_module
