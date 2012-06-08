! If flag = .true. then put the state into perturbational form.  That
! is, construct S' = S - S_0 for each of the desired state variables.
! Here, the base state, S_0, is left unchanged.
!
! If flag = .false. then convert from perturbational form to full 
! state form (S = S' + S_0)
!

module pert_form_module

  use multifab_module

  implicit none

  private

  public :: put_in_pert_form

  
contains

  subroutine put_in_pert_form(mla,s,base,dx,comp,bc_comp,flag,the_bc_level)

    use geometry, only: spherical
    use ml_layout_module
    use define_bc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module
    use multifab_physbc_module
    use bl_prof_module

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: comp,bc_comp
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: base(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    logical        , intent(in   ) :: flag
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: i,ng,n,dm,nlevs
    real(kind=dp_t) :: dx_temp(mla%nlevel,mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_in_pert_form")

    dm = mla%dim
    nlevs = mla%nlevel

    ng = nghost(s(1))

    ! the reason I go through this dx_temp thing is because when you
    ! restart into finer, we haven't initialized dx when we add the base state
    ! back to the perturbational state.  This should be fixed.
    dx_temp(1,:) = dx(1,:)
    do n=2,nlevs
       dx_temp(n,:) = dx_temp(n-1,:)/2.d0
    end do

    do n=1,nlevs
       do i = 1, nboxes(s(n))
          if ( multifab_remote(s(n),i) ) cycle
          sp => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call pert_form_2d(sp(:,:,1,:),base(n,:),lo,hi,ng,comp,flag)
          case (3)
             if (spherical .eq. 1) then
                call pert_form_3d_sphr(sp(:,:,:,:),base(1,:),lo,hi,ng,dx_temp(n,:),comp,flag)
             else
                call pert_form_3d_cart(sp(:,:,:,:),base(n,:),lo,hi,ng,comp,flag)
             end if
          end select
       end do
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s(nlevs),comp,1)
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),comp,bc_comp,1,the_bc_level(nlevs))
       
    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(s(n-1),comp,s(n),comp,mla%mba%rr(n-1,:),1)
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1), &
                                         nghost(s(n)),mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         comp,bc_comp,1,fill_crse_input=.false.)
       end do

    end if

    call destroy(bpt)

  end subroutine put_in_pert_form

  subroutine pert_form_2d(s,s0,lo,hi,ng,comp,flag)

    integer        , intent(in   ) ::  lo(:),hi(:),ng,comp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,mult

    if (flag) then
       mult = -1
    else
       mult = +1
    end if

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          s(i,j,comp) = s(i,j,comp) + mult*s0(j)
       end do
    end do

  end subroutine pert_form_2d

  subroutine pert_form_3d_cart(s,s0,lo,hi,ng,comp,flag)

    integer        , intent(in   ) ::  lo(:),hi(:),ng,comp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,k,mult

    if (flag) then
       mult = -1
    else
       mult = +1
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             s(i,j,k,comp) = s(i,j,k,comp) + mult*s0(k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine pert_form_3d_cart

  subroutine pert_form_3d_sphr(s,s0,lo,hi,ng,dx,comp,flag)

    use fill_3d_module

    integer        , intent(in   ) :: lo(:),hi(:),ng,comp
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    logical        , intent(in   ) :: flag

    ! local
    integer :: i,j,k,mult

    real(kind=dp_t), allocatable :: s0_cart(:,:,:,:)

    allocate(s0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0,s0_cart,lo,hi,dx,0)

    if (flag) then
       mult = -1
    else
       mult = +1
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             s(i,j,k,comp) = s(i,j,k,comp) + mult*s0_cart(i,j,k,1)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(s0_cart)

  end subroutine pert_form_3d_sphr

end module pert_form_module
