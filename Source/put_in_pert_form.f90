module pert_form_module

  ! Put the state into perturbational form

  use multifab_module

  implicit none

  private

  public :: put_in_pert_form

  
contains

  subroutine put_in_pert_form(nlevs,s,base,dx,startcomp,numcomp,flag,mla,the_bc_level)

    use geometry, only: spherical
    use variables, only: foextrap_comp, nscal, is_pert_form
    use ml_layout_module
    use define_bc_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module

    integer        , intent(in   ) :: nlevs,startcomp,numcomp
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: base(:,0:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    logical        , intent(in   ) :: flag
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    integer :: i,lo(s(1)%dim),hi(s(1)%dim),ng,dm,n,comp

    ng = s(1)%ng
    dm = s(1)%dim

    ! sanity checks -- if we are attempting to go into pert form and any of those
    ! variables are already in pert form, abort.  Likewise, if we are going out of
    ! pert form, make sure that all of the variable are in pert form to begin with.
    if ( (      flag .and.       any(is_pert_form(startcomp:startcomp+numcomp-1)) ) .or. &
         (.not. flag .and. .not. all(is_pert_form(startcomp:startcomp+numcomp-1)) ) ) then
       call bl_error('Illegal call to put_in_pert_form')
    endif

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sp => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call pert_form_2d(sp(:,:,1,:),base(n,:,:),lo,hi,ng,startcomp,numcomp,flag)
          case (3)
             if (spherical .eq. 1) then
                call pert_form_3d_sphr(n,sp(:,:,:,:),base(n,:,:),lo,hi,ng,dx(n,:), &
                                       startcomp,numcomp,flag)
             else
                call pert_form_3d_cart(sp(:,:,:,:),base(n,:,:),lo,hi,ng,startcomp,numcomp,flag)
             end if
          end select
       end do

       call multifab_fill_boundary_c(s(n),startcomp,numcomp)

    end do

    ! is_pert_form is .true. is the state is now in perturbational form
    if (flag) then
       is_pert_form(startcomp:startcomp+numcomp-1) = .true.
    else
       is_pert_form(startcomp:startcomp+numcomp-1) = .false.
    endif


    do n=nlevs,2,-1
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))

       if(flag) then
          do comp=startcomp,startcomp+numcomp-1
             call multifab_fill_ghost_cells(s(n),s(n-1), &
                                            s(n)%ng,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1),the_bc_level(n), &
                                            comp,foextrap_comp,1)
          end do
       else
          call multifab_fill_ghost_cells(s(n),s(n-1), &
                                         s(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         startcomp,dm+startcomp,numcomp)
       end if

    end do

  end subroutine put_in_pert_form

  subroutine pert_form_2d(s,s0,lo,hi,ng,startcomp,numcomp,flag)

    integer        , intent(in   ) ::  lo(:),hi(:),ng,startcomp,numcomp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,comp

    if (flag) then
       do comp = startcomp, startcomp+numcomp-1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                s(i,j,comp) = s(i,j,comp) - s0(j,comp)
             end do
          end do
       end do
    else
       do comp = startcomp, startcomp+numcomp-1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                s(i,j,comp) = s(i,j,comp) + s0(j,comp)
             end do
          end do
       end do
    end if

  end subroutine pert_form_2d

  subroutine pert_form_3d_cart(s,s0,lo,hi,ng,startcomp,numcomp,flag)

    integer        , intent(in   ) ::  lo(:),hi(:),ng,startcomp,numcomp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,k,comp

    if (flag) then
       do comp = startcomp, startcomp+numcomp-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp) - s0(k,comp)
                end do
             end do
          end do
       end do
    else
       do comp = startcomp, startcomp+numcomp-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp) + s0(k,comp)
                end do
             end do
          end do
       end do
    end if

  end subroutine pert_form_3d_cart

  subroutine pert_form_3d_sphr(n,s,s0,lo,hi,ng,dx,startcomp,numcomp,flag)

    use fill_3d_module

    integer        , intent(in   ) :: n,lo(:),hi(:),ng,startcomp,numcomp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    logical        , intent(in   ) :: flag

    real(kind=dp_t), allocatable :: s0_cart(:,:,:)
    integer                      :: i,j,k,comp

    allocate(s0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

    if (flag) then
       do comp = startcomp,startcomp+numcomp-1
          call fill_3d_data(n,s0_cart,s0(0:,comp),lo,hi,dx,0)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp) - s0_cart(i,j,k)
                end do
             end do
          end do
       end do
    else
       do comp = startcomp,startcomp+numcomp-1
          call fill_3d_data(n,s0_cart,s0(0:,comp),lo,hi,dx,0)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp) + s0_cart(i,j,k)
                end do
             end do
          end do
       end do
    end if

    deallocate(s0_cart)

  end subroutine pert_form_3d_sphr

end module pert_form_module
