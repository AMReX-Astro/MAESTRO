module fill_force_module

  use bl_types
  use bc_module
  use multifab_module
  use make_grav_module
  use network
  use variables

  implicit none
  
contains

  subroutine fill_force_with_grav (force,rho0,dr,spherical)

    type(multifab) , intent(inout) :: force
    real(kind=dp_t), intent(in   ) ::  rho0(:)
    real(kind=dp_t), intent(in   ) ::  dr
    integer        , intent(in   ) :: spherical

    real(kind=dp_t), allocatable:: grav_cell(:)
    real(kind=dp_t), pointer:: fp(:,:,:,:)
    integer :: lo(force%dim),hi(force%dim)
    integer :: i,nz,ng,dm

    ng = force%ng
    dm = force%dim

    nz = size(rho0,dim=1)
    allocate(grav_cell(nz))
    call make_grav_cell(grav_cell,rho0,dr,spherical)

    do i = 1, force%nboxes
       if ( multifab_remote(force, i) ) cycle
       fp => dataptr(force, i)
       lo =  lwb(get_box(force, i))
       hi =  upb(get_box(force, i))
       select case (dm)
       case (2)
          call fill_force_2d(fp(:,:,1,dm), grav_cell, lo, hi, ng)
       case (3)
          call fill_force_3d(fp(:,:,:,dm), grav_cell, lo, hi, ng)
       end select
    end do
    deallocate(grav_cell)

  end subroutine fill_force_with_grav

  subroutine fill_force_2d(force,grav,lo,hi,ng)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: force(lo(1)-ng:,lo(2)-ng:)
    real (kind = dp_t), intent(in   ) ::  grav(lo(2):)

    ! Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          force(i,j) = grav(j)
       enddo
    enddo

  end subroutine fill_force_2d

  subroutine fill_force_3d(force,grav,lo,hi,ng)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: force(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in   ) ::  grav(lo(3):)

    ! Local variables
    integer :: i, j, k

    do k = lo(3), hi(3)
     do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          force(i,j,k) = grav(k)
       enddo
     enddo
    enddo

  end subroutine fill_force_3d

end module fill_force_module
