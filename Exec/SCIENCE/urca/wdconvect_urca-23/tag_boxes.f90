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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MUST set this to .true. if tagging uses ghost cells (e.g., tagging on gradient). !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical, save :: tagging_needs_ghost_cells = .false.

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    use variables, ONLY: rho_comp

    type( multifab)          , intent(in   ) :: mf
    type(lmultifab)          , intent(inout) :: tagboxes
    real(dp_t)               , intent(in   ) :: dx
    integer                  , intent(in   ) :: lev
    type( multifab), optional, intent(in   ) :: aux_tag_mf

    real(kind = dp_t), pointer :: sp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(get_dim(mf)), ng

    ng = nghost(mf)

    do i = 1, nfabs(mf)
       sp => dataptr(mf, i)
       tp => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))

       select case (get_dim(mf))
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),sp(:,:,1,rho_comp),lo,ng,lev)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),sp(:,:,:,:),lo,ng,lev)
       end select
    end do

  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,mf,lo,ng,lev)

    use probin_module, ONLY : base_cutoff_density, &
         tag_density_1, tag_density_2, tag_density_3

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:)
    integer, optional, intent(in   ) :: lev
    integer :: i,j,nx,ny,llev

    llev = 1; if (present(lev)) llev = lev
    nx = size(tagbox,dim=1)
    ny = size(tagbox,dim=2)

    tagbox = .false.

    select case(llev)
    case (1)
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             if (mf(i,j) .gt. tag_density_1) then
                tagbox(i,j) = .true.
             end if
          end do
       enddo
    case (2)
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             if (mf(i,j) .gt. tag_density_2) then
                tagbox(i,j) = .true.
             end if
          end do
       end do
    case (3)
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             if (mf(i,j) .gt. tag_density_3) then
                tagbox(i,j) = .true.
             end if
          end do
       end do
    end select

  end subroutine tag_boxes_2d

  subroutine tag_boxes_3d(tagbox,mf,lo,ng,lev)

    use variables, ONLY: rho_comp, spec_comp, temp_comp
    use network, ONLY: get_electron_fraction, nspec
    use probin_module, ONLY : base_cutoff_density, &
         tag_density_1, &
         tag_rhoye_lo_2, &
         tag_rhoye_lo_3

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    integer, optional, intent(in   ) :: lev

    real(kind = dp_t) :: ye, rhoye

    integer :: i,j,k,nx,ny,nz,llev

    llev = 1; if (present(lev)) llev = lev
    nx = size(tagbox,dim=1)
    ny = size(tagbox,dim=2)
    nz = size(tagbox,dim=3)

    tagbox = .false.

    select case(llev)
    case (1)
!$omp parallel do private(i,j,k)
       do k = lo(3),lo(3)+nz-1
          do j = lo(2),lo(2)+ny-1
             do i = lo(1),lo(1)+nx-1
                if (mf(i,j,k, rho_comp) .gt. tag_density_1) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          enddo
       end do
!$omp end parallel do
    case (2)
!$omp parallel do private(i,j,k)
       do k = lo(3),lo(3)+nz-1
          do j = lo(2),lo(2)+ny-1
             do i = lo(1),lo(1)+nx-1
                call get_electron_fraction(ye, mf(i, j, k, spec_comp:spec_comp + nspec - 1))
                rhoye = mf(i, j, k, rho_comp) * ye
                if (rhoye .gt. tag_rhoye_lo_2) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          end do
       end do
!$omp end parallel do
    case (3)
!$omp parallel do private(i,j,k)
       do k = lo(3),lo(3)+nz-1
          do j = lo(2),lo(2)+ny-1
             do i = lo(1),lo(1)+nx-1
                call get_electron_fraction(ye, mf(i, j, k, spec_comp:spec_comp + nspec - 1))
                rhoye = mf(i, j, k, rho_comp) * ye
                if (rhoye .gt. tag_rhoye_lo_3) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          end do
       end do
!$omp end parallel do
    end select

  end subroutine tag_boxes_3d

end module tag_boxes_module
