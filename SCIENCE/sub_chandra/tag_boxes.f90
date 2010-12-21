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

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    use variables, ONLY: rho_comp, spec_comp

    use network, ONLY: network_species_index

    type( multifab)          , intent(in   ) :: mf
    type(lmultifab)          , intent(inout) :: tagboxes
    real(dp_t)               , intent(in   ) :: dx
    integer                  , intent(in   ) :: lev
    type( multifab), optional, intent(in   ) :: aux_tag_mf

    real(kind = dp_t), pointer :: sp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(mf%dim), ng

    integer           :: ihe4

    ihe4 = network_species_index("helium-4")

    ng = mf%ng

    do i = 1, mf%nboxes
       if ( multifab_remote(mf, i) ) cycle
       sp => dataptr(mf, i)
       tp => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))

       select case (mf%dim)
       case (2)
          call bl_error("ERROR: 2-d tag_boxes not implemented")
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),sp(:,:,:,rho_comp),sp(:,:,:,spec_comp-1+ihe4),lo,ng,lev)
       end select
    end do

  end subroutine tag_boxes


  subroutine tag_boxes_3d(tagbox,rho,rho_Xhe,lo,ng,lev)

    use bl_constants_module, ONLY: HALF
    use probin_module, ONLY: base_cutoff_density, prob_lo, prob_hi

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    real(kind = dp_t), intent(in   ) :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind = dp_t), intent(in   ) :: rho_Xhe(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer, optional, intent(in   ) :: lev

    real(kind = dp_t) :: Xhe, x, y, z

    integer :: i,j,k,nx,ny,nz,llev

    llev = 1; if (present(lev)) llev = lev
    nx = size(tagbox,dim=1)
    ny = size(tagbox,dim=2)
    nz = size(tagbox,dim=3)

    tagbox = .false.

!$omp parallel do private(i,j,k,Xhe)
    do k = lo(3),lo(3)+nz-1
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             Xhe = rho_Xhe(i,j,k)/rho(i,j,k)
             if (Xhe > 0.01_dp_t .and. rho(i,j,k) >= base_cutoff_density) then
                tagbox(i,j,k) = .true.
             end if
          end do
       enddo
    end do
!$omp end parallel do

    ! refine the very center of the star, for average
!$omp parallel do private(i,j,k,x,y,z)
    do k = lo(3),lo(3)+nz-1
       !z = prob_lo(3) + (dble(k)+HALF)*dx(3)

       do j = lo(2),lo(2)+ny-1
          !y = prob_lo(2) + (dble(j)+HALF)*dx(2)

          do i = lo(1),lo(1)+nx-1
             !x = prob_lo(1) + (dble(i)+HALF)*dx(1)

             !if (sqrt(x**2 + y**2 + z**2) < 0.1*prob_hi(1)) then
             if (sqrt(dble(i)**2 + dble(j)**2 + dble(k)**2) < 10.d0) then
                tagbox(i,j,k) = .true.
             end if
          end do
       enddo
    end do
!$omp end parallel do

  end subroutine tag_boxes_3d

end module tag_boxes_module
