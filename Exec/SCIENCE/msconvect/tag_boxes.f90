module tag_boxes_module

  use multifab_module
  use bl_error_module

  implicit none 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MUST set this to .true. if tagging uses ghost cells (e.g., tagging on gradient). !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical, save :: tagging_needs_ghost_cells = .false.

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)
    
    use variables, only: rho_comp
    use probin_module, only: lo_dens_tag, hi_dens_tag, steep_tag
                                 
    type( multifab)         , intent(in   ) :: mf
    type(lmultifab)         , intent(inout) :: tagboxes
    real(dp_t)              , intent(in   ) :: dx
    integer                 , intent(in   ) :: lev
    type(multifab), optional, intent(in   ) :: aux_tag_mf
    ! aux_tag_mf allows user to pass in additional multifabs for tagging logic
    
    ! local variables
    real(kind = dp_t), pointer :: mfp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i
    integer           :: lo(get_dim(mf)), hi(get_dim(mf))
    integer           :: tlo(4), mflo(4)
    type(mfiter)      :: mfi
    type(box)         :: bx
    real(dp_t)        :: lo_tag, hi_tag
   ! if (present(aux_tag_mf)) do nothing. Not implemented
    
    
    
    !$omp parallel private(mfp,tp,i,lo,hi,lo_tag,hi_tag,mfi,bx,tlo,mflo)
    call mfiter_build(mfi,tagboxes,.true.)
    do while(next_tile(mfi,i))
       bx = get_tilebox(mfi)
       lo =  lwb(bx)
       hi =  upb(bx)
       
       lo_tag = (1. + (steep_tag ** (lev-1.))) * lo_dens_tag
       hi_tag = (1. - (steep_tag ** (lev-1.))) * hi_dens_tag
       
       mfp => dataptr(mf, i)
       tp  => dataptr(tagboxes, i)


       mflo = lbound(mfp)
       tlo = lbound(tp)

       select case (get_dim(mf))
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),tlo(1:2),mfp(:,:,1,rho_comp),mflo(1:2),lo,hi,dx,lev)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),tlo(1:3),mfp(:,:,:,rho_comp),mflo(1:3),lo,hi,lo_tag,hi_tag,dx,lev)
       end select
    end do
    !$omp end parallel

  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,tlo,mf,mflo,lo,hi,dx,lev)

    integer          , intent(in   ) :: lo(2),hi(2),tlo(2), mflo(2)
    logical          , intent(inout) :: tagbox( tlo(1):, tlo(2):)
    real(kind = dp_t), intent(in   ) ::     mf(mflo(1):,mflo(2):)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j

    ! initially say that we do not want to tag any cells for refinement
    tagbox(lo(1):hi(1),lo(2):hi(2)) = .false.

    select case(lev)
    case (1)
       ! tag all boxes where the first component of mf >= 1.01
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.01_dp_t) then
                tagbox(i,j) = .true.
             end if
          end do
       enddo
    case (2)
       ! for level 2 tag all boxes where the first component of mf >= 1.1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.1_dp_t) then
                tagbox(i,j) = .true.
             end if
          end do
       end do
    case default
       ! for level 3 or greater tag all boxes where the first component of mf >= 1.5
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.5_dp_t) then
                tagbox(i,j) = .true.
             end if
          end do
       end do
    end select

  end subroutine tag_boxes_2d

  subroutine tag_boxes_3d(tagbox,tlo,mf,mflo,lo,hi,lo_tag,hi_tag,dx,lev)


    integer          , intent(in   ) :: lo(3),hi(3),tlo(3),mflo(3)
    logical          , intent(inout) :: tagbox( tlo(1):, tlo(2):, tlo(3):)
    real(kind = dp_t), intent(in   ) ::     mf(mflo(1):,mflo(2):,mflo(3):)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev
    real(kind = dp_t), intent(in   ) :: lo_tag, hi_tag

    ! local variables
    integer :: i,j,k

    ! initially say that we do not want to tag any cells for refinement
    tagbox(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = .false.

    do k = lo(3),hi(3)
      do j = lo(2),hi(2)
	do i = lo(1),hi(1)
	  if (mf(i,j,k) .gt. lo_tag .and. mf(i,j,k) .lt. hi_tag) then
	  tagbox(i,j,k) = .true.
	  end if
	end do
      end do
    end do
    
  end subroutine tag_boxes_3d

end module tag_boxes_module
