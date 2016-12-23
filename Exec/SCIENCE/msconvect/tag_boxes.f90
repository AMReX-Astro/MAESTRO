module tag_boxes_module

  use multifab_module
  use bl_error_module
  use bl_constants_module

  implicit none 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MUST set this to .true. if tagging uses ghost cells (e.g., tagging on gradient). !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical, save :: tagging_needs_ghost_cells = .true.

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    use bl_constants_module
    use variables, only: spec_comp, rho_comp
    use probin_module, only: radiative_X, X_grad_min
                                 
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
   ! if (present(aux_tag_mf)) do nothing. Not implemented
    
    
    call mfiter_build(mfi,tagboxes,.true.)
    

    do while(next_tile(mfi,i))
       bx = get_tilebox(mfi)
       lo =  lwb(bx)
       hi =  upb(bx)
       mfp => dataptr(mf, i)
       tp  => dataptr(tagboxes, i)


       mflo = lbound(mfp)
       tlo = lbound(tp)

       select case (get_dim(mf))
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),tlo(1:2),mfp(:,:,1,spec_comp),mfp(:,:,1,rho_comp),mflo(1:2),lo,hi,radiative_X,X_grad_min,dx,lev)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),tlo(1:3),mfp(:,:,:,spec_comp),mfp(:,:,:,rho_comp),mflo(1:3),lo,hi,radiative_X,X_grad_min,dx,lev)
       end select
    end do
    
  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,tlo,mf_spec,mf_rho,mflo,lo,hi,radiative_X,X_grad_min,dx,lev)

    integer          , intent(in   ) :: lo(2),hi(2),tlo(2), mflo(2)
    logical          , intent(inout) :: tagbox( tlo(1):, tlo(2):)
    real(kind = dp_t), intent(in   ) ::     mf_spec(mflo(1):,mflo(2):)
    real(kind = dp_t), intent(in   ) ::     mf_rho(mflo(1):,mflo(2):)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev
    real(kind = dp_t), intent(in   ) :: radiative_X, X_grad_min

    ! local variables
    integer :: i,j
    real(kind = dp_t)  :: dspec

    ! initially say that we do not want to tag any cells for refinement
    tagbox(lo(1):hi(1),lo(2):hi(2)) = .false.

    select case (lev)
    case(1)
    !$OMP PARALLEL DO PRIVATE(i,j)  
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
	  if (mf_spec(i,j)/mf_rho(i,j) .lt. radiative_X) then
	    tagbox(i,j) = .true.
	  endif 
      enddo
    enddo
    !$OMP END PARALLEL DO
    case(2)

    !$OMP PARALLEL DO PRIVATE(i,j,dspec)
    do j = lo(2), hi(2)
	do i = lo(1), hi(1)
	    
	    dspec = ZERO             
	    
	    if (j == lo(2)) then
	      dspec = dspec + abs(mf_spec(i,j+1)/mf_rho(i,j+1) - mf_spec(i,j)/mf_rho(i,j))
	      ! backward difference
	    else if (j == hi(2)) then
	      dspec = dspec + abs(mf_spec(i,j)/mf_rho(i,j) - mf_spec(i,j-1)/mf_rho(i,j-1))
	      ! centered difference
	    else
	      dspec = dspec + abs(mf_spec(i,j+1)/mf_rho(i,j+1) - mf_spec(i,j-1)/mf_rho(i,j-1))
	    endif
	    
	    
	    if (i == lo(1)) then
	      dspec = dspec + abs(mf_spec(i+1,j)/mf_rho(i+1,j) - mf_spec(i,j)/mf_rho(i,j))
	      ! backward difference
	    else if (i == hi(1)) then
	      dspec = dspec + abs(mf_spec(i,j)/mf_rho(i,j) - mf_spec(i-1,j)/mf_rho(i-1,j))
	      ! centered difference
	    else
	      dspec = dspec + abs(mf_spec(i+1,j)/mf_rho(i+1,j) - mf_spec(i-1,j)/mf_rho(i-1,j))
	    endif
	    
	    
	    if (dspec .gt. X_grad_min .and. mf_spec(i,j)/mf_rho(i,j) .lt. radiative_X) then
	      tagbox(i,j) = .true.
	    endif

	enddo
      enddo
    !$OMP END PARALLEL DO
    endselect
    
  end subroutine tag_boxes_2d

  subroutine tag_boxes_3d(tagbox,tlo,mf_spec,mf_rho,mflo,lo,hi,radiative_X,X_grad_min,dx,lev)


    integer          , intent(in   ) :: lo(3),hi(3),tlo(3),mflo(3)
    logical          , intent(inout) :: tagbox( tlo(1):, tlo(2):, tlo(3):)
    real(kind = dp_t), intent(in   ) ::     mf_spec(mflo(1):,mflo(2):,mflo(3):)
    real(kind = dp_t), intent(in   ) ::     mf_rho(mflo(1):,mflo(2):,mflo(3):)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev
    real(kind = dp_t), intent(in   ) :: radiative_X, X_grad_min 
   

    ! local variables
    integer :: i,j,k
    real(kind = dp_t)  :: dspec
     ! initially say that we do not want to tag any cells for refinement
    tagbox(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = .false.


    select case (lev)
    case(1)
    !$OMP PARALLEL DO PRIVATE(i,j,k)    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (mf_spec(i,j,k)/mf_rho(i,j,k) .lt. radiative_X) then
               tagbox(i,j,k) = .true.
             endif 
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    case(2)

    !$OMP PARALLEL DO PRIVATE(i,j,k,dspec)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
	     
	     dspec = ZERO             
             if (k == lo(3)) then
                dspec = dspec + abs(mf_spec(i,j,k+1)/mf_rho(i,j,k+1) - mf_spec(i,j,k)/mf_rho(i,j,k))
                ! backward difference
             else if (k == hi(3)) then
                dspec = dspec + abs(mf_spec(i,j,k)/mf_rho(i,j,k) - mf_spec(i,j,k-1)/mf_rho(i,j,k-1))
                ! centered difference
             else
                dspec = dspec + abs(mf_spec(i,j,k+1)/mf_rho(i,j,k+1) - mf_spec(i,j,k-1)/mf_rho(i,j,k-1))
             endif

             
             if (j == lo(2)) then
                dspec = dspec + abs(mf_spec(i,j+1,k)/mf_rho(i,j+1,k) - mf_spec(i,j,k)/mf_rho(i,j,k))
                ! backward difference
             else if (j == hi(2)) then
                dspec = dspec + abs(mf_spec(i,j,k)/mf_rho(i,j,k) - mf_spec(i,j-1,k)/mf_rho(i,j-1,k))
                ! centered difference
             else
                dspec = dspec + abs(mf_spec(i,j+1,k)/mf_rho(i,j+1,k) - mf_spec(i,j-1,k)/mf_rho(i,j-1,k))
             endif
             
             
             if (i == lo(1)) then
                dspec = dspec + abs(mf_spec(i+1,j,k)/mf_rho(i+1,j,k) - mf_spec(i,j,k)/mf_rho(i,j,k))
                ! backward difference
             else if (i == hi(1)) then
                dspec = dspec + abs(mf_spec(i,j,k)/mf_rho(i,j,k) - mf_spec(i-1,j,k)/mf_rho(i-1,j,k))
                ! centered difference
             else
                dspec = dspec + abs(mf_spec(i+1,j,k)/mf_rho(i+1,j,k) - mf_spec(i-1,j,k)/mf_rho(i-1,j,k))
             endif
             
             
             if (dspec .gt. X_grad_min .and. mf_spec(i,j,k)/mf_rho(i,j,k) .lt. radiative_X) then
                tagbox(i,j,k) = .true.
             endif

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    endselect
    
  end subroutine tag_boxes_3d

end module tag_boxes_module
