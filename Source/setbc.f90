! set the boundary conditions

module setbc_module

  use bl_types
  use bl_error_module

  implicit none

  private

  public :: setbc_1d, setbc_2d, setbc_3d

contains

  subroutine setbc_1d(s,lo,hi,ng,bc)

    use bc_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:)
    integer        , intent(in   ) :: bc(:,:)

    !     Local variables
    integer :: i

    if(ng == 0) return

    !--------------------------------------------------------------------------
    ! lower 
    !--------------------------------------------------------------------------
    if (bc(1,1) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(lo(1)-ng:lo(1)-1) = 0.d0
    else if (bc(1,1) .eq. FOEXTRAP) then
       s(lo(1)-ng:lo(1)-1) = s(lo(1))
    else if (bc(1,1) .eq. HOEXTRAP) then
       s(lo(1)-ng:lo(1)-1) = &
            ( 15.d0 * s(lo(1)  ) &
             -10.d0 * s(lo(1)+1) &
             + 3.d0 * s(lo(1)+2) ) * EIGHTH
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       do i = 1,ng
          s(lo(1)-i) = s(lo(1)+i-1)
       end do
    else if (bc(1,1) .eq. REFLECT_ODD) then
       do i = 1,ng
          s(lo(1)-i) = -s(lo(1)+i-1)
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(1,1) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_1d: bc(1,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper 
    !--------------------------------------------------------------------------
    if (bc(1,2) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(hi(1)+1:hi(1)+ng) = 0.d0
    else if (bc(1,2) .eq. FOEXTRAP) then
       s(hi(1)+1:hi(1)+ng) = s(hi(1))
    else if (bc(1,2) .eq. HOEXTRAP) then
       s(hi(1)+1:hi(1)+ng) = &
            ( 15.d0 * s(hi(1)  ) &
             -10.d0 * s(hi(1)-1) &
             + 3.d0 * s(hi(1)-2) ) * EIGHTH
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       do i = 1,ng
          s(hi(1)+i) = s(hi(1)-i+1)
       end do
    else if (bc(1,2) .eq. REFLECT_ODD) then
       do i = 1,ng
          s(hi(1)+i) = -s(hi(1)-i+1)
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(1,2) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_1d: bc(1,2) not yet supported")
    end if

  end subroutine setbc_1d

  subroutine setbc_2d(s,lo,hi,ng,bc)

    use bc_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)

    !     Local variables
    integer :: i,j

    if(ng == 0) return

    !-------------------------------------------------------------------------- 
    ! lower X                                                                  
    !--------------------------------------------------------------------------
    if (bc(1,1) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(lo(1)-ng:lo(1)-1,:) = 0.d0
    else if (bc(1,1) .eq. FOEXTRAP) then
       do j = lo(2)-1,hi(2)+1
          s(lo(1)-ng:lo(1)-1,j) = s(lo(1),j)
       end do
    else if (bc(1,1) .eq. HOEXTRAP) then
       do j = lo(2)-1,hi(2)+1
          s(lo(1)-ng:lo(1)-1,j) = &
               ( 15.d0 * s(lo(1)  ,j) &
                -10.d0 * s(lo(1)+1,j) &
                + 3.d0 * s(lo(1)+2,j) ) * EIGHTH
       end do
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       do i = 1,ng
          s(lo(1)-i,lo(2)-1:hi(2)+1) = s(lo(1)+i-1,lo(2)-1:hi(2)+1)
       end do
    else if (bc(1,1) .eq. REFLECT_ODD) then
       do i = 1,ng
          s(lo(1)-i,lo(2)-1:hi(2)+1) = -s(lo(1)+i-1,lo(2)-1:hi(2)+1)
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(1,1) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_2d: bc(1,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper X
    !--------------------------------------------------------------------------
    if (bc(1,2) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(hi(1)+1:hi(1)+ng,:) = 0.d0
    else if (bc(1,2) .eq. FOEXTRAP) then
       do j = lo(2)-1,hi(2)+1
          s(hi(1)+1:hi(1)+ng,j) = s(hi(1),j)
       end do
    else if (bc(1,2) .eq. HOEXTRAP) then
       do j = lo(2)-1,hi(2)+1
          s(hi(1)+1:hi(1)+ng,j) = &
               ( 15.d0 * s(hi(1)  ,j) &
                -10.d0 * s(hi(1)-1,j) &
                + 3.d0 * s(hi(1)-2,j) ) * EIGHTH
       end do
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       do i = 1,ng
          s(hi(1)+i,lo(2)-ng:hi(2)+ng) = s(hi(1)-i+1,lo(2)-ng:hi(2)+ng)
       end do
    else if (bc(1,2) .eq. REFLECT_ODD) then
       do i = 1,ng
          s(hi(1)+i,lo(2)-ng:hi(2)+ng) = -s(hi(1)-i+1,lo(2)-ng:hi(2)+ng)
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(1,2) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_2d: bc(1,2) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! lower Y
    !--------------------------------------------------------------------------
    if (bc(2,1) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(:,lo(2)-ng:lo(2)-1) = 0.d0
    else if (bc(2,1) .eq. FOEXTRAP) then
       do i = lo(1)-ng,hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = s(i,lo(2))
       end do
    else if (bc(2,1) .eq. HOEXTRAP) then
       do i = lo(1)-ng,hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = &
               ( 15.d0 * s(i,lo(2)  ) &
                -10.d0 * s(i,lo(2)+1) &
                + 3.d0 * s(i,lo(2)+2) ) * EIGHTH
       end do
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       do j = 1,ng
          s(lo(1)-ng:hi(1)+ng,lo(2)-j) = s(lo(1)-ng:hi(1)+ng,lo(2)+j-1)
       end do
    else if (bc(2,1) .eq. REFLECT_ODD) then
       do j = 1,ng
          s(lo(1)-ng:hi(1)+ng,lo(2)-j) = -s(lo(1)-ng:hi(1)+ng,lo(2)+j-1)
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(2,1) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_2d: bc(2,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper Y
    !--------------------------------------------------------------------------
    if (bc(2,2) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(:,hi(2)+1:hi(2)+ng) = 0.d0
    else if (bc(2,2) .eq. FOEXTRAP) then
       do i = lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2))
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       do i = lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = &
               ( 15.d0 * s(i,hi(2)  ) &
                -10.d0 * s(i,hi(2)-1) &
                + 3.d0 * s(i,hi(2)-2) ) * EIGHTH
       end do
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       do j = 1,ng
          s(lo(1)-ng:hi(1)+ng,hi(2)+j) = s(lo(1)-ng:hi(1)+ng,hi(2)-j+1)
       end do
    else if (bc(2,2) .eq. REFLECT_ODD) then
       do j = 1,ng
          s(lo(1)-ng:hi(1)+ng,hi(2)+j) = -s(lo(1)-ng:hi(1)+ng,hi(2)-j+1)
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(2,2) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_2d: bc(2,2) not yet supported")
    end if

  end subroutine setbc_2d

  subroutine setbc_3d(s,lo,hi,ng,bc)

    use bc_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)

    !     Local variables
    integer :: i,j,k

    if (ng == 0) return

    !--------------------------------------------------------------------------
    ! lower X
    !--------------------------------------------------------------------------
    if (bc(1,1) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+1
             s(lo(1)-ng:lo(1)-1,j,k) = ZERO
          end do
       end do
    else if (bc(1,1) .eq. FOEXTRAP) then
       do i = 1,ng
          s(lo(1)-i,:,:) = s(lo(1),:,:)
       end do
    else if (bc(1,1) .eq. HOEXTRAP) then
       s(lo(1)-1,:,:) = &
         ( 15.d0 * s(lo(1)  ,:,:) &
          -10.d0 * s(lo(1)+1,:,:) &
          + 3.d0 * s(lo(1)+2,:,:) ) * EIGHTH
       do i = 2,ng
          s(lo(1)-i,:,:) = s(lo(1)-1,:,:)
       end do
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       do i = 1,ng
          s(lo(1)-i  ,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = &
          s(lo(1)+i-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
       end do
    else if (bc(1,1) .eq. REFLECT_ODD) then
       do i = 1,ng
          s(lo(1)-i  ,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = &
         -s(lo(1)+i-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(1,1) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_3d: bc(1,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper X
    !--------------------------------------------------------------------------
    if (bc(1,2) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(hi(1)+1:hi(1)+ng,:,:) = ZERO
    else if (bc(1,2) .eq. FOEXTRAP) then
       do i = 1,ng
          s(hi(1)+i,:,:) = s(hi(1),:,:)
       end do
    else if (bc(1,2) .eq. HOEXTRAP) then
       s(hi(1)+1,:,:) = &
         ( 15.d0 * s(hi(1)  ,:,:) &
          -10.d0 * s(hi(1)-1,:,:) &
          + 3.d0 * s(hi(1)-2,:,:) ) * EIGHTH
       do i = 2,ng
          s(hi(1)+i,:,:) = s(hi(1)+1,:,:)
       end do
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       do i = 1,ng
          s(hi(1)+i,:,:) =  s(hi(1)-i+1,:,:)
       end do
    else if (bc(1,2) .eq. REFLECT_ODD) then
       do i = 1,ng
          s(hi(1)+i,:,:) = -s(hi(1)-i+1,:,:)
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(1,2) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_3d: bc(1,2) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! lower Y
    !--------------------------------------------------------------------------
    if (bc(2,1) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(:,lo(2)-ng:lo(2)-1,:) = ZERO
    else if (bc(2,1) .eq. FOEXTRAP) then
       do j = 1,ng
          s(:,lo(2)-j,:) = s(:,lo(2),:)
       end do
    else if (bc(2,1) .eq. HOEXTRAP) then
       s(:,lo(2)-1,:) = &
            ( 15.d0 * s(:,lo(2)  ,:) &
             -10.d0 * s(:,lo(2)+1,:) &
             + 3.d0 * s(:,lo(2)+2,:) ) * EIGHTH
       do j = 2,ng
          s(:,lo(2)-j,:) = s(:,lo(2)-1,:)
       end do
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       do j = 1,ng
          s(:,lo(2)-j  ,:) =  s(:,lo(2)+j-1,:)
       end do
    else if (bc(2,1) .eq. REFLECT_ODD) then
       do j = 1,ng
          s(:,lo(2)-j  ,:) = -s(:,lo(2)+j-1,:)
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(2,1) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_3d: bc(2,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper Y
    !--------------------------------------------------------------------------
    if (bc(2,2) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(:,hi(2)+1:hi(2)+ng,:) = ZERO
    else if (bc(2,2) .eq. FOEXTRAP) then
       do j = 1,ng
          s(:,hi(2)+j,:) = s(:,hi(2),:)
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       s(:,hi(2)+1,:) = &
            ( 15.d0 * s(:,hi(2)  ,:) &
             -10.d0 * s(:,hi(2)-1,:) &
             + 3.d0 * s(:,hi(2)-2,:) ) * EIGHTH
       do j = 2,ng
          s(:,hi(2)+j,:) = s(:,hi(2)+1,:)
       end do
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       do j = 1,ng
          s(:,hi(2)+j,:) =  s(:,hi(2)-j+1,:)
       end do
    else if (bc(2,2) .eq. REFLECT_ODD) then
       do j = 1,ng
          s(:,hi(2)+j,:) = -s(:,hi(2)-j+1,:)
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(2,2) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_3d: bc(2,2) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! lower Z
    !--------------------------------------------------------------------------   
    if (bc(3,1) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(:,:,lo(3)-ng:lo(3)-1) = ZERO
    else if (bc(3,1) .eq. FOEXTRAP) then
       do k = 1,ng
          s(:,:,lo(3)-k) = s(:,:,lo(3))
       end do
    else if (bc(3,1) .eq. HOEXTRAP) then
       s(:,:,lo(3)-1) = &
            ( 15.d0 * s(:,:,lo(3)  ) &
             -10.d0 * s(:,:,lo(3)+1) &
             + 3.d0 * s(:,:,lo(3)+2) ) * EIGHTH
       do k = 2,ng
          s(:,:,lo(3)-k) = s(:,:,lo(3)-1)
       end do
    else if (bc(3,1) .eq. REFLECT_EVEN) then
       do k = 1,ng
          s(:,:,lo(3)-k) =  s(:,:,lo(3)+k-1)
       end do
    else if (bc(3,1) .eq. REFLECT_ODD) then
       do k = 1,ng
          s(:,:,lo(3)-k) = -s(:,:,lo(3)+k-1)
       end do
    else if (bc(3,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(3,1) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_3d: bc(3,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper Z
    !--------------------------------------------------------------------------
    if (bc(3,2) .eq. EXT_DIR) then
       ! NOTE: if you want an inhomogeneous dirichlet condition, you need to write it
       s(:,:,hi(3)+1:hi(3)+ng) = ZERO
    else if (bc(3,2) .eq. FOEXTRAP) then
       do k = 1,ng
          s(:,:,hi(3)+k) = s(:,:,hi(3))
       end do
    else if (bc(3,2) .eq. HOEXTRAP) then
       s(:,:,hi(3)+1) = &
            ( 15.d0 * s(:,:,hi(3)  ) &
             -10.d0 * s(:,:,hi(3)-1) &
             + 3.d0 * s(:,:,hi(3)-2) ) * EIGHTH
       do k = 2,ng
          s(:,:,hi(3)+k) = s(:,:,hi(3)+1)
       end do
    else if (bc(3,2) .eq. REFLECT_EVEN) then
       do k = 1,ng
          s(:,:,hi(3)+k) =  s(:,:,hi(3)-k+1)
       end do
    else if (bc(3,2) .eq. REFLECT_ODD) then
       do k = 1,ng
          s(:,:,hi(3)+k) = -s(:,:,hi(3)-k+1)
       end do
    else if (bc(3,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(3,2) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       call bl_error("setbc_3d: bc(3,2) not yet supported")
    end if

  end subroutine setbc_3d

end module setbc_module
