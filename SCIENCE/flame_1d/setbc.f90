! set the boundary conditions for a laminar flame.  On one end we
! will have inflow and on the other end we will have outflow.  
!
! In 1-d the flame is assumed to move in the X direction
! In 2-d the flame is assumed to move in the Y direction
! In 3-d the flame is assumed to move in the Z direction
!
! In the inputs file, these are specified via
!
! bcy_lo = 11   ! 2-d
! bcy_hi = 12
!
! where boxlib/bc.f90 defines the integer parameters 
!
! INLET        = 11
! OUTLET       = 12
!
! note: INLET and OUTLET are "macros" that define_bc_tower.f90
! interprets and uses to set the BC type for each of the components
!

module setbc_module

  use bl_types
  use bl_error_module

  implicit none

  private

  public :: setbc_1d, setbc_2d, setbc_3d

contains

  !============================================================================
  ! setbc_1d
  !============================================================================
  subroutine setbc_1d(s,lo,hi,ng,bc,icomp)

    use bc_module
    use bl_constants_module
    use inlet_bc_module
    use network
    use variables

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    character (len=256) err_string
    integer :: dm

    if(ng == 0) return

    if ( (bc(1,1) == EXT_DIR .or. &
          bc(1,2) == EXT_DIR) .and. .NOT. inlet_bc_initialized) then
       call bl_error("ERROR: in setbc, but inlet BCs not initialized")
    endif

    dm = 1
    

    !--------------------------------------------------------------------------
    ! lower X
    !--------------------------------------------------------------------------
    if (bc(1,1) == EXT_DIR) then

       ! velocity components
       if (icomp == 1) s(lo(1)-ng:lo(1)-1) = INLET_VEL 

       ! density
       if (icomp == dm+rho_comp) s(lo(1)-ng:lo(1)-1) = INLET_RHO

       ! rho * h
       if (icomp == dm+rhoh_comp) s(lo(1)-ng:lo(1)-1) = INLET_RHOH
       
       ! species
       if (icomp >= dm+spec_comp .and. icomp < dm+spec_comp+nspec) then
          s(lo(1)-ng:lo(1)-1) = &
               INLET_RHOX(icomp-(spec_comp+dm)+1)
       endif

       ! others
       if (icomp == dm+temp_comp) s(lo(1)-ng:lo(1)-1) = INLET_TEMP
       if (icomp == dm+trac_comp) s(lo(1)-ng:lo(1)-1) = INLET_TRA

    else if (bc(1,1) == FOEXTRAP) then
       s(lo(1)-ng:lo(1)-1) = s(lo(1))

    else if (bc(1,1) == INTERIOR) then
       ! nothing to do 

    else 
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(1,1) = ", bc(1,1)
       call bl_error(trim(err_string))
    end if

    !--------------------------------------------------------------------------
    ! upper X
    !--------------------------------------------------------------------------
    if (bc(1,2) == FOEXTRAP) then
       s(hi(1)+1:hi(1)+ng) = s(hi(1))

    else if (bc(1,2) == INTERIOR) then
       ! nothing to do

    else
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(2,2) = ", bc(1,2)
       call bl_error(trim(err_string))
    end if

  end subroutine setbc_1d

  !============================================================================
  ! setbc_2d
  !============================================================================
  subroutine setbc_2d(s,lo,hi,ng,bc,icomp)

    use bc_module
    use bl_constants_module
    use inlet_bc_module
    use network
    use variables
    
    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    !     Local variables
    integer         :: i, dm
    
    character (len=256) err_string

    if(ng == 0) return

    if ( (bc(1,1) == EXT_DIR .or. &
          bc(1,2) == EXT_DIR .or. &
          bc(2,1) == EXT_DIR .or. &
          bc(2,2) == EXT_DIR) .and. .NOT. inlet_bc_initialized) then
       call bl_error("ERROR: in setbc, but inlet BCs not initialized")
    endif

    dm = 2

    !--------------------------------------------------------------------------
    ! lower X 
    !--------------------------------------------------------------------------

    ! we are assuming that the domain is periodic in X
    if (bc(1,1) /= INTERIOR) then
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(1,1) = ", bc(1,1)
       call bl_error(trim(err_string))
    end if


    !--------------------------------------------------------------------------
    ! upper X 
    !--------------------------------------------------------------------------

    ! we are assuming that the domain is periodic in X
    if (bc(1,2) /= INTERIOR) then
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(1,2) = ", bc(1,2)
       call bl_error(trim(err_string))
    end if


    !--------------------------------------------------------------------------
    ! lower Y
    !--------------------------------------------------------------------------
    if (bc(2,1) == EXT_DIR) then

       ! velocity components
       if (icomp == 1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = ZERO
       if (icomp == 2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_VEL

       ! density
       if (icomp == dm+rho_comp) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_RHO

       ! rho * h
       if (icomp == dm+rhoh_comp) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_RHOH
       
       ! species
       if (icomp >= dm+spec_comp .and. icomp < dm+spec_comp+nspec) then
          s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = &
               INLET_RHOX(icomp-(spec_comp+dm)+1)
       endif

       ! others
       if (icomp == dm+temp_comp) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_TEMP
       if (icomp == dm+trac_comp) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_TRA

    else if (bc(2,1) == FOEXTRAP) then
       do i = lo(1)-ng,hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = s(i,lo(2))
       end do

    else if (bc(2,1) == INTERIOR) then
       ! nothing to do 

    else 
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(2,1) = ", bc(2,1)
       call bl_error(trim(err_string))
    end if


    !--------------------------------------------------------------------------
    ! upper Y
    !--------------------------------------------------------------------------
    if (bc(2,2) == FOEXTRAP) then
       do i = lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2))
       end do

    else if (bc(2,2) == INTERIOR) then
       ! nothing to do

    else
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(2,2) = ", bc(2,2)
       call bl_error(trim(err_string))
    end if

  end subroutine setbc_2d


  !============================================================================
  ! setbc_3d
  !============================================================================
  subroutine setbc_3d(s,lo,hi,ng,bc,icomp)

    use bc_module
    use bl_constants_module
    use inlet_bc_module
    use network
    use variables

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    !     Local variables
    integer :: i,j, dm

    character (len=256) err_string

    if (ng == 0) return

    if ( (bc(1,1) == EXT_DIR .or. &
          bc(1,2) == EXT_DIR .or. &
          bc(2,1) == EXT_DIR .or. &
          bc(2,2) == EXT_DIR .or. &
          bc(3,1) == EXT_DIR .or. &
          bc(3,2) == EXT_DIR) .and. .NOT. inlet_bc_initialized) then
       call bl_error("ERROR: in setbc, but inlet BCs not initialized")
    endif

    dm = 3

    !--------------------------------------------------------------------------
    ! lower X 
    !--------------------------------------------------------------------------

    ! we are assuming that the domain is periodic in X
    if (bc(1,1) /= INTERIOR) then
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(1,1) = ", bc(1,1)
       call bl_error(trim(err_string))
    end if


    !--------------------------------------------------------------------------
    ! upper X 
    !--------------------------------------------------------------------------

    ! we are assuming that the domain is periodic in X
    if (bc(1,2) /= INTERIOR) then
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(1,2) = ", bc(1,2)
       call bl_error(trim(err_string))
    end if


    !--------------------------------------------------------------------------
    ! lower Y
    !--------------------------------------------------------------------------

    ! we are assuming that the domain is periodic in X
    if (bc(2,1) /= INTERIOR) then
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(2,1) = ", bc(2,1)
       call bl_error(trim(err_string))
    end if


    !--------------------------------------------------------------------------
    ! upper Y 
    !--------------------------------------------------------------------------

    ! we are assuming that the domain is periodic in X
    if (bc(2,2) /= INTERIOR) then
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(2,2) = ", bc(2,2)
       call bl_error(trim(err_string))
    end if


    !--------------------------------------------------------------------------
    ! lower Z                                                                  
    !--------------------------------------------------------------------------
    if (bc(3,1) == EXT_DIR) then

       ! velocity components
       if (icomp == 1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = ZERO
       if (icomp == 2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = ZERO
       if (icomp == 3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_VEL

       ! density
       if (icomp == dm+rho_comp) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_RHO

       ! rho * h
       if (icomp == dm+rhoh_comp) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_RHOH

       ! species
       if (icomp >= dm+spec_comp .and. icomp < dm+spec_comp+nspec) then
          s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = &
               INLET_RHOX(icomp-(spec_comp+dm)+1)
       endif

       ! other
       if (icomp == dm+temp_comp) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_TEMP
       if (icomp == dm+trac_comp) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_TRA

    else if (bc(3,1) == FOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = s(i,j,lo(3))
          end do
       end do

    else if (bc(3,1) == INTERIOR) then
       ! nothing to do
       
    else 
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(3,1) = ", bc(3,1)
       call bl_error(trim(err_string))
    end if


    !--------------------------------------------------------------------------
    ! upper Z
    !--------------------------------------------------------------------------
    if (bc(3,2) == FOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = s(i,j,hi(3))
          end do
       end do

    else if (bc(3,2) == INTERIOR) then
       ! nothing to do

    else 
       write (unit=err_string, fmt="(a, i3)") &
            "ERROR: setbc does not implement bc(3,2) = ", bc(3,2)
       call bl_error(trim(err_string))
    end if

  end subroutine setbc_3d

end module setbc_module
