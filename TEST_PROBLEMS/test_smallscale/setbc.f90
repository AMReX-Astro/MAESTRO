module setbc_module

  use bl_types
  use bl_error_module

  implicit none

  private

  public :: setbc_2d, setbc_3d

contains

  subroutine setbc_2d(s,lo,ng,bc,icomp)

    use bc_module
    use bl_constants_module
    use inlet_bc_module
    use network

    integer        , intent(in   ) :: lo(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    !     Local variables
    integer         :: i,j,hi(2)
    real(kind=dp_t) :: x,dir_val

    if (ng == 0) return

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

    if (bc(1,1) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,1) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,1) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,1) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,1) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,1) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. INTERIOR) then

    else
       print *,''
       print *,'*******************************************'
       print *,'In setbc.f90: bc(1,1) = ',bc(1,1),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

    if (bc(1,2) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,2) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''
    else if (bc(1,2) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,2) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,2) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,2) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,2) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(1,2) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,2) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''
    else if (bc(1,2) .eq. INTERIOR) then

    else 
       print *,''
       print *,'*******************************************'
       print *,'In setbc.f90: bc(1,2) = ',bc(1,2),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

    if (bc(2,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_VT
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_VN
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_RHO
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_RHOH
       if (icomp .eq. 5 .or. icomp .eq. 6 .or. icomp .eq. 7) then
          if(spec_names(icomp-4) .eq. "carbon-12") then
             s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_RHOC12
          else if(spec_names(icomp-4) .eq. "magnesium-24") then
             s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_RHOMG24
          else if(spec_names(icomp-4) .eq. "oxygen-16") then
             s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_RHOO16
          endif
       endif
       if (icomp.eq.8) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_TEMP
       if (icomp.eq.9) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_TRA
    else if (bc(2,1) .eq. FOEXTRAP) then
       do i = lo(1)-ng,hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = s(i,lo(2))
       end do
    else if (bc(2,1) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,1) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,1) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(2,1) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,1) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''
    else if (bc(2,1) .eq. INTERIOR) then

    else 
       print *,''
       print *,'*******************************************'
       print *,'In setbc.f90: bc(2,1) = ',bc(2,1),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

    if (bc(2,2) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,2) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''
    else if (bc(2,2) .eq. FOEXTRAP) then
       do i = lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2))
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,2) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,2) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(2,2) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,2) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''
    else if (bc(2,2) .eq. INTERIOR) then

    else 
       print *,''
       print *,'*******************************************'
       print *,'In setbc.f90: bc(2,2) = ',bc(2,2),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

  end subroutine setbc_2d

  subroutine setbc_3d(s,lo,ng,bc,icomp)

    use bc_module
    use bl_constants_module
    use inlet_bc_module
    use network

    integer        , intent(in   ) :: lo(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    !     Local variables
    integer :: i,j,k,hi(3)

    if (ng == 0) return

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
    hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

    if (bc(1,1) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,1) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,1) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,1) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,1) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,1) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,'' 
    else if (bc(1,1) .eq. INTERIOR) then

    else
       print *,''
       print *,'*******************************************'
       print *,'In setbc.f90: bc(1,1) = ',bc(1,1),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

    if (bc(1,2) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,2) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''         
    else if (bc(1,2) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,2) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''         
    else if (bc(1,2) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,2) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''        
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,2) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''         
    else if (bc(1,2) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(1,2) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''         
    else if (bc(1,2) .eq. INTERIOR) then

    else 
       print *,''
       print *,'*******************************************'
       print *,'In setbc.f90: bc(1,2) = ',bc(1,2),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''         
    end if

    if (bc(2,1) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,1) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,1) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,1) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,1) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,1) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,1) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,1) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,1) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,1) .eq. INTERIOR) then

    else
       print *,''
       print *,'*******************************************'
       print *,'In setbc.f90: bc(2,1) = ',bc(2,1),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''         
    end if

    if (bc(2,2) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,2) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,2) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,2) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,2) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,2) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,2) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,2) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(2,2) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,2) .eq. INTERIOR) then

    else 
       print *,''
       print *,'*******************************************'
       print *,'In setbc.f90: bc(2,2) = ',bc(2,2),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''        
    end if

    if (bc(3,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_VT
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_VT
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_VN
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_RHO
       if (icomp.eq.5) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_RHOH
       if (icomp .eq. 6 .or. icomp .eq. 7 .or. icomp .eq. 8) then
          if(spec_names(icomp-5) .eq. "carbon-12") then
             s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_RHOC12
          else if(spec_names(icomp-5) .eq. "magnesium-24") then
             s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_RHOMG24
          else if(spec_names(icomp-5) .eq. "oxygen-16") then
             s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_RHOO16
          endif
       endif
       if (icomp.eq.9) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_TEMP
       if(icomp.eq.10) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_TRA
    else if (bc(3,1) .eq. FOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = s(i,j,lo(3))
          end do
       end do
    else if(bc(3,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(3,1) .eq. FOEXTRAP .or. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''         
    else if (bc(3,1) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(3,1) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''        
    else if (bc(3,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(3,1) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''         
    else if (bc(3,1) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(3,1) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''         
    else if (bc(3,1) .eq. INTERIOR) then

    else 
       print *,''
       print *,'*******************************************'
       print *,'In setbc.f90: bc(3,1) = ',bc(3,1),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''         
    end if

    if (bc(3,2) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(3,2) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''
    else if (bc(3,2) .eq. FOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = s(i,j,hi(3))
          end do
       end do
    else if (bc(3,2) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(3,2) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(3,2) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(3,2) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(3,2) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In setbc.f90: bc(3,2) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''
    else if (bc(3,2) .eq. INTERIOR) then

    else 
       print *,''
       print *,'*******************************************'
       print *,'In setbc.f90: bc(3,2) = ',bc(3,2),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

  end subroutine setbc_3d

end module setbc_module
