module multifab_physbc_module

  use multifab_module
  use define_bc_module
  use bl_types
  use bl_error_module

  implicit none

  private

  public :: multifab_physbc

contains

  subroutine multifab_physbc(s,start_scomp,start_bccomp,num_comp,the_bc_level)

    use bl_prof_module

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level

    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s)),dm
    integer                  :: i,ng,scomp,bccomp
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    type(bl_prof_timer), save :: bpt
    
    dm = get_dim(s)

    ng = nghost(s)

    if (ng == 0) return

    call build(bpt, "multifab_physbc")
    
    do i=1,nboxes(s)
       if ( multifab_remote(s,i) ) cycle
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (1)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_1d(sp(:,1,1,scomp), lo, hi, ng, &
                           the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
          end do
       case (2)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_2d(sp(:,:,1,scomp), lo, hi, ng, &
                           the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
          end do
       case (3)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_3d(sp(:,:,:,scomp), lo, hi, ng, &
                           the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
          end do
       end select
    end do

    call destroy(bpt)
 
  end subroutine multifab_physbc

  subroutine physbc_1d(s,lo,hi,ng,bc,icomp)

    use bc_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    print *,''
    print *,'*******************************************'
    print *,'WARNING: physbc_1d not supported'
    print *,'*******************************************'
    print *,''

  end subroutine physbc_1d

  subroutine physbc_2d(s,lo,hi,ng,bc,icomp)

    use bc_module
    use bl_constants_module
    use inlet_bc_module
    use network

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    !     Local variables
    integer         :: i,j
    real(kind=dp_t) :: dir_val

    if(ng == 0) return

    if (bc(1,1) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,1) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,1) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,1) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,1) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,1) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. INTERIOR) then
    else
       print *,''
       print *,'*******************************************'
       print *,'In physbc.f90: bc(1,1) = ',bc(1,1),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

    if (bc(1,2) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,2) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''
    else if (bc(1,2) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,2) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,2) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,2) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,2) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(1,2) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,2) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''
    else if (bc(1,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       print *,''
       print *,'*******************************************'
       print *,'In physbc.f90: bc(1,2) = ',bc(1,2),' NOT YET SUPPORTED '
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
       print *,'WARNING: In physbc.f90: bc(2,1) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,1) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(2,1) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,1) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''
    else if (bc(2,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       print *,''
       print *,'*******************************************'
       print *,'In physbc.f90: bc(2,1) = ',bc(2,1),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

    if (bc(2,2) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,2) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''
    else if (bc(2,2) .eq. FOEXTRAP) then
       do i = lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2))
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,2) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,2) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(2,2) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,2) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''
    else if (bc(2,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       print *,''
       print *,'*******************************************'
       print *,'In physbc.f90: bc(2,2) = ',bc(2,2),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

  end subroutine physbc_2d

  subroutine physbc_3d(s,lo,hi,ng,bc,icomp)

    use bc_module
    use bl_constants_module
    use inlet_bc_module
    use network

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    !     Local variables
    integer :: i,j,k

    if (ng == 0) return

    if (bc(1,1) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,1) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,1) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,1) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,1) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(1,1) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,1) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,'' 
    else if (bc(1,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       print *,''
       print *,'*******************************************'
       print *,'In physbc.f90: bc(1,1) = ',bc(1,1),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

    if (bc(1,2) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,2) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''         
    else if (bc(1,2) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,2) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''         
    else if (bc(1,2) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,2) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''        
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,2) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''         
    else if (bc(1,2) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(1,2) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''         
    else if (bc(1,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       print *,''
       print *,'*******************************************'
       print *,'In physbc.f90: bc(1,2) = ',bc(1,2),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''         
    end if

    if (bc(2,1) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,1) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,1) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,1) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,1) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,1) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,1) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,1) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,1) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       print *,''
       print *,'*******************************************'
       print *,'In physbc.f90: bc(2,1) = ',bc(2,1),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''         
    end if

    if (bc(2,2) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,2) .eq. EXT_DIR'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,2) .eq. FOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,2) .eq. FOEXTRAP'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,2) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,2) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,2) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,2) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(2,2) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''         
    else if (bc(2,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       print *,''
       print *,'*******************************************'
       print *,'In physbc.f90: bc(2,2) = ',bc(2,2),' NOT YET SUPPORTED '
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
       print *,'WARNING: In physbc.f90: bc(3,1) .eq. FOEXTRAP .or. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''         
    else if (bc(3,1) .eq. HOEXTRAP) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(3,1) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''        
    else if (bc(3,1) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(3,1) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''         
    else if (bc(3,1) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(3,1) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''         
    else if (bc(3,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       print *,''
       print *,'*******************************************'
       print *,'In physbc.f90: bc(3,1) = ',bc(3,1),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''         
    end if

    if (bc(3,2) .eq. EXT_DIR) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(3,2) .eq. EXT_DIR'
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
       print *,'WARNING: In physbc.f90: bc(3,2) .eq. HOEXTRAP'
       print *,'*******************************************'
       print *,''
    else if (bc(3,2) .eq. REFLECT_EVEN) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(3,2) .eq. REFLECT_EVEN'
       print *,'*******************************************'
       print *,''
    else if (bc(3,2) .eq. REFLECT_ODD) then
       print *,''
       print *,'*******************************************'
       print *,'WARNING: In physbc.f90: bc(3,2) .eq. REFLECT_ODD'
       print *,'*******************************************'
       print *,''
    else if (bc(3,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else 
       print *,''
       print *,'*******************************************'
       print *,'In physbc.f90: bc(3,2) = ',bc(3,2),' NOT YET SUPPORTED '
       print *,'*******************************************'
       print *,''
    end if

  end subroutine physbc_3d

end module multifab_physbc_module
