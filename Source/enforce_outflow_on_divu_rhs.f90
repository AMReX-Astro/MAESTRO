module enforce_outflow_on_divu_module

  use bl_types
  use bc_module
  use define_bc_module
  use multifab_module
  use bl_constants_module
 
  implicit none
 
  private
 
  public :: enforce_outflow_on_divu_rhs

contains

    subroutine enforce_outflow_on_divu_rhs(divu_rhs,the_bc_tower)

      use geometry, only : dm,nlevs

      type(multifab) , intent(inout) :: divu_rhs(:)
      type(bc_tower) , intent(in   ) :: the_bc_tower

      integer        :: i,n,ng_d
      type(bc_level) :: bc
      real(kind=dp_t), pointer :: divp(:,:,:,:) 

      ng_d = nghost(divu_rhs(1))

      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, nboxes(divu_rhs(n))
            if ( multifab_remote(divu_rhs(n), i) ) cycle
            divp => dataptr(divu_rhs(n)     , i)
            select case (dm)
            case (1)
               call enforce_outflow_1d(divp(:,1,1,1), ng_d, bc%phys_bc_level_array(i,:,:))
            case (2)
               call enforce_outflow_2d(divp(:,:,1,1), ng_d, bc%phys_bc_level_array(i,:,:))
            case (3)
               call enforce_outflow_3d(divp(:,:,:,1), ng_d, bc%phys_bc_level_array(i,:,:))
            end select
         end do
      end do

    end subroutine enforce_outflow_on_divu_rhs

    ! ******************************************************************************** !

    subroutine enforce_outflow_1d(divu_rhs,ng_d,phys_bc)

      integer        , intent(in   ) :: ng_d
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_d:)
      integer        , intent(in   ) :: phys_bc(:,:)

      integer :: nx
      nx = size(divu_rhs,dim=1)-1

      if (phys_bc(1,1) .eq. OUTLET) divu_rhs(0 ) = ZERO
      if (phys_bc(1,2) .eq. OUTLET) divu_rhs(nx) = ZERO

    end subroutine enforce_outflow_1d


    ! ******************************************************************************** !

    subroutine enforce_outflow_2d(divu_rhs,ng_d,phys_bc)

      integer        , intent(in   ) :: ng_d
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_d:,-ng_d:)
      integer        , intent(in   ) :: phys_bc(:,:)

      integer :: nx,ny
      nx = size(divu_rhs,dim=1)-1
      ny = size(divu_rhs,dim=2)-1

      if (phys_bc(1,1) .eq. OUTLET) divu_rhs(0,  :) = ZERO
      if (phys_bc(1,2) .eq. OUTLET) divu_rhs(nx, :) = ZERO
      if (phys_bc(2,1) .eq. OUTLET) divu_rhs(: , 0) = ZERO
      if (phys_bc(2,2) .eq. OUTLET) divu_rhs(: ,ny) = ZERO

    end subroutine enforce_outflow_2d

    ! ******************************************************************************** !

    subroutine enforce_outflow_3d(divu_rhs,ng_d,phys_bc)

      integer        , intent(in   ) :: ng_d
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_d:,-ng_d:,-ng_d:)
      integer        , intent(in   ) :: phys_bc(:,:)

      integer :: nx,ny,nz
      nx = size(divu_rhs,dim=1)-1
      ny = size(divu_rhs,dim=2)-1
      nz = size(divu_rhs,dim=3)-1

      if (phys_bc(1,1) .eq. OUTLET) divu_rhs(0,  :, :) = ZERO
      if (phys_bc(1,2) .eq. OUTLET) divu_rhs(nx, :, :) = ZERO
      if (phys_bc(2,1) .eq. OUTLET) divu_rhs( :, 0, :) = ZERO
      if (phys_bc(2,2) .eq. OUTLET) divu_rhs( :,ny, :) = ZERO
      if (phys_bc(3,1) .eq. OUTLET) divu_rhs( :,: , 0) = ZERO
      if (phys_bc(3,2) .eq. OUTLET) divu_rhs( :,: ,nz) = ZERO

    end subroutine enforce_outflow_3d

end module enforce_outflow_on_divu_module
