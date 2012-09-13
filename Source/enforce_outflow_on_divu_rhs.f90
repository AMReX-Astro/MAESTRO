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

      type(multifab) , intent(inout) :: divu_rhs(:)
      type(bc_tower) , intent(in   ) :: the_bc_tower

      integer        :: i,n,ng_d,dm,nlevs,gid
      integer        :: lo(get_dim(divu_rhs(1))),hi(get_dim(divu_rhs(1)))
      type(bc_level) :: bc
      real(kind=dp_t), pointer :: divp(:,:,:,:) 

      dm = get_dim(divu_rhs(1))
      nlevs = size(divu_rhs)

      ng_d = nghost(divu_rhs(1))

      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, nfabs(divu_rhs(n))
            gid =   global_index(divu_rhs(n),i)
            divp => dataptr(divu_rhs(n) , i)
            lo =    lwb(get_box(divu_rhs(n),i))
            hi =    upb(get_box(divu_rhs(n),i))
            select case (dm)
            case (1)
               call enforce_outflow_1d(divp(:,1,1,1), lo, hi, ng_d, &
                                       bc%phys_bc_level_array(gid,:,:))
            case (2)
               call enforce_outflow_2d(divp(:,:,1,1), lo, hi, ng_d, &
                                       bc%phys_bc_level_array(gid,:,:))
            case (3)
               call enforce_outflow_3d(divp(:,:,:,1), lo, hi, ng_d, &
                                       bc%phys_bc_level_array(gid,:,:))
            end select
         end do
      end do

    end subroutine enforce_outflow_on_divu_rhs

    ! ******************************************************************************** !

    subroutine enforce_outflow_1d(divu_rhs,lo,hi,ng_d,phys_bc)

      integer        , intent(in   ) :: ng_d,lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: divu_rhs(lo(1)-ng_d:)
      integer        , intent(in   ) :: phys_bc(:,:)

      if (phys_bc(1,1) .eq. OUTLET) divu_rhs(lo(1)  ) = ZERO
      if (phys_bc(1,2) .eq. OUTLET) divu_rhs(hi(1)+1) = ZERO

    end subroutine enforce_outflow_1d


    ! ******************************************************************************** !

    subroutine enforce_outflow_2d(divu_rhs,lo,hi,ng_d,phys_bc)

      integer        , intent(in   ) :: ng_d,lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: divu_rhs(lo(1)-ng_d:,lo(2)-ng_d:)
      integer        , intent(in   ) :: phys_bc(:,:)

      if (phys_bc(1,1) .eq. OUTLET) divu_rhs(lo(1)  ,:) = ZERO
      if (phys_bc(1,2) .eq. OUTLET) divu_rhs(hi(1)+1,:) = ZERO
      if (phys_bc(2,1) .eq. OUTLET) divu_rhs(:,lo(2)  ) = ZERO
      if (phys_bc(2,2) .eq. OUTLET) divu_rhs(:,hi(2)+1) = ZERO

    end subroutine enforce_outflow_2d

    ! ******************************************************************************** !

    subroutine enforce_outflow_3d(divu_rhs,lo,hi,ng_d,phys_bc)

      integer        , intent(in   ) :: ng_d,lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: divu_rhs(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:)
      integer        , intent(in   ) :: phys_bc(:,:)

      if (phys_bc(1,1) .eq. OUTLET) divu_rhs(lo(1)  ,:,:) = ZERO
      if (phys_bc(1,2) .eq. OUTLET) divu_rhs(hi(1)+1,:,:) = ZERO
      if (phys_bc(2,1) .eq. OUTLET) divu_rhs(:,lo(2)  ,:) = ZERO
      if (phys_bc(2,2) .eq. OUTLET) divu_rhs(:,hi(2)+1,:) = ZERO
      if (phys_bc(3,1) .eq. OUTLET) divu_rhs(:,:,lo(3)  ) = ZERO
      if (phys_bc(3,2) .eq. OUTLET) divu_rhs(:,:,hi(3)+1) = ZERO

    end subroutine enforce_outflow_3d

end module enforce_outflow_on_divu_module
