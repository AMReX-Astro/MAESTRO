module ppm_module

  use bl_types

  implicit none

  private

  public :: ppm_2d, ppm_fpu_2d, ppm_3d, ppm_fpu_3d

contains

  ! characteristics based on u
  subroutine ppm_2d(s,ng_s,u,ng_u,slx,sly,lo,hi,bc)

    use bc_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_u
    real(kind=dp_t), intent(in   ) ::   s(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(in   ) ::   u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real(kind=dp_t), intent(  out) :: slx(lo(1)-1   :,lo(2)-1   :,:) 
    real(kind=dp_t), intent(  out) :: sly(lo(1)-1   :,lo(2)-1   :,:) 
    integer        , intent(in   ) :: bc(:,:,:)

    ! local



  end subroutine ppm_2d

  ! characteristics based on umac
  subroutine ppm_fpu_2d(s,ng_s,umac,ng_u,slx,sly,lo,hi,bc)

    use bc_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_u
    real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real(kind=dp_t), intent(  out) ::  slx(lo(1)-1   :,lo(2)-1   :,:) 
    real(kind=dp_t), intent(  out) ::  sly(lo(1)-1   :,lo(2)-1   :,:) 
    integer        , intent(in   ) :: bc(:,:,:)

    ! local



  end subroutine ppm_fpu_2d

  ! characteristics based on u
  subroutine ppm_3d(s,ng_s,u,ng_u,slx,sly,slz,lo,hi,bc)

    use bc_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_u
    real(kind=dp_t), intent(in   ) ::   s(lo(1)-ng_s:,lo(2)-ng_s:,hi(3)-ng_s:)
    real(kind=dp_t), intent(in   ) ::   u(lo(1)-ng_u:,lo(2)-ng_u:,hi(3)-ng_u:,:)
    real(kind=dp_t), intent(  out) :: slx(lo(1)-1   :,lo(2)-1   :,hi(3)-1   :,:) 
    real(kind=dp_t), intent(  out) :: sly(lo(1)-1   :,lo(2)-1   :,hi(3)-1   :,:)  
    real(kind=dp_t), intent(  out) :: slz(lo(1)-1   :,lo(2)-1   :,hi(3)-1   :,:)
    integer        , intent(in   ) :: bc(:,:,:)

    ! local



  end subroutine ppm_3d

  ! characteristics based on umac
  subroutine ppm_fpu_3d(s,ng_s,umac,ng_u,slx,sly,slz,lo,hi,bc)

    use bc_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_u
    real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng_s:,lo(2)-ng_s:,hi(3)-ng_s:)
    real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:,lo(2)-ng_u:,hi(3)-ng_u:,:)
    real(kind=dp_t), intent(  out) ::  slx(lo(1)-1   :,lo(2)-1   :,hi(3)-1   :,:) 
    real(kind=dp_t), intent(  out) ::  sly(lo(1)-1   :,lo(2)-1   :,hi(3)-1   :,:)  
    real(kind=dp_t), intent(  out) ::  slz(lo(1)-1   :,lo(2)-1   :,hi(3)-1   :,:)
    integer        , intent(in   ) :: bc(:,:,:)

    ! local



  end subroutine ppm_fpu_3d


end module ppm_module
