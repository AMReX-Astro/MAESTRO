module heating_module

  use bl_types

  implicit none
  
  private
  public :: get_rho_Hext_2d, get_rho_Hext_3d
  
contains
  
  subroutine get_rho_Hext_2d(rho_Hext,s,lo,hi,ng,dx,time)
    
    use bl_constants_module
    
    integer, intent(in) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1):,lo(2):)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    rho_Hext = 0.0_dp_t
    
  end subroutine get_rho_Hext_2d
  
  subroutine get_rho_Hext_3d(rho_Hext,s,lo,hi,ng,dx,time)
    
    use bl_constants_module
    
    integer, intent(in) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1):,lo(2):,lo(3):)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    rho_Hext = 0.0_dp_t
    
  end subroutine get_rho_Hext_3d
  
end module heating_module
