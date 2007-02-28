module burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network

contains

  subroutine burner(dens, temp, Xin, hin, dt, Xout, hout, rho_omegadot)


    implicit none
    
    real(kind=dp_t), intent(in) :: dens, temp, Xin(nspec), hin, dt
    real(kind=dp_t), intent(out) :: Xout(nspec), hout, rho_omegadot(nspec)
    
    integer :: n
    real(kind=dp_t) :: enuc, dX
    
    Xout(:) = Xin(:)
    
    enuc = 0.0_dp_t
    do n = 1, nspec
       dX = Xout(n)-Xin(n) 
       
       enuc = enuc - ebin(n) * dX
       
       rho_omegadot(n) = dens * dX / dt
    enddo
  
    hout = hin + enuc
  
  end subroutine burner

end module burner_module
