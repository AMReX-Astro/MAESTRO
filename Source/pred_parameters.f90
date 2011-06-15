module pred_parameters

   implicit none

   ! species prediction
   integer, parameter :: predict_rhoprime_and_X = 1
   integer, parameter :: predict_rhoX           = 2
   integer, parameter :: predict_rho_and_X      = 3

   ! enthalpy prediction
   integer, parameter :: predict_rhohprime        = 1
   integer, parameter :: predict_h                = 2
   integer, parameter :: predict_T_then_rhohprime = 3
   integer, parameter :: predict_T_then_h         = 4
   integer, parameter :: predict_hprime           = 5
   integer, parameter :: predict_Tprime_then_h    = 6

end module pred_parameters
