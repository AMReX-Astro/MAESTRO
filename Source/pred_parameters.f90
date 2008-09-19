module pred_parameters

   implicit none

   integer, parameter :: predict_rhohprime        = 1
   integer, parameter :: predict_h                = 2
   integer, parameter :: predict_T_then_rhohprime = 3
   integer, parameter :: predict_T_then_h         = 4
   integer, parameter :: predict_hprime           = 5

end module pred_parameters
