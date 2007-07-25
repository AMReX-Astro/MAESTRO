module inlet_bc

   use bl_types,  only: dp_t
   implicit none

   real(kind=dp_t), parameter :: INLET_VX  =  1.0d0
   real(kind=dp_t), parameter :: INLET_VY  =  0.0d0
   real(kind=dp_t), parameter :: INLET_VZ  =  0.0d0
   real(kind=dp_t), parameter :: INLET_DEN =  1.0d0
   real(kind=dp_t), parameter :: INLET_X1  =  0.5d0
   real(kind=dp_t), parameter :: INLET_X2  =  0.0d0
   real(kind=dp_t), parameter :: INLET_X3  =  0.5d0
   real(kind=dp_t), parameter :: INLET_TRA =  1.0d0

end module inlet_bc
