module inlet_bc

   use bl_types,  only: dp_t
   implicit none

   real(kind=dp_t), parameter :: INLET_VN      =  0.0d0
   real(kind=dp_t), parameter :: INLET_VT      =  0.0d0
   real(kind=dp_t), parameter :: INLET_RHO     =  4.0e7
   real(kind=dp_t), parameter :: INLET_RHOH    =  1.8954332646e25
   real(kind=dp_t), parameter :: INLET_RHOBIGH = -2.8317120266e26
   real(kind=dp_t), parameter :: INLET_RHOC12  =  2.0e7
   real(kind=dp_t), parameter :: INLET_RHOO16  =  2.0e7
   real(kind=dp_t), parameter :: INLET_RHOMG24 =  0.0d0
   real(kind=dp_t), parameter :: INLET_TRA     =  0.0d0

end module inlet_bc
