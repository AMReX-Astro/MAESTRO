module inlet_bc

   use bl_types,  only: dp_t
   implicit none

   real(kind=dp_t), parameter :: INLET_VN   =  0.0d0
   real(kind=dp_t), parameter :: INLET_VT   =  0.0d0
   real(kind=dp_t), parameter :: INLET_RHO  =  4.0e7
   real(kind=dp_t), parameter :: INLET_RHOH =  3.009433262473e25
   real(kind=dp_t), parameter :: INLET_C12  =  0.5d0
   real(kind=dp_t), parameter :: INLET_O16  =  0.5d0
   real(kind=dp_t), parameter :: INLET_MG24 =  0.0d0
   real(kind=dp_t), parameter :: INLET_TRA  =  0.0d0

end module inlet_bc
