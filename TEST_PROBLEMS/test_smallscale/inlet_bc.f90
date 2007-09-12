module inlet_bc

   use bl_types,  only: dp_t
   implicit none

   real(kind=dp_t), parameter :: INLET_VN      =  0.0d0
   real(kind=dp_t), parameter :: INLET_VT      =  0.0d0
   real(kind=dp_t), parameter :: INLET_RHO     =  3.9999999999798551d7
   real(kind=dp_t), parameter :: INLET_RHOH    =  1.89543326464019d25
!  This is for the ebin formulation
   real(kind=dp_t), parameter :: INLET_RHOBIGH = -2.83171202661189d26
!  This is for the qs formulation
!   real(kind=dp_t), parameter :: INLET_RHOBIGH =  3.00943325470894d25
   real(kind=dp_t), parameter :: INLET_RHOC12  =  0.5d0*3.9999999999798551d7
   real(kind=dp_t), parameter :: INLET_RHOO16  =  0.5d0*3.9999999999798551d7
   real(kind=dp_t), parameter :: INLET_RHOMG24 =  0.0d0
   real(kind=dp_t), parameter :: INLET_TEMP    =  1.00000000047000d8
   real(kind=dp_t), parameter :: INLET_TRA     =  0.0d0

end module inlet_bc
