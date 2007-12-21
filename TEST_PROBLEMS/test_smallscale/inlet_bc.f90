module inlet_bc_module

  use bl_types
  use bl_space

  implicit none

  real(dp_t), save    :: INLET_VN
  real(dp_t), save    :: INLET_VT
  real(dp_t), save    :: INLET_RHO
  real(dp_t), save    :: INLET_RHOH
  real(dp_t), save    :: INLET_RHOC12
  real(dp_t), save    :: INLET_RHOO16
  real(dp_t), save    :: INLET_RHOMG24
  real(dp_t), save    :: INLET_TEMP
  real(dp_t), save    :: INLET_TRA

  namelist /inlet_bc/ INLET_VN
  namelist /inlet_bc/ INLET_VT
  namelist /inlet_bc/ INLET_RHO
  namelist /inlet_bc/ INLET_RHOH
  namelist /inlet_bc/ INLET_RHOC12
  namelist /inlet_bc/ INLET_RHOO16
  namelist /inlet_bc/ INLET_RHOMG24
  namelist /inlet_bc/ INLET_TEMP
  namelist /inlet_bc/ INLET_TRA

end module inlet_bc_module
