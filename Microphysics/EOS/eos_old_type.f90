! The old EOS interface had the variables passed through the eos()
! interface declared in a module for convenience.  Now if we want to
! use that old interface, we need to include this module.

module eos_old_interface

  use bl_types
  use network, only: nspec
  use bl_space

  real(kind=dp_t), public :: xn_eos(nspec)
  real(kind=dp_t), public :: temp_eos
  real(kind=dp_t), public :: den_eos
  real(kind=dp_t), public :: abar_eos
  real(kind=dp_t), public :: zbar_eos
  real(kind=dp_t), public :: e_eos
  real(kind=dp_t), public :: p_eos
  real(kind=dp_t), public :: h_eos
  real(kind=dp_t), public :: cv_eos
  real(kind=dp_t), public :: cp_eos
  real(kind=dp_t), public :: xne_eos
  real(kind=dp_t), public :: eta_eos
  real(kind=dp_t), public :: pele_eos
  real(kind=dp_t), public :: dpdt_eos
  real(kind=dp_t), public :: dpdr_eos
  real(kind=dp_t), public :: dedr_eos
  real(kind=dp_t), public :: dedt_eos
  real(kind=dp_t), public :: gam1_eos
  real(kind=dp_t), public ::   cs_eos
  real(kind=dp_t), public ::    s_eos
  real(kind=dp_t), public :: dsdt_eos
  real(kind=dp_t), public :: dsdr_eos
  real(kind=dp_t), public :: dpdX_eos(nspec)
  real(kind=dp_t), public :: dhdX_eos(nspec)
  real(kind=dp_t), public :: conduct_eos

  integer, public         :: pt_index_eos(MAX_SPACEDIM)

  common /eos_common/ xn_eos,temp_eos,den_eos,abar_eos,zbar_eos,e_eos,p_eos,h_eos
  common /eos_common/ cv_eos,cp_eos,xne_eos,eta_eos,pele_eos,dpdt_eos,dpdr_eos,dedr_eos
  common /eos_common/ dedt_eos,gam1_eos,cs_eos,s_eos,dsdt_eos,dsdr_eos,dpdX_eos,dhdX_eos
  common /eos_common/ conduct_eos,pt_index_eos
  SAVE /eos_common/
!$omp threadprivate(/eos_common/)

end module eos_old_interface

