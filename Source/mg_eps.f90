! set the tolerances used in the various multigrid solves throughout
! the algorithm.  This is described in MAESTRO/docs/mg/.

module mg_eps_module

  use bl_types
  implicit none

  ! tolerances for the initial projection
  real (kind=dp_t) :: eps_init_proj_cart = 1.d-12
  real (kind=dp_t) :: eps_init_proj_sph  = 1.d-10


  ! tolerances for the divu iterations
  real (kind=dp_t) :: eps_divu_cart = 1.d-12
  real (kind=dp_t) :: eps_divu_sph  = 1.d-10

  real (kind=dp_t) :: divu_iter_factor = 100.d0
  real (kind=dp_t) :: divu_level_factor = 10.d0


  ! tolerances for the MAC projection
  real (kind=dp_t) :: eps_mac = 1.d-10
  real (kind=dp_t) :: eps_mac_bottom = 1.d-3


  ! tolerances for the HG projection
  real (kind=dp_t) :: eps_hg = 1.d-12

  real (kind=dp_t) :: eps_hg_min = 1.d-10
  real (kind=dp_t) :: hg_level_factor = 10.d0

end module mg_eps_module
