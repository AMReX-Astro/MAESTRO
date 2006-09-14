module eos_module

      use bl_types
      use network

      integer, parameter :: NP = 1
      integer, parameter :: npts = 1

      real(kind=dp_t) :: xmass(nspec)= (/  0.3_dp_t,  0.7_dp_t, 0.0_dp_t /)

      real(kind=dp_t) :: xn_zone(nspec)

      real(kind=dp_t) :: temp_row(NP)
      real(kind=dp_t) :: den_row(NP)
      real(kind=dp_t) :: abar_row(NP)
      real(kind=dp_t) :: zbar_row(NP)
      real(kind=dp_t) :: e_row(NP)
      real(kind=dp_t) :: p_row(NP)
      real(kind=dp_t) :: h_row(NP)
      real(kind=dp_t) :: cv_row(NP)
      real(kind=dp_t) :: cp_row(NP)
      real(kind=dp_t) :: xne_row(NP)
      real(kind=dp_t) :: eta_row(NP)
      real(kind=dp_t) :: pele_row(NP)
      real(kind=dp_t) :: dpdt_row(NP)
      real(kind=dp_t) :: dpdr_row(NP)
      real(kind=dp_t) :: dedr_row(NP)
      real(kind=dp_t) :: dedt_row(NP)
      real(kind=dp_t) :: gam1_row(NP)
      real(kind=dp_t) ::   cs_row(NP)
      real(kind=dp_t) ::    s_row(NP)
      integer :: input_flag
      logical :: do_diag

end module eos_module
