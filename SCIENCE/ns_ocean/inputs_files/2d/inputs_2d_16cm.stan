&PROBIN
 model_file = "Xtable_lo_16cm.hse"
 spherical_in = 0
 plot_base = T

 max_levs = 1
 n_cellx = 1024
 n_celly = 1024

! 4x16 = 64
 max_grid_size = 128

 ppm_type = 1

 do_burning = F

 mg_bottom_solver = 4
 hg_bottom_solver = 4
 max_mg_bottom_nlevels = 4

 burner_threshold_species = "He4"
 
 use_eos_coulomb = T

 prob_lo_x = 0.0d0
 prob_lo_y = 7616.0d0

 prob_hi_x = 16384.0d0
 prob_hi_y = 24000.0d0

 grav_const = -2.0e14

 max_step  = 1000
 init_iter = 3

 stop_time = 2.5

 plot_int  = -1
 plot_deltat = 1.e-4
 chk_int   = 200

 cflfac = 0.7
 init_shrink = 0.1
 max_dt_growth = 1.1
 use_soundspeed_firstdt = T
 use_divu_firstdt = F

 bcx_lo = -1
 bcx_hi = -1
 bcy_lo = 14
 bcy_hi = 12


    verbose = 1
 mg_verbose = 0
 cg_verbose = 0

 do_initial_projection = T
 init_divu_iter = 3

 use_thermal_diffusion = T
 dpdt_factor = 0.3d0

 base_cutoff_density = 1.d4
 anelastic_cutoff = 1.d4

 perturb_model = F
 xrb_pert_size = 50
 xrb_pert_factor = 1.0d-5
 xrb_pert_type = 1

 do_sponge = F
 xrb_use_bottom_sponge = F

!restart = 50

 apply_vel_field    = F
 velpert_scale      = 1.0d2
 velpert_amplitude  = 1.0d2
 velpert_height_loc = 6.5d3
 num_vortices       = 1

 limit_conductivity = F

 xrb_pert_height = 6160.d0
/
