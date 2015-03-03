f90sources += main.f90
f90sources += pred_parameters.f90
f90sources += probin.f90
f90sources += varden.f90
# |
# --> 
 f90sources += varden_aux.f90
f90sources += variables.f90
f90sources += geometry.f90

f90sources += average.f90
# |
# --> 
 f90sources += restrict_base.f90
 

f90sources += react_state.f90
# |
# --> 
f90sources += define_bc_tower.f90
f90sources += rhoh_vs_t.f90
f90sources += make_heating.f90
  # |
  # ----> 
  f90sources += fill_3d_data.f90
f90sources += multifab_fill_ghost_cells.f90
  # |
  # ----> 
  f90sources += fillpatch.f90
    # |
    # ------> 
    f90sources += multifab_physbc.f90
f90sources += build_info.f90
f90sources += write_job_info.f90
f90sources += mg_eps.f90
f90sources += cputime.f90
f90sources += constants_cgs.f90 
f90sources += simple_log.f90

# To be organized
f90sources += convert_rhoX_to_X.f90
f90sources += put_in_pert_form.f90

