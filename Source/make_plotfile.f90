module make_plotfile_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_boxarray_module

  implicit none

  private
  public :: get_plot_names, make_plotfile, make_filename

  integer, parameter :: MAX_FILENAME_LEN = 256
  public :: MAX_FILENAME_LEN

contains

  function make_filename(base_name, index) result (filename)

    character (len=*), intent(in) :: base_name
    integer, intent(in) :: index

    character (len=MAX_FILENAME_LEN) :: filename

    character(len=5) :: str_index
    character(len=6) :: str_index6
    character(len=7) :: str_index7
    
    if (index <= 99999) then
       write(unit=str_index,fmt='(i5.5)') index
       filename = trim(base_name) // str_index
    else if (index <= 999999) then
       write(unit=str_index6,fmt='(i6.6)') index
       filename = trim(base_name) // str_index6
    else
       write(unit=str_index7,fmt='(i7.7)') index
       filename = trim(base_name) // str_index7
    endif
  end function make_filename


  subroutine get_plot_names(p)

    use plot_variables_module
    use variables
    use network, only: nspec, short_spec_names
    use probin_module, only: plot_spec, plot_trac, plot_base, &
                             use_thermal_diffusion, plot_omegadot, plot_Hnuc, &
                             plot_Hext, plot_eta, plot_ad_excess, &
                             use_tfromp, plot_h_with_use_tfromp, plot_gpi, plot_cs, &
                             plot_sponge_fdamp, dm_in, use_particles, &
                             plot_processors, plot_pidivu
    use geometry, only: spherical

    type(plot_t), intent(inout) :: p

    ! Local variables
    integer :: comp

    allocate(p%names(p%n_plot_comps))
    
    if (p%icomp_vel > 0) then
       p%names(p%icomp_vel  ) = "x_vel"
       if (dm_in > 1) then
          p%names(p%icomp_vel+1) = "y_vel"
       end if
       if (dm_in > 2) then
          p%names(p%icomp_vel+2) = "z_vel"
       end if
    endif

    if (p%icomp_rho > 0) p%names(p%icomp_rho)  = "density"
    if (p%icomp_rhoh > 0) p%names(p%icomp_rhoh) = "rhoh"
    if (p%icomp_h > 0) p%names(p%icomp_h)    = "h"

    if (p%icomp_spec > 0) then
       do comp = 1, nspec
          p%names(p%icomp_spec+comp-1) = "X(" // trim(short_spec_names(comp)) // ")"
       end do
    end if

    if (p%n_species > 0) then
       do comp = 1, p%n_species
          p%names(p%ipf_spec(comp)) = "X(" // trim(short_spec_names(p%spec(comp))) // ")"
       enddo
    endif

    if (p%icomp_trac > 0) then
       do comp = 1, ntrac
          p%names(p%icomp_trac+comp-1) = "tracer"
       end do
    end if

    if (p%icomp_w0 > 0) then
       p%names(p%icomp_w0)   = "w0_x"
       if (dm_in > 1) p%names(p%icomp_w0+1) = "w0_y"
       if (dm_in > 2) p%names(p%icomp_w0+2) = "w0_z"
    endif

    if (p%icomp_divw0 > 0) p%names(p%icomp_divw0) = "divw0"
    if (p%icomp_rho0 > 0) p%names(p%icomp_rho0)  = "rho0"
    if (p%icomp_rhoh0 > 0) p%names(p%icomp_rhoh0) = "rhoh0"
    if (p%icomp_h0 > 0) p%names(p%icomp_h0)    = "h0"
    if (p%icomp_p0 > 0) p%names(p%icomp_p0)    = "p0"

    if (p%icomp_velr > 0) p%names(p%icomp_velr) = "radial_velocity"
    if (p%icomp_velc > 0) p%names(p%icomp_velc) = "circum_velocity"

    if (p%icomp_magvel > 0) p%names(p%icomp_magvel)      = "magvel"
    if (p%icomp_mom > 0) p%names(p%icomp_mom)         = "momentum"
    if (p%icomp_vort > 0) p%names(p%icomp_vort)        = "vort"
    if (p%icomp_s_cc > 0) p%names(p%icomp_s_cc)         = "S"
    if (p%icomp_rhopert > 0) p%names(p%icomp_rhopert)     = "rhopert"
    if (p%icomp_rhohpert > 0) p%names(p%icomp_rhohpert)    = "rhohpert"
  
    if (p%icomp_tfromp > 0) p%names(p%icomp_tfromp)      = "tfromp"
    if (p%icomp_tfromH > 0) p%names(p%icomp_tfromH)      = "tfromh"
    if (p%icomp_dT > 0) p%names(p%icomp_dT)          = "deltaT"
    if (p%icomp_dp > 0) p%names(p%icomp_dp)          = "deltap"

    if (p%icomp_tpert > 0) p%names(p%icomp_tpert)       = "tpert"
    if (p%icomp_machno > 0) p%names(p%icomp_machno)      = "Machnumber"
    if (p%icomp_cs > 0) p%names(p%icomp_cs)          = "soundspeed"

    if (p%icomp_dg > 0) p%names(p%icomp_dg)          = "deltagamma"
    if (p%icomp_entropy > 0) p%names(p%icomp_entropy)     = "entropy"
    if (p%icomp_entropypert > 0) p%names(p%icomp_entropypert) = "entropypert"
    if (p%icomp_sponge > 0) then
       if (plot_sponge_fdamp) then
          p%names(p%icomp_sponge)      = "sponge_fdamp"
       else
          p%names(p%icomp_sponge)      = "sponge"
       end if
    endif

    if (p%icomp_pi > 0) p%names(p%icomp_pi)          = "pi"
    if (p%icomp_gpi > 0) then
       p%names(p%icomp_gpi)         = "gpi_x"
       if (dm_in > 1) p%names(p%icomp_gpi+1) = "gpi_y"
       if (dm_in > 2) p%names(p%icomp_gpi+2) = "gpi_z"
    endif

    if (p%icomp_pioverp0 > 0) p%names(p%icomp_pioverp0)    = "pioverp0"
    if (p%icomp_p0pluspi > 0) p%names(p%icomp_p0pluspi)    = "p0pluspi"

    if (p%icomp_omegadot > 0) then
       do comp = 1, nspec
          p%names(p%icomp_omegadot+comp-1) = &
               "omegadot(" // trim(short_spec_names(comp)) // ")"
       end do
    end if

    if (p%icomp_enuc > 0) p%names(p%icomp_enuc) = "enucdot"
    if (p%icomp_Hext > 0) p%names(p%icomp_Hext) = "Hext"

    if (p%icomp_eta > 0) p%names(p%icomp_eta) = "eta_rho"

    if (p%icomp_thermal > 0) p%names(p%icomp_thermal) = "thermal"
    if (p%icomp_conductivity > 0) p%names(p%icomp_conductivity) = "conductivity"

    if (p%icomp_ad_excess > 0) p%names(p%icomp_ad_excess) = "ad_excess"

    if (p%icomp_part > 0) p%names(p%icomp_part) = "particle_count"

    if (p%icomp_proc > 0) p%names(p%icomp_proc) = "processor_number"

    if (p%icomp_pidivu > 0) p%names(p%icomp_pidivu) = "pi_divu"

  end subroutine get_plot_names

  subroutine make_plotfile(p,dirname,mla,u,s,pi,gpi,rho_omegadot, &
                           rho_Hnuc,rho_Hext, &
                           thermal,S_cc,sponge,mba,dx, &
                           the_bc_tower,w0,rho0,rhoh0,p0, &
                           tempbar,gamma1bar,etarho_cc, &
                           normal,dt,particles,write_pf_time)

    use bl_prof_module
    use fabio_module
    use variables
    use plot_variables_module
    use fill_3d_module
    use probin_module, only: nOutFiles, lUsingNFiles, plot_spec, plot_trac, & 
                             plot_base, plot_omegadot, plot_Hnuc, plot_Hext, &
                             plot_eta, plot_ad_excess, &
                             single_prec_plotfiles, &
                             do_smallscale, use_thermal_diffusion, &
                             evolve_base_state, prob_lo, prob_hi, &
                             use_tfromp, plot_h_with_use_tfromp, plot_gpi, &
                             plot_cs, sponge_kappa, plot_sponge_fdamp, use_particles, &
                             plot_processors, plot_pidivu, use_alt_energy_fix
    use geometry, only: spherical, nr_fine, nlevs_radial, numdisjointchunks, &
         r_start_coord, r_end_coord
    use average_module
    use ml_restrict_fill_module
    use bl_constants_module
    use network, only: nspec
    use time_module, only: time
    use particle_module, only: particle_container, make_particle_count
    use make_grav_module
    use make_div_coeff_module
    use make_pi_cc_module

    type(plot_t)     , intent(in   ) :: p
    character(len=*) , intent(in   ) :: dirname
    type(ml_layout)  , intent(in   ) :: mla
    type(multifab)   , intent(in   ) :: u(:)
    type(multifab)   , intent(in   ) :: s(:)
    type(multifab)   , intent(in   ) :: pi(:)
    type(multifab)   , intent(in   ) :: gpi(:)
    type(multifab)   , intent(in   ) :: rho_omegadot(:)
    type(multifab)   , intent(in   ) :: rho_Hnuc(:)
    type(multifab)   , intent(in   ) :: rho_Hext(:)
    type(multifab)   , intent(in   ) :: thermal(:)
    type(multifab)   , intent(in   ) :: S_cc(:)
    type(multifab)   , intent(in   ) :: sponge(:)
    type(ml_boxarray), intent(in   ) :: mba
    real(dp_t)       , intent(in   ) :: dt,dx(:,:)
    type(bc_tower)   , intent(in   ) :: the_bc_tower
    real(dp_t)       , intent(in   ) :: w0(:,0:)
    real(dp_t)       , intent(in   ) :: rho0(:,0:)
    real(dp_t)       , intent(in   ) :: rhoh0(:,0:)
    real(dp_t)       , intent(in   ) :: p0(:,0:)
    real(dp_t)       , intent(in   ) :: tempbar(:,0:)
    real(dp_t)       , intent(in   ) :: gamma1bar(:,0:)
    real(dp_t)       , intent(in   ) :: etarho_cc(:,0:)
    type(multifab)   , intent(in   ) :: normal(:)
    type(particle_container), intent(inout) :: particles
    real(dp_t)       , intent(  out) :: write_pf_time

    type(multifab) :: plotdata(mla%nlevel)
    type(multifab) ::  tempfab(mla%nlevel)
    type(multifab) ::    w0mac(mla%nlevel,mla%dim)
    type(multifab) :: w0r_cart(mla%nlevel)
    type(multifab) ::    pi_cc(mla%nlevel)

    real(dp_t)  :: div_coeff(nlevs_radial,0:nr_fine-1)
    real(dp_t)  :: grav_cell(nlevs_radial,0:nr_fine-1)

    real(dp_t) :: entropybar(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::         h0(nlevs_radial,0:nr_fine-1)

    real(dp_t) :: tempval

    integer :: n,r,j,n_1d,prec,comp,dm,nlevs

    type(bl_prof_timer), save :: bpt

    real(dp_t) :: writetime1, writetime2

    call build(bpt, "make_plotfile")

    dm = mla%dim
    nlevs = mla%nlevel

    if (single_prec_plotfiles) then
       prec = FABIO_SINGLE
    else
       prec = FABIO_DOUBLE
    endif

    do n = 1,nlevs

       call multifab_build(plotdata(n), mla%la(n), p%n_plot_comps, 0)

       ! for temporary storage, as needed
       call multifab_build(tempfab(n),  mla%la(n), dm,           1)
              
       ! velocity
       if (p%icomp_vel > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_vel,u(n),1,dm)
       endif

       ! density, (rho h), and h
       if (p%icomp_rho > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_rho,s(n),rho_comp,1)
       endif

       if (p%icomp_rhoh > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_rhoh,s(n),rhoh_comp,1)
       endif

       if (p%icomp_h > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_h,   s(n),rhoh_comp,1)
          call multifab_div_div_c(plotdata(n),p%icomp_h,s(n),rho_comp,1)
       end if


       ! rhopert and rhohpert
       if (spherical .eq. 1) then
          if (p%icomp_rhopert > 0) then
             call make_rhopert( plotdata(n),p%icomp_rhopert, s(n), rho0(1,:),dx(n,:))
          endif
          if (p%icomp_rhohpert > 0) then
             call make_rhohpert(plotdata(n),p%icomp_rhohpert,s(n),rhoh0(1,:),dx(n,:))
          endif
       else
          if (p%icomp_rhopert > 0) then
             call make_rhopert( plotdata(n),p%icomp_rhopert, s(n), rho0(n,:),dx(n,:))
          endif
          if (p%icomp_rhohpert > 0) then
             call make_rhohpert(plotdata(n),p%icomp_rhohpert,s(n),rhoh0(n,:),dx(n,:))
          endif
       endif


       ! species
       if (p%icomp_spec > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_spec,s(n),spec_comp,nspec)
          do comp=1,nspec
             call multifab_div_div_c(plotdata(n),p%icomp_spec+comp-1,s(n),rho_comp,1)
          end do
       endif

       ! individual species -- if specified by name for a mini plotfile
       if (p%n_species > 0) then
          do comp = 1, p%n_species
             call multifab_copy_c(plotdata(n),p%ipf_spec(comp), &
                                  s(n),spec_comp+p%spec(comp)-1,1)
             call multifab_div_div_c(plotdata(n),p%ipf_spec(comp),s(n),rho_comp,1)
          enddo
       endif

       ! omegadot
       if (p%icomp_omegadot > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_omegadot,rho_omegadot(n),1,nspec)
          do comp=1,nspec
             call multifab_div_div_c(plotdata(n),p%icomp_omegadot+comp-1,s(n),rho_comp,1)
          end do
       end if


       ! enucdot
       if (p%icomp_enuc > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_enuc,rho_Hnuc(n),1)
          call multifab_div_div_c(plotdata(n),p%icomp_enuc,s(n),rho_comp,1)
       end if

       if (p%icomp_Hext > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_Hext,rho_Hext(n),1)
          call multifab_div_div_c(plotdata(n),p%icomp_Hext,s(n),rho_comp,1)
       end if

       ! thermal = del dot kappa grad T
       if (p%icomp_thermal > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_thermal,thermal(n),1)
       endif

       ! tracer
       if (p%icomp_trac > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_trac,s(n),trac_comp,ntrac)
       end if

    end do

    if (spherical .eq. 1) then

       do n=1,nlevs

          do comp=1,dm
             ! w0mac will contain an edge-centered w0 on a Cartesian grid,
             ! for use in computing divergences.
             call multifab_build_edge(w0mac(n,comp), mla%la(n),1,1,comp)
             call setval(w0mac(n,comp), ZERO, all=.true.)
          enddo

          ! w0r_cart is w0 but onto a Cartesian grid in cell-centered as
          ! a scalar.  Since w0 is the radial expansion velocity, w0r_cart
          ! is the radial w0 in a zone
          call multifab_build(w0r_cart(n), mla%la(n),1,0)
          call setval(w0r_cart(n), ZERO, all=.true.)

       end do

       if (evolve_base_state) then
          ! put w0 on Cartesian edges as a vector
          call make_w0mac(mla,w0,w0mac,dx,the_bc_tower%bc_tower_array)

          ! put w0 in Cartesian cell-centers as a scalar (the radial
          ! expansion velocity)
          call put_1d_array_on_cart(w0,w0r_cart,1,.true.,.false.,dx, &
                                    the_bc_tower%bc_tower_array,mla)
       end if

    end if  ! spherical

    ! w0
    if (p%icomp_w0 > 0) then
       ! this puts w0 onto Cartsian cell-centers as a vector, so we
       ! have all components here
       if (evolve_base_state) then
          call put_1d_array_on_cart(w0,tempfab,1,.true.,.true.,dx, &
                                    the_bc_tower%bc_tower_array,mla)
       else
          do n=1,nlevs
             call setval(tempfab(n), ZERO, all=.true.)
          end do
       end if

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),p%icomp_w0,tempfab(n),1,dm)
       end do
    endif

    ! divw0
    if (p%icomp_divw0 > 0) then
       do n=1,nlevs
          if (spherical .eq. 1) then
             n_1d = 1
          else
             n_1d = n
          end if
          call make_divw0(plotdata(n),p%icomp_divw0,w0(n_1d,:),w0mac(n,:),dx(n,:))
       end do
    endif

    ! rho0
    if (p%icomp_rho0 > 0) then
       call put_1d_array_on_cart(rho0,tempfab,dm+rho_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),p%icomp_rho0,tempfab(n),1,1)
       end do
    endif

    ! rhoh0
    if (p%icomp_rhoh0 > 0) then
       call put_1d_array_on_cart(rhoh0,tempfab,dm+rhoh_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),p%icomp_rhoh0,tempfab(n),1,1)
       end do
    endif

    ! h0
    if (p%icomp_h0 > 0) then
       if (do_smallscale) then
          h0 = ZERO
       else
          if (spherical .eq. 1) then
             ! only one level of base state so rho0 is defined everywhere
             h0 = rhoh0 / rho0
          else
             ! for Cartesian, rho0 will be zero in index locations where 
             ! there is no grid at that resolution.  Prevent dividing
             ! by zero
             h0(:,:) = ZERO

             do n = 1, nlevs
                do j = 1,numdisjointchunks(n)
                   do r = r_start_coord(n,j), r_end_coord(n,j)
                      h0(n,r) = rhoh0(n,r)/rho0(n,r)
                   enddo
                enddo
             enddo
          endif

       end if

       call put_1d_array_on_cart(h0,tempfab,foextrap_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),p%icomp_h0,tempfab(n),1,1)
       end do
    endif

    ! p0
    if (p%icomp_p0 > 0) then
       call put_1d_array_on_cart(p0,tempfab,foextrap_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)
       do n=1,nlevs
          call multifab_copy_c(plotdata(n),p%icomp_p0,tempfab(n),1,1)
       end do

    end if


    if (p%icomp_eta > 0) then
       call put_1d_array_on_cart(etarho_cc,tempfab,foextrap_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)
       do n=1,nlevs
          call multifab_copy_c(plotdata(n),p%icomp_eta,tempfab(n),1,1)
       end do
    endif


    do n = 1,nlevs

       if (spherical .eq. 1) then
          n_1d = 1
       else
          n_1d = n
       end if

       ! RADIAL AND CIRCUMFERENTIAL VELOCITY (spherical only)
       if (spherical .eq. 1) then
          if (p%icomp_velr > 0 .or. p%icomp_velc > 0) then
             call make_velrc(plotdata(n),p%icomp_velr,p%icomp_velc, &
                             u(n),w0r_cart(n),normal(n))
          endif
       endif

       ! MAGVEL = |U + w0|
       if (p%icomp_magvel > 0 .or. p%icomp_mom > 0) then
          call make_magvel(plotdata(n),p%icomp_magvel,p%icomp_mom, &
                           s(n),u(n),w0(n_1d,:),w0mac(n,:))
       endif

       ! VORTICITY
       if (p%icomp_vort > 0) then
          call make_vorticity(plotdata(n),p%icomp_vort,u(n),dx(n,:), &
                              the_bc_tower%bc_tower_array(n))
       endif

       ! DIVU
       if (p%icomp_s_cc > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_s_cc,S_cc(n),1,1)
       endif

    end do

    do n=1,nlevs

       ! make_tfromp -> TFROMP, TPERT, MACHNUMBER, CS, DELTAGAMMA, and ENTROPY
       ! make_tfromH -> TFROMP AND DELTA_P
       if (spherical .eq. 1) then

          if ( p%icomp_tfromp > 0 .or. p%icomp_tpert > 0 .or. &
               p%icomp_machno > 0 .or. p%icomp_cs > 0 .or. p%icomp_dg > 0 .or. &
               p%icomp_entropy > 0 .or. p%icomp_magvel > 0 ) then
             
             call make_tfromp(p, plotdata(n), s(n), &
                              tempbar(1,:),gamma1bar(1,:),p0(1,:),dx(n,:))
          endif

          if ( p%icomp_tfromH > 0 .or. p%icomp_dp > 0 ) then
             call make_tfromH(plotdata(n),p%icomp_tfromH,p%icomp_tpert, &
                              p%icomp_dp,s(n),p0(1,:), &
                              tempbar(1,:),dx(n,:))
          endif

       else

          if ( p%icomp_tfromp > 0 .or. p%icomp_tpert > 0 .or. &
               p%icomp_machno > 0 .or. p%icomp_cs > 0 .or. p%icomp_dg > 0 .or. &
               p%icomp_entropy > 0 .or. p%icomp_magvel > 0 ) then

             call make_tfromp(p, plotdata(n), s(n), &
                              tempbar(n,:),gamma1bar(n,:),p0(n,:),dx(n,:))
          endif

          if ( p%icomp_tfromH > 0 .or. p%icomp_dp > 0 ) then
             call make_tfromH(plotdata(n),p%icomp_tfromH,p%icomp_tpert, &
                              p%icomp_dp,s(n),p0(n,:), &
                              tempbar(n,:),dx(n,:))
          endif

       end if
       
    end do

    ! CONDUCTIVITY
    if (p%icomp_conductivity > 0) then
       do n=1,nlevs
          ! this just uses (rho, T, X_k) ---> conductivity
          ! it doesn't need to do anything fancy for spherical
          call make_conductivity(plotdata(n),p%icomp_conductivity,s(n))
       end do
    end if

    ! ADIABATIC EXCESS
    if (p%icomp_ad_excess > 0) then
       do n = 1, nlevs
          call make_ad_excess(plotdata(n),p%icomp_ad_excess,s(n),normal(n))
       enddo
    endif


    ! PARTICLES
    if (p%icomp_part > 0) then
       do n = 1, nlevs
          call multifab_setval_c(plotdata(n),ZERO,p%icomp_part,1)
       enddo
       call make_particle_count(mla,plotdata,p%icomp_part,particles)
    endif


    ! processor number
    if (p%icomp_proc > 0) then
       do n = 1, nlevs
          call make_processor_number(plotdata(n),p%icomp_proc)
       enddo
    endif
    

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       if (p%icomp_tfromp > 0) then
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_tfromp,plotdata(n),p%icomp_tfromp, &
                                   mla%mba%rr(n-1,:),1)
       endif
       if (p%icomp_tpert > 0) then
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_tpert,plotdata(n),p%icomp_tpert, &
                                   mla%mba%rr(n-1,:),1)
       endif
       if (p%icomp_rhopert > 0) then
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_rhopert,plotdata(n),p%icomp_rhopert, &
                                   mla%mba%rr(n-1,:),1)
       endif
       if (p%icomp_rhohpert > 0) then
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_rhohpert,plotdata(n),p%icomp_rhohpert, &
                                   mla%mba%rr(n-1,:),1)
       endif
       if (p%icomp_tfromH > 0) then
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_tfromH,plotdata(n),p%icomp_tfromH, &
                                   mla%mba%rr(n-1,:),1)
       endif
       if (p%icomp_dp > 0) then
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_dp,plotdata(n),p%icomp_dp, &
                                   mla%mba%rr(n-1,:),1)
       endif
       if (p%icomp_cs > 0) then
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_cs,plotdata(n),p%icomp_cs, &
                                   mla%mba%rr(n-1,:),1)
       end if
       if (p%icomp_machno > 0) then
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_machno,plotdata(n),p%icomp_machno, &
                                   mla%mba%rr(n-1,:),1)
       endif
       if (p%icomp_dg > 0) then
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_dg,plotdata(n),p%icomp_dg, &
                                   mla%mba%rr(n-1,:),1)
       endif
       if (p%icomp_entropy > 0) then
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_entropy,plotdata(n),p%icomp_entropy, &
                                   mla%mba%rr(n-1,:),1)
       endif
    end do

    ! build a cell-centered multifab to hold pi
    do n=1,nlevs
       call multifab_build(pi_cc(n), mla%la(n), 1, 0)
       call setval(pi_cc(n), ZERO, all=.true.)
    end do


    if (use_alt_energy_fix) then
       ! make beta_0 on a multifab
       call make_grav_cell(grav_cell,rho0)
       call make_div_coeff(div_coeff,rho0,p0,gamma1bar,grav_cell)
       
       call put_1d_array_on_cart(div_coeff,tempfab,foextrap_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)
    endif

    ! new function that average the nodal pi to cell-centers, then
    ! normalized the entire signal to sum to 0.  If we are doing
    ! use_alt_energy_fix = T, then we first want to convert
    ! (pi/beta_0) to pi.
    call make_pi_cc(mla,pi,pi_cc,1,the_bc_tower%bc_tower_array,tempfab)

    do n=1,nlevs

       ! DELTA_T
       if (p%icomp_dT > 0 .and. p%icomp_tfromp > 0 .and. p%icomp_tfromH > 0) then
          call make_deltaT(plotdata(n),p%icomp_dT,p%icomp_tfromp,p%icomp_tfromH)
       endif

       ! PI
       if (p%icomp_pi > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_pi,pi_cc(n),1,1)
       endif

       ! GRAD PI
       if (p%icomp_gpi > 0) then
          call multifab_copy_c(plotdata(n),p%icomp_gpi,gpi(n),1,dm)
       endif

       ! pi * div(U)
       if (p%icomp_pidivu > 0) then
          call make_pidivu(plotdata(n),p%icomp_pidivu,pi_cc(n),u(n),dx(n,:))
       endif

    end do


    ! SPONGE
    if (p%icomp_sponge > 0) then
       if (plot_sponge_fdamp) then
          tempval = dt*sponge_kappa
          do n=1,nlevs
             ! compute f_damp assuming sponge=1/(1+dt*kappa*fdamp)
             ! therefore fdamp = (1/sponge-1)/(dt*kappa)

             ! plotdata = 1
             call multifab_setval_c(plotdata(n),ONE,p%icomp_sponge,1)
             ! tempfab = 1
             call multifab_setval(tempfab(n),ONE)
             ! plotdata = 1/sponge
             call multifab_div_div_c(plotdata(n),p%icomp_sponge,sponge(n),1,1)
             ! plotdata = 1/sponge-1
             call multifab_sub_sub_c(plotdata(n),p%icomp_sponge,tempfab(n),1,1)
             ! plotdata = (1/sponge-1)/(dt*kappa)
             call multifab_div_div_s_c(plotdata(n),p%icomp_sponge,tempval,1)
          end do
       else
          do n=1,nlevs
             call multifab_copy_c(plotdata(n),p%icomp_sponge,sponge(n),1,1)
          end do
       end if
    endif

    !PIOVERP0 and P0PLUSPI
    if (plot_base) then
       do n=1,nlevs
          if (p%icomp_pioverp0 > 0 .and. p%icomp_p0 > 0) then
             call multifab_copy_c(plotdata(n),p%icomp_pioverp0,pi_cc(n),1,1)
             call multifab_div_div_c(plotdata(n),p%icomp_pioverp0,plotdata(n),p%icomp_p0,1)
          endif
          if (p%icomp_p0pluspi > 0 .and. p%icomp_p0 > 0) then
             call multifab_copy_c(plotdata(n),p%icomp_p0pluspi,pi_cc(n),1,1)
             call multifab_plus_plus_c(plotdata(n),p%icomp_p0pluspi,plotdata(n),p%icomp_p0,1)
          endif
       end do
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! we just made the entropy above.  To compute s - sbar, we need to average
    ! the entropy first, and then compute that.
    ! an average quantity needs ghostcells, so copy entropy into
    ! tempfab
    if (p%icomp_entropy > 0 .and. p%icomp_entropypert > 0) then
       do n=1,nlevs
          call multifab_copy_c(tempfab(n),1,plotdata(n),p%icomp_entropy,1)
       end do

       ! restrict data and fill all ghost cells
       call ml_restrict_and_fill(nlevs,tempfab,mla%mba%rr,the_bc_tower%bc_tower_array, &
                                 icomp=1, &
                                 bcomp=foextrap_comp, &
                                 nc=1, &
                                 ng=tempfab(1)%ng)
    
       call average(mla,tempfab,entropybar,dx,1)

       call make_entropypert(plotdata,p%icomp_entropy,p%icomp_entropypert,entropybar,dx)

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(plotdata(n-1),p%icomp_entropypert,plotdata(n), &
                                   p%icomp_entropypert,mla%mba%rr(n-1,:),1)
       end do
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (parallel_IOProcessor()) then
       write(6,*) 'Writing state to plotfile ',trim(dirname)
    end if

    writetime1 = parallel_wtime()

    call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), dirname, p%names, &
                                   mba%pd(1), prob_lo, prob_hi, time, dx(1,:), &
                                   nOutFiles = nOutFiles, &
                                   lUsingNFiles = lUsingNFiles, prec = prec)

    writetime2 = parallel_wtime() - writetime1
    call parallel_reduce(writetime1, writetime2, MPI_MAX, proc=parallel_IOProcessorNode())
    if (parallel_IOProcessor()) then
       print*,'Time to write plotfile: ',writetime1,' seconds'
       print*,''
    end if

    write_pf_time = writetime1

    do n = 1,nlevs
       call destroy(plotdata(n))
       call destroy(tempfab(n))
       call destroy(pi_cc(n))
    end do

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(w0r_cart(n))
          do comp=1,dm
             call destroy(w0mac(n,comp))
          end do
       end do
    end if

    call destroy(bpt)

  end subroutine make_plotfile

end module make_plotfile_module
