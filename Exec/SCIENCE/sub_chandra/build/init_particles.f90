! these routines are used to initialize and update the particles.
!
! init_mode = 1 means that we are doing the t = 0 initialization.
! init_mode = 2 is called during the normal evolution

module init_particles_module
  
  use bl_types
  use bl_constants_module
  use multifab_module
  use particle_module

  implicit none

  private

  public :: init_particles

contains

  subroutine init_particles(particles,s,rho0,rhoh0,p0,tempbar,mla,dx,init_mode)

     use bl_prof_module
     use geometry, only: spherical, center, nr_fine, r_cc_loc, anelastic_cutoff_coord
     use variables, only: spec_comp
     use network, only: network_species_index
     use probin_module, only: octant, prob_lo
     use fill_3d_module
     use average_module, only: average

     type(particle_container), intent(inout) :: particles
     type(multifab) ,          intent(in   ) :: s(:)
     real(kind=dp_t),          intent(in   ) :: rho0(:,0:)
     real(kind=dp_t),          intent(in   ) :: rhoh0(:,0:)
     real(kind=dp_t),          intent(in   ) :: p0(:,0:)
     real(kind=dp_t),          intent(in   ) :: tempbar(:,0:)
     type(ml_layout),          intent(inout) :: mla
     real(kind=dp_t),          intent(in   ) :: dx(:,:)
     integer        ,          intent(in   ) :: init_mode
    
     integer, parameter :: NRAYS   = 1   !Number of rays of particles
     integer, parameter :: NPERRAY = 10  !Number of particles per ray
     
     type(bl_prof_timer), save :: bpt

     real(kind=dp_t) :: theta, theta_max, theta_delta
     real(kind=dp_t) :: phi,   phi_max,   phi_delta
     real(kind=dp_t) :: r, r_delta, r_lo, r_hi, x, y, z
     real(kind=dp_t) :: xhe4bar(1,nr_fine), point(mla%dim)

     integer :: i, j, ipxhe4, xhe4_comp
     
     call build(bpt, "init_particles")

     !Do some consistency/error checking
     if (spherical /= 1 .or. mla%dim /= 3) then
        call bl_error('ERROR: Sub-Chandra is meant to be run in 3D spherical')
     endif
    
     !Average out helium mass fraction
     xhe4_comp = spec_comp - 1 + network_species_index('He4')
     call average(mla, s, xhe4bar, dx, xhe4_comp)

     !Determine index of radial bin with X_He = 0.1
     do i=0, nr_fine
        if (xhe4bar(1,i) > 0.1) then
           ipxhe4 = i-1
           exit
        endif
     enddo

     !Determine maximum angle based on geometry
     if (octant) then
        theta_max = M_PI/2.0
        phi_max   = M_PI/2.0
     else
        theta_max = 2.0*M_PI
        phi_max   = M_PI
     endif
   
     !Set lower radial bound to be where X_He4 ~= 0.1
     r_lo = r_cc_loc(1,ipxhe4)
     
     !Set upper radial bound to be the anelastic cutoff
     r_hi = r_cc_loc(1,anelastic_cutoff_coord(1))

     !Distribute particles along rays 
     theta_delta = theta_max / dble(NRAYS+1)
     phi_delta   = phi_max   / dble(NRAYS+1)
     r_delta = (r_hi - r_lo)/NPERRAY
     do i=1, NRAYS
        theta = i * theta_delta
        phi   = i * phi_delta
        do j=0, NPERRAY-1
           r = r_lo + j*r_delta

           x = r*cos(theta)*sin(phi)
           y = r*sin(theta)*sin(phi)
           z = r*cos(phi)

           point(1) = center(1) + x
           point(2) = center(2) + y
           point(3) = center(3) + z

           call add(particles,point,mla,dx,prob_lo)
        enddo
     enddo

     call destroy(bpt)

  end subroutine init_particles

end module init_particles_module
