module init_particles_module
  
  use bl_types
  use bl_constants_module
  use multifab_module
  use particle_module

  implicit none

  private

  public :: init_particles

contains

  subroutine init_particles(particles,s,rho0,rhoh0,p0,tempbar,mla,dx)

    use bl_prof_module
    use geometry, only: spherical, nlevs_radial, nr_fine
    use variables, only: rho_comp, temp_comp
    use fill_3d_module

    type(particle_container), intent(inout) :: particles
    type(multifab) ,          intent(in   ) :: s(:)
    real(kind=dp_t),          intent(in   ) :: rho0(:,0:)
    real(kind=dp_t),          intent(in   ) :: rhoh0(:,0:)
    real(kind=dp_t),          intent(in   ) :: p0(:,0:)
    real(kind=dp_t),          intent(in   ) :: tempbar(:,0:)
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)

    real(kind=dp_t), pointer::   sop(:,:,:,:)
    
    integer :: lo(mla%dim),hi(mla%dim),i,n,comp,dm,nlevs
    integer :: ng_s
    
    integer :: iseed

    type(bl_prof_timer), save :: bpt
    
    call build(bpt, "init_particles")
    
    dm = mla%dim
    nlevs = mla%nlevel
    
    ng_s  = nghost(s(1))

    do n=1,nlevs
       
       do i = 1, nboxes(s(n))
          if ( multifab_remote(s(n), i) ) cycle
          sop   => dataptr(s(n), i)
          
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          
          select case (dm)
          case (1)
             call init_particles_1d(n, sop(:,1,1,:), ng_s, &
                                    rho0(n,:), rhoh0(n,:), p0(n,:), tempbar(n,:), lo, hi, &
                                    particles, dx, mla)

          case (2)
             ! the seed for the random number generator -- we'll set
             ! it to the box number to be reproducable.
             iseed = i
             call init_particles_2d(n, sop(:,:,1,:), ng_s, &
                                    rho0(n,:), rhoh0(n,:), p0(n,:), tempbar(n,:), lo, hi, &
                                    particles, dx, mla, iseed)

          case (3)
             if (spherical .eq. 1) then
                call bl_error("ERROR: spherical init_particles not written")

             else
                call init_particles_3d_cart(n, sop(:,:,:,:), ng_s, &
                                            rho0(n,:), rhoh0(n,:), p0(n,:), tempbar(n,:), lo, hi, &
                                            particles, dx, mla)
             end if
          end select

       end do

    enddo

    call destroy(bpt)

  end subroutine init_particles

  subroutine init_particles_1d(n, s, ng_s, &
                               rho0, rhoh0, p0, tempbar, lo, hi, &
                               particles, dx, mla)

    use probin_module, only: prob_lo
    use variables, only: rho_comp

    integer,                  intent(in) :: n, lo(:), hi(:), ng_s
    real (kind = dp_t),       intent(in   ) ::     s(lo(1)-ng_s :,:)  
    real (kind = dp_t),       intent(in   ) :: rho0(0:), rhoh0(0:), p0(0:), tempbar(0:)
    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)    

    ! local variables
    integer :: i
    real(kind=dp_t) :: x
    real(kind=dp_t) :: point(1)


    do i = lo(1), hi(1)
       x = prob_lo(1) + (dble(i) + HALF) * dx(n,1)
       
       ! point(1) = x
       
       ! call add(particles,point,mla,dx,prob_lo)

    enddo

  end subroutine init_particles_1d
  
  subroutine init_particles_2d(n, s, ng_s, &
                               rho0, rhoh0, p0, tempbar, lo, hi, &
                               particles, dx, mla, iseed)

    use probin_module, only : prob_lo, &
                              particle_temp_threshold, particle_probability
    use variables, only: temp_comp
    use mt19937_module

    integer,                  intent(in) :: n, lo(:), hi(:), ng_s
    real (kind = dp_t),       intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,:)  
    real (kind = dp_t),       intent(in   ) :: rho0(0:), rhoh0(0:), p0(0:), tempbar(0:)
    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)    
    integer,                  intent(in   ) :: iseed

    ! local variables
    integer :: i, j
    real(kind=dp_t) :: x, y
    real(kind=dp_t) :: point(2)
    real(kind=dp_t) :: rnd

    call init_genrand(iseed)

    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j) + HALF) * dx(n,2)

       do i = lo(1), hi(1)
          x = prob_lo(1) + (dble(i) + HALF) * dx(n,1)
       
          ! add particles randomly if the temperature is above
          ! a threshold
          if (s(i,j,temp_comp) > particle_temp_threshold) then

             rnd = genrand_real3()
             if (rnd <= particle_probability) then
                point(1) = x
                point(2) = y

                call add(particles,point,mla,dx,prob_lo)
             endif
          endif

       enddo
    enddo

  end subroutine init_particles_2d
  
  subroutine init_particles_3d_cart(n, s, ng_s, &
                                    rho0, rhoh0, p0, tempbar, lo, hi, &
                                    particles, dx, mla)

    use probin_module, only: prob_lo
    use variables, only: rho_comp

    integer,                  intent(in) :: n, lo(:), hi(:), ng_s
    real (kind = dp_t),       intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)  
    real (kind = dp_t),       intent(in   ) :: rho0(0:), rhoh0(0:), p0(0:), tempbar(0:)
    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)    

    ! local variables
    integer :: i, j, k
    real(kind=dp_t) :: x, y, z
    real(kind=dp_t) :: point(3)

    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k) + HALF) * dx(n,3)

       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j) + HALF) * dx(n,2)

          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i) + HALF) * dx(n,1)

             ! point(1) = x
             ! point(2) = y
             
             ! call add(particles,point,mla,dx,prob_lo)
    
          enddo
       enddo
    enddo

  end subroutine init_particles_3d_cart
  
end module init_particles_module
