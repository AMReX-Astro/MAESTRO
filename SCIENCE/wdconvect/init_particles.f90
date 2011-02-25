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
    use geometry, only: spherical, nlevs_radial, nr_fine
    use variables, only: rho_comp, temp_comp

    type(particle_container), intent(inout) :: particles
    type(multifab) ,          intent(in   ) :: s(:)
    real(kind=dp_t),          intent(in   ) :: rho0(:,0:)
    real(kind=dp_t),          intent(in   ) :: rhoh0(:,0:)
    real(kind=dp_t),          intent(in   ) :: p0(:,0:)
    real(kind=dp_t),          intent(in   ) :: tempbar(:,0:)
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)
    integer,                  intent(in   ) :: init_mode

    real(kind=dp_t), pointer::   sop(:,:,:,:)
    
    integer :: lo(mla%dim),hi(mla%dim),i,n,comp,dm,nlevs
    integer :: ng_s
    
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
                                    rho0(n,:), rhoh0(n,:), &
                                    p0(n,:), tempbar(n,:), lo, hi, &
                                    particles, dx, mla, &
                                    init_mode)

          case (2)
             call init_particles_2d(n, sop(:,:,1,:), ng_s, &
                                    rho0(n,:), rhoh0(n,:), &
                                    p0(n,:), tempbar(n,:), lo, hi, &
                                    particles, dx, mla, &
                                    init_mode)

          case (3)
             if (spherical .eq. 1) then
                call init_particles_3d_sph(n, sop(:,:,:,:), ng_s, &
                                           rho0(n,:), rhoh0(n,:), &
                                           p0(n,:), tempbar(n,:), lo, hi, &
                                           particles, dx, mla, &
                                           init_mode)
             else
                call init_particles_3d_cart(n, sop(:,:,:,:), ng_s, &
                                            rho0(n,:), rhoh0(n,:), &
                                            p0(n,:), tempbar(n,:), lo, hi, &
                                            particles, dx, mla, &
                                            init_mode)
             end if
          end select

       end do

    enddo

    call destroy(bpt)

  end subroutine init_particles

  subroutine init_particles_1d(n, s, ng_s, &
                               rho0, rhoh0, &
                               p0, tempbar, lo, hi, &
                               particles, dx, mla, &
                               init_mode)

    use variables, only: rho_comp

    integer,                  intent(in) :: n, lo(:), hi(:), ng_s
    real (kind = dp_t),       intent(in   ) ::     s(lo(1)-ng_s :,:)  
    real (kind = dp_t),       intent(in   ) :: rho0(0:), rhoh0(0:), p0(0:), tempbar(0:)
    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)    
    integer,                  intent(in   ) :: init_mode

    ! local variables
    integer :: i


    do i = lo(1), hi(1)

    enddo

  end subroutine init_particles_1d
  
  subroutine init_particles_2d(n, s, ng_s, &
                               rho0, rhoh0, &
                               p0, tempbar, lo, hi, &
                               particles, dx, mla, &
                               init_mode)

    integer,                  intent(in) :: n, lo(:), hi(:), ng_s
    real (kind = dp_t),       intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,:)  
    real (kind = dp_t),       intent(in   ) :: rho0(0:), rhoh0(0:), p0(0:), tempbar(0:)
    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)    
    integer,                  intent(in   ) :: init_mode

    ! local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

       enddo
    enddo

  end subroutine init_particles_2d

  
  subroutine init_particles_3d_cart(n, s, ng_s, &
                                    rho0, rhoh0, &
                                    p0, tempbar, lo, hi, &
                                    particles, dx, mla, &
                                    init_mode)

    use variables, only: rho_comp

    integer,                  intent(in) :: n, lo(:), hi(:), ng_s
    real (kind = dp_t),       intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)  
    real (kind = dp_t),       intent(in   ) :: rho0(0:), rhoh0(0:), p0(0:), tempbar(0:)
    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)    
    integer,                  intent(in   ) :: init_mode

    ! local variables
    integer :: i, j, k


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
    
          enddo
       enddo
    enddo

  end subroutine init_particles_3d_cart


  subroutine init_particles_3d_sph(n, s, ng_s, &
                                   rho0, rhoh0, &
                                   p0, tempbar, lo, hi, &
                                   particles, dx, mla, &
                                   init_mode)

    use probin_module, only : prob_lo, particle_temp_cutoff, particle_tpert_threshold
    use variables, only: temp_comp
    use fill_3d_module
    
    integer,                  intent(in) :: n, lo(:), hi(:), ng_s
    real (kind = dp_t),       intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)  
    real (kind = dp_t),       intent(in   ) :: rho0(0:), rhoh0(0:), p0(0:), tempbar(0:)
    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)    
    integer,                  intent(in   ) :: init_mode

    ! local variables
    integer :: i, j, k
    real(kind=dp_t) :: x, y, z
    real(kind=dp_t) :: point(3)
    real(kind=dp_t) :: Tpert
    logical :: status

    real(kind=dp_t), allocatable :: tempbar_cart(:,:,:,:)


    ! put the basestate temperature onto a Cart grid
    allocate(tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,tempbar,tempbar_cart,lo,hi,dx(n,:),0)

    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k) + HALF) * dx(n,3)

       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j) + HALF) * dx(n,2)

          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i) + HALF) * dx(n,1)
          
             Tpert = s(i,j,k,temp_comp) - tempbar_cart(i,j,k,1)
             
             ! only add a particle if the basestate temperature at
             ! this radius is high enough and the temperature
             ! perturbation is large
             if (tempbar_cart(i,j,k,1) > particle_temp_cutoff .and. &
                 Tpert > particle_tpert_threshold) then
                
                point(1) = x
                point(2) = y
                point(3) = z

                ! do a conditional add -- don't want to add a particle
                ! if there is already one nearby
                call add(particles,point,mla,dx,prob_lo, conditional=.true.)
             endif

          enddo
       enddo
    enddo

    deallocate(tempbar_cart)


  end subroutine init_particles_3d_sph
  
end module init_particles_module
