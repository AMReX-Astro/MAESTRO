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
    use geometry, only: spherical
    use fill_3d_module

    type(particle_container), intent(inout) :: particles
    type(multifab) ,          intent(in   ) :: s(:)
    real(kind=dp_t),          intent(in   ) :: rho0(:,0:)
    real(kind=dp_t),          intent(in   ) :: rhoh0(:,0:)
    real(kind=dp_t),          intent(in   ) :: p0(:,0:)
    real(kind=dp_t),          intent(in   ) :: tempbar(:,0:)
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)
    integer        ,          intent(in   ) :: init_mode
    
    real(kind=dp_t), pointer::   sop(:,:,:,:)
    
    integer :: lo(mla%dim),hi(mla%dim),i,n,nlevs
    integer :: ng_s
    
    type(bl_prof_timer), save :: bpt
    
    call build(bpt, "init_particles")
    
    nlevs = mla%nlevel
    ng_s  = nghost(s(1))

    do n=1,nlevs
       
       do i = 1, nfabs(s(n))
          sop   => dataptr(s(n), i)
          
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))

          if (spherical == 1) then
             call init_particles_sph(n, sop(:,:,:,:), ng_s, &
                                     rho0(1,:), rhoh0(1,:), &
                                     p0(1,:), tempbar(1,:), lo, hi, &
                                     particles, dx, mla, &
                                     init_mode)
          else
             call init_particles_cart(n, sop(:,:,:,:), ng_s, &
                                      rho0(n,:), rhoh0(n,:), &
                                      p0(n,:), tempbar(n,:), lo, hi, &
                                      particles, dx, mla, &
                                      init_mode)
          endif
       end do
    enddo

    call destroy(bpt)

  end subroutine init_particles

  subroutine init_particles_sph(n, s, ng_s, &
                                 rho0, rhoh0, &
                                 p0, tempbar, lo, hi, &
                                 particles, dx, mla, &
                                 init_mode)

    use probin_module, only: prob_lo

    integer,                  intent(in   ) :: n, lo(:), hi(:), ng_s
    real (kind = dp_t),       intent(in   ) :: s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real (kind = dp_t),       intent(in   ) :: rho0(0:), rhoh0(0:), p0(0:), tempbar(0:)
    type(particle_container), intent(inout) :: particles
    type(ml_layout),          intent(inout) :: mla
    real(kind=dp_t),          intent(in   ) :: dx(:,:)    
    integer,                  intent(in   ) :: init_mode

    ! local variables
    integer :: i, j, k, dm
    real(kind=dp_t) :: x, y, z
    real(kind=dp_t) :: point(mla%dim)

    !NOTE: For this stub code, init_particles_sph and init_particles_cart are
    !identical.  We maintain the two because for particular implementations
    !these may differ and we want the stub code to be a useful template.
     
    dm = mla%dim

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i) + HALF) * dx(n,1)

             select case (dm)
             case (1)
                ! point(1) = x
                
                ! call add(particles,point,mla,dx,prob_lo)
             case (2)
                y = prob_lo(2) + (dble(j) + HALF) * dx(n,2)
                
                ! point(1) = x
                ! point(2) = y
                
                ! call add(particles,point,mla,dx,prob_lo)
             case (3)
                z = prob_lo(3) + (dble(k) + HALF) * dx(n,3)
                y = prob_lo(2) + (dble(j) + HALF) * dx(n,2)
                
                ! point(1) = x
                ! point(2) = y
                ! point(3) = z
                
                ! call add(particles,point,mla,dx,prob_lo)
             end select
    
          enddo
       enddo
    enddo

  end subroutine init_particles_sph

  subroutine init_particles_cart(n, s, ng_s, &
                                 rho0, rhoh0, &
                                 p0, tempbar, lo, hi, &
                                 particles, dx, mla, &
                                 init_mode)

     use probin_module, only: prob_lo

     integer,                  intent(in   ) :: n, lo(:), hi(:), ng_s
     real (kind = dp_t),       intent(in   ) :: s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
     real (kind = dp_t),       intent(in   ) :: rho0(0:), rhoh0(0:), p0(0:), tempbar(0:)
     type(particle_container), intent(inout) :: particles
     type(ml_layout),          intent(inout) :: mla
     real(kind=dp_t),          intent(in   ) :: dx(:,:)    
     integer,                  intent(in   ) :: init_mode

     ! local variables
     integer :: i, j, k, dm
     real(kind=dp_t) :: x, y, z
     real(kind=dp_t) :: point(mla%dim)

     !NOTE: For this stub code, init_particles_sph and init_particles_cart are
     !identical.  We maintain the two because for particular implementations
     !these may differ and we want the stub code to be a useful template.
      
     dm = mla%dim

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              x = prob_lo(1) + (dble(i) + HALF) * dx(n,1)

              select case (dm)
              case (1)
                 ! point(1) = x
                 
                 ! call add(particles,point,mla,dx,prob_lo)
              case (2)
                 y = prob_lo(2) + (dble(j) + HALF) * dx(n,2)
                 
                 ! point(1) = x
                 ! point(2) = y
                 
                 ! call add(particles,point,mla,dx,prob_lo)
              case (3)
                 z = prob_lo(3) + (dble(k) + HALF) * dx(n,3)
                 y = prob_lo(2) + (dble(j) + HALF) * dx(n,2)
                 
                 ! point(1) = x
                 ! point(2) = y
                 ! point(3) = z
                 
                 ! call add(particles,point,mla,dx,prob_lo)
              end select
     
           enddo
        enddo
     enddo

  end subroutine init_particles_cart
  
end module init_particles_module
