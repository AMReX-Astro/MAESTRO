! compute w0 -- the base state velocity.  This is based on the average 
! heating in a layer (Sbar) and the mixing (the eta quantities).  The
! computation of w0 for plane-parallel atmospheres was first described
! in paper II, with modifications due to mixing in paper III.  For 
! spherical geometry, it was first described in paper III.

module make_w0_module

  use bl_types

  implicit none

  private

  public :: make_w0

contains

  subroutine make_w0(w0,w0_old,w0_force,Sbar_in, &
                     rho0_old,rho0_new,p0_old,p0_new, &
                     gamma1bar_old,gamma1bar_new,p0_minus_pthermbar, &
                     psi,etarho_ec,etarho_cc,dt,dtold)

    use parallel
    use bl_prof_module
    use geometry, only: spherical, dr, r_start_coord, r_end_coord, nlevs_radial
    use bl_constants_module
    use probin_module, only: verbose, do_planar_invsq_grav

    real(kind=dp_t), intent(  out) ::                 w0(:,0:)
    real(kind=dp_t), intent(in   ) ::             w0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::                psi(:,0:)
    real(kind=dp_t), intent(in   ) ::          etarho_ec(:,0:)
    real(kind=dp_t), intent(in   ) ::          etarho_cc(:,0:)
    real(kind=dp_t), intent(inout) ::           w0_force(:,0:)
    real(kind=dp_t), intent(in   ) ::           rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::           rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) ::             p0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::             p0_new(:,0:)
    real(kind=dp_t), intent(in   ) ::      gamma1bar_old(:,0:)
    real(kind=dp_t), intent(in   ) ::      gamma1bar_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_minus_pthermbar(:,0:)
    real(kind=dp_t), intent(in   ) ::            Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    integer         :: r,n
    real(kind=dp_t) :: max_w0

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_w0")

    w0_force = ZERO

    if (spherical .eq. 0) then

       if (.not. do_planar_invsq_grav) then
          call make_w0_planar(w0,w0_old,Sbar_in, &
                              p0_old,p0_new,gamma1bar_old,gamma1bar_new, &
                              p0_minus_pthermbar,psi,w0_force, &
                              dt,dtold)
       else
          
          call make_w0_planar_invsq(w0,w0_old,Sbar_in, &
                                    rho0_old,rho0_new,p0_old,p0_new, &
                                    gamma1bar_old,gamma1bar_new, &
                                    p0_minus_pthermbar, &
                                    etarho_cc,w0_force, &
                                    dt,dtold)
       endif


    else

       call make_w0_spherical(w0(1,:),w0_old(1,:),Sbar_in(1,:), &
                              rho0_old(1,:),rho0_new(1,:), &
                              p0_old(1,:),p0_new(1,:), &
                              gamma1bar_old(1,:),gamma1bar_new(1,:), &
                              p0_minus_pthermbar(1,:), &
                              etarho_ec(1,:),etarho_cc(1,:),w0_force(1,:), &
                              dt,dtold)

    end if

    if (verbose .ge. 1) then
       do n=1,nlevs_radial
          max_w0 = zero
          do r=r_start_coord(n,1),r_end_coord(n,1)+1
             max_w0 = max(max_w0, abs(w0(n,r)))
          end do
          if (parallel_IOProcessor()) then
             write(6,*) '... max CFL of w0: ',max_w0 * dt / dr(n)
          end if
       end do
       if (parallel_IOProcessor()) then
          write(6,*) ''
       end if
    end if

    call destroy(bpt)

  end subroutine make_w0



  subroutine make_w0_planar(w0,w0_old,Sbar_in,p0_old,p0_new, &
                            gamma1bar_old,gamma1bar_new,p0_minus_pthermbar, &
                            psi,w0_force,dt,dtold)

    use geometry, only: nr_fine, r_start_coord, r_end_coord, dr, &
         base_cutoff_density_coord,numdisjointchunks, nr
    use variables, only: rho_comp
    use bl_constants_module
    use probin_module, only: grav_const, dpdt_factor, base_cutoff_density
    use restrict_base_module

    real(kind=dp_t), intent(  out) ::                 w0(:,0:)
    real(kind=dp_t), intent(in   ) ::             w0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::            Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) ::             p0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::             p0_new(:,0:)
    real(kind=dp_t), intent(in   ) ::      gamma1bar_old(:,0:)
    real(kind=dp_t), intent(in   ) ::      gamma1bar_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_minus_pthermbar(:,0:)
    real(kind=dp_t), intent(in   ) ::                psi(:,0:)
    real(kind=dp_t), intent(  out) ::           w0_force(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! Local variables
    integer         :: r, n, i, j, refrat, nlevs
    real(kind=dp_t) :: w0_old_cen(size(w0,dim=1),0:nr_fine-1)
    real(kind=dp_t) :: w0_new_cen(size(w0,dim=1),0:nr_fine-1)
    real(kind=dp_t) :: w0_avg, div_avg, dt_avg, gamma1bar_p0_avg
    real(kind=dp_t) :: volume_discrepancy, offset

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Multilevel Outline
    !
    ! Compute w0 at level 1 only
    ! Initialize new w0 at bottom of coarse base array to zero.
    ! do n=2,nlevs
    !   Compute w0 on edges at level n
    !   Obtain the starting value of w0 from the coarser grid
    !   if n>1, compare the difference between w0 at top of level n to the 
    !           corresponding point on level n-1
    !   do i=n-1,1,-1
    !     Restrict w0 from level n to level i
    !     Offset the w0 on level i above the top of level n
    !   end do
    ! end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nlevs = size(w0,dim=1)

    w0 = ZERO
    
    ! Compute w0 on edges at level n
    do n=1,nlevs

       do j=1,numdisjointchunks(n)
          
          if (n .eq. 1) then
             ! Initialize new w0 at bottom of coarse base array to zero.
             w0(1,0) = ZERO
          else
             ! Obtain the starting value of w0 from the coarser grid
             w0(n,r_start_coord(n,j)) = w0(n-1,r_start_coord(n,j)/2)
          end if

          do r=r_start_coord(n,j)+1,r_end_coord(n,j)+1
             gamma1bar_p0_avg = (gamma1bar_old(n,r-1)+gamma1bar_new(n,r-1)) * &
                  (p0_old(n,r-1)+p0_new(n,r-1))/4.0d0

             if (r-1 .lt. base_cutoff_density_coord(n)) then
                volume_discrepancy = dpdt_factor * p0_minus_pthermbar(n,r-1)/dt
             else
                volume_discrepancy = 0.0d0
             end if

             w0(n,r) = w0(n,r-1) + Sbar_in(n,r-1) * dr(n) &
                  - ( (psi(n,r-1)+volume_discrepancy) / gamma1bar_p0_avg ) * dr(n)
          end do

          if (n .gt. 1) then

             ! Compare the difference between w0 at top of level n to
             ! the corresponding point on level n-1
             offset = w0(n,r_end_coord(n,j)+1) - w0(n-1,(r_end_coord(n,j)+1)/2)

             do i=n-1,1,-1

                refrat = 2**(n-i)
                
                ! Restrict w0 from level n to level i
                do r=r_start_coord(n,j),r_end_coord(n,j)+1
                   if (mod(r,refrat) .eq. 0) then
                      w0(i,r/refrat) = w0(n,r)
                   end if
                end do
                
                ! Offset the w0 on level i above the top of level n
                do r=(r_end_coord(n,j)+1)/refrat+1,nr(i)
                   w0(i,r) = w0(i,r) + offset
                end do
                
             end do
             
          end if

       end do

    end do

       ! zero w0 where there is no corresponding full state array
       do n=2,nlevs
          do j=1,numdisjointchunks(n)
             if (j .eq. numdisjointchunks(n)) then
                do r=r_end_coord(n,j)+2,nr(n)
                   w0(n,r) = ZERO
                end do
             else
                do r=r_end_coord(n,j)+2,r_start_coord(n,j+1)-1
                   w0(n,r) = ZERO
                end do
             end if
          end do
       end do

    call restrict_base(w0,.false.)
    call fill_ghost_base(w0,.false.)

    do n=1,nlevs
       do j=1,numdisjointchunks(n)

          ! Compute the forcing term in the base state velocity
          ! equation, - 1/rho0 grad pi0
          dt_avg = HALF * (dt + dtold)
          do r=r_start_coord(n,j),r_end_coord(n,j)
             w0_old_cen(n,r) = HALF * (w0_old(n,r) + w0_old(n,r+1))
             w0_new_cen(n,r) = HALF * (w0(n,r) + w0(n,r+1))
             w0_avg = HALF * (dt * w0_old_cen(n,r) + dtold *  w0_new_cen(n,r)) / dt_avg
             div_avg = HALF * (dt * (w0_old(n,r+1)-w0_old(n,r)) + &
                  dtold * (w0(n,r+1)-w0(n,r))) / dt_avg
             w0_force(n,r) = (w0_new_cen(n,r)-w0_old_cen(n,r))/dt_avg + w0_avg*div_avg/dr(n)
          end do
          
       end do
    end do

    call restrict_base(w0_force,.true.)
    call fill_ghost_base(w0_force,.true.)

  end subroutine make_w0_planar



  subroutine make_w0_planar_invsq(w0,w0_old,Sbar_in, &
                                  rho0_old,rho0_new,p0_old,p0_new, &
                                  gamma1bar_old,gamma1bar_new, &
                                  p0_minus_pthermbar, &
                                  etarho_cc,w0_force, &
                                  dt,dtold)

    use geometry, only: r_cc_loc, nr_fine, nlevs_radial, r_edge_loc, dr, &
         base_cutoff_density_coord, r_start_coord, r_end_coord, numdisjointchunks, nr
    use make_grav_module
    use bl_constants_module
    use fundamental_constants_module, only: Gconst
    use probin_module, only: dpdt_factor, base_cutoff_density
    use prolong_base_to_uniform_module
    use restrict_base_module
    use make_grav_module

    real(kind=dp_t), intent(  out) ::                 w0(:,0:)
    real(kind=dp_t), intent(in   ) ::             w0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::            Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) ::           rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::           rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) ::             p0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::             p0_new(:,0:)
    real(kind=dp_t), intent(in   ) ::      gamma1bar_old(:,0:)
    real(kind=dp_t), intent(in   ) ::      gamma1bar_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_minus_pthermbar(:,0:)
    real(kind=dp_t), intent(in   ) ::          etarho_cc(:,0:)
    real(kind=dp_t), intent(  out) ::           w0_force(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! local variables
    real(kind=dp_t), allocatable :: w0_fine(:)
    real(kind=dp_t), allocatable :: w0bar_fine(:)
    real(kind=dp_t), allocatable :: deltaw0_fine(:)
    real(kind=dp_t), allocatable :: p0_old_fine(:)
    real(kind=dp_t), allocatable :: p0_new_fine(:)
    real(kind=dp_t), allocatable :: p0_nph_fine(:)
    real(kind=dp_t), allocatable :: rho0_old_fine(:)
    real(kind=dp_t), allocatable :: rho0_new_fine(:)
    real(kind=dp_t), allocatable :: rho0_nph_fine(:)
    real(kind=dp_t), allocatable :: gamma1bar_old_fine(:)
    real(kind=dp_t), allocatable :: gamma1bar_new_fine(:)
    real(kind=dp_t), allocatable :: gamma1bar_nph_fine(:)
    real(kind=dp_t), allocatable :: p0_minus_pthermbar_fine(:)
    real(kind=dp_t), allocatable :: etarho_cc_fine(:)
    real(kind=dp_t), allocatable :: Sbar_in_fine(:)
    real(kind=dp_t), allocatable :: grav_edge_fine(:)

    integer :: n, r, j

    real(kind=dp_t), allocatable :: A(:), B(:), C(:), u(:), F(:)

    real(kind=dp_t) :: gamma1bar_p0_avg, volume_discrepancy, dpdr
    real(kind=dp_t) :: dt_avg, w0_avg, div_avg
    real(kind=dp_t) :: w0_old_cen(nlevs_radial,0:nr_fine-1)
    real(kind=dp_t) :: w0_new_cen(nlevs_radial,0:nr_fine-1)


    ! The planar 1/r**2 gravity constraint equation is solved 
    ! by calling the tridiagonal solver, just like spherical.  
    ! This is accomplished by putting all the requisite data
    ! on the finest basestate grid, solving for w0, and then
    ! restricting w0 back down to the coarse grid.
    

    ! 1) allocate the finely-gridded temporary basestate arrays
    allocate(                w0_fine(0:nr_fine))
    allocate(             w0bar_fine(0:nr_fine))
    allocate(           deltaw0_fine(0:nr_fine))
    allocate(            p0_old_fine(0:nr_fine-1))
    allocate(            p0_new_fine(0:nr_fine-1))
    allocate(            p0_nph_fine(0:nr_fine-1))
    allocate(          rho0_old_fine(0:nr_fine-1))
    allocate(          rho0_new_fine(0:nr_fine-1))
    allocate(          rho0_nph_fine(0:nr_fine-1))
    allocate(     gamma1bar_old_fine(0:nr_fine-1))
    allocate(     gamma1bar_new_fine(0:nr_fine-1))
    allocate(     gamma1bar_nph_fine(0:nr_fine-1))
    allocate(p0_minus_pthermbar_fine(0:nr_fine-1))
    allocate(         etarho_cc_fine(0:nr_fine-1))
    allocate(           Sbar_in_fine(0:nr_fine-1))
    allocate(         grav_edge_fine(0:nr_fine))


    ! 2) copy the data into the temp, uniformly-gridded basestate arrays.
    call prolong_base_to_uniform(p0_old,p0_old_fine)
    call prolong_base_to_uniform(p0_new,p0_new_fine)
    call prolong_base_to_uniform(rho0_old,rho0_old_fine)
    call prolong_base_to_uniform(rho0_new,rho0_new_fine)
    call prolong_base_to_uniform(gamma1bar_old,gamma1bar_old_fine)
    call prolong_base_to_uniform(gamma1bar_new,gamma1bar_new_fine)
    call prolong_base_to_uniform(p0_minus_pthermbar,p0_minus_pthermbar_fine)
    call prolong_base_to_uniform(etarho_cc,etarho_cc_fine)
    call prolong_base_to_uniform(Sbar_in,Sbar_in_fine)

    ! create time-centered base-state quantities
    do r=0,nr_fine-1
       p0_nph_fine(r)        = HALF*(p0_old_fine(r)        + p0_new_fine(r))
       rho0_nph_fine(r)      = HALF*(rho0_old_fine(r)      + rho0_new_fine(r))
       gamma1bar_nph_fine(r) = HALF*(gamma1bar_old_fine(r) + gamma1bar_new_fine(r))       
    enddo


    ! 3) solve to w0bar -- here we just take into account the Sbar and
    !    volume discrepancy terms
    w0bar_fine(:) = ZERO

    ! lower boundary condition
    w0bar_fine(0) = ZERO

    do r = 1, nr_fine
       gamma1bar_p0_avg = gamma1bar_nph_fine(r-1) * p0_nph_fine(r-1)

       if (r-1 .lt. base_cutoff_density_coord(nlevs_radial)) then
          volume_discrepancy = dpdt_factor * p0_minus_pthermbar_fine(r-1)/dt
       else
          volume_discrepancy = ZERO
       end if

       w0bar_fine(r) =  w0bar_fine(r-1) + Sbar_in_fine(r-1) * dr(nlevs_radial) &
            - (volume_discrepancy / gamma1bar_p0_avg ) * dr(nlevs_radial)
       
    enddo



    ! 4) get the edge-centered gravity on the uniformly-gridded
    ! basestate arrays
    call make_grav_edge_uniform(grav_edge_fine, rho0_nph_fine)


    ! 5) solve for delta w0
    deltaw0_fine(:) = ZERO

    ! this takes the form of a tri-diagonal matrix:
    ! A_j (dw_0)_{j-3/2} + 
    ! B_j (dw_0)_{j-1/2} + 
    ! C_j (dw_0)_{j+1/2} = F_j

    allocate(A(0:nr_fine))
    allocate(B(0:nr_fine))
    allocate(C(0:nr_fine))
    allocate(u(0:nr_fine))
    allocate(F(0:nr_fine))

    A   = ZERO
    B   = ZERO
    C   = ZERO
    F   = ZERO
    u   = ZERO
   
    do r=1,base_cutoff_density_coord(nlevs_radial)
       A(r) = gamma1bar_nph_fine(r-1) * p0_nph_fine(r-1) 
       A(r) = A(r) / dr(nlevs_radial)**2
    end do

    do r=1,base_cutoff_density_coord(nlevs_radial)
       B(r) = -(gamma1bar_nph_fine(r-1) * p0_nph_fine(r-1) + &
                gamma1bar_nph_fine(r  ) * p0_nph_fine(r  )) / dr(nlevs_radial)**2 

       dpdr = (p0_nph_fine(r)-p0_nph_fine(r-1))/dr(nlevs_radial)
       B(r) = B(r) - TWO * dpdr / (r_edge_loc(nlevs_radial,r))
    end do

    do r=1,base_cutoff_density_coord(nlevs_radial)
       C(r) = gamma1bar_nph_fine(r) * p0_nph_fine(r) 
       C(r) = C(r) / dr(nlevs_radial)**2
    end do

    do r=1,base_cutoff_density_coord(nlevs_radial)
       dpdr = (p0_nph_fine(r)-p0_nph_fine(r-1))/dr(nlevs_radial)
       F(r) = TWO * dpdr * w0bar_fine(r) / r_edge_loc(nlevs_radial,r) - &
              grav_edge_fine(r) * (etarho_cc_fine(r) - etarho_cc_fine(r-1)) / &
              dr(nlevs_radial)
    end do

    ! Lower boundary
    A(0) = zero
    B(0) = one
    C(0) = zero
    F(0) = zero

    ! Upper boundary
    A(base_cutoff_density_coord(nlevs_radial)+1) = -one
    B(base_cutoff_density_coord(nlevs_radial)+1) = one
    C(base_cutoff_density_coord(nlevs_radial)+1) = zero
    F(base_cutoff_density_coord(nlevs_radial)+1) = zero

    ! Call the tridiagonal solver
    call tridiag(A, B, C, F, u, base_cutoff_density_coord(nlevs_radial)+2)

    do r=1,base_cutoff_density_coord(nlevs_radial)+1
       deltaw0_fine(r) = u(r)
    end do

    do r=base_cutoff_density_coord(nlevs_radial)+2,nr_fine
       deltaw0_fine(r) = deltaw0_fine(base_cutoff_density_coord(nlevs_radial)+1)
    end do


    ! 6) compute w0 = w0bar + deltaw0
    w0_fine = w0bar_fine + deltaw0_fine


    ! 7) fill the multilevel w0 array from the uniformly-gridded w0 we
    ! just solved for.  Here, we make the coarse edge underneath equal
    ! to the fine edge value.
    w0(nlevs_radial,:) = w0_fine(:)
    do n = nlevs_radial, 2, -1
       do r = 0, nr(n), 2
          w0(n-1,r/2) = w0(n,r)
       enddo
    enddo


    ! 8) zero w0 where there is no corresponding full state array
    do n=2,nlevs_radial
       do j=1,numdisjointchunks(n)
          if (j .eq. numdisjointchunks(n)) then
             do r=r_end_coord(n,j)+2,nr(n)
                w0(n,r) = ZERO
             end do
          else
             do r=r_end_coord(n,j)+2,r_start_coord(n,j+1)-1
                w0(n,r) = ZERO
             end do
          end if
       end do
    end do

    call restrict_base(w0,.false.)
    call fill_ghost_base(w0,.false.)


    ! compute the forcing terms
    do n=1,nlevs_radial
       do j=1,numdisjointchunks(n)

          ! Compute the forcing term in the base state velocity
          ! equation, - 1/rho0 grad pi0
          dt_avg = HALF * (dt + dtold)
          do r=r_start_coord(n,j),r_end_coord(n,j)
             w0_old_cen(n,r) = HALF * (w0_old(n,r) + w0_old(n,r+1))
             w0_new_cen(n,r) = HALF * (w0(n,r) + w0(n,r+1))
             w0_avg = HALF * (dt * w0_old_cen(n,r) + dtold *  w0_new_cen(n,r)) / dt_avg
             div_avg = HALF * (dt * (w0_old(n,r+1)-w0_old(n,r)) + &
                  dtold * (w0(n,r+1)-w0(n,r))) / dt_avg
             w0_force(n,r) = (w0_new_cen(n,r)-w0_old_cen(n,r))/dt_avg + w0_avg*div_avg/dr(n)
          end do
          
       end do
    end do

    call restrict_base(w0_force,.true.)
    call fill_ghost_base(w0_force,.true.)


  end subroutine make_w0_planar_invsq



  subroutine make_w0_spherical(w0,w0_old,Sbar_in, &
                               rho0_old,rho0_new,p0_old,p0_new, &
                               gamma1bar_old,gamma1bar_new, &
                               p0_minus_pthermbar, &
                               etarho_ec,etarho_cc,w0_force, &
                               dt,dtold)

    use geometry, only: r_cc_loc, nr_fine, r_edge_loc, dr, &
         base_cutoff_density_coord
    use make_grav_module
    use bl_constants_module
    use fundamental_constants_module, only: Gconst
    use probin_module, only: dpdt_factor, base_cutoff_density

    real(kind=dp_t), intent(  out) ::                 w0(0:)
    real(kind=dp_t), intent(in   ) ::             w0_old(0:)
    real(kind=dp_t), intent(in   ) ::            Sbar_in(0:)
    real(kind=dp_t), intent(in   ) ::           rho0_old(0:)
    real(kind=dp_t), intent(in   ) ::           rho0_new(0:)
    real(kind=dp_t), intent(in   ) ::             p0_old(0:)
    real(kind=dp_t), intent(in   ) ::             p0_new(0:)
    real(kind=dp_t), intent(in   ) ::      gamma1bar_old(0:)
    real(kind=dp_t), intent(in   ) ::      gamma1bar_new(0:)
    real(kind=dp_t), intent(in   ) :: p0_minus_pthermbar(0:)
    real(kind=dp_t), intent(in   ) ::          etarho_ec(0:)
    real(kind=dp_t), intent(in   ) ::          etarho_cc(0:)
    real(kind=dp_t), intent(  out) ::           w0_force(0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! Local variables
    integer                    :: r
    real(kind=dp_t)            :: dpdr, volume_discrepancy, w0_avg, div_avg, dt_avg

    real(kind=dp_t) ::    w0_old_cen(0:nr_fine-1)
    real(kind=dp_t) ::    w0_new_cen(0:nr_fine-1)
    real(kind=dp_t) :: gamma1bar_nph(0:nr_fine-1)
    real(kind=dp_t) ::        p0_nph(0:nr_fine-1)
    real(kind=dp_t) ::             A(0:nr_fine)
    real(kind=dp_t) ::             B(0:nr_fine)
    real(kind=dp_t) ::             C(0:nr_fine)
    real(kind=dp_t) ::             u(0:nr_fine)
    real(kind=dp_t) ::             F(0:nr_fine)
    real(kind=dp_t) ::  w0_from_Sbar(0:nr_fine)

    ! These need the extra dimension so we can call make_grav_edge
    real(kind=dp_t) ::  rho0_nph(1,0:nr_fine-1)
    real(kind=dp_t) :: grav_edge(1,0:nr_fine-1)

    ! create time-centered base-state quantities
    do r=0,nr_fine-1
       p0_nph(r)        = HALF*(p0_old(r)        + p0_new(r))
       rho0_nph(1,r)    = HALF*(rho0_old(r)      + rho0_new(r))
       gamma1bar_nph(r) = HALF*(gamma1bar_old(r) + gamma1bar_new(r))       
    enddo

    ! NOTE: We first solve for the w0 resulting only from Sbar,
    !      w0_from_sbar by integrating d/dr (r^2 w0_from_sbar) =
    !      (r^2 Sbar).  Then we will solve for the update, delta w0.

    w0_from_Sbar = ZERO
    do r=1,nr_fine

       if (rho0_old(r-1) .gt. base_cutoff_density) then
          volume_discrepancy = dpdt_factor * p0_minus_pthermbar(r-1)/dt
       else
          volume_discrepancy = ZERO
       endif

       w0_from_Sbar(r) = w0_from_Sbar(r-1) + dr(1) * Sbar_in(r-1) * r_cc_loc(1,r-1)**2 - &
            dr(1)* volume_discrepancy * r_cc_loc(1,r-1)**2 / (gamma1bar_nph(r-1)*p0_nph(r-1))

    end do

    do r=1,nr_fine
       w0_from_Sbar(r) = w0_from_Sbar(r) / r_edge_loc(1,r)**2
    end do

    ! make the edge-centered gravity
    call make_grav_edge(grav_edge,rho0_nph)

    ! NOTE:  now we solve for the remainder, (r^2 * delta w0)
    ! this takes the form of a tri-diagonal matrix:
    ! A_j (r^2 dw_0)_{j-3/2} + 
    ! B_j (r^2 dw_0)_{j-1/2} + 
    ! C_j (r^2 dw_0)_{j+1/2} = F_j

    A   = ZERO
    B   = ZERO
    C   = ZERO
    F   = ZERO
    u   = ZERO
   
    ! Note that we are solving for (r^2 delta w0), not just w0. 

    do r=1,base_cutoff_density_coord(1)
       A(r) = gamma1bar_nph(r-1) * p0_nph(r-1) / r_cc_loc(1,r-1)**2
       A(r) = A(r) / dr(1)**2
    end do

    do r=1,base_cutoff_density_coord(1)
       B(r) = -( gamma1bar_nph(r-1) * p0_nph(r-1) / r_cc_loc(1,r-1)**2 &
                +gamma1bar_nph(r  ) * p0_nph(r  ) / r_cc_loc(1,r  )**2 ) / dr(1)**2 

       dpdr = (p0_nph(r)-p0_nph(r-1))/dr(1)
       B(r) = B(r) - four * dpdr / (r_edge_loc(1,r))**3
    end do

    do r=1,base_cutoff_density_coord(1)
       C(r) = gamma1bar_nph(r) * p0_nph(r) / r_cc_loc(1,r)**2
       C(r) = C(r) / dr(1)**2
    end do

    do r=1,base_cutoff_density_coord(1)
       dpdr = (p0_nph(r)-p0_nph(r-1))/dr(1)
       F(r) = four * dpdr * w0_from_Sbar(r) / r_edge_loc(1,r) - &
              grav_edge(1,r) * (r_cc_loc(1,r  )**2 * etarho_cc(r  ) - &
              r_cc_loc(1,r-1)**2 * etarho_cc(r-1)) / &
              (dr(1) * r_edge_loc(1,r)**2) - &
              four * M_PI * Gconst * HALF * &
              (rho0_nph(1,r) + rho0_nph(1,r-1)) * etarho_ec(r)
    end do

    ! Lower boundary
    A(0) = zero
    B(0) = one
    C(0) = zero
    F(0) = zero

    ! Upper boundary
    A(base_cutoff_density_coord(1)+1) = -one
    B(base_cutoff_density_coord(1)+1) = one
    C(base_cutoff_density_coord(1)+1) = zero
    F(base_cutoff_density_coord(1)+1) = zero

    ! Call the tridiagonal solver
    call tridiag(A, B, C, F, u, base_cutoff_density_coord(1)+2)

    w0(0) = ZERO
    do r=1,base_cutoff_density_coord(1)+1
       w0(r) = u(r) / r_edge_loc(1,r)**2
    end do

    do r=0,base_cutoff_density_coord(1)+1
       w0(r) = w0(r) + w0_from_Sbar(r)
    end do

    do r=base_cutoff_density_coord(1)+2,nr_fine
       w0(r) = w0(base_cutoff_density_coord(1)+1)&
            *r_edge_loc(1,base_cutoff_density_coord(1)+1)**2/r_edge_loc(1,r)**2
    end do

    ! Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0 
    dt_avg = HALF * (dt + dtold)
    do r = 0,nr_fine-1
       w0_old_cen(r) = HALF * (w0_old(r) + w0_old(r+1))
       w0_new_cen(r) = HALF * (w0    (r) + w0    (r+1))
       w0_avg = HALF * (dt *  w0_old_cen(r)           + dtold *  w0_new_cen(r)  ) / dt_avg
       div_avg = HALF * (dt * (w0_old(r+1)-w0_old(r)) + dtold * (w0(r+1)-w0(r))) / dt_avg
       w0_force(r) = (w0_new_cen(r)-w0_old_cen(r)) / dt_avg + w0_avg * div_avg / dr(1)
    end do

  end subroutine make_w0_spherical

  subroutine tridiag(a,b,c,r,u,n)

    use bl_error_module

    real(kind=dp_t), intent(in   ) :: a(:), b(:), c(:), r(:)
    real(kind=dp_t), intent(  out) :: u(:)
    integer, intent(in)            :: n

    ! local
    real(kind=dp_t) :: bet
    real(kind=dp_t) :: gam(n)
    integer         :: j

    if ( b(1) .eq. 0 ) call bl_error('tridiag: CANT HAVE B(1) = ZERO')
    
    bet = b(1)
    u(1) = r(1)/bet
    
    do j = 2,n
       gam(j) = c(j-1)/bet
       bet = b(j) - a(j)*gam(j)
       if ( bet .eq. 0 ) call bl_error('tridiag: TRIDIAG FAILED')
       u(j) = (r(j)-a(j)*u(j-1))/bet
    end do
    
    do j = n-1,1,-1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do
    
  end subroutine tridiag
  
end module make_w0_module
