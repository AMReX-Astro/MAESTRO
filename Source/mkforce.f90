module mk_vel_force_module

  ! Compute the force that appears in the velocity (or momentum)
  ! equations.  This is used both when predicting the interface
  ! states and in the final, conservative update.

  ! for the final conservative update of the velocity, we need to
  ! time-center the Coriolis term ( -2 omega x U ), which means we
  ! should use umac.  This is selected by setting is_final_update = T

  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private
  public :: mk_vel_force, add_utilde_force

contains

  subroutine mk_vel_force(vel_force,is_final_update, &
                          uold,umac,w0,w0mac,gpi,s,index_rho, &
                          rho0,grav,dx,w0_force,w0_force_cart,the_bc_level,mla)

    ! index_rho refers to the index into s where the density lives.
    ! Usually s will be the full state array, and index_rho would
    ! be rho_comp, but sometimes we pass in only a single-variable
    ! multifab array, so index_rho may be different.

    use bl_prof_module
    use geometry, only: spherical, dm, nlevs
    use variables, only: foextrap_comp
    use bl_constants_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module
    use multifab_physbc_module
    use probin_module, only: evolve_base_state
    use fill_3d_module, only : make_w0mac, put_1d_array_on_cart

    type(multifab) , intent(inout) :: vel_force(:)
    logical        , intent(in   ) :: is_final_update
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    type(multifab) , intent(in   ) :: gpi(:)
    type(multifab) , intent(in   ) :: s(:)
    integer                        :: index_rho  
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: w0_force(:,0:)
    type(multifab) , intent(in   ) :: w0_force_cart(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    ! Local variables
    real(kind=dp_t), pointer ::  uop(:,:,:,:)
    real(kind=dp_t), pointer ::  ump(:,:,:,:)
    real(kind=dp_t), pointer ::  vmp(:,:,:,:)
    real(kind=dp_t), pointer ::  wmp(:,:,:,:)
    real(kind=dp_t), pointer :: w0cp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer ::  gpp(:,:,:,:)
    real(kind=dp_t), pointer ::   fp(:,:,:,:)
    real(kind=dp_t), pointer ::   rp(:,:,:,:)
    real(kind=dp_t), pointer ::  w0p(:,:,:,:)
    integer                  :: i,comp,lo(dm),hi(dm)
    integer                  :: ng_s,ng_f,ng_gp,n,ng_uo,ng_um
   
    type(multifab) :: w0_cart(mla%nlevel)
    integer :: ng_wc, ng_wm, ng_w

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mk_vel_force")

    ng_s  = nghost(s(1))
    ng_f  = nghost(vel_force(1))
    ng_gp = nghost(gpi(1))

    ! put w0 on cell centers
    if (spherical .eq. 1) then
       do n=1,nlevs
          ! w0_cart will contain the cell-centered Cartesian components of w0, 
          ! for use in computing the Coriolis term in the prediction
          ! w0mac is passed in and is used to compute the Coriolis term
          ! in the cell update
          call build(w0_cart(n),mla%la(n),dm,0)
          call setval(w0_cart(n), ZERO, all=.true.)
       enddo       

       if (evolve_base_state) then
          ! fill the all dm components of the cell-centered w0_cart
          call put_1d_array_on_cart(w0,w0_cart,foextrap_comp,.true.,.true.,dx, &
                                    the_bc_level,mla)
       end if

    endif

    do n=1,nlevs
       do i=1,nboxes(s(n))
          if ( multifab_remote(s(n), i) ) cycle
          fp  => dataptr(vel_force(n),i)
          gpp => dataptr(gpi(n),i)
          rp  => dataptr(s(n),i)
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (1)
             call mk_vel_force_1d(fp(:,1,1,1),ng_f,gpp(:,1,1,1),ng_gp,rp(:,1,1,index_rho), &
                                  ng_s,rho0(n,:),grav(n,:),w0_force(n,:),lo,hi)
          case (2)
             call mk_vel_force_2d(fp(:,:,1,:),ng_f,gpp(:,:,1,:),ng_gp,rp(:,:,1,index_rho), &
                                  ng_s,rho0(n,:),grav(n,:),w0_force(n,:),lo,hi)
          case (3)
             uop => dataptr(uold(n),i)
             ump => dataptr(umac(n,1),i)
             vmp => dataptr(umac(n,2),i)
             wmp => dataptr(umac(n,3),i)
             ng_uo = nghost(uold(1))
             ng_um = nghost(umac(1,1))

             if (spherical .eq. 1) then
                w0cp  => dataptr(w0_cart(n), i)
                w0xp  => dataptr(w0mac(n,1),i)
                w0yp  => dataptr(w0mac(n,2),i)
                w0p   => dataptr(w0_force_cart(n), i)
                ng_wm = nghost(w0mac(1,1))
                ng_wc = nghost(w0_cart(1))
                ng_w = nghost(w0_force_cart(1))
                call mk_vel_force_3d_sphr(fp(:,:,:,:),ng_f,is_final_update, &
                                          uop(:,:,:,:),ng_uo, &
                                          ump(:,:,:,1),vmp(:,:,:,1),ng_um, &
                                          w0cp(:,:,:,:),ng_wc, &
                                          w0xp(:,:,:,1),w0yp(:,:,:,1),ng_wm, &
                                          gpp(:,:,:,:),ng_gp,rp(:,:,:,index_rho),ng_s, &
                                          rho0(1,:),grav(1,:),w0p(:,:,:,:),ng_w,lo,hi,dx(n,:))
             else
                call mk_vel_force_3d_cart(fp(:,:,:,:),ng_f,is_final_update, &
                                          uop(:,:,:,:),ng_uo, &
                                          ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),ng_um, &
                                          w0(n,:), &
                                          gpp(:,:,:,:),ng_gp,rp(:,:,:,index_rho),ng_s, &
                                          rho0(n,:),grav(n,:),w0_force(n,:),lo,hi)
             end if
          end select
       end do
    enddo

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(w0_cart(n))
       end do
    end if

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(vel_force(nlevs))

       ! fill non-periodic domain boundary ghost cells
       do i=1,dm
          call multifab_physbc(vel_force(nlevs),i,foextrap_comp,1,the_bc_level(nlevs))
       end do

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(vel_force(n-1),vel_force(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          do i=1,dm
             call multifab_fill_ghost_cells(vel_force(n),vel_force(n-1), &
                                            ng_f,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1),the_bc_level(n), &
                                            i,foextrap_comp,1,fill_crse_input=.false.)
          end do

       enddo

    end if

    call destroy(bpt)

  end subroutine mk_vel_force

  subroutine mk_vel_force_1d(vel_force,ng_f,gpi,ng_gp,rho,ng_s,rho0,grav,w0_force,lo,hi)

    use bl_constants_module
    use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     gpi(lo(1)-ng_gp:)
    real(kind=dp_t), intent(in   ) ::       rho(lo(1)-ng_s :)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: w0_force(0:)

    integer         :: i
    real(kind=dp_t) :: rhopert

    vel_force = ZERO

    do i = lo(1),hi(1)

       rhopert = rho(i) - rho0(i)
       
       ! cutoff the buoyancy term if we are outside of the star
       if (rho(i) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
          rhopert = 0.d0
       end if

       vel_force(i) =  rhopert / rho(i) * grav(i) - gpi(i) / rho(i) - w0_force(i)

    end do

  end subroutine mk_vel_force_1d

  subroutine mk_vel_force_2d(vel_force,ng_f,gpi,ng_gp,rho,ng_s,rho0,grav,w0_force,lo,hi)

    use bl_constants_module
    use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,:)
    real(kind=dp_t), intent(in   ) ::       rho(lo(1)-ng_s :,lo(2)-ng_s :)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: w0_force(0:)

    integer         :: i,j
    real(kind=dp_t) :: rhopert

    vel_force = ZERO

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          rhopert = rho(i,j) - rho0(j)
          
          ! cutoff the buoyancy term if we are outside of the star
          if (rho(i,j) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
             rhopert = 0.d0
          end if

          vel_force(i,j,1) = - gpi(i,j,1) / rho(i,j)
          vel_force(i,j,2) =  rhopert / rho(i,j) * grav(j) &
               - gpi(i,j,2) / rho(i,j) - w0_force(j)
       end do
    end do

  end subroutine mk_vel_force_2d

  subroutine mk_vel_force_3d_cart(vel_force,ng_f,is_final_update, &
                                  uold,ng_uo, &
                                  umac,vmac,wmac,ng_um, &
                                  w0, &
                                  gpi,ng_gp,rho,ng_s, &
                                  rho0,grav,w0_force,lo,hi)

    use geometry,  only: sin_theta, cos_theta, omega
    use bl_constants_module
    use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor, &
                             rotation_radius

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s, ng_uo, ng_um
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    logical        , intent(in   ) :: is_final_update
    real(kind=dp_t), intent(in   ) ::      uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
    real(kind=dp_t), intent(in   ) ::      umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::      vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::      wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
    real(kind=dp_t), intent(in   ) ::       rho(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: w0_force(0:)

    integer         :: i,j,k
    real(kind=dp_t) :: rhopert

    real(kind=dp_t) :: coriolis_term(3), centrifugal_term(3)

    vel_force = ZERO

    ! CURRENTLY for rotation in plane-parallel, we make the (bad) assumption 
    ! that all points within the patch have the same centrifugal forcing terms.
    !
    ! We assume the centrifugal term applies at a constant radius, 
    ! rotation_radius, for the patch.  In otherwords, the patch lives on the
    ! surface of a sphere of radius rotation_radius.
    !
    ! Furthermore, we assume the patch lives at longitude = 0.
    !
    ! Then the orientation of the patch is such that e_z is in the 
    ! outward radial direction of the star, e_x is in the co_latitude (polar) 
    ! angle direction and e_y is in the global y-direction.
    !
    ! centrifugal_term = omega x (omega x r) = (omega dot r) * omega
    !                                          - omega^2 * r
    ! where omega = (-|omega| sin_theta) e_x + (|omega| cos_theta) e_z
    !           r = rotation_radius e_z
    !
    ! See docs/rotation for derivation and figures.
    ! 

    centrifugal_term(1) = - omega**2 * rotation_radius * sin_theta * sin_theta
    centrifugal_term(2) = ZERO
    centrifugal_term(3) = omega**2 * rotation_radius * cos_theta * sin_theta &
                          - omega**2 * rotation_radius

    !$OMP PARALLEL DO PRIVATE(i,j,k,rhopert,coriolis_term)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             rhopert = rho(i,j,k) - rho0(k)
             
             ! cutoff the buoyancy term if we are outside of the star
             if (rho(i,j,k) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
                rhopert = 0.d0
             end if

             ! the coriolis term is:
             !    TWO * omega x U
             ! where omega is given above and U = (u, v, w) is the velocity

             if (is_final_update) then

                ! use umac so we are time-centered
                coriolis_term(1) = -TWO * omega * &
                     HALF*(vmac(i,j,k) + vmac(i,j+1,k)) * cos_theta

                coriolis_term(2) =  TWO * omega * &
                     (HALF*(wmac(i,j,k)   + w0(k) + &
                            wmac(i,j,k+1) + w0(k+1)) * sin_theta + &
                      HALF*(umac(i,j,k) + umac(i+1,j,k)) * cos_theta)

                coriolis_term(3) = -TWO * omega * &
                     HALF*(vmac(i,j,k) + vmac(i,j+1,k)) * sin_theta

             else
                coriolis_term(1) = -TWO * omega * uold(i,j,k,2) * cos_theta

                coriolis_term(2) =  TWO * omega * ((uold(i,j,k,3) + HALF*(w0(k) + w0(k+1))) * sin_theta + &
                                                   uold(i,j,k,1) * cos_theta)

                coriolis_term(3) = -TWO * omega * uold(i,j,k,2) * sin_theta
             endif

             vel_force(i,j,k,1) = -coriolis_term(1) - centrifugal_term(1) - &
                  gpi(i,j,k,1) / rho(i,j,k) 

             vel_force(i,j,k,2) = -coriolis_term(2) - centrifugal_term(2) - &
                  gpi(i,j,k,2) / rho(i,j,k) 

             vel_force(i,j,k,3) = -coriolis_term(3) - centrifugal_term(3) + &
                  ( rhopert * grav(k) - gpi(i,j,k,3) ) / rho(i,j,k) &
                  - w0_force(k)

          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine mk_vel_force_3d_cart

  subroutine mk_vel_force_3d_sphr(vel_force,ng_f,is_final_update, &
                                  uold,ng_uo, &
                                  umac,vmac,ng_um, &
                                  w0_cart,ng_wc, &
                                  w0macx,w0macy,ng_wm, &
                                  gpi,ng_gp,rho,ng_s, &
                                  rho0,grav,w0_force_cart,ng_w,lo,hi,dx)

    use fill_3d_module
    use bl_constants_module
    use geometry,  only: omega, center
    use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor, prob_lo

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_gp,ng_s,ng_uo,ng_um,ng_wc,ng_wm,ng_w
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    logical        , intent(in   ) :: is_final_update
    real(kind=dp_t), intent(in   ) ::      uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
    real(kind=dp_t), intent(in   ) ::      umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::      vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   w0_cart(lo(1)-ng_wc:,lo(2)-ng_wc:,lo(3)-ng_wc:,:)
    real(kind=dp_t), intent(in   ) ::    w0macx(lo(1)-ng_wm:,lo(2)-ng_wm:,lo(3)-ng_wm:)
    real(kind=dp_t), intent(in   ) ::    w0macy(lo(1)-ng_wm:,lo(2)-ng_wm:,lo(3)-ng_wm:)
    real(kind=dp_t), intent(in   ) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
    real(kind=dp_t), intent(in   ) ::       rho(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) :: w0_force_cart(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:,:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) ::   dx(:)

    integer         :: i,j,k

    real(kind=dp_t), allocatable :: rho0_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: grav_cart(:,:,:,:)

    real(kind=dp_t) :: rhopert
    real(kind=dp_t) :: xx, yy, zz
    real(kind=dp_t) :: centrifugal_term(3), coriolis_term(3)

    allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(grav_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3))

    vel_force = ZERO

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,rho0,rho0_cart,lo,hi,dx,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.true.,grav,grav_cart,lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k,xx,yy,zz,rhopert,centrifugal_term,coriolis_term)
    do k = lo(3),hi(3)
       zz = prob_lo(3) + (dble(k) + HALF)*dx(3) - center(3)
       do j = lo(2),hi(2)
          yy = prob_lo(2) + (dble(j) + HALF)*dx(2) - center(2)
          do i = lo(1),hi(1)
             xx = prob_lo(1) + (dble(i) + HALF)*dx(1) - center(1)

             rhopert = rho(i,j,k) - rho0_cart(i,j,k,1)

             ! cutoff the buoyancy term if we are outside of the star
             if (rho(i,j,k) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
                rhopert = 0.d0
             end if


             ! Coriolis and centrifugal forces.  We assume that the
             ! rotation axis is the z direction, with angular velocity
             ! omega

             ! omega x (omega x r ) = - omega^2 x e_x  - omega^2 y e_y    
             ! (with omega = omega e_z)
             centrifugal_term(1) = -omega * omega * xx
             centrifugal_term(2) = -omega * omega * yy
             centrifugal_term(3) = ZERO

             ! cutoff the centrifugal term if we are outside the star
             if (rho(i,j,k) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
                centrifugal_term(:) = 0.d0
             end if


             ! 2 omega x U = - 2 omega v e_x  + 2 omega u e_y
             ! (with omega = omega e_z)
             if (is_final_update) then

                ! use umac so we are time-centered
                coriolis_term(1) = -TWO * omega * &
                     HALF*(vmac(i,j,k)   + w0macy(i,j,k) + &
                           vmac(i,j+1,k) + w0macy(i,j+1,k))

                coriolis_term(2) =  TWO * omega * &
                     HALF*(umac(i,j,k)   + w0macx(i,j,k) + &
                           umac(i+1,j,k) + w0macx(i+1,j,k))

                coriolis_term(3) = ZERO

             else
                coriolis_term(1) = -TWO * omega * (uold(i,j,k,2) + w0_cart(i,j,k,2))
                coriolis_term(2) =  TWO * omega * (uold(i,j,k,1) + w0_cart(i,j,k,1))
                coriolis_term(3) = ZERO
             endif


             ! F_Coriolis = -2 omega x U  
             ! F_centrifugal = - omega x (omega x r)

             ! we just computed the absolute value of the forces above, so use
             ! the right sign here
             vel_force(i,j,k,1) = -coriolis_term(1) - centrifugal_term(1) + &
                  ( rhopert * grav_cart(i,j,k,1) - gpi(i,j,k,1) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,1)

             vel_force(i,j,k,2) = -coriolis_term(2) - centrifugal_term(2) + &
                  ( rhopert * grav_cart(i,j,k,2) - gpi(i,j,k,2) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,2)

             vel_force(i,j,k,3) = -coriolis_term(3) - centrifugal_term(3) + &
                  ( rhopert * grav_cart(i,j,k,3) - gpi(i,j,k,3) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,3)

          end do
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(rho0_cart,grav_cart)

  end subroutine mk_vel_force_3d_sphr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_utilde_force(vel_force,normal,umac,w0,dx,the_bc_level,mla)

    use geometry, only: spherical, nr_fine, dr, dm, nlevs
    use fill_3d_module, only: put_1d_array_on_cart
    use variables, only: foextrap_comp

    type(multifab) , intent(inout) :: vel_force(:)
    type(multifab) , intent(in   ) ::    normal(:)
    type(multifab) , intent(in   ) ::      umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(in   ) :: mla
    
    ! local
    integer :: n,r,i
    integer :: ng_f, ng_u, ng_w0, ng_n

    real(kind=dp_t), pointer ::  fp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: gwp(:,:,:,:)
    real(kind=dp_t), pointer ::  np(:,:,:,:)

    ! stuff for spherical only
    real(kind=dp_t) :: gradw0_rad(1,0:nr_fine-1)
    type(multifab)  :: gradw0_cart(nlevs)
    integer         :: lo(dm), hi(dm), n_1d

    type(bl_prof_timer), save :: bpt

    call build(bpt, "add_utilde_force")

    do n=1,nlevs
       call multifab_build(gradw0_cart(n),get_layout(vel_force(n)),1,1)
    end do

    if (spherical .eq. 1) then
       do r=0,nr_fine-1
          gradw0_rad(1,r) = (w0(1,r+1) - w0(1,r)) / dr(1)
       enddo
       call put_1d_array_on_cart(gradw0_rad,gradw0_cart,foextrap_comp, &
                                 .false.,.false.,dx,the_bc_level,mla)
    end if

    ng_f = nghost(vel_force(1))
    ng_n = nghost(normal(1))
    ng_u = nghost(umac(1,1))
    ng_w0 = nghost(gradw0_cart(1))

    do n=1,nlevs
       do i=1, nboxes(vel_force(n))
          if ( multifab_remote(vel_force(n),i) ) cycle
          fp   => dataptr(vel_force(n),i)
          ump  => dataptr(umac(n,1),i)
          lo   =  lwb(get_box(vel_force(n),i))
          hi   =  upb(get_box(vel_force(n),i))
          select case (dm)
          case (1)
          call add_utilde_force_1d(n, lo, hi, fp(:,1,1,:), ng_f, &
                                   ump(:,1,1,1), ng_u, w0(n,:))
          case (2)
          vmp  => dataptr(umac(n,2),i)
          call add_utilde_force_2d(n, lo, hi, fp(:,:,1,:), ng_f, &
                                   ump(:,:,1,1), vmp(:,:,1,1), ng_u, w0(n,:))
          case (3)
          vmp  => dataptr(umac(n,2),i)
          wmp  => dataptr(umac(n,3),i)
          gwp  => dataptr(gradw0_cart(n),i)
          np   => dataptr(normal(n),i)
          if (spherical .eq. 1) then
             n_1d = 1
          else
             n_1d = n
          end if
          call add_utilde_force_3d(n, lo, hi, fp(:,:,:,:), ng_f, np(:,:,:,:), ng_n, &
                                   ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_u, &
                                   gwp(:,:,:,1), ng_w0, w0(n_1d,:))
          end select
       end do
    end do

    do n=1,nlevs
       call destroy(gradw0_cart(n))
    end do

    call destroy(bpt)

  end subroutine add_utilde_force

  subroutine add_utilde_force_1d(n,lo,hi,force,ng_f,umac,ng_u,w0)

    use geometry, only: spherical, nr, dr
    use bl_constants_module

    integer,         intent(in   ) :: n,lo(:),hi(:),ng_f,ng_u
    real(kind=dp_t), intent(inout) ::       force(lo(1)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::        umac(lo(1)-ng_u :)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    
    integer :: i

    do i=lo(1)-1,hi(1)+1

       if (i .eq. -1) then
          ! do not modify force since dw0/dr=0
       else if (i .eq. nr(n)) then
          ! do not modify force since dw0/dr=0
       else
          force(i,1) = force(i,1) &
               - (umac(i+1)+umac(i))*(w0(i+1)-w0(i)) / (TWO*dr(n))
       end if

    end do

  end subroutine add_utilde_force_1d

  subroutine add_utilde_force_2d(n,lo,hi,force,ng_f,umac,vmac,ng_u,w0)

    use geometry, only: spherical, nr, dr
    use bl_constants_module

    integer,         intent(in   ) :: n,lo(:),hi(:),ng_f,ng_u
    real(kind=dp_t), intent(inout) ::       force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::        umac(lo(1)-ng_u :,lo(2)-ng_u :)
    real(kind=dp_t), intent(in   ) ::        vmac(lo(1)-ng_u :,lo(2)-ng_u :)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    
    integer :: i,j


    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1

          if (j .eq. -1) then
             ! do not modify force since dw0/dr=0
          else if (j .eq. nr(n)) then
             ! do not modify force since dw0/dr=0
          else
             force(i,j,2) = force(i,j,2) &
                  - (vmac(i,j+1)+vmac(i,j))*(w0(j+1)-w0(j)) / (TWO*dr(n))
          end if

       end do
    end do

  end subroutine add_utilde_force_2d

  subroutine add_utilde_force_3d(n,lo,hi,force,ng_f,normal,ng_n,umac,vmac,wmac,ng_u, &
                                 gradw0_cart,ng_w0,w0)

    use geometry, only: spherical, nr, dr
    use bl_constants_module

    integer,         intent(in   ) :: n,lo(:),hi(:),ng_f,ng_n,ng_u,ng_w0
    real(kind=dp_t), intent(inout) ::       force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::      normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:)
    real(kind=dp_t), intent(in   ) ::        umac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
    real(kind=dp_t), intent(in   ) ::        vmac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
    real(kind=dp_t), intent(in   ) ::        wmac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
    real(kind=dp_t), intent(in   ) :: gradw0_cart(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)

    integer :: i,j,k

    real(kind=dp_t) :: Ut_dot_er

    if (spherical .eq. 1) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,Ut_dot_er)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                Ut_dot_er = &
                     HALF*(umac(i,j,k)+umac(i+1,j  ,k  ))*normal(i,j,k,1) + &
                     HALF*(vmac(i,j,k)+vmac(i  ,j+1,k  ))*normal(i,j,k,2) + &
                     HALF*(wmac(i,j,k)+wmac(i  ,j,  k+1))*normal(i,j,k,3)

                force(i,j,k,1) = force(i,j,k,1) - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,1)
                force(i,j,k,2) = force(i,j,k,2) - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,2)
                force(i,j,k,3) = force(i,j,k,3) - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,3)

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                if (k .eq. -1) then
                   ! do not modify force since dw0/dr=0
                else if (k .eq. nr(n)) then
                   ! do not modify force since dw0/dr=0
                else
                   force(i,j,k,3) = force(i,j,k,3) &
                        - (wmac(i,j,k+1)+wmac(i,j,k))*(w0(k+1)-w0(k)) / (TWO*dr(n))
                end if

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    end if

  end subroutine add_utilde_force_3d

end module mk_vel_force_module
