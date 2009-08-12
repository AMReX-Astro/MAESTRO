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
  public :: mk_vel_force, add_w0_force

contains

  subroutine mk_vel_force(vel_force,is_final_update, &
                          uold,umac,w0,gpres,s,index_rho,normal, &
                          rho0,grav,dx,the_bc_level,mla)

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
    use probin_module, only: edge_nodal_flag
    use fill_3d_module, only : make_w0mac, put_1d_array_on_cart

    type(multifab) , intent(inout) :: vel_force(:)
    logical        , intent(in   ) :: is_final_update
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: gpres(:)
    type(multifab) , intent(in   ) :: s(:)
    integer                        :: index_rho  
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
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
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    real(kind=dp_t), pointer ::  gpp(:,:,:,:)
    real(kind=dp_t), pointer ::   fp(:,:,:,:)
    real(kind=dp_t), pointer ::   rp(:,:,:,:)
    real(kind=dp_t), pointer ::   np(:,:,:,:)
    integer                  :: i,comp,lo(dm),hi(dm)
    integer                  :: ng_s,ng_f,ng_n,ng_gp,n,ng_uo,ng_um
   
    type(multifab) :: w0_cart(mla%nlevel)
    type(multifab) :: w0mac(mla%nlevel,dm)
    integer :: ng_wc, ng_wm

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mk_vel_force")

    ng_s = s(1)%ng
    ng_f = vel_force(1)%ng
    ng_gp = gpres(1)%ng


    ! put w0 on cart both cell-centered and on edges
    if (spherical .eq. 1) then
       do n=1,nlevs
          do comp=1,dm
             ! w0mac will contain an edge-centered w0 on a Cartesian grid,                                          
             ! for use in computing the Coriolis term in the final update
             call multifab_build(w0mac(n,comp), mla%la(n),1,1,nodal=edge_nodal_flag(comp,:))
             call setval(w0mac(n,comp), ZERO, all=.true.)
          enddo

          ! w0_cart will contain the cell-centered Cartesian components
          ! of w0, for use in computing the Coriolis term in the prediction
          call build(w0_cart(n),mla%la(n),dm,0)
          call setval(w0_cart(n), ZERO, all=.true.)
          
       enddo

       ! fill the edge-centered w0mac
       call make_w0mac(mla,w0,w0mac,dx,the_bc_level)          
       ng_wm = w0mac(1,1)%ng

       ! fill the all dm components of the cell-centered w0_cart
       call put_1d_array_on_cart(w0,w0_cart,foextrap_comp,.true.,.true.,dx,the_bc_level,mla)
       ng_wc = w0_cart(1)%ng

    endif

    do n=1,nlevs
       do i=1,s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          fp  => dataptr(vel_force(n),i)
          gpp => dataptr(gpres(n),i)
          rp  => dataptr(s(n),i)
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call mk_vel_force_2d(fp(:,:,1,:),ng_f,gpp(:,:,1,:),ng_gp,rp(:,:,1,index_rho), &
                                  ng_s,rho0(n,:),grav(n,:),lo,hi)
          case (3)
             ng_uo = uold(1)%ng
             uop => dataptr(uold(n),i)

             ng_um = umac(1,1)%ng
             ump => dataptr(umac(n,1),i)
             vmp => dataptr(umac(n,2),i)
             wmp => dataptr(umac(n,3),i)

             if (spherical .eq. 1) then
                ng_n = normal(1)%ng
                np => dataptr(normal(n), i)
                w0cp  => dataptr(w0_cart(n), i)
                w0xp  => dataptr(w0mac(n,1),i)
                w0yp  => dataptr(w0mac(n,2),i)

                call mk_vel_force_3d_sphr(fp(:,:,:,:),ng_f,is_final_update, &
                                          uop(:,:,:,:),ng_uo, &
                                          ump(:,:,:,1),vmp(:,:,:,1),ng_um, &
                                          w0cp(:,:,:,:),ng_wc, &
                                          w0xp(:,:,:,1),w0yp(:,:,:,1),ng_wm, &
                                          gpp(:,:,:,:),ng_gp,rp(:,:,:,index_rho),ng_s, &
                                          np(:,:,:,:),ng_n,rho0(1,:),grav(1,:),lo,hi,dx(n,:))
             else
                call mk_vel_force_3d_cart(fp(:,:,:,:),ng_f,is_final_update, &
                                          uop(:,:,:,:),ng_uo, &
                                          ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),ng_um, &
                                          w0(n,:), &
                                          gpp(:,:,:,:),ng_gp,rp(:,:,:,index_rho),ng_s, &
                                          rho0(n,:),grav(n,:),lo,hi)
             end if
          end select
       end do
    enddo

    if (spherical .eq. 1) then
       do n=1,nlevs
          do comp=1,dm
             call destroy(w0mac(n,comp))
          end do
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

  subroutine mk_vel_force_2d(vel_force,ng_f,gpres,ng_gp,rho,ng_s,rho0,grav,lo,hi)

    use bl_constants_module

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::     gpres(lo(1)-ng_gp:,lo(2)-ng_gp:,:)
    real(kind=dp_t), intent(in   ) ::       rho(lo(1)-ng_s :,lo(2)-ng_s :)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)

    integer         :: i,j
    real(kind=dp_t) :: rhopert

    vel_force = ZERO

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          rhopert = rho(i,j) - rho0(j)

          vel_force(i,j,1) = - gpres(i,j,1) / rho(i,j)
          vel_force(i,j,2) =  rhopert / rho(i,j) * grav(j) &
               - gpres(i,j,2) / rho(i,j)
       end do
    end do

  end subroutine mk_vel_force_2d

  subroutine mk_vel_force_3d_cart(vel_force,ng_f,is_final_update, &
                                  uold,ng_uo, &
                                  umac,vmac,wmac,ng_um, &
                                  w0, &
                                  gpres,ng_gp,rho,ng_s, &
                                  rho0,grav,lo,hi)

    use geometry,  only: sin_theta, cos_theta, omega, centrifugal_term
    use bl_constants_module

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s, ng_uo, ng_um
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    logical        , intent(in   ) :: is_final_update
    real(kind=dp_t), intent(in   ) ::      uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
    real(kind=dp_t), intent(in   ) ::      umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::      vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::      wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) ::     gpres(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
    real(kind=dp_t), intent(in   ) ::       rho(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)

    integer         :: i,j,k
    real(kind=dp_t) :: rhopert

    real(kind=dp_t) :: coriolis_term(3)

    vel_force = ZERO

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             rhopert = rho(i,j,k) - rho0(k)
             
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
                  gpres(i,j,k,1) / rho(i,j,k) 

             vel_force(i,j,k,2) = -coriolis_term(2) - centrifugal_term(2) - &
                  gpres(i,j,k,2) / rho(i,j,k) 

             vel_force(i,j,k,3) = -coriolis_term(3) - centrifugal_term(3) + &
                  ( rhopert * grav(k) - gpres(i,j,k,3) ) / rho(i,j,k)

          end do
       end do
    end do

  end subroutine mk_vel_force_3d_cart

  subroutine mk_vel_force_3d_sphr(vel_force,ng_f,is_final_update, &
                                  uold,ng_uo, &
                                  umac,vmac,ng_um, &
                                  w0_cart,ng_wc, &
                                  w0macx,w0macy,ng_wm, &
                                  gpres,ng_gp,rho,ng_s, &
                                  normal,ng_n,rho0,grav,lo,hi,dx)

    use fill_3d_module
    use bl_constants_module
    use geometry,  only: omega, center
    use probin_module, only: base_cutoff_density

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_gp,ng_s,ng_n,ng_uo,ng_um,ng_wc,ng_wm
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    logical        , intent(in   ) :: is_final_update
    real(kind=dp_t), intent(in   ) ::      uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
    real(kind=dp_t), intent(in   ) ::      umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::      vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   w0_cart(lo(1)-ng_wc:,lo(2)-ng_wc:,lo(3)-ng_wc:,:)
    real(kind=dp_t), intent(in   ) ::    w0macx(lo(1)-ng_wm:,lo(2)-ng_wm:,lo(3)-ng_wm:)
    real(kind=dp_t), intent(in   ) ::    w0macy(lo(1)-ng_wm:,lo(2)-ng_wm:,lo(3)-ng_wm:)
    real(kind=dp_t), intent(in   ) ::     gpres(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
    real(kind=dp_t), intent(in   ) ::       rho(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::    normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) ::   dx(:)

    integer         :: i,j,k

    real(kind=dp_t), allocatable :: rho0_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: grav_cart(:,:,:,:)

    real(kind=dp_t) :: rhopert
    real(kind=dp_t) :: xx, yy, zz, distance, cos_theta
    real(kind=dp_t) :: centrifugal_term(3), coriolis_term(3)

    allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(grav_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3))

    vel_force = ZERO

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,rho0,rho0_cart,lo,hi,dx,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.true.,grav,grav_cart,lo,hi,dx,0)

    do k = lo(3),hi(3)
       zz = (dble(k) + HALF)*dx(3) - center(3)
       do j = lo(2),hi(2)
          yy = (dble(j) + HALF)*dx(2) - center(2)
          do i = lo(1),hi(1)
             xx = (dble(i) + HALF)*dx(1) - center(1)

             rhopert = rho(i,j,k) - rho0_cart(i,j,k,1)

             if (rho(i,j,k) .lt. 5.d0*base_cutoff_density) then
                rhopert = 0.d0
             end if

             distance = sqrt(xx**2 + yy**2 + zz**2)
             cos_theta = normal(i,j,k,3)

             ! omega x (omega x r ) = - omega^2 x e_x  - omega^2 y e_y    
             ! (with omega = omega e_z)
             centrifugal_term(1) = -omega * omega * distance * normal(i,j,k,1)
             centrifugal_term(2) = -omega * omega * distance * normal(i,j,k,2)
             centrifugal_term(3) = ZERO

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

             vel_force(i,j,k,1) = -coriolis_term(1) - centrifugal_term(1) + &
                  ( rhopert * grav_cart(i,j,k,1) - gpres(i,j,k,1) ) / rho(i,j,k)

             vel_force(i,j,k,2) = -coriolis_term(2) - centrifugal_term(2) + &
                  ( rhopert * grav_cart(i,j,k,2) - gpres(i,j,k,2) ) / rho(i,j,k)

             vel_force(i,j,k,3) = -coriolis_term(3) - centrifugal_term(3) + &
                  ( rhopert * grav_cart(i,j,k,3) - gpres(i,j,k,3) ) / rho(i,j,k)

          end do
       end do
    end do

    deallocate(rho0_cart,grav_cart)

  end subroutine mk_vel_force_3d_sphr

  subroutine add_w0_force(vel_force,w0_force,w0_force_cart,the_bc_level,mla)

    use bl_prof_module
    use geometry, only: spherical, dm, nlevs
    use variables, only: foextrap_comp
    use bl_constants_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module
    use multifab_physbc_module

    type(multifab) , intent(inout) :: vel_force(:)
    real(kind=dp_t), intent(in   ) :: w0_force(:,0:)
    type(multifab) , intent(in   ) :: w0_force_cart(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    ! Local variables
    real(kind=dp_t), pointer:: fp(:,:,:,:)
    real(kind=dp_t), pointer:: w0p(:,:,:,:)
    integer                 :: i,n,ng_f,ng_w
    integer                 :: lo(dm),hi(dm)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "add_w0_force")

    ng_f = vel_force(1)%ng

    do n=1,nlevs
       do i=1,vel_force(n)%nboxes
          if ( multifab_remote(vel_force(n), i) ) cycle
          fp => dataptr(vel_force(n),i)
          lo = lwb(get_box(vel_force(n),i))
          hi = upb(get_box(vel_force(n),i))
          select case (dm)
          case (2)
             call add_w0_force_2d(fp(:,:,1,:),ng_f,w0_force(n,:),lo,hi)
          case (3)
             if (spherical .eq. 1) then
                ng_w = w0_force_cart(1)%ng
                w0p => dataptr(w0_force_cart(n), i)
                call add_w0_force_3d_sphr(fp(:,:,:,:),ng_f,w0p(:,:,:,:),ng_w,lo,hi)
             else
                call add_w0_force_3d_cart(fp(:,:,:,:),ng_f,w0_force(n,:),lo,hi)
             end if
          end select
       end do
    enddo

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
                                            1,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1),the_bc_level(n), &
                                            i,foextrap_comp,1,fill_crse_input=.false.)
          end do

       enddo

    end if

    call destroy(bpt)

  end subroutine add_w0_force

  subroutine add_w0_force_2d(vel_force,ng_f,w0_force,lo,hi)

    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_f
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f:,lo(2)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: w0_force(0:)

    integer         :: i,j

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          vel_force(i,j,2) = vel_force(i,j,2) - w0_force(j)
       end do
    end do

  end subroutine add_w0_force_2d

  subroutine add_w0_force_3d_cart(vel_force,ng_f,w0_force,lo,hi)

    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_f
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: w0_force(0:)

    integer         :: i,j,k

    do k=lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             vel_force(i,j,k,3) = vel_force(i,j,k,3) - w0_force(k)
          end do
       end do
    end do

  end subroutine add_w0_force_3d_cart

  subroutine add_w0_force_3d_sphr(vel_force,ng_f,w0_force_cart,ng_w,lo,hi)

    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_w
    real(kind=dp_t), intent(inout) :: vel_force    (lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(inout) :: w0_force_cart(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:,:)

    integer         :: i,j,k

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             vel_force(i,j,k,:) = vel_force(i,j,k,:) - w0_force_cart(i,j,k,:)
          end do
       end do
    end do

  end subroutine add_w0_force_3d_sphr

end module mk_vel_force_module
