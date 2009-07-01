module mk_vel_force_module

  ! Compute the force that appears in the velocity (or momentum)
  ! equations.  This is used both when predicting the interface
  ! states and in the final, conservative update.

  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private
  public :: mk_vel_force, add_w0_force

contains

  subroutine mk_vel_force(vel_force,uold,gpres,s,normal,rho0,grav,dx,the_bc_level,mla)

    use bl_prof_module
    use geometry, only: spherical, dm, nlevs
    use variables, only: foextrap_comp, rho_comp
    use bl_constants_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module
    use multifab_physbc_module

    type(multifab) , intent(inout) :: vel_force(:)
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: gpres(:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    ! Local variables
    real(kind=dp_t), pointer:: uop(:,:,:,:)
    real(kind=dp_t), pointer:: gpp(:,:,:,:)
    real(kind=dp_t), pointer:: fp(:,:,:,:)
    real(kind=dp_t), pointer:: rp(:,:,:,:)
    real(kind=dp_t), pointer:: np(:,:,:,:)
    integer                 :: i,lo(dm),hi(dm),ng_s,ng_f,ng_n,ng_gp,n,ng_uo

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mk_vel_force")

    ng_s = s(1)%ng
    ng_f = vel_force(1)%ng
    ng_gp = gpres(1)%ng

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
             call mk_vel_force_2d(fp(:,:,1,:),ng_f,gpp(:,:,1,:),ng_gp,rp(:,:,1,rho_comp), &
                                  ng_s,rho0(n,:),grav(n,:),lo,hi)
          case (3)
             ng_uo = uold(1)%ng
             uop => dataptr(uold(n),i)

             if (spherical .eq. 1) then
                ng_n = normal(1)%ng
                np => dataptr(normal(n), i)
                call mk_vel_force_3d_sphr(fp(:,:,:,:),ng_f,uop(:,:,:,:),ng_uo, &
                                          gpp(:,:,:,:),ng_gp,rp(:,:,:,rho_comp),ng_s, &
                                          np(:,:,:,:),ng_n,rho0(1,:),grav(1,:),lo,hi,dx(n,:))
             else
                call mk_vel_force_3d_cart(fp(:,:,:,:),ng_f,uop(:,:,:,:),ng_uo, &
                                          gpp(:,:,:,:),ng_gp,rp(:,:,:,rho_comp),ng_s, &
                                          rho0(n,:),grav(n,:),lo,hi)
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
                                            ng_f,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1),the_bc_level(n), &
                                            i,foextrap_comp,1,fill_crse_input=.false.)
          end do

       enddo

    end if

    call destroy(bpt)

  end subroutine mk_vel_force

  subroutine mk_vel_force_2d(vel_force,ng_f,gpres,ng_gp,rho,ng_s,rho0,grav,lo,hi)

    use variables, only: rho_comp
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

  subroutine mk_vel_force_3d_cart(vel_force,ng_f,uold,ng_uo,gpres,ng_gp,rho,ng_s,rho0,grav,lo,hi)

    use variables, only: rho_comp
    use geometry,  only: sin_theta, cos_theta, omega, centrifugal_term
    use bl_constants_module

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s, ng_uo
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::      uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
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

             coriolis_term(1) = -TWO * omega * uold(i,j,k,2) * cos_theta
             coriolis_term(2) =  TWO * omega * (uold(i,j,k,3) * sin_theta + uold(i,j,k,1) * cos_theta)
             coriolis_term(3) = -TWO * omega * uold(i,j,k,2) * sin_theta

             vel_force(i,j,k,1) = -coriolis_term(1) - centrifugal_term(1) - gpres(i,j,k,1) / rho(i,j,k) 
             vel_force(i,j,k,2) = -coriolis_term(2) - centrifugal_term(2) - gpres(i,j,k,2) / rho(i,j,k) 
             vel_force(i,j,k,3) = -coriolis_term(3) - centrifugal_term(3) + ( rhopert * grav(k) - gpres(i,j,k,3) ) / rho(i,j,k)

          end do
       end do
    end do

  end subroutine mk_vel_force_3d_cart

  subroutine mk_vel_force_3d_sphr(vel_force,ng_f,uold,ng_uo,gpres,ng_gp,rho,ng_s,normal, &
                                  ng_n,rho0,grav,lo,hi,dx)

    use variables, only: rho_comp
    use fill_3d_module
    use bl_constants_module
    use geometry,  only: omega, center

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_gp,ng_s,ng_n,ng_uo
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::      uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
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

             distance = sqrt(xx**2 + yy**2 + zz**2)
             cos_theta = normal(i,j,k,3)

             ! omega x (omega x r ) = - omega^2 x e_x  - omega^2 y e_y
             centrifugal_term(1) = -omega * omega * distance * normal(i,j,k,1)
             centrifugal_term(2) = -omega * omega * distance * normal(i,j,k,2)
             centrifugal_term(3) = ZERO

             ! 2 omega x U = - 2 omega v e_x  + 2 omega u e_y
             coriolis_term(1) = -TWO * omega * uold(i,j,k,2)
             coriolis_term(2) =  TWO * omega * uold(i,j,k,1)
             coriolis_term(3) = ZERO

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
    use variables, only: foextrap_comp, rho_comp
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

    use variables, only: rho_comp
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

    use variables, only: rho_comp
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

    use variables, only: rho_comp
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
