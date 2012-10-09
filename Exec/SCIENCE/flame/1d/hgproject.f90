module hgproject_module

  use bl_types
  use mg_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private
  public :: hgproject

contains 

  subroutine hgproject(proj_type,mla,unew,uold,rhohalf,pres,gpres,dx,dt,the_bc_tower, &
                       press_comp, divu_rhs,div_coeff_1d,div_coeff_3d,eps_in)

    use bc_module
    use bl_prof_module
    use proj_parameters
    use nodal_divu_module
    use stencil_module
    use ml_solve_module
    use ml_restriction_module
    use multifab_fill_ghost_module
    use probin_module, only: verbose, mg_verbose, cg_verbose, hg_dense_stencil, nodal
    use geometry, only: dm, nlevs, spherical

    integer        , intent(in   ) :: proj_type
    type(ml_layout), intent(inout) :: mla
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: uold(:)
    type(multifab ), intent(inout) :: rhohalf(:)
    type(multifab ), intent(inout) :: pres(:)
    type(multifab ), intent(inout) :: gpres(:)
    real(dp_t)     , intent(in   ) :: dx(:,:),dt
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: press_comp

    type(multifab ), intent(inout), optional :: divu_rhs(:)
    real(dp_t)     , intent(in   ), optional :: div_coeff_1d(:,:)
    type(multifab ), intent(in   ), optional :: div_coeff_3d(:)
    real(dp_t)     , intent(in   ), optional :: eps_in

    ! Local  
    type(multifab) :: phi(mla%nlevel)
    type(multifab) :: gphi(mla%nlevel)

    integer                   :: n,stencil_type
    real(dp_t)                :: umin,umax,vmin,vmax,wmin,wmax
    logical                   :: use_div_coeff_1d, use_div_coeff_3d
    type(bl_prof_timer), save :: bpt

    call build(bpt, "hgproject")

    if (hg_dense_stencil) then
       stencil_type = ST_DENSE
    else
       stencil_type = ST_CROSS
    end if

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       print *,'PROJ_TYPE IN HGPROJECT:',proj_type
       select case (stencil_type)
       case (ST_CROSS)
          print *,'PROJ_TYPE IN HGPROJECT: ST_CROSS'
       case (ST_DENSE)
          print *,'PROJ_TYPE IN HGPROJECT: ST_DENSE'
       end select
    end if

    use_div_coeff_1d = .false.
    if (present(div_coeff_1d)) use_div_coeff_1d = .true.

    use_div_coeff_3d = .false.
    if (present(div_coeff_3d)) use_div_coeff_3d = .true.

    if (use_div_coeff_1d .and. use_div_coeff_3d) then
       call bl_error('CANT HAVE 1D and 3D DIV_COEFF IN HGPROJECT')
    end if

    if (spherical .eq. 1 .and. use_div_coeff_1d) then
       call bl_error('CANT HAVE SPHERICAL .eq. 1 and 1D DIV_COEFF IN HGPROJECT')
    end if

    if (spherical .eq. 0 .and. use_div_coeff_3d) then
       call bl_error('CANT HAVE SPHERICAL .eq. 0 and 3D DIV_COEFF IN HGPROJECT')
    end if

    if (verbose .ge. 1) then
       umin =  HUGE(umin)
       vmin =  HUGE(vmin)
       wmin =  HUGE(wmin)
       umax = -HUGE(umax)
       vmax = -HUGE(vmax)
       wmax = -HUGE(wmax)
       do n = 1, nlevs
          umin = min(umin,multifab_min_c(unew(n),1))
          umax = max(umax,multifab_max_c(unew(n),1))
          vmin = min(vmin,multifab_min_c(unew(n),2))
          vmax = max(vmax,multifab_max_c(unew(n),2))
          if (dm .eq. 3) then
             wmin = min(wmin,multifab_min_c(unew(n),3))
             wmax = max(wmax,multifab_max_c(unew(n),3))
          end if
       end do
       if (parallel_IOProcessor()) then
          write(6,1001) umin,umax
          write(6,1002) vmin,vmax
          if (dm .eq. 3) write(6,1003) wmin,wmax
          write(6,1004)
       end if
    end if

1001 format('... x-velocity before projection ',e17.10,2x,e17.10)
1002 format('... y-velocity before projection ',e17.10,2x,e17.10)
1003 format('... z-velocity before projection ',e17.10,2x,e17.10)
1004 format(' ')

    call create_uvec_for_projection(unew,uold,rhohalf,gpres,dt,the_bc_tower,proj_type)

    if (use_div_coeff_1d) then
       call mult_by_1d_coeff(unew,div_coeff_1d,.true.)
       call mult_by_1d_coeff(rhohalf,div_coeff_1d,.false.)
    else if (use_div_coeff_3d) then
       call mult_by_3d_coeff(unew,div_coeff_3d,.true.)
       call mult_by_3d_coeff(rhohalf,div_coeff_3d,.false.)
    end if

    if (present(divu_rhs)) then
       call enforce_outflow_on_divu_rhs(divu_rhs,the_bc_tower)
    end if

    do n = 1, nlevs
       call multifab_build(phi(n), mla%la(n), 1, 1, nodal)
       call setval(phi(n),ZERO,all=.true.)
    end do

    if (present(eps_in)) then
       call hg_multigrid(mla,unew,rhohalf,phi,dx,the_bc_tower, &
                         press_comp,stencil_type,divu_rhs,eps_in)
    else
       call hg_multigrid(mla,unew,rhohalf,phi,dx,the_bc_tower, &
                         press_comp,stencil_type,divu_rhs)
    end if

    if (use_div_coeff_1d) then
       call mult_by_1d_coeff(unew,div_coeff_1d,.false.)
       call mult_by_1d_coeff(rhohalf,div_coeff_1d,.true.)
    else if (use_div_coeff_3d) then
       call mult_by_3d_coeff(unew,div_coeff_3d,.false.)
       call mult_by_3d_coeff(rhohalf,div_coeff_3d,.true.)
    end if

    do n = 1, nlevs
       call multifab_build(gphi(n), mla%la(n), dm, 0)
    end do

    call mkgphi(gphi,phi,dx)

    call hg_update(proj_type,unew,uold,gpres,gphi,rhohalf,  &
                   pres,phi,dt,mla,the_bc_tower%bc_tower_array)

    do n = 1,nlevs
       call destroy(phi(n))
       call destroy(gphi(n))
    end do

    if (verbose .ge. 1) then
       umin =  HUGE(umin)
       vmin =  HUGE(vmin)
       wmin =  HUGE(wmin)
       umax = -HUGE(umax)
       vmax = -HUGE(vmax)
       wmax = -HUGE(wmax)
       do n = 1, nlevs
          umin = min(umin,multifab_min_c(unew(n),1))
          umax = max(umax,multifab_max_c(unew(n),1))
          vmin = min(vmin,multifab_min_c(unew(n),2))
          vmax = max(vmax,multifab_max_c(unew(n),2))
          if (dm .eq. 3) then
             wmin = min(wmin,multifab_min_c(unew(n),3))
             wmax = max(wmax,multifab_max_c(unew(n),3))
          end if
       end do
       if (parallel_IOProcessor()) then
          write(6,1101) umin,umax
          write(6,1102) vmin,vmax
          if (dm .eq. 3) write(6,1103) wmin,wmax
          write(6,1104)
       end if
    end if

1101 format('... x-velocity  after projection ',e17.10,2x,e17.10)
1102 format('... y-velocity  after projection ',e17.10,2x,e17.10)
1103 format('... z-velocity  after projection ',e17.10,2x,e17.10)
1104 format(' ')

    call destroy(bpt)

  contains

    ! ******************************************************************************* !

    subroutine create_uvec_for_projection(unew,uold,rhohalf,gpres,dt,the_bc_tower, &
                                          proj_type)

      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(in   ) :: uold(:)
      type(multifab) , intent(in   ) :: rhohalf(:)
      type(multifab) , intent(inout) :: gpres(:)
      real(kind=dp_t), intent(in   ) :: dt
      type(bc_tower) , intent(in   ) :: the_bc_tower
      integer        , intent(in   ) :: proj_type
  
      type(bc_level) :: bc
  
      real(kind=dp_t), pointer :: unp(:,:,:,:)
      real(kind=dp_t), pointer :: uop(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:)
  
      integer :: i,n
      integer :: ng_un,ng_uo,ng_rh,ng_gp

      type(bl_prof_timer), save :: bpt

      call build(bpt, "create_uvec_for_projection")
  
      ng_un = unew(1)%ng
      ng_uo = uold(1)%ng
      ng_rh = rhohalf(1)%ng
      ng_gp = gpres(1)%ng
  
      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, unew(n)%nboxes
            if ( multifab_remote(unew(n), i) ) cycle
            unp => dataptr(unew(n)     , i) 
            uop => dataptr(uold(n)     , i) 
            gpp => dataptr(gpres(n)       , i)
             rp => dataptr(  rhohalf(n), i)
            select case (dm)
               case (2)
                 call create_uvec_2d(unp(:,:,1,:), ng_un, uop(:,:,1,:), ng_uo, &
                                     rp(:,:,1,1), ng_rh, gpp(:,:,1,:), ng_gp, &
                                     dt,bc%phys_bc_level_array(i,:,:), proj_type)
               case (3)
                 call create_uvec_3d(unp(:,:,:,:), ng_un, uop(:,:,:,:), ng_uo, &
                                     rp(:,:,:,1), ng_rh, gpp(:,:,:,:), ng_gp, &
                                     dt, bc%phys_bc_level_array(i,:,:), proj_type)
            end select
         end do
         call multifab_fill_boundary(unew(n))
      end do 

      call destroy(bpt)

    end subroutine create_uvec_for_projection

    !   ******************************************************************************** !

    subroutine create_uvec_2d(unew,ng_un,uold,ng_uo,rhohalf,ng_rh,gpres,ng_gp, &
                              dt,phys_bc,proj_type)

      use proj_parameters

      integer        , intent(in   ) :: ng_un,ng_uo,ng_rh,ng_gp
      real(kind=dp_t), intent(inout) ::    unew(-ng_un:,-ng_un:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng_uo:,-ng_uo:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-ng_rh:,-ng_rh:)
      real(kind=dp_t), intent(inout) ::   gpres(-ng_gp:,-ng_gp:,:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      integer :: nx,ny
      nx = size(gpres,dim=1) - 2
      ny = size(gpres,dim=2) - 2

      if (phys_bc(1,1) .eq. INLET) gpres(-1,-1:ny,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gpres(nx,-1:ny,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gpres(-1:nx,-1,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gpres(-1:nx,ny,:) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection_comp) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters_comp) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters_comp) then

         unew(-1:nx,-1:ny,1) = ( unew(-1:nx,-1:ny,1) - uold(-1:nx,-1:ny,1) ) / dt
         unew(-1:nx,-1:ny,2) = ( unew(-1:nx,-1:ny,2) - uold(-1:nx,-1:ny,2) ) / dt
     
      ! quantity projected is Ustar + dt * (1/rho) gpres
      else if (proj_type .eq. regular_timestep_comp) then

         unew(-1:nx,-1:ny,1) = &
              unew(-1:nx,-1:ny,1) + dt*gpres(-1:nx,-1:ny,1)/rhohalf(-1:nx,-1:ny)
         unew(-1:nx,-1:ny,2) = &
              unew(-1:nx,-1:ny,2) + dt*gpres(-1:nx,-1:ny,2)/rhohalf(-1:nx,-1:ny)

       else
     
          call bl_error('No proj_type by this number ')

      end if

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL) unew(-1,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL) unew(nx,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL) unew(:,-1,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL) unew(:,ny,:) = ZERO

    end subroutine create_uvec_2d

    !  *********************************************************************************** !

    subroutine create_uvec_3d(unew,ng_un,uold,ng_uo,rhohalf,ng_rh,gpres,ng_gp, &
                              dt,phys_bc,proj_type)

      use proj_parameters

      integer        , intent(in   ) :: ng_un,ng_uo,ng_rh,ng_gp
      real(kind=dp_t), intent(inout) ::    unew(-ng_un:,-ng_un:,-ng_un:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng_uo:,-ng_uo:,-ng_uo:,:)
      real(kind=dp_t), intent(inout) ::   gpres(-ng_gp:,-ng_gp:,-ng_gp:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-ng_rh:,-ng_rh:,-ng_rh:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      integer :: nx,ny,nz

      nx = size(gpres,dim=1) - 2
      ny = size(gpres,dim=2) - 2
      nz = size(gpres,dim=3) - 2

      if (phys_bc(1,1) .eq. INLET) gpres(-1,-1:ny,-1:nz,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gpres(nx,-1:ny,-1:nz,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gpres(-1:nx,-1,-1:nz,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gpres(-1:nx,ny,-1:nz,:) = ZERO
      if (phys_bc(3,1) .eq. INLET) gpres(-1:nx,-1:ny,-1,:) = ZERO
      if (phys_bc(3,2) .eq. INLET) gpres(-1:nx,-1:ny,nz,:) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection_comp) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters_comp) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters_comp) then

         unew(-1:nx,-1:ny,-1:nz,:) =&
              ( unew(-1:nx,-1:ny,-1:nz,:) - uold(-1:nx,-1:ny,-1:nz,:) ) / dt

      ! quantity projected is Ustar + dt * (1/rho) gpres
      else if (proj_type .eq. regular_timestep_comp) then

         unew(-1:nx,-1:ny,-1:nz,1) = unew(-1:nx,-1:ny,-1:nz,1) + &
                                    dt*gpres(-1:nx,-1:ny,-1:nz,1)/rhohalf(-1:nx,-1:ny,-1:nz)
         unew(-1:nx,-1:ny,-1:nz,2) = unew(-1:nx,-1:ny,-1:nz,2) + &
                                    dt*gpres(-1:nx,-1:ny,-1:nz,2)/rhohalf(-1:nx,-1:ny,-1:nz)
         unew(-1:nx,-1:ny,-1:nz,3) = unew(-1:nx,-1:ny,-1:nz,3) + &
                                    dt*gpres(-1:nx,-1:ny,-1:nz,3)/rhohalf(-1:nx,-1:ny,-1:nz)

      else

          call bl_error('No proj_type by this number ')

      end if

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL) unew(-1,:,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL) unew(nx,:,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL) unew(:,-1,:,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL) unew(:,ny,:,:) = ZERO
      if (phys_bc(3,1)==SLIP_WALL .or. phys_bc(3,1)==NO_SLIP_WALL) unew(:,:,-1,:) = ZERO
      if (phys_bc(3,2)==SLIP_WALL .or. phys_bc(3,2)==NO_SLIP_WALL) unew(:,:,nz,:) = ZERO

    end subroutine create_uvec_3d

    ! ******************************************************************************* !

    subroutine mkgphi(gphi,phi,dx)

      type(multifab), intent(inout) :: gphi(:)
      type(multifab), intent(in   ) :: phi(:)
      real(dp_t) :: dx(:,:)

      integer :: i,n,ng_p,ng_gp

      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: pp(:,:,:,:) 

      type(bl_prof_timer), save :: bpt

      call build(bpt, "mkgphi")

      ng_p = phi(1)%ng
      ng_gp = gphi(1)%ng

      do n = 1, nlevs

         do i = 1, phi(n)%nboxes
            if ( multifab_remote(phi(n),i) ) cycle
            gph => dataptr(gphi(n),i)
            pp  => dataptr(phi(n),i)
            select case (dm)
            case (2)
               call mkgphi_2d(gph(:,:,1,:), ng_gp, pp(:,:,1,1), ng_p, dx(n,:))
            case (3)
               call mkgphi_3d(gph(:,:,:,:), ng_gp, pp(:,:,:,1), ng_p, dx(n,:))
            end select
         end do

      end do

      call destroy(bpt)

    end subroutine mkgphi

    !   ********************************************************************************* !

    subroutine mkgphi_2d(gphi,ng_gp,phi,ng_p,dx)

      integer        , intent(in   ) :: ng_gp,ng_p
      real(kind=dp_t), intent(inout) :: gphi(-ng_gp:,-ng_gp:,:)
      real(kind=dp_t), intent(inout) ::  phi(-ng_p :,-ng_p :)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,nx,ny

      nx = size(gphi,dim=1)
      ny = size(gphi,dim=2)

      do j = 0,ny-1
         do i = 0,nx-1
            gphi(i,j,1) = HALF*(phi(i+1,j) + phi(i+1,j+1) - &
                 phi(i  ,j) - phi(i  ,j+1) ) /dx(1)
            gphi(i,j,2) = HALF*(phi(i,j+1) + phi(i+1,j+1) - &
                 phi(i,j  ) - phi(i+1,j  ) ) /dx(2)
         end do
      end do

    end subroutine mkgphi_2d

    !   ******************************************************************************** !

    subroutine mkgphi_3d(gphi,ng_gp,phi,ng_p,dx)

      integer        , intent(in   ) :: ng_gp,ng_p
      real(kind=dp_t), intent(inout) :: gphi(-ng_gp:,-ng_gp:,-ng_gp:,:)
      real(kind=dp_t), intent(inout) ::  phi(-ng_p :,-ng_p :,-ng_p :)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k,nx,ny,nz

      nx = size(gphi,dim=1)
      ny = size(gphi,dim=2)
      nz = size(gphi,dim=3)

      do k = 0,nz-1
         do j = 0,ny-1
            do i = 0,nx-1
               gphi(i,j,k,1) = FOURTH*(phi(i+1,j,k  ) + phi(i+1,j+1,k  ) &
                    +phi(i+1,j,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i  ,j,k  ) - phi(i  ,j+1,k  ) &
                    -phi(i  ,j,k+1) - phi(i  ,j+1,k+1) ) /dx(1)
               gphi(i,j,k,2) = FOURTH*(phi(i,j+1,k  ) + phi(i+1,j+1,k  ) &
                    +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                    -phi(i,j  ,k+1) - phi(i+1,j  ,k+1) ) /dx(2)
               gphi(i,j,k,3) = FOURTH*(phi(i,j  ,k+1) + phi(i+1,j  ,k+1) &
                    +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                    -phi(i,j+1,k  ) - phi(i+1,j+1,k  ) ) /dx(3)
            end do
         end do
      end do

    end subroutine mkgphi_3d

    ! ******************************************************************************* !

    subroutine hg_update(proj_type,unew,uold,gpres,gphi,rhohalf,pres,phi,dt, &
                         mla,the_bc_level)

      use multifab_physbc_module
      use variables, only: foextrap_comp

      integer        , intent(in   ) :: proj_type
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(in   ) :: uold(:)
      type(multifab) , intent(inout) :: gpres(:)
      type(multifab) , intent(in   ) :: gphi(:)
      type(multifab) , intent(in   ) :: rhohalf(:)
      type(multifab) , intent(inout) :: pres(:)
      type(multifab) , intent(in   ) :: phi(:)
      real(kind=dp_t), intent(in   ) :: dt
      type(ml_layout), intent(inout) :: mla
      type(bc_level) , intent(in   ) :: the_bc_level(:)

      ! local
      integer :: i,n
      integer :: ng_un,ng_uo,ng_gp,ng_gh,ng_rh,ng_p,ng_h

      real(kind=dp_t), pointer :: upn(:,:,:,:) 
      real(kind=dp_t), pointer :: uon(:,:,:,:) 
      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:) 
      real(kind=dp_t), pointer ::  ph(:,:,:,:) 
      real(kind=dp_t), pointer ::  pp(:,:,:,:) 

      type(bl_prof_timer), save :: bpt

      call build(bpt, "hg_update")

      ng_un = unew(1)%ng
      ng_uo = uold(1)%ng
      ng_gp = gpres(1)%ng
      ng_gh = gphi(1)%ng
      ng_rh = rhohalf(1)%ng
      ng_p = pres(1)%ng
      ng_h = phi(1)%ng

      do n = 1, nlevs

         do i = 1, unew(n)%nboxes
            if ( multifab_remote(unew(n),i) ) cycle
            upn => dataptr(unew(n),i)
            uon => dataptr(uold(n),i)
            gpp => dataptr(gpres(n),i)
            gph => dataptr(gphi(n),i)
            rp  => dataptr(rhohalf(n),i)
            pp  => dataptr(pres(n),i)
            ph  => dataptr(phi(n),i)
            select case (dm)
            case (2)
               call hg_update_2d(proj_type, upn(:,:,1,:), ng_un, uon(:,:,1,:), ng_uo, &
                                 gpp(:,:,1,:), ng_gp, gph(:,:,1,:), ng_gh, rp(:,:,1,1), &
                                 ng_rh, pp(:,:,1,1), ng_p, ph(:,:,1,1), ng_h, dt)
            case (3)
               call hg_update_3d(proj_type, upn(:,:,:,:), ng_un, uon(:,:,:,:), ng_uo, &
                                 gpp(:,:,:,:), ng_gp, gph(:,:,:,:), ng_gh, rp(:,:,:,1), &
                                 ng_rh, pp(:,:,:,1), ng_p, ph(:,:,:,1), ng_h, dt)
            end select
         end do

         ! fill ghost cells for two adjacent grids at the same level
         ! this includes periodic domain boundary ghost cells
         call multifab_fill_boundary(pres(n))

      end do

      if (nlevs .eq. 1) then

         ! fill ghost cells for two adjacent grids at the same level
         ! this includes periodic domain boundary ghost cells
         call multifab_fill_boundary(unew(nlevs))
         call multifab_fill_boundary(gpres(nlevs))

         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(unew(nlevs),1,1,dm,the_bc_level(nlevs))
         do i=1,dm
            call multifab_physbc(gpres(nlevs),i,foextrap_comp,1,the_bc_level(nlevs))
         end do

      else

         ! the loop over nlevs must count backwards to make sure the finer grids are 
         ! done first
         do n=nlevs,2,-1

            ! set level n-1 data to be the average of the level n data covering it
            call ml_cc_restriction(unew(n-1),unew(n),mla%mba%rr(n-1,:)) 
            call ml_cc_restriction(gpres(n-1),gpres(n),mla%mba%rr(n-1,:))

            ! fill level n ghost cells using interpolation from level n-1 data
            ! note that multifab_fill_boundary and multifab_physbc are called for
            ! both levels n-1 and n
            call multifab_fill_ghost_cells(unew(n),unew(n-1),unew(n)%ng,mla%mba%rr(n-1,:), &
                                           the_bc_level(n-1),the_bc_level(n),1,1,dm, &
                                           fill_crse_input=.false.)
            do i=1,dm
               call multifab_fill_ghost_cells(gpres(n),gpres(n-1),1,mla%mba%rr(n-1,:), &
                                              the_bc_level(n-1),the_bc_level(n),i, &
                                              foextrap_comp,1,fill_crse_input=.false.)
            end do

         end do

      end if

      call destroy(bpt)

    end subroutine hg_update

    !   ****************************************************************************** !

    subroutine hg_update_2d(proj_type,unew,ng_un,uold,ng_uo,gpres,ng_gp,gphi,ng_gh, &
                            rhohalf,ng_rh,pres,ng_p,phi,ng_h,dt)

      use proj_parameters

      integer        , intent(in   ) :: ng_un,ng_uo,ng_gp,ng_gh,ng_rh,ng_p,ng_h
      integer        , intent(in   ) :: proj_type
      real(kind=dp_t), intent(inout) ::    unew(-ng_un:,-ng_un:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng_uo:,-ng_uo:,:)
      real(kind=dp_t), intent(inout) ::   gpres(-ng_gp:,-ng_gp:,:)
      real(kind=dp_t), intent(in   ) ::    gphi(-ng_gh:,-ng_gh:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-ng_rh:,-ng_rh:)
      real(kind=dp_t), intent(inout) ::    pres(-ng_p :,-ng_p :)
      real(kind=dp_t), intent(in   ) ::     phi(-ng_h :,-ng_h :)
      real(kind=dp_t), intent(in   ) :: dt

      integer         :: nx,ny
      integer j

      nx = size(gphi,dim=1)-1
      ny = size(gphi,dim=2)-1

      !     Subtract off the density-weighted gradient.
      unew(0:nx,0:ny,1) = unew(0:nx,0:ny,1) - gphi(0:nx,0:ny,1)/rhohalf(0:nx,0:ny) 
      unew(0:nx,0:ny,2) = unew(0:nx,0:ny,2) - gphi(0:nx,0:ny,2)/rhohalf(0:nx,0:ny) 

      unew(:,:,1) = ZERO
      do j = 0, ny
         unew(:,j,2) = unew(0,j,2)
      enddo

      if (proj_type .eq. pressure_iters_comp) &    ! unew held the projection of (ustar-uold)
           unew(0:nx,0:ny,:) = uold(0:nx,0:ny,:) + dt * unew(0:nx,0:ny,:)

      if ( (proj_type .eq. initial_projection_comp) .or. &
           (proj_type .eq. divu_iters_comp) ) then

         gpres = ZERO
         pres = ZERO

      else if (proj_type .eq. pressure_iters_comp) then

         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         gpres(0:nx,0:ny,:) = gpres(0:nx,0:ny,:) + gphi(0:nx,0:ny,:)
         pres(0:nx,0:ny  ) =  pres(0:nx,0:ny  )  +  phi(0:nx,0:ny  )

      else if (proj_type .eq. regular_timestep_comp) then

         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         gpres(0:nx,0:ny,:) = (ONE/dt) * gphi(0:nx,0:ny,:)
         pres(0:nx,0:ny  ) = (ONE/dt) * phi(0:nx,0:ny  )

      end if

    end subroutine hg_update_2d

    !   ******************************************************************************* !

    subroutine hg_update_3d(proj_type,unew,ng_un,uold,ng_uo,gpres,ng_gp,gphi,ng_gh, &
                            rhohalf,ng_rh,pres,ng_p,phi,ng_h,dt)

      use proj_parameters

      integer        , intent(in   ) :: ng_un,ng_uo,ng_gp,ng_gh,ng_rh,ng_p,ng_h
      integer        , intent(in   ) :: proj_type
      real(kind=dp_t), intent(inout) ::    unew(-ng_un:,-ng_un:,-ng_un:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng_uo:,-ng_uo:,-ng_uo:,:)
      real(kind=dp_t), intent(inout) ::   gpres(-ng_gp:,-ng_gp:,-ng_gp:,:)
      real(kind=dp_t), intent(in   ) ::    gphi(-ng_gh:,-ng_gh:,-ng_gh:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-ng_rh:,-ng_rh:,-ng_rh:)
      real(kind=dp_t), intent(inout) ::    pres(-ng_p :,-ng_p :,-ng_p :)
      real(kind=dp_t), intent(in   ) ::     phi(-ng_h :,-ng_h :,-ng_h :)
      real(kind=dp_t), intent(in   ) :: dt

      integer         :: nx,ny,nz

      nx = size(gphi,dim=1)-1
      ny = size(gphi,dim=2)-1
      nz = size(gphi,dim=3)-1

      !     Subtract off the density-weighted gradient.
      unew(0:nx,0:ny,0:nz,1) = &
           unew(0:nx,0:ny,0:nz,1) - gphi(0:nx,0:ny,0:nz,1)/rhohalf(0:nx,0:ny,0:nz) 
      unew(0:nx,0:ny,0:nz,2) = &
           unew(0:nx,0:ny,0:nz,2) - gphi(0:nx,0:ny,0:nz,2)/rhohalf(0:nx,0:ny,0:nz) 
      unew(0:nx,0:ny,0:nz,3) = &
           unew(0:nx,0:ny,0:nz,3) - gphi(0:nx,0:ny,0:nz,3)/rhohalf(0:nx,0:ny,0:nz) 

      if (proj_type .eq. pressure_iters_comp) &    ! unew held the projection of (ustar-uold)
           unew(0:nx,0:ny,0:nz,:) = uold(0:nx,0:ny,0:nz,:) + dt * unew(0:nx,0:ny,0:nz,:)

      if ( (proj_type .eq. initial_projection_comp) .or. &
           (proj_type .eq. divu_iters_comp) ) then

         gpres = ZERO
         pres = ZERO

      else if (proj_type .eq. pressure_iters_comp) then

         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         gpres(0:nx,0:ny,0:nz,:) = gpres(0:nx,0:ny,0:nz,:) + gphi(0:nx,0:ny,0:nz,:)
         pres(0:nx,0:ny,0:nz  ) =  pres(0:nx,0:ny,0:nz  )  +  phi(0:nx,0:ny,0:nz  )

      else if (proj_type .eq. regular_timestep_comp) then

         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         gpres(0:nx,0:ny,0:nz,:) = (ONE/dt) * gphi(0:nx,0:ny,0:nz,:)
         pres(0:nx,0:ny,0:nz  ) = (ONE/dt) *  phi(0:nx,0:ny,0:nz  )

      end if

    end subroutine hg_update_3d

    ! ********************************************************************************** !

    subroutine enforce_outflow_on_divu_rhs(divu_rhs,the_bc_tower)

      type(multifab) , intent(inout) :: divu_rhs(:)
      type(bc_tower) , intent(in   ) :: the_bc_tower

      integer        :: i,n,ng,ng_d
      type(bc_level) :: bc
      real(kind=dp_t), pointer :: divp(:,:,:,:) 

      ng_d = divu_rhs(1)%ng

      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, divu_rhs(n)%nboxes
            if ( multifab_remote(divu_rhs(n), i) ) cycle
            divp => dataptr(divu_rhs(n)     , i)
            select case (dm)
            case (2)
               call enforce_outflow_2d(divp(:,:,1,1), ng_d, bc%phys_bc_level_array(i,:,:))
            case (3)
               call enforce_outflow_3d(divp(:,:,:,1), ng_d, bc%phys_bc_level_array(i,:,:))
            end select
         end do
      end do

    end subroutine enforce_outflow_on_divu_rhs

    ! ******************************************************************************** !

    subroutine enforce_outflow_2d(divu_rhs,ng_d,phys_bc)

      integer        , intent(in   ) :: ng_d
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_d:,-ng_d:)
      integer        , intent(in   ) :: phys_bc(:,:)

      integer :: nx,ny
      nx = size(divu_rhs,dim=1)-1
      ny = size(divu_rhs,dim=2)-1

      if (phys_bc(1,1) .eq. OUTLET) divu_rhs(0,  :) = ZERO
      if (phys_bc(1,2) .eq. OUTLET) divu_rhs(nx, :) = ZERO
      if (phys_bc(2,1) .eq. OUTLET) divu_rhs(: , 0) = ZERO
      if (phys_bc(2,2) .eq. OUTLET) divu_rhs(: ,ny) = ZERO

    end subroutine enforce_outflow_2d

    ! ******************************************************************************** !

    subroutine enforce_outflow_3d(divu_rhs,ng_d,phys_bc)

      integer        , intent(in   ) :: ng_d
      real(kind=dp_t), intent(inout) :: divu_rhs(-ng_d:,-ng_d:,-ng_d:)
      integer        , intent(in   ) :: phys_bc(:,:)

      integer :: nx,ny,nz
      nx = size(divu_rhs,dim=1)-1
      ny = size(divu_rhs,dim=2)-1
      nz = size(divu_rhs,dim=3)-1

      if (phys_bc(1,1) .eq. OUTLET) divu_rhs(0,  :, :) = ZERO
      if (phys_bc(1,2) .eq. OUTLET) divu_rhs(nx, :, :) = ZERO
      if (phys_bc(2,1) .eq. OUTLET) divu_rhs( :, 0, :) = ZERO
      if (phys_bc(2,2) .eq. OUTLET) divu_rhs( :,ny, :) = ZERO
      if (phys_bc(3,1) .eq. OUTLET) divu_rhs( :,: , 0) = ZERO
      if (phys_bc(3,2) .eq. OUTLET) divu_rhs( :,: ,nz) = ZERO

    end subroutine enforce_outflow_3d

    ! ******************************************************************************** !

  end subroutine hgproject

  ! ******************************************************************************** !

  subroutine hg_multigrid(mla,unew,rhohalf,phi,dx,the_bc_tower, &
                          press_comp,stencil_type,divu_rhs,eps_in)

    use bl_prof_module
    use bl_constants_module
    use stencil_module
    use coeffs_module
    use ml_solve_module
    use nodal_divu_module
    use probin_module, only : hg_bottom_solver, max_mg_bottom_nlevels, verbose, mg_verbose, cg_verbose, nodal
    use geometry, only: dm, nlevs

    type(ml_layout), intent(inout) :: mla
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: rhohalf(:)
    type(multifab ), intent(inout) :: phi(:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: press_comp
    integer        , intent(in   ) :: stencil_type

    type(multifab ), intent(in   ), optional :: divu_rhs(:)
    real(dp_t)     , intent(in)   , optional :: eps_in 

    ! Local variables
    type(box     )  :: pd
    type(  layout)  :: la

    type(mg_tower) :: mgt(mla%nlevel)

    type(multifab) :: rh(mla%nlevel)
    type(multifab) :: one_sided_ss(2:mla%nlevel)

    type(multifab), allocatable :: coeffs(:)

    ! Bottom MGT stuff
    type(mg_tower)  :: bottom_mgt
    type(layout)    :: old_coarse_la,new_coarse_la
    type(layout)    :: old_la_grown, new_la_grown
    type(box)       :: coarse_pd,bxs
    type(boxarray)  :: ba_cc,new_coarse_ba
    type(multifab)  :: stored_coeffs, stored_coeffs_grown
    type(multifab)  :: new_coeffs_grown
    type(multifab), allocatable :: coarse_coeffs(:)
    integer         :: j,nx,mglev,bottom_box_size
    real(dp_t), pointer :: sc_orig(:,:,:,:), sc_grown(:,:,:,:)
    real(dp_t)      :: coarse_dx(dm)

    real(dp_t) :: bottom_solver_eps
    real(dp_t) :: eps
    real(dp_t) :: omega

    integer :: i, ns
    integer :: bottom_solver, bottom_max_iter
    integer :: max_iter
    integer :: min_width
    integer :: max_nlevel
    integer :: nu1, nu2, gamma, cycle_type, smoother
    integer :: n
    integer :: max_nlevel_in
    integer :: do_diagnostics

    type(bl_prof_timer), save :: bpt

    call build(bpt, "hg_multigrid")

    !! Defaults:
    max_nlevel        = mgt(nlevs)%max_nlevel
    max_iter          = mgt(nlevs)%max_iter
    eps               = mgt(nlevs)%eps
    smoother          = mgt(nlevs)%smoother
    nu1               = mgt(nlevs)%nu1
    nu2               = mgt(nlevs)%nu2
    gamma             = mgt(nlevs)%gamma
    omega             = mgt(nlevs)%omega
    cycle_type        = mgt(nlevs)%cycle_type
    bottom_solver     = mgt(nlevs)%bottom_solver
    bottom_solver_eps = mgt(nlevs)%bottom_solver_eps
    bottom_max_iter   = mgt(nlevs)%bottom_max_iter
    min_width         = mgt(nlevs)%min_width

    ! Note: put this here to minimize asymmetries - ASA
    eps = 1.d-12

    if ( hg_bottom_solver >= 0 ) then
        if (hg_bottom_solver == 4 .and. phi(1)%nboxes == 1) then
           if (parallel_IOProcessor()) then
              print *,'Dont use hg_bottom_solver == 4 with only one grid -- '
              print *,'  Reverting to default bottom solver ',bottom_solver
           end if
        else if (hg_bottom_solver == 4 .and. max_mg_bottom_nlevels < 2) then
           if (parallel_IOProcessor()) then
              print *,'Dont use hg_bottom_solver == 4 with max_mg_bottom_nlevels < 2'
              print *,'  Reverting to default bottom solver ',bottom_solver
           end if
        else
           bottom_solver = hg_bottom_solver
        end if
    end if

    min_width = 2

    ! Note: put this here for robustness
    max_iter = 100

    if (stencil_type .eq. ST_DENSE) then
       if (dm .eq. 3) then
          i = mgt(nlevs)%nlevels
          if ( (dx(nlevs,1) .eq. dx(nlevs,2)) .and. &
               (dx(nlevs,1) .eq. dx(nlevs,3)) ) then
             ns = 21
          else
             ns = 27
          end if
       else if (dm .eq. 2) then
          ns = 9
       end if
    else
       ns = 2*dm+1
       do n = nlevs, 2, -1
          la = mla%la(n)
          call multifab_build(one_sided_ss(n), la, ns, 0, nodal, stencil=.true.)
          call setval(one_sided_ss(n), ZERO,all=.true.)
       end do
    end if

    do n = nlevs, 1, -1

       if (n == 1) then
          max_nlevel_in = max_nlevel
       else
          if ( all(mla%mba%rr(n-1,:) == 2) ) then
             max_nlevel_in = 1
          else if ( all(mla%mba%rr(n-1,:) == 4) ) then
             max_nlevel_in = 2
          else 
             call bl_error("HG_MULTIGRID: confused about ref_ratio")
          end if
       end if

       pd = layout_get_pd(mla%la(n))

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                       the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,press_comp), &
                           dh = dx(n,:), &
                           ns = ns, &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           nub = nu2, &
                           gamma = gamma, &
                           cycle_type = cycle_type, &
                           omega = omega, &
                           bottom_solver = bottom_solver, &
                           bottom_max_iter = bottom_max_iter, &
                           bottom_solver_eps = bottom_solver_eps, &
                           max_iter = max_iter, &
                           max_nlevel = max_nlevel_in, &
                           min_width = min_width, &
                           eps = eps, &
                           verbose = mg_verbose, &
                           cg_verbose = cg_verbose, &
                           nodal = nodal)
       
    end do

    !! Fill coefficient array
    do n = nlevs,1,-1

       allocate(coeffs(mgt(n)%nlevels))

       la = mla%la(n)
       pd = layout_get_pd(la)

       call multifab_build(coeffs(mgt(n)%nlevels), la, 1, 1)
       call setval(coeffs(mgt(n)%nlevels), 0.0_dp_t, 1, all=.true.)

       call mkcoeffs(rhohalf(n),coeffs(mgt(n)%nlevels))
       call multifab_fill_boundary(coeffs(mgt(n)%nlevels))

       do i = mgt(n)%nlevels-1, 1, -1
          call multifab_build(coeffs(i), mgt(n)%ss(i)%la, 1, 1)
          call setval(coeffs(i), 0.0_dp_t, 1, all=.true.)
          call coarsen_coeffs(coeffs(i+1),coeffs(i))
          call multifab_fill_boundary(coeffs(i))
       end do

       do i = mgt(n)%nlevels, 1, -1
          call stencil_fill_nodal(mgt(n)%ss(i), coeffs(i), mgt(n)%dh(:,i), &
                                  mgt(n)%mm(i), mgt(n)%face_type, stencil_type)
          pd  = coarsen(pd,2)
       end do

       if (stencil_type .eq. ST_CROSS .and. n .gt. 1) then
          i = mgt(n)%nlevels
          call stencil_fill_one_sided(one_sided_ss(n), coeffs(i), &
                                      mgt(n    )%dh(:,i), &
                                      mgt(n)%mm(i), mgt(n)%face_type)
       end if

       ! Need to hang on to these coeffs at the bottom level of this mg
       if ( n == 1 .and. bottom_solver == 4 ) then
          call multifab_build(stored_coeffs, mgt(n)%ss(1)%la, 1, 1)
          call multifab_copy_c(stored_coeffs,1,coeffs(n),1,1,ng=coeffs(1)%ng)
       end if

       do i = mgt(n)%nlevels, 1, -1
          call destroy(coeffs(i))
       end do
       deallocate(coeffs)

    end do

    ! START OF BOTTOM_SOLVER == 4
    if (bottom_solver == 4) then

       ! Get the old/new coarse problem domain
       old_coarse_la = mgt(1)%ss(1)%la
       coarse_pd = layout_get_pd(old_coarse_la)

       ! Get the new coarse boxarray and layout
       call box_build_2(bxs,coarse_pd%lo(1:dm),coarse_pd%hi(1:dm))
       call boxarray_build_bx(new_coarse_ba,bxs)

       ! This is how many levels could be built if we made just one grid
       n = max_mg_levels(new_coarse_ba,min_width)

       ! This is the user-imposed limit
       n = min(n,max_mg_bottom_nlevels)

       if ( n .eq. 1) then
          call bl_error("DONT USE HG_BOTTOM_SOLVER == 4 WHEN BOTTOM GRID NOT PROPERLY DIVISIBLE : n = 1 ")
       end if

       bottom_box_size = 2**n
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          print *,'N ',n
          print *,'BOTTOM_BOX SIZE ',bottom_box_size
       end if

       do j = 1,dm
          nx = extent(bxs,j)
          if ( (bottom_box_size * (nx/bottom_box_size)) .ne. nx ) then
             call bl_error("DONT USE MG_BOTTOM_SOLVER == 4 WHEN BOTTOM GRID NOT PROPERLY DIVISIBLE ")
          end if
       end do

       call boxarray_maxsize(new_coarse_ba,bottom_box_size)
       call layout_build_ba(new_coarse_la,new_coarse_ba,coarse_pd,pmask=old_coarse_la%lap%pmask)

       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          call print(layout_get_pd(old_coarse_la),"COARSE PD")
          print *,'ORIG HG NBOXES ',old_coarse_la%lap%nboxes
          print *,'NEW  HG NBOXES ',new_coarse_la%lap%nboxes
       end if

       coarse_dx(:) = dx(1,:) * 2**(mgt(1)%nlevels-1)

       call mg_tower_build(bottom_mgt, new_coarse_la, coarse_pd, &
                           the_bc_tower%bc_tower_array(1)%ell_bc_level_array(0,:,:,press_comp),&
                           dh = coarse_dx, &
                           ns = ns, &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           gamma = gamma, &
                           cycle_type = cycle_type, &
                           omega = omega, &
                           bottom_solver = 1, &
                           bottom_max_iter = bottom_max_iter, &
                           bottom_solver_eps = bottom_solver_eps, &
                           max_iter = max_iter, &
                           max_nlevel = max_nlevel, &
                           min_width = min_width, &
                           eps = eps, &
                           verbose = mg_verbose, &
                           cg_verbose = cg_verbose, &
                           nodal = nodal)

       ! START SPECIAL COPY
       ! Here we do special stuff to be able to copy the ghost cells of stored_coeffs into
       !   the ghost cells of coarse_coeffs(bottom)

       ! Make sure to do this before the copy so we get all the data
       call multifab_fill_boundary(stored_coeffs)

       mglev = bottom_mgt%nlevels

       allocate(coarse_coeffs(mglev))
       call multifab_build(coarse_coeffs(mglev),new_coarse_la,1,1)
       call setval(coarse_coeffs(mglev),ZERO,all=.true.)

       call boxarray_build_copy(ba_cc,get_boxarray(stored_coeffs))
       call boxarray_grow(ba_cc,1)
       call layout_build_ba(old_la_grown,ba_cc,pmask=old_coarse_la%lap%pmask, &
                            explicit_mapping=get_proc(old_coarse_la))
       call destroy(ba_cc)
       call multifab_build(stored_coeffs_grown,old_la_grown,1,ng=0)

       do i = 1, stored_coeffs_grown%nboxes
          if (remote(stored_coeffs_grown,i)) cycle 
          sc_orig  => dataptr(stored_coeffs      ,i,get_pbox(stored_coeffs_grown,i),1,1)
          sc_grown => dataptr(stored_coeffs_grown,i,get_pbox(stored_coeffs_grown,i),1,1)
          sc_grown = sc_orig
       end do

       call boxarray_build_copy(ba_cc,new_coarse_ba)
       call boxarray_grow(ba_cc,1)
       call layout_build_ba(new_la_grown,ba_cc,pmask=old_coarse_la%lap%pmask, &
                            explicit_mapping=get_proc(new_coarse_la))
       call destroy(ba_cc)
       call multifab_build(new_coeffs_grown,new_la_grown,1,ng=0)
       call multifab_copy_c(new_coeffs_grown,1,stored_coeffs_grown,1,1)

       do i = 1, new_coeffs_grown%nboxes
          if (remote(new_coeffs_grown,i)) cycle 
          sc_orig  => dataptr(coarse_coeffs(mglev),i,get_pbox(new_coeffs_grown,i),1,1)
          sc_grown => dataptr(new_coeffs_grown    ,i,get_pbox(new_coeffs_grown,i),1,1)
          sc_orig = sc_grown
       end do

       call destroy(new_coeffs_grown)
       call destroy(new_coarse_ba)
       !   END SPECIAL COPY

       do i = mglev-1, 1, -1
          call multifab_build(coarse_coeffs(i), bottom_mgt%ss(i)%la, 1, 1)
          call setval(coarse_coeffs(i), ZERO, 1, 1, all=.true.)
          call coarsen_coeffs(coarse_coeffs(i+1),coarse_coeffs(i))
          call multifab_fill_boundary(coarse_coeffs(i))
       end do

       do i = mglev, 1, -1
          call stencil_fill_nodal(bottom_mgt%ss(i), coarse_coeffs(i), bottom_mgt%dh(:,i), &
                                  bottom_mgt%mm(i), bottom_mgt%face_type, stencil_type)
       end do

       do i = mglev, 1, -1
          call destroy(coarse_coeffs(i))
       end do
       deallocate(coarse_coeffs)
       call destroy(stored_coeffs)
       call destroy(stored_coeffs_grown)
       call destroy(old_la_grown)
       call destroy(new_la_grown)
    end if
    ! END   OF BOTTOM_SOLVER == 4

    do n = 1, nlevs
       call multifab_build(rh(n),mla%la(n),1,1,nodal)
       call setval(rh(n),ZERO,all=.true.)
    end do

    call divu(nlevs,mgt,unew,rh,mla%mba%rr,verbose,nodal)

    ! Do rh = rh - divu_rhs (this routine preserves rh=0 on
    !  nodes which have bc_dirichlet = true.
    if (present(divu_rhs)) &
       call subtract_divu_from_rh(nlevs,mgt,rh,divu_rhs)

    if ( mg_verbose >= 3 ) then
       do_diagnostics = 1
    else
       do_diagnostics = 0
    end if

    if (present(eps_in)) then
       if (bottom_solver == 4) then
          call ml_nd_solve(mla,mgt,rh,phi,one_sided_ss,mla%mba%rr,do_diagnostics,&
                           eps_in=eps_in,bottom_mgt=bottom_mgt)
       else
          call ml_nd_solve(mla,mgt,rh,phi,one_sided_ss,mla%mba%rr,do_diagnostics,eps_in=eps_in)
       end if
    else
       if (bottom_solver == 4) then
          call ml_nd_solve(mla,mgt,rh,phi,one_sided_ss,mla%mba%rr,do_diagnostics,bottom_mgt=bottom_mgt)
       else
          call ml_nd_solve(mla,mgt,rh,phi,one_sided_ss,mla%mba%rr,do_diagnostics)
       end if
    end if

    do n = nlevs,1,-1
       call multifab_fill_boundary(phi(n))
    end do

    do n = 1, nlevs
       call mg_tower_destroy(mgt(n))
       call destroy(rh(n))
    end do

    if (bottom_solver == 4) then
       call destroy(new_coarse_la)
       call mg_tower_destroy(bottom_mgt)
    end if

    if (stencil_type .ne. ST_DENSE) then
       do n = nlevs, 2, -1
          call destroy(one_sided_ss(n))
       end do
    end if

    call destroy(bpt)

  end subroutine hg_multigrid

  !   ********************************************************************************* !

  subroutine mkcoeffs(rho,coeffs)

    use geometry, only: dm

    type(multifab) , intent(in   ) :: rho
    type(multifab) , intent(inout) :: coeffs

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    integer :: i,ng_r,ng_c

    ng_r = rho%ng
    ng_c = coeffs%ng

    do i = 1, rho%nboxes
       if ( multifab_remote(rho, i) ) cycle
       rp => dataptr(rho   , i)
       cp => dataptr(coeffs, i)
       select case (dm)
       case (2)
          call mkcoeffs_2d(cp(:,:,1,1), ng_c, rp(:,:,1,1), ng_r)
       case (3)
          call mkcoeffs_3d(cp(:,:,:,1), ng_c, rp(:,:,:,1), ng_r)
       end select
    end do

  end subroutine mkcoeffs

  !   *********************************************************************************** !

  subroutine mkcoeffs_2d(coeffs,ng_c,rho,ng_r)

      use bl_constants_module

    integer                        :: ng_c,ng_r
    real(kind=dp_t), intent(inout) :: coeffs(1-ng_c:,1-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(1-ng_r:,1-ng_r:)

    integer :: i,j
    integer :: nx,ny

    nx = size(coeffs,dim=1) - 2
    ny = size(coeffs,dim=2) - 2

    do j = 1,ny
       do i = 1,nx
          coeffs(i,j) = ONE / rho(i,j)
       end do
    end do

  end subroutine mkcoeffs_2d

  !   ********************************************************************************** !

  subroutine mkcoeffs_3d(coeffs,ng_c,rho,ng_r)

      use bl_constants_module

    integer                        :: ng_c,ng_r
    real(kind=dp_t), intent(inout) :: coeffs(1-ng_c:,1-ng_c:,1-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(1-ng_r:,1-ng_r:,1-ng_r:)

    integer :: i,j,k
    integer :: nx,ny,nz

    nx = size(coeffs,dim=1) - 2
    ny = size(coeffs,dim=2) - 2
    nz = size(coeffs,dim=3) - 2

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             coeffs(i,j,k) = ONE / rho(i,j,k)
          end do
       end do
    end do

  end subroutine mkcoeffs_3d

  !   ********************************************************************************** !

  subroutine mult_by_1d_coeff(u,div_coeff,do_mult)

    use geometry, only: dm, nlevs

    type(multifab), intent(inout)           :: u(:)
    real(dp_t)    , intent(in   )           :: div_coeff(:,:)
    logical       , intent(in   ), optional :: do_mult

    ! local
    real(kind=dp_t), pointer :: ump(:,:,:,:) 
    integer :: i,ng_u,n
    integer :: lo(dm),hi(dm)
    logical :: local_do_mult

    local_do_mult = .true.
    if (present(do_mult)) local_do_mult = do_mult

    ng_u = u(1)%ng

    do n = 1, nlevs

       ! Multiply u by div coeff
       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          ump => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call mult_by_1d_coeff_2d(ump(:,:,1,:),ng_u,div_coeff(n,:),lo,hi,local_do_mult)
          case (3)
             call mult_by_1d_coeff_3d(ump(:,:,:,:),ng_u,div_coeff(n,:),lo,hi,local_do_mult)
          end select
       end do

    end do

  end subroutine mult_by_1d_coeff

  !   ********************************************************************************** !

  subroutine mult_by_1d_coeff_2d(u,ng_u,div_coeff,lo,hi,do_mult)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u
    real(kind=dp_t), intent(inout) :: u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real(dp_t)     , intent(in   ) :: div_coeff(0:)
    logical        , intent(in   ) :: do_mult

    integer :: j

    if (do_mult) then
       do j = lo(2),hi(2)
          u(:,j,:) = u(:,j,:) * div_coeff(j)
       end do
    else
       do j = lo(2),hi(2)
          u(:,j,:) = u(:,j,:) / div_coeff(j)
       end do
    end if

  end subroutine mult_by_1d_coeff_2d

  !   ********************************************************************************** !

  subroutine mult_by_1d_coeff_3d(u,ng_u,div_coeff,lo,hi,do_mult)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u
    real(kind=dp_t), intent(inout) :: u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real(dp_t)     , intent(in   ) :: div_coeff(0:)
    logical        , intent(in   ) :: do_mult

    integer :: k

    if (do_mult) then
       do k = lo(3),hi(3)
          u(:,:,k,:) = u(:,:,k,:) * div_coeff(k)
       end do
    else
       do k = lo(3),hi(3)
          u(:,:,k,:) = u(:,:,k,:) / div_coeff(k)
       end do
    end if

  end subroutine mult_by_1d_coeff_3d

  !   *********************************************************************************** !

  subroutine mult_by_3d_coeff(u,div_coeff,do_mult)

    use geometry, only: dm, nlevs

    type(multifab) , intent(inout) :: u(:)
    type(multifab) , intent(in   ) :: div_coeff(:)
    logical        , intent(in   ) :: do_mult

    ! local
    real(kind=dp_t), pointer :: ump(:,:,:,:) 
    real(kind=dp_t), pointer ::  dp(:,:,:,:) 
    integer :: i,ng_u,ng_d,n

    ng_u = u(1)%ng
    ng_d = div_coeff(1)%ng

    do n = 1, nlevs
       ! Multiply u by div coeff
       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          ump => dataptr(u(n),i)
          dp => dataptr(div_coeff(n),i)
          select case (dm)
          case (3)
             call mult_by_3d_coeff_3d(ump(:,:,:,:), ng_u, dp(:,:,:,1), ng_d, do_mult)
          end select
       end do

    end do

  end subroutine mult_by_3d_coeff

  !   ********************************************************************************** !
  
  subroutine mult_by_3d_coeff_3d(u,ng_u,div_coeff,ng_d,do_mult)

    integer        , intent(in   ) :: ng_u,ng_d
    real(kind=dp_t), intent(inout) ::         u(-ng_u:,-ng_u:,-ng_u:,:)
    real(dp_t)     , intent(in   ) :: div_coeff(-ng_d:,-ng_d:,-ng_d:)
    logical        , intent(in   ) :: do_mult

    integer :: i,j,k,m,nx,ny,nz
    nx = size(u,dim=1) - 2*ng_u
    ny = size(u,dim=2) - 2*ng_u
    nz = size(u,dim=3) - 2*ng_u

    if (do_mult) then
       do m = 1, size(u,dim=4)
          do k = 0,nz-1 
             do j = 0,ny-1 
                do i = 0,nx-1 
                   u(i,j,k,m) = u(i,j,k,m) * div_coeff(i,j,k)
                end do
             end do
          end do
       end do
    else
       do m = 1, size(u,dim=4)
          do k = 0,nz-1 
             do j = 0,ny-1 
                do i = 0,nx-1 
                   u(i,j,k,m) = u(i,j,k,m) / div_coeff(i,j,k)
                end do
             end do
          end do
       end do
    end if

  end subroutine mult_by_3d_coeff_3d

  !   ********************************************************************************** !

end module hgproject_module
