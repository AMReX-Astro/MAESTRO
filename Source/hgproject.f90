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

  subroutine hgproject(proj_type,mla,unew,uold,rhohalf,p,gp,dx,dt,the_bc_tower, &
                       verbose,mg_verbose,cg_verbose,press_comp, &
                       divu_rhs,div_coeff_1d,div_coeff_3d,eps_in)

    use bc_module
    use proj_parameters
    use nodal_divu_module
    use stencil_module
    use ml_solve_module
    use ml_restriction_module
    use multifab_fill_ghost_module

    integer        , intent(in   ) :: proj_type
    type(ml_layout), intent(inout) :: mla
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: uold(:)
    type(multifab ), intent(inout) :: rhohalf(:)
    type(multifab ), intent(inout) :: p(:)
    type(multifab ), intent(inout) :: gp(:)
    real(dp_t)     , intent(in   ) :: dx(:,:),dt
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: verbose,mg_verbose,cg_verbose
    integer        , intent(in   ) :: press_comp

    type(multifab ), intent(inout), optional :: divu_rhs(:)
    real(dp_t)     , intent(in   ), optional :: div_coeff_1d(:,:)
    type(multifab ), intent(in   ), optional :: div_coeff_3d(:)
    real(dp_t)     , intent(in   ), optional :: eps_in

    ! Local  
    type(multifab), allocatable :: phi(:),gphi(:)
    logical,        allocatable :: nodal(:)
    integer                     :: n,nlevs,dm,ng
    real(dp_t)                  :: umin,umax,vmin,vmax,wmin,wmax
    integer                     :: stencil_type
    logical                     :: use_div_coeff_1d, use_div_coeff_3d

    ! stencil_type = ST_DENSE
    stencil_type = ST_CROSS

    nlevs = mla%nlevel
    dm = mla%dim
    ng = unew(nlevs)%ng

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       print *,'PROJ_TYPE IN HGPROJECT:',proj_type
    endif

    allocate(phi(nlevs), gphi(nlevs), nodal(dm))
    nodal = .true.

    use_div_coeff_1d = .false.
    if (present(div_coeff_1d)) use_div_coeff_1d = .true.

    use_div_coeff_3d = .false.
    if (present(div_coeff_3d)) use_div_coeff_3d = .true.

    if (use_div_coeff_1d .and. use_div_coeff_3d) then
       print *,'CANT HAVE 1D and 3D DIV_COEFF IN HGPROJECT '
       stop
    end if

    do n = 1, nlevs
       call multifab_build( phi(n), mla%la(n), 1, 1, nodal)
       call multifab_build(gphi(n), mla%la(n), dm, 0) 
       call multifab_copy(phi(n),p(n))
       call multifab_mult_mult_s(phi(n),dt,phi(n)%ng)
    end do

    if (verbose .ge. 1) then
       umin = 1.d30
       vmin = 1.d30
       wmin = 1.d30
       umax = -1.d30
       vmax = -1.d30
       wmax = -1.d30
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

    call create_uvec_for_projection(nlevs,unew,uold,rhohalf,gp,dt,the_bc_tower,proj_type)

    if (use_div_coeff_1d) then
       call mult_by_1d_coeff(nlevs,unew,div_coeff_1d,.true.)
       call mult_by_1d_coeff(nlevs,rhohalf,div_coeff_1d,.false.)
    else if (use_div_coeff_3d) then
       call mult_by_3d_coeff(nlevs,unew,div_coeff_3d,.true.)
       call mult_by_3d_coeff(nlevs,rhohalf,div_coeff_3d,.false.)
    end if

    do n = 1, nlevs
       call setval(phi(n),ZERO,all=.true.)
    end do

    if (present(divu_rhs)) then
       call enforce_outflow_on_divu_rhs(divu_rhs,the_bc_tower)
    endif

    if (present(eps_in)) then
       call hg_multigrid(mla,unew,rhohalf,phi,dx,the_bc_tower,verbose, &
                         mg_verbose,cg_verbose,press_comp,stencil_type,divu_rhs,eps_in)
    else
       call hg_multigrid(mla,unew,rhohalf,phi,dx,the_bc_tower,verbose, &
                         mg_verbose,cg_verbose,press_comp,stencil_type,divu_rhs)
    end if

    if (use_div_coeff_1d) then
       call mult_by_1d_coeff(nlevs,unew,div_coeff_1d,.false.)
       call mult_by_1d_coeff(nlevs,rhohalf,div_coeff_1d,.true.)
    else if (use_div_coeff_3d) then
       call mult_by_3d_coeff(nlevs,unew,div_coeff_3d,.false.)
       call mult_by_3d_coeff(nlevs,rhohalf,div_coeff_3d,.true.)
    end if

    call mkgphi(nlevs,gphi,phi,dx)

    call hg_update(nlevs,proj_type,unew,uold,gp,gphi,rhohalf,  &
                   p,phi,ng,dt,mla,the_bc_tower%bc_tower_array)

    if (verbose .ge. 1) then
       umin = 1.d30
       vmin = 1.d30
       wmin = 1.d30
       umax = -1.d30
       vmax = -1.d30
       wmax = -1.d30
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
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
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

    do n = 1,nlevs
       call multifab_destroy(phi(n))
       call multifab_destroy(gphi(n))
    end do

    deallocate(phi)
    deallocate(gphi)
    deallocate(nodal)

  contains

    subroutine create_uvec_for_projection(nlevs,unew,uold,rhohalf,gp,dt,the_bc_tower,proj_type)

      integer        , intent(in   ) :: nlevs
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(in   ) :: uold(:)
      type(multifab) , intent(in   ) :: rhohalf(:)
      type(multifab) , intent(inout) :: gp(:)
      real(kind=dp_t), intent(in   ) :: dt
      type(bc_tower) , intent(in   ) :: the_bc_tower
      integer        , intent(in   ) :: proj_type
  
      type(bc_level) :: bc
  
      real(kind=dp_t), pointer :: unp(:,:,:,:)
      real(kind=dp_t), pointer :: uop(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:)
  
      integer :: i,n,dm,ng
  
      dm = unew(nlevs)%dim
      ng = unew(nlevs)%ng
  
      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, unew(n)%nboxes
            if ( multifab_remote(unew(n), i) ) cycle
            unp => dataptr(unew(n)     , i) 
            uop => dataptr(uold(n)     , i) 
            gpp => dataptr(gp(n)       , i)
             rp => dataptr(  rhohalf(n), i)
            select case (dm)
               case (2)
                 call create_uvec_2d(unp(:,:,1,:), uop(:,:,1,:), rp(:,:,1,1), gpp(:,:,1,:), dt, &
                                     bc%phys_bc_level_array(i,:,:), ng, proj_type)
               case (3)
                 call create_uvec_3d(unp(:,:,:,:), uop(:,:,:,:), rp(:,:,:,1), gpp(:,:,:,:), dt, &
                                     bc%phys_bc_level_array(i,:,:), ng, proj_type)
            end select
         end do
         call multifab_fill_boundary(unew(n))
      end do 

    end subroutine create_uvec_for_projection


    !   ******************************************************************************* !

    subroutine mkgphi(nlevs,gphi,phi,dx)

      integer       , intent(in   ) :: nlevs
      type(multifab), intent(inout) :: gphi(:)
      type(multifab), intent(in   ) :: phi(:)
      real(dp_t) :: dx(:,:)

      integer :: i,dm,n

      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: pp(:,:,:,:) 

      dm = phi(1)%dim

      do n = 1, nlevs

         do i = 1, phi(n)%nboxes
            if ( multifab_remote(phi(n),i) ) cycle
            gph => dataptr(gphi(n),i)
            pp  => dataptr(phi(n),i)
            select case (dm)
            case (2)
               call mkgphi_2d(gph(:,:,1,:), pp(:,:,1,1), dx(n,:))
            case (3)
               call mkgphi_3d(gph(:,:,:,:), pp(:,:,:,1), dx(n,:))
            end select
         end do

         call multifab_fill_boundary(gphi(n))

      end do

    end subroutine mkgphi

    !   ********************************************************************************** !

    subroutine hg_update(nlevs,proj_type,unew,uold,gp,gphi,rhohalf,p,phi,ng,dt, &
                         mla,the_bc_level)

      integer        , intent(in   ) :: nlevs
      integer        , intent(in   ) :: proj_type
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(in   ) :: uold(:)
      type(multifab) , intent(inout) :: gp(:)
      type(multifab) , intent(in   ) :: gphi(:)
      type(multifab) , intent(in   ) :: rhohalf(:)
      type(multifab) , intent(inout) :: p(:)
      type(multifab) , intent(in   ) :: phi(:)
      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(in   ) :: dt
      type(ml_layout), intent(inout) :: mla
      type(bc_level) , intent(in   ) :: the_bc_level(:)

      ! local
      integer :: i,dm,n

      real(kind=dp_t), pointer :: upn(:,:,:,:) 
      real(kind=dp_t), pointer :: uon(:,:,:,:) 
      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:) 
      real(kind=dp_t), pointer ::  ph(:,:,:,:) 
      real(kind=dp_t), pointer ::  pp(:,:,:,:) 

      dm = unew(1)%dim

      do n = 1, nlevs

         do i = 1, unew(n)%nboxes
            if ( multifab_remote(unew(n),i) ) cycle
            upn => dataptr(unew(n),i)
            uon => dataptr(uold(n),i)
            gpp => dataptr(gp(n),i)
            gph => dataptr(gphi(n),i)
            rp  => dataptr(rhohalf(n),i)
            pp  => dataptr(p(n),i)
            ph  => dataptr(phi(n),i)
            select case (dm)
            case (2)
               call hg_update_2d(proj_type, upn(:,:,1,:), uon(:,:,1,:), gpp(:,:,1,:), &
                                 gph(:,:,1,:),rp(:,:,1,1),pp(:,:,1,1), ph(:,:,1,1), ng, dt)
            case (3)
               call hg_update_3d(proj_type, upn(:,:,:,:), uon(:,:,:,:), gpp(:,:,:,:), &
                                 gph(:,:,:,:),rp(:,:,:,1),pp(:,:,:,1), ph(:,:,:,1), ng, dt)
            end select
         end do

         call multifab_fill_boundary(unew(n))
         call multifab_fill_boundary(gp(n))
         call multifab_fill_boundary(p(n))

      end do

      do n = nlevs, 2, -1
         call ml_cc_restriction(unew(n-1),unew(n),mla%mba%rr(n-1,:)) 
         call ml_cc_restriction(  gp(n-1),  gp(n),mla%mba%rr(n-1,:))

         call multifab_fill_ghost_cells(unew(n),unew(n-1),ng,mla%mba%rr(n-1,:), &
                                        the_bc_level(n-1),the_bc_level(n),1,1,dm)
      end do

    end subroutine hg_update

    !   ********************************************************************************** !

    subroutine enforce_outflow_on_divu_rhs(divu_rhs,the_bc_tower)

      type(multifab) , intent(inout) :: divu_rhs(:)
      type(bc_tower) , intent(in   ) :: the_bc_tower

      integer        :: i,n,dm,ng,nlevs
      type(bc_level) :: bc
      real(kind=dp_t), pointer :: divp(:,:,:,:) 

      nlevs = size(divu_rhs,dim=1)
      dm = divu_rhs(1)%dim

      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, divu_rhs(n)%nboxes
            if ( multifab_remote(divu_rhs(n), i) ) cycle
            divp => dataptr(divu_rhs(n)     , i)
            select case (dm)
            case (2)
               call enforce_outflow_2d(divp(:,:,1,1), bc%phys_bc_level_array(i,:,:))
            case (3)
               call enforce_outflow_3d(divp(:,:,:,1), bc%phys_bc_level_array(i,:,:))
            end select
         end do
      end do

    end subroutine enforce_outflow_on_divu_rhs

    !   ******************************************************************************** !

    subroutine enforce_outflow_2d(divu_rhs,phys_bc)

      real(kind=dp_t), intent(inout) :: divu_rhs(0:,0:)
      integer        , intent(in   ) :: phys_bc(:,:)

      integer :: nx,ny
      nx = size(divu_rhs,dim=1)-1
      ny = size(divu_rhs,dim=2)-1

      if (phys_bc(1,1) .eq. OUTLET) divu_rhs(0,  :) = ZERO
      if (phys_bc(1,2) .eq. OUTLET) divu_rhs(nx, :) = ZERO
      if (phys_bc(2,1) .eq. OUTLET) divu_rhs(: , 0) = ZERO
      if (phys_bc(2,2) .eq. OUTLET) divu_rhs(: ,ny) = ZERO

    end subroutine enforce_outflow_2d

    !   ******************************************************************************** !

    subroutine enforce_outflow_3d(divu_rhs,phys_bc)

      real(kind=dp_t), intent(inout) :: divu_rhs(0:,0:,0:)
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

    !   ******************************************************************************** !

    subroutine create_uvec_2d(unew,uold,rhohalf,gp,dt,phys_bc,ng,proj_type)

      use proj_parameters

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::    unew(-ng:,-ng:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng:,-ng:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf( -1:, -1:)
      real(kind=dp_t), intent(inout) ::      gp( -1:, -1:,:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      integer :: nx,ny
      nx = size(gp,dim=1) - 2
      ny = size(gp,dim=2) - 2

      if (phys_bc(1,1) .eq. INLET) gp(-1,-1:ny,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gp(nx,-1:ny,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gp(-1:nx,-1,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gp(-1:nx,ny,:) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters) then

         unew(-1:nx,-1:ny,1) = ( unew(-1:nx,-1:ny,1) - uold(-1:nx,-1:ny,1) ) / dt
         unew(-1:nx,-1:ny,2) = ( unew(-1:nx,-1:ny,2) - uold(-1:nx,-1:ny,2) ) / dt
     
      ! quantity projected is Ustar + dt * (1/rho) Gp
      else if (proj_type .eq. regular_timestep) then

         unew(-1:nx,-1:ny,1) = unew(-1:nx,-1:ny,1) + dt*gp(-1:nx,-1:ny,1)/rhohalf(-1:nx,-1:ny)
         unew(-1:nx,-1:ny,2) = unew(-1:nx,-1:ny,2) + dt*gp(-1:nx,-1:ny,2)/rhohalf(-1:nx,-1:ny)

       else
     
          call bl_error('No proj_type by this number ')

      end if

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL) unew(-1,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL) unew(nx,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL) unew(:,-1,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL) unew(:,ny,:) = ZERO

    end subroutine create_uvec_2d

!   *************************************************************************************** !

    subroutine create_uvec_3d(unew,uold,rhohalf,gp,dt,phys_bc,ng,proj_type)

      use proj_parameters

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::    unew(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::      gp( -1:, -1:, -1:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf( -1:, -1:, -1:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      integer :: nx,ny,nz

      nx = size(gp,dim=1) - 2
      ny = size(gp,dim=2) - 2
      nz = size(gp,dim=3) - 2

      if (phys_bc(1,1) .eq. INLET) gp(-1,-1:ny,-1:nz,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gp(nx,-1:ny,-1:nz,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gp(-1:nx,-1,-1:nz,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gp(-1:nx,ny,-1:nz,:) = ZERO
      if (phys_bc(3,1) .eq. INLET) gp(-1:nx,-1:ny,-1,:) = ZERO
      if (phys_bc(3,2) .eq. INLET) gp(-1:nx,-1:ny,nz,:) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters) then

         unew(-1:nx,-1:ny,-1:nz,:) = ( unew(-1:nx,-1:ny,-1:nz,:) - uold(-1:nx,-1:ny,-1:nz,:) ) / dt

      ! quantity projected is Ustar + dt * (1/rho) Gp
      else if (proj_type .eq. regular_timestep) then

         unew(-1:nx,-1:ny,-1:nz,1) = unew(-1:nx,-1:ny,-1:nz,1) + &
                                    dt*gp(-1:nx,-1:ny,-1:nz,1)/rhohalf(-1:nx,-1:ny,-1:nz)
         unew(-1:nx,-1:ny,-1:nz,2) = unew(-1:nx,-1:ny,-1:nz,2) + &
                                    dt*gp(-1:nx,-1:ny,-1:nz,2)/rhohalf(-1:nx,-1:ny,-1:nz)
         unew(-1:nx,-1:ny,-1:nz,3) = unew(-1:nx,-1:ny,-1:nz,3) + &
                                    dt*gp(-1:nx,-1:ny,-1:nz,3)/rhohalf(-1:nx,-1:ny,-1:nz)

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

    !   ********************************************************************************* !

    subroutine mkgphi_2d(gp,phi,dx)

      real(kind=dp_t), intent(inout) ::  gp(0:,0:,:)
      real(kind=dp_t), intent(inout) :: phi(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,nx,ny

      nx = size(gp,dim=1)
      ny = size(gp,dim=2)

      do j = 0,ny-1
         do i = 0,nx-1
            gp(i,j,1) = HALF*(phi(i+1,j) + phi(i+1,j+1) - &
                 phi(i  ,j) - phi(i  ,j+1) ) /dx(1)
            gp(i,j,2) = HALF*(phi(i,j+1) + phi(i+1,j+1) - &
                 phi(i,j  ) - phi(i+1,j  ) ) /dx(2)
         end do
      end do

    end subroutine mkgphi_2d

    !   ******************************************************************************** !

    subroutine mkgphi_3d(gp,phi,dx)

      real(kind=dp_t), intent(inout) ::  gp(0:,0:,0:,1:)
      real(kind=dp_t), intent(inout) :: phi(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k,nx,ny,nz

      nx = size(gp,dim=1)
      ny = size(gp,dim=2)
      nz = size(gp,dim=3)

      do k = 0,nz-1
         do j = 0,ny-1
            do i = 0,nx-1
               gp(i,j,k,1) = FOURTH*(phi(i+1,j,k  ) + phi(i+1,j+1,k  ) &
                    +phi(i+1,j,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i  ,j,k  ) - phi(i  ,j+1,k  ) &
                    -phi(i  ,j,k+1) - phi(i  ,j+1,k+1) ) /dx(1)
               gp(i,j,k,2) = FOURTH*(phi(i,j+1,k  ) + phi(i+1,j+1,k  ) &
                    +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                    -phi(i,j  ,k+1) - phi(i+1,j  ,k+1) ) /dx(2)
               gp(i,j,k,3) = FOURTH*(phi(i,j  ,k+1) + phi(i+1,j  ,k+1) &
                    +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                    -phi(i,j+1,k  ) - phi(i+1,j+1,k  ) ) /dx(3)
            end do
         end do
      end do

    end subroutine mkgphi_3d

    !   ****************************************************************************** !

    subroutine hg_update_2d(proj_type,unew,uold,gp,gphi,rhohalf,p,phi,ng,dt)

      use proj_parameters

      integer        , intent(in   ) :: proj_type
      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::    unew(-ng:,-ng:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::      gp( -1:, -1:,:)
      real(kind=dp_t), intent(in   ) ::    gphi(  0:,  0:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf( -1:, -1:)
      real(kind=dp_t), intent(inout) ::       p( -1:, -1:)
      real(kind=dp_t), intent(in   ) ::     phi( -1:, -1:)
      real(kind=dp_t), intent(in   ) :: dt

      integer         :: nx,ny

      nx = size(gphi,dim=1)-1
      ny = size(gphi,dim=2)-1

      !     Subtract off the density-weighted gradient.
      unew(0:nx,0:ny,1) = unew(0:nx,0:ny,1) - gphi(0:nx,0:ny,1)/rhohalf(0:nx,0:ny) 
      unew(0:nx,0:ny,2) = unew(0:nx,0:ny,2) - gphi(0:nx,0:ny,2)/rhohalf(0:nx,0:ny) 

      if (proj_type .eq. pressure_iters) &    ! unew held the projection of (ustar-uold)
           unew(0:nx,0:ny,:) = uold(0:nx,0:ny,:) + dt * unew(0:nx,0:ny,:)

      if ( (proj_type .eq. initial_projection) .or. (proj_type .eq. divu_iters) ) then

         gp = ZERO
         p = ZERO

      else if (proj_type .eq. pressure_iters) then

         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         gp(0:nx,0:ny,:) = gp(0:nx,0:ny,:) +            gphi(0:nx,0:ny,:)
         p(0:nx,0:ny  ) =  p(0:nx,0:ny  ) +             phi(0:nx,0:ny  )

      else if (proj_type .eq. regular_timestep) then

         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         gp(0:nx,0:ny,:) = (ONE/dt) * gphi(0:nx,0:ny,:)
         p(0:nx,0:ny  ) = (ONE/dt) *  phi(0:nx,0:ny  )

      end if

    end subroutine hg_update_2d

    !   ******************************************************************************* !

    subroutine hg_update_3d(proj_type,unew,uold,gp,gphi,rhohalf,p,phi,ng,dt)

      use proj_parameters

      integer        , intent(in   ) :: proj_type
      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::    unew(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::      gp(-1:,-1:,-1:,:)
      real(kind=dp_t), intent(in   ) ::    gphi( 0:, 0:, 0:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) ::       p( -1:, -1:,-1:)
      real(kind=dp_t), intent(in   ) ::     phi( -1:, -1:,-1:)
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

      if (proj_type .eq. pressure_iters) &    ! unew held the projection of (ustar-uold)
           unew(0:nx,0:ny,0:nz,:) = uold(0:nx,0:ny,0:nz,:) + dt * unew(0:nx,0:ny,0:nz,:)

      if ( (proj_type .eq. initial_projection) .or. (proj_type .eq. divu_iters) ) then

         gp = ZERO
         p = ZERO

      else if (proj_type .eq. pressure_iters) then

         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         gp(0:nx,0:ny,0:nz,:) = gp(0:nx,0:ny,0:nz,:) +            gphi(0:nx,0:ny,0:nz,:)
         p(0:nx,0:ny,0:nz  ) =  p(0:nx,0:ny,0:nz  ) +             phi(0:nx,0:ny,0:nz  )

      else if (proj_type .eq. regular_timestep) then

         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         gp(0:nx,0:ny,0:nz,:) = (ONE/dt) * gphi(0:nx,0:ny,0:nz,:)
         p(0:nx,0:ny,0:nz  ) = (ONE/dt) *  phi(0:nx,0:ny,0:nz  )

      end if

    end subroutine hg_update_3d

  end subroutine hgproject

  subroutine hg_multigrid(mla,unew,rhohalf,phi,dx,the_bc_tower,divu_verbose,mg_verbose, &
                          cg_verbose,press_comp,stencil_type,divu_rhs,eps_in)

    use bl_constants_module
    use stencil_module
    use coeffs_module
    use ml_solve_module
    use nodal_divu_module

    type(ml_layout), intent(inout) :: mla
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: rhohalf(:)
    type(multifab ), intent(inout) :: phi(:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: divu_verbose,mg_verbose,cg_verbose
    integer        , intent(in   ) :: press_comp
    integer        , intent(in   ) :: stencil_type

    type(multifab ), intent(in   ), optional :: divu_rhs(:)
    real(dp_t)     , intent(in)   , optional :: eps_in 
    !
    ! Local variables
    !

    type(box     )  :: pd
    type(  layout)  :: la

    type(mg_tower), allocatable :: mgt(:)
    type(multifab), allocatable :: coeffs(:),rh(:)
    type(multifab), allocatable :: one_sided_ss(:)

    real(dp_t) :: bottom_solver_eps
    real(dp_t) :: eps
    real(dp_t) :: omega

    integer :: i, dm, nlevs, ns
    integer :: bottom_solver, bottom_max_iter
    integer :: max_iter
    integer :: min_width
    integer :: max_nlevel
    integer :: nu1, nu2, gamma, cycle, smoother
    integer :: n
    integer :: max_nlevel_in
    integer :: verbose
    integer :: do_diagnostics
    logical, allocatable :: nodal(:)

    !! Defaults:

    dm    = mla%dim
    nlevs = mla%nlevel

    allocate(mgt(nlevs), nodal(dm))
    allocate(one_sided_ss(2:nlevs))
    nodal = .true.

    max_nlevel        = mgt(nlevs)%max_nlevel
    max_iter          = mgt(nlevs)%max_iter
    eps               = mgt(nlevs)%eps
    smoother          = mgt(nlevs)%smoother
    nu1               = mgt(nlevs)%nu1
    nu2               = mgt(nlevs)%nu2
    gamma             = mgt(nlevs)%gamma
    omega             = mgt(nlevs)%omega
    cycle             = mgt(nlevs)%cycle
    bottom_solver     = mgt(nlevs)%bottom_solver
    bottom_solver_eps = mgt(nlevs)%bottom_solver_eps
    bottom_max_iter   = mgt(nlevs)%bottom_max_iter
    min_width         = mgt(nlevs)%min_width
    verbose           = mgt(nlevs)%verbose

    verbose = mg_verbose

    ! Note: put this here to minimize asymmetries - ASA
    eps = 1.d-12

    bottom_solver = 2
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
          call multifab_build(one_sided_ss(n), la, ns, 0, nodal)
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
                           cycle = cycle, &
                           omega = omega, &
                           bottom_solver = bottom_solver, &
                           bottom_max_iter = bottom_max_iter, &
                           bottom_solver_eps = bottom_solver_eps, &
                           max_iter = max_iter, &
                           max_nlevel = max_nlevel_in, &
                           min_width = min_width, &
                           eps = eps, &
                           verbose = verbose, &
                           cg_verbose = cg_verbose, &
                           nodal = nodal)
       
    end do

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

       !    NOTE: we define the stencils with the finest dx.
       do i = mgt(n)%nlevels, 1, -1
          call stencil_fill_nodal(mgt(n)%ss(i), coeffs(i), mgt(n)%dh(:,i), &
                                  mgt(n)%mm(i), mgt(n)%face_type,stencil_type)
          pd  = coarsen(pd,2)
       end do
       if (stencil_type .eq. ST_CROSS .and. n .gt. 1) then
          i = mgt(n)%nlevels
          call stencil_fill_one_sided(one_sided_ss(n), coeffs(i), &
                                      mgt(n    )%dh(:,i), &
                                      mgt(n)%mm(i), mgt(n)%face_type)
       end if

       do i = mgt(n)%nlevels, 1, -1
          call multifab_destroy(coeffs(i))
       end do
       deallocate(coeffs)

    end do

    allocate(rh(nlevs))
    do n = 1, nlevs
       call multifab_build(rh(n),mla%la(n),1,1,nodal)
       call setval(rh(n),ZERO,all=.true.)
    end do

    call divu(nlevs,mgt,unew,rh,mla%mba%rr,divu_verbose,nodal)

    ! Do rh = rh - divu_rhs
    if (present(divu_rhs)) then
       do n = 1, nlevs
          call multifab_sub_sub(rh(n),divu_rhs(n))
       end do
    end if

    if ( mg_verbose >= 3 ) then
       do_diagnostics = 1
    else
       do_diagnostics = 0
    end if

    if (present(eps_in)) then
       call ml_nd_solve(mla,mgt,rh,phi,one_sided_ss,mla%mba%rr,do_diagnostics,eps_in)
    else
       call ml_nd_solve(mla,mgt,rh,phi,one_sided_ss,mla%mba%rr,do_diagnostics)
    end if

    do n = nlevs,1,-1
       call multifab_fill_boundary(phi(n))
    end do

    do n = 1, nlevs
       call mg_tower_destroy(mgt(n))
       call multifab_destroy(rh(n))
    end do

    if (stencil_type .ne. ST_DENSE) then
       do n = nlevs, 2, -1
          call multifab_destroy(one_sided_ss(n))
       end do
    endif

    deallocate(mgt)
    deallocate(rh)
    deallocate(nodal)
    deallocate(one_sided_ss)

  end subroutine hg_multigrid

  !   ********************************************************************************* !

  subroutine mkcoeffs(rho,coeffs)

    type(multifab) , intent(in   ) :: rho
    type(multifab) , intent(inout) :: coeffs

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    integer :: i,dm,ng

    dm = rho%dim
    ng = rho%ng

    do i = 1, rho%nboxes
       if ( multifab_remote(rho, i) ) cycle
       rp => dataptr(rho   , i)
       cp => dataptr(coeffs, i)
       select case (dm)
       case (2)
          call mkcoeffs_2d(cp(:,:,1,1), rp(:,:,1,1), ng)
       case (3)
          call mkcoeffs_3d(cp(:,:,:,1), rp(:,:,:,1), ng)
       end select
    end do

  end subroutine mkcoeffs

  !   *********************************************************************************** !

  subroutine mkcoeffs_2d(coeffs,rho,ng)

      use bl_constants_module

    integer :: ng
    real(kind=dp_t), intent(inout) :: coeffs(   0:,   0:)
    real(kind=dp_t), intent(in   ) ::  rho(1-ng:,1-ng:)

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

  subroutine mkcoeffs_3d(coeffs,rho,ng)

      use bl_constants_module

    integer :: ng
    real(kind=dp_t), intent(inout) :: coeffs(   0:,   0:, 0:)
    real(kind=dp_t), intent(in   ) ::  rho(1-ng:,1-ng:,1-ng:)

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

  subroutine mult_by_1d_coeff(nlevs,u,div_coeff,do_mult)

    integer       , intent(in   )           :: nlevs
    type(multifab), intent(inout)           :: u(:)
    real(dp_t)    , intent(in   )           :: div_coeff(:,:)
    logical       , intent(in   ), optional :: do_mult

    ! local
    real(kind=dp_t), pointer :: ump(:,:,:,:) 
    integer :: i,ng,n,dm
    integer :: lo(u(1)%dim),hi(u(1)%dim)
    logical :: local_do_mult

    local_do_mult = .true.
    if (present(do_mult)) local_do_mult = do_mult

    ng = u(1)%ng
    dm = u(1)%dim

    do n = 1, nlevs

       ! Multiply u by div coeff
       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          ump => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call mult_by_1d_coeff_2d(ump(:,:,1,:),div_coeff(n,:),lo,hi,ng,local_do_mult)
          case (3)
             call mult_by_1d_coeff_3d(ump(:,:,:,:),div_coeff(n,:),lo,hi,ng,local_do_mult)
          end select
       end do

    end do

  end subroutine mult_by_1d_coeff

  subroutine mult_by_1d_coeff_2d(u,div_coeff,lo,hi,ng,do_mult)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: u(lo(1)-ng:,lo(2)-ng:,:)
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

  subroutine mult_by_1d_coeff_3d(u,div_coeff,lo,hi,ng,do_mult)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
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

  subroutine mult_by_3d_coeff(nlevs,u,div_coeff,do_mult)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:)
    type(multifab) , intent(in   ) :: div_coeff(:)
    logical        , intent(in   ) :: do_mult

    ! local
    real(kind=dp_t), pointer :: ump(:,:,:,:) 
    real(kind=dp_t), pointer ::  dp(:,:,:,:) 
    integer :: i,ngu,ngd,n,dm

    ngu = u(1)%ng
    ngd = div_coeff(1)%ng
    dm = u(1)%dim

    do n = 1, nlevs
       ! Multiply u by div coeff
       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          ump => dataptr(u(n),i)
          dp => dataptr(div_coeff(n),i)
          select case (dm)
          case (3)
             call mult_by_3d_coeff_3d(ump(:,:,:,:), ngu, dp(:,:,:,1), ngd, do_mult)
          end select
       end do

    end do

  end subroutine mult_by_3d_coeff

  subroutine mult_by_3d_coeff_3d(u,ngu,div_coeff,ngd,do_mult)

    integer        , intent(in   ) :: ngu,ngd
    real(kind=dp_t), intent(inout) ::         u(-ngu:,-ngu:,-ngu:,:)
    real(dp_t)     , intent(in   ) :: div_coeff(-ngd:,-ngd:,-ngd:)
    logical        , intent(in   ) :: do_mult

    integer :: i,j,k,nx,ny,nz
    nx = size(u,dim=1) - 2*ngu
    ny = size(u,dim=2) - 2*ngu
    nz = size(u,dim=3) - 2*ngu

    if (do_mult) then
       do k = 0,nz-1 
          do j = 0,ny-1 
             do i = 0,nx-1 
                u(i,j,k,:) = u(i,j,k,:) * div_coeff(i,j,k)
             end do
          end do
       end do
    else
       do k = 0,nz-1 
          do j = 0,ny-1 
             do i = 0,nx-1 
                u(i,j,k,:) = u(i,j,k,:) / div_coeff(i,j,k)
             end do
          end do
       end do
    end if

  end subroutine mult_by_3d_coeff_3d

end module hgproject_module
