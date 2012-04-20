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

  subroutine hgproject(proj_type,mla,unew,uold,rhohalf,pi,gpi,dx,dt,the_bc_tower, &
                       div_coeff_3d,divu_rhs,eps_in)

    use bc_module
    use bl_constants_module
    use bl_prof_module
    use proj_parameters
    use multifab_fill_ghost_module , only : multifab_fill_ghost_cells
    use ml_restriction_module      , only : ml_cc_restriction
    use hg_multigrid_module        , only : hg_multigrid
    use hg_hypre_module            , only : hg_hypre

    use mg_eps_module              , only : eps_hg, eps_hg_max, hg_level_factor
    use probin_module              , only: verbose, mg_verbose, cg_verbose, hg_dense_stencil, nodal, &
                                           use_hypre
    use enforce_outflow_on_divu_module, only : enforce_outflow_on_divu_rhs


    integer        , intent(in   ) :: proj_type
    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: uold(:)
    type(multifab ), intent(inout) :: rhohalf(:)
    type(multifab ), intent(inout) :: pi(:)
    type(multifab ), intent(inout) :: gpi(:)
    real(dp_t)     , intent(in   ) :: dx(:,:),dt
    type(bc_tower ), intent(in   ) :: the_bc_tower
    type(multifab ), intent(in   ) :: div_coeff_3d(:)

    type(multifab ), intent(inout), optional :: divu_rhs(:)
    real(dp_t)     , intent(in   ), optional :: eps_in

    ! Local  
    type(multifab) :: phi(mla%nlevel)
    type(multifab) :: gphi(mla%nlevel)
    type(multifab) ::   rh(mla%nlevel)

    integer                   :: n,d,stencil_type,dm,nlevs
    real(dp_t)                :: umin,umax,vmin,vmax,wmin,wmax
    real(dp_t)                :: rel_solver_eps
    real(dp_t)                :: abs_solver_eps
    type(bl_prof_timer), save :: bpt

    call build(bpt, "hgproject")

    dm = mla%dim
    nlevs = mla%nlevel

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
          if (dm .ge. 2) then
             vmin = min(vmin,multifab_min_c(unew(n),2))
             vmax = max(vmax,multifab_max_c(unew(n),2))
          end if
          if (dm .eq. 3) then
             wmin = min(wmin,multifab_min_c(unew(n),3))
             wmax = max(wmax,multifab_max_c(unew(n),3))
          end if
       end do
       if (parallel_IOProcessor()) then
          write(6,1001) umin,umax
          if (dm .ge. 2) write(6,1002) vmin,vmax
          if (dm .eq. 3) write(6,1003) wmin,wmax
          write(6,1004)
       end if
    end if

1001 format('... x-velocity before projection ',e17.10,2x,e17.10)
1002 format('... y-velocity before projection ',e17.10,2x,e17.10)
1003 format('... z-velocity before projection ',e17.10,2x,e17.10)
1004 format(' ')

    call create_uvec_for_projection(unew,uold,rhohalf,gpi,dt,the_bc_tower,proj_type)

    do n=1,nlevs
       do d=1,dm
          call multifab_mult_mult_c(unew(n),d,div_coeff_3d(n),1,1,nghost(div_coeff_3d(n)))
       end do
       call multifab_div_div_c(rhohalf(n),1,div_coeff_3d(n),1,1,nghost(div_coeff_3d(n)))
    end do

    if (present(divu_rhs)) then
       call enforce_outflow_on_divu_rhs(divu_rhs,the_bc_tower)
    end if

    do n = 1, nlevs
       call multifab_build(phi(n), mla%la(n), 1, 1, nodal)
       call multifab_build( rh(n), mla%la(n), 1, 1, nodal)
       call setval(phi(n),ZERO,all=.true.)
       call setval( rh(n),ZERO,all=.true.)
    end do

!   if (dm .eq. 1) then
!      call hg_1d_solver(mla,unew,rhohalf,phi,dx,the_bc_tower,divu_rhs)
!   else 

    if (present(eps_in)) then
       rel_solver_eps = eps_in
    else
       rel_solver_eps = min( eps_hg_max, eps_hg*hg_level_factor**(nlevs-1) )
    end if

    abs_solver_eps = -1.d0

    if (use_hypre) then 
       if (present(divu_rhs)) then
          call hg_hypre(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                        stencil_type,rel_solver_eps,abs_solver_eps,divu_rhs)
       else
          call hg_hypre(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                        stencil_type,rel_solver_eps,abs_solver_eps)
       end if
    else
       if (present(divu_rhs)) then
          call hg_multigrid(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                            stencil_type,rel_solver_eps,abs_solver_eps,divu_rhs)
       else
          call hg_multigrid(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                            stencil_type,rel_solver_eps,abs_solver_eps)
       end if
    endif

    do n = 1,nlevs
       call destroy(rh(n))
    end do

    do n=1,nlevs
       do d=1,dm
          call multifab_div_div_c(unew(n),d,div_coeff_3d(n),1,1,nghost(div_coeff_3d(n)))
       end do
       call multifab_mult_mult_c(rhohalf(n),1,div_coeff_3d(n),1,1,nghost(div_coeff_3d(n)))
    end do

    do n = 1, nlevs
       call multifab_build(gphi(n), mla%la(n), dm, 0)
    end do

    call mkgphi(gphi,phi,dx)

    call hg_update(proj_type,unew,uold,gpi,gphi,rhohalf,  &
                   pi,phi,dt,mla,the_bc_tower%bc_tower_array)

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
          if (dm .ge. 2) then
             vmin = min(vmin,multifab_min_c(unew(n),2))
             vmax = max(vmax,multifab_max_c(unew(n),2))
          end if
          if (dm .eq. 3) then
             wmin = min(wmin,multifab_min_c(unew(n),3))
             wmax = max(wmax,multifab_max_c(unew(n),3))
          end if
       end do
       if (parallel_IOProcessor()) then
          write(6,1101) umin,umax
          if (dm .ge. 2) write(6,1102) vmin,vmax
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

    subroutine create_uvec_for_projection(unew,uold,rhohalf,gpi,dt,the_bc_tower, &
                                          proj_type)

      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(in   ) :: uold(:)
      type(multifab) , intent(in   ) :: rhohalf(:)
      type(multifab) , intent(inout) :: gpi(:)
      real(kind=dp_t), intent(in   ) :: dt
      type(bc_tower) , intent(in   ) :: the_bc_tower
      integer        , intent(in   ) :: proj_type
  
      type(bc_level) :: bc
  
      real(kind=dp_t), pointer :: unp(:,:,:,:)
      real(kind=dp_t), pointer :: uop(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:)
  
      integer :: i,n
      integer :: lo(unew(1)%dim),hi(unew(1)%dim)
      integer :: ng_un,ng_uo,ng_rh,ng_gp

      type(bl_prof_timer), save :: bpt

      call build(bpt, "create_uvec_for_projection")
  
      ng_un = nghost(unew(1))
      ng_uo = nghost(uold(1))
      ng_rh = nghost(rhohalf(1))
      ng_gp = nghost(gpi(1))
  
      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, nboxes(unew(n))
            if ( multifab_remote(unew(n), i) ) cycle
            lo = lwb(get_box(unew(n),i))
            hi = lwb(get_box(unew(n),i))
            unp => dataptr(unew(n)     , i) 
            uop => dataptr(uold(n)     , i) 
            gpp => dataptr(gpi(n)       , i)
             rp => dataptr(  rhohalf(n), i)
            select case (dm)
               case (1)
                 call create_uvec_1d(unp(:,1,1,1), ng_un, uop(:,1,1,1), ng_uo, &
                                      rp(:,1,1,1), ng_rh, gpp(:,1,1,1), ng_gp, &
                                     lo,hi,dt,bc%phys_bc_level_array(i,:,:), proj_type)
               case (2)
                 call create_uvec_2d(unp(:,:,1,:), ng_un, uop(:,:,1,:), ng_uo, &
                                      rp(:,:,1,1), ng_rh, gpp(:,:,1,:), ng_gp, &
                                     lo,hi,dt,bc%phys_bc_level_array(i,:,:), proj_type)
               case (3)
                 call create_uvec_3d(unp(:,:,:,:), ng_un, uop(:,:,:,:), ng_uo, &
                                      rp(:,:,:,1), ng_rh, gpp(:,:,:,:), ng_gp, &
                                     lo,hi,dt, bc%phys_bc_level_array(i,:,:), proj_type)
            end select
         end do
         call multifab_fill_boundary(unew(n))
      end do 

      call destroy(bpt)

    end subroutine create_uvec_for_projection

    !   ******************************************************************************** !

    subroutine create_uvec_1d(unew,ng_un,uold,ng_uo,rhohalf,ng_rh,gpi,ng_gp, &
                              lo,hi,dt,phys_bc,proj_type)

      use proj_parameters

      integer        , intent(in   ) :: ng_un,ng_uo,ng_rh,ng_gp
      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::    unew(lo(1)-ng_un:)
      real(kind=dp_t), intent(in   ) ::    uold(lo(1)-ng_uo:)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-ng_rh:)
      real(kind=dp_t), intent(inout) ::     gpi(lo(1)-ng_gp:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      if (phys_bc(1,1) .eq. INLET) gpi(lo(1)-1) = ZERO
      if (phys_bc(1,2) .eq. INLET) gpi(hi(1)  ) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection_comp) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters_comp) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters_comp) then

         unew(lo(1)-1:hi(1)+1) = ( unew(lo(1)-1:hi(1)+1) - uold(lo(1)-1:hi(1)+1) ) / dt

      ! quantity projected is Ustar + dt * (1/rho) gpi
      else if (proj_type .eq. regular_timestep_comp) then

         unew(lo(1)-1:hi(1)+1) = unew(lo(1)-1:hi(1)+1) + dt*gpi(lo(1)-1:hi(1)+1)/rhohalf(lo(1)-1:hi(1)+1)

       else

          call bl_error('No proj_type by this number ')

      end if

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL .or. &
          phys_bc(1,1)==SYMMETRY ) unew(lo(1)-1) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL .or. &
          phys_bc(1,2)==SYMMETRY ) unew(hi(1)+1) = ZERO

    end subroutine create_uvec_1d


    !   ******************************************************************************** !

    subroutine create_uvec_2d(unew,ng_un,uold,ng_uo,rhohalf,ng_rh,gpi,ng_gp, &
                              lo,hi,dt,phys_bc,proj_type)

      use proj_parameters

      integer        , intent(in   ) :: ng_un,ng_uo,ng_rh,ng_gp
      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::    unew(lo(1)-ng_un:,lo(2)-ng_un:,:)
      real(kind=dp_t), intent(in   ) ::    uold(lo(1)-ng_uo:,lo(2)-ng_uo:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-ng_rh:,lo(2)-ng_rh:)
      real(kind=dp_t), intent(inout) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      if (phys_bc(1,1) .eq. INLET) gpi(lo(1)-1,lo(2)-1:hi(2)+2,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gpi(hi(1)+1,lo(2)-1:hi(2)+2,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gpi(lo(1)-1:hi(1)+1,lo(2)-1,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gpi(lo(1)-1:hi(1)+1,hi(2)+2,:) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection_comp) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters_comp) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters_comp) then

         unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1) = ( &
            unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1) - uold(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1) ) / dt
         unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2) = ( &
            unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2) - uold(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2) ) / dt
     
      ! quantity projected is Ustar + dt * (1/rho) gpi
      else if (proj_type .eq. regular_timestep_comp) then

         unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1) = unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1) + &
            dt*gpi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)/rhohalf(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)
         unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2) = unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2) + &
            dt*gpi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2)/rhohalf(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)

       else
     
          call bl_error('No proj_type by this number ')

      end if

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL .or. &
          phys_bc(1,1)==SYMMETRY ) unew(lo(1)-1,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL .or. &
          phys_bc(1,2)==SYMMETRY ) unew(hi(1)+1,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL .or. &
          phys_bc(2,1)==SYMMETRY ) unew(:,lo(2)-1,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL .or. &
          phys_bc(2,2)==SYMMETRY ) unew(:,hi(2)+1,:) = ZERO

    end subroutine create_uvec_2d

    !  *********************************************************************************** !

    subroutine create_uvec_3d(unew,ng_un,uold,ng_uo,rhohalf,ng_rh,gpi,ng_gp, &
                              lo,hi,dt,phys_bc,proj_type)

      use proj_parameters

      integer        , intent(in   ) :: ng_un,ng_uo,ng_rh,ng_gp
      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::    unew(lo(1)-ng_un:,lo(2)-ng_un:,lo(3)-ng_un:,:)
      real(kind=dp_t), intent(in   ) ::    uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
      real(kind=dp_t), intent(inout) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-ng_rh:,lo(2)-ng_rh:,lo(3)-ng_rh:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      integer :: i,j,k,m

      if (phys_bc(1,1) .eq. INLET) gpi(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gpi(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gpi(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gpi(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1,:) = ZERO
      if (phys_bc(3,1) .eq. INLET) gpi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1,:) = ZERO
      if (phys_bc(3,2) .eq. INLET) gpi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1,:) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection_comp) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters_comp) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters_comp) then

         !$OMP PARALLEL PRIVATE(i,j,k,m)
         do m=1,3
            !$OMP DO
            do k=lo(3)-1,hi(3)+1
               do j=lo(2)-1,hi(2)+1
                  do i=lo(1)-1,hi(1)+1
                     unew(i,j,k,m) = (unew(i,j,k,m) - uold(i,j,k,m)) / dt
                  end do
               end do
            end do
            !$OMP END DO NOWAIT
         end do
         !$OMP END PARALLEL

      ! quantity projected is Ustar + dt * (1/rho) gpi
      else if (proj_type .eq. regular_timestep_comp) then

         !$OMP PARALLEL PRIVATE(i,j,k,m)
         do m=1,3
            !$OMP DO
            do k=lo(3)-1,hi(3)+1
               do j=lo(2)-1,hi(2)+1
                  do i=lo(1)-1,hi(1)+1
                     unew(i,j,k,m) = unew(i,j,k,m) + dt*gpi(i,j,k,m)/rhohalf(i,j,k)
                  end do
               end do
            end do
            !$OMP END DO NOWAIT
         end do
         !$OMP END PARALLEL

      else

          call bl_error('No proj_type by this number ')

      end if

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL .or. &
          phys_bc(1,1)==SYMMETRY) unew(lo(1)-1,:,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL .or. &
          phys_bc(1,2)==SYMMETRY) unew(hi(1)+1,:,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL .or. &
          phys_bc(2,1)==SYMMETRY) unew(:,lo(2)-1,:,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL .or. &
          phys_bc(2,2)==SYMMETRY) unew(:,hi(2)+1,:,:) = ZERO
      if (phys_bc(3,1)==SLIP_WALL .or. phys_bc(3,1)==NO_SLIP_WALL .or. &
          phys_bc(3,1)==SYMMETRY) unew(:,:,lo(3)-1,:) = ZERO
      if (phys_bc(3,2)==SLIP_WALL .or. phys_bc(3,2)==NO_SLIP_WALL .or. &
          phys_bc(3,2)==SYMMETRY) unew(:,:,hi(3)+1,:) = ZERO

    end subroutine create_uvec_3d

    ! ******************************************************************************* !

    subroutine mkgphi(gphi,phi,dx)

      type(multifab), intent(inout) :: gphi(:)
      type(multifab), intent(in   ) :: phi(:)
      real(dp_t) :: dx(:,:)

      integer :: lo(phi(1)%dim),hi(phi(1)%dim)
      integer :: i,n,ng_p,ng_gp

      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: pp(:,:,:,:) 

      type(bl_prof_timer), save :: bpt

      call build(bpt, "mkgphi")

      ng_p  = nghost(phi(1))
      ng_gp = nghost(gphi(1))

      do n = 1, nlevs

         do i = 1, nboxes(phi(n))
            if ( multifab_remote(phi(n),i) ) cycle
            lo = lwb(get_box(gphi(n),i))
            hi = lwb(get_box(gphi(n),i))
            gph => dataptr(gphi(n),i)
            pp  => dataptr(phi(n),i)
            select case (dm)
            case (1)
               call mkgphi_1d(gph(:,1,1,1), ng_gp, pp(:,1,1,1), ng_p, lo, hi, dx(n,:))
            case (2)
               call mkgphi_2d(gph(:,:,1,:), ng_gp, pp(:,:,1,1), ng_p, lo, hi, dx(n,:))
            case (3)
               call mkgphi_3d(gph(:,:,:,:), ng_gp, pp(:,:,:,1), ng_p, lo, hi, dx(n,:))
            end select
         end do

      end do

      call destroy(bpt)

    end subroutine mkgphi

    !   ********************************************************************************* !

    subroutine mkgphi_1d(gphi,ng_gp,phi,ng_p,lo,hi,dx)

      integer        , intent(in   ) :: ng_gp,ng_p
      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: gphi(lo(1)-ng_gp:)
      real(kind=dp_t), intent(inout) ::  phi(lo(1)-ng_p :)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i

      do i = lo(1),hi(1)
         gphi(i) = ( phi(i+1) - phi(i) ) / dx(1)
      end do

    end subroutine mkgphi_1d

    !   ********************************************************************************* !

    subroutine mkgphi_2d(gphi,ng_gp,phi,ng_p,lo,hi,dx)

      integer        , intent(in   ) :: ng_gp,ng_p
      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: gphi(lo(1)-ng_gp:,lo(2)-ng_gp:,:)
      real(kind=dp_t), intent(inout) ::  phi(lo(1)-ng_p :,lo(2)-ng_p :)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            gphi(i,j,1) = HALF*(phi(i+1,j) + phi(i+1,j+1) - &
                 phi(i  ,j) - phi(i  ,j+1) ) /dx(1)
            gphi(i,j,2) = HALF*(phi(i,j+1) + phi(i+1,j+1) - &
                 phi(i,j  ) - phi(i+1,j  ) ) /dx(2)
         end do
      end do

    end subroutine mkgphi_2d

    !   ******************************************************************************** !

    subroutine mkgphi_3d(gphi,ng_gp,phi,ng_p,lo,hi,dx)

      integer        , intent(in   ) :: ng_gp,ng_p
      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: gphi(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
      real(kind=dp_t), intent(inout) ::  phi(lo(1)-ng_p :,lo(2)-ng_p :,lo(3)-ng_p :)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
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
      !$OMP END PARALLEL DO

    end subroutine mkgphi_3d

    ! ******************************************************************************* !

    subroutine hg_update(proj_type,unew,uold,gpi,gphi,rhohalf,pi,phi,dt, &
                         mla,the_bc_level)

      use multifab_physbc_module
      use variables, only: foextrap_comp

      integer        , intent(in   ) :: proj_type
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(in   ) :: uold(:)
      type(multifab) , intent(inout) :: gpi(:)
      type(multifab) , intent(in   ) :: gphi(:)
      type(multifab) , intent(in   ) :: rhohalf(:)
      type(multifab) , intent(inout) :: pi(:)
      type(multifab) , intent(in   ) :: phi(:)
      real(kind=dp_t), intent(in   ) :: dt
      type(ml_layout), intent(in   ) :: mla
      type(bc_level) , intent(in   ) :: the_bc_level(:)

      ! local
      integer :: i,n
      integer :: ng_un,ng_uo,ng_gp,ng_gh,ng_rh,ng_p,ng_h
      integer :: lo(unew(1)%dim),hi(unew(1)%dim)

      real(kind=dp_t), pointer :: upn(:,:,:,:) 
      real(kind=dp_t), pointer :: uon(:,:,:,:) 
      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:) 
      real(kind=dp_t), pointer ::  ph(:,:,:,:) 
      real(kind=dp_t), pointer ::  pp(:,:,:,:) 

      type(bl_prof_timer), save :: bpt

      call build(bpt, "hg_update")

      ng_un = nghost(unew(1))
      ng_uo = nghost(uold(1))
      ng_gp = nghost(gpi(1))
      ng_gh = nghost(gphi(1))
      ng_rh = nghost(rhohalf(1))
      ng_p  = nghost(pi(1))
      ng_h  = nghost(phi(1))

      do n = 1, nlevs

         do i = 1, nboxes(unew(n))
            if ( multifab_remote(unew(n),i) ) cycle
            lo = lwb(get_box(unew(n),i))
            hi = lwb(get_box(unew(n),i))
            upn => dataptr(unew(n),i)
            uon => dataptr(uold(n),i)
            gpp => dataptr(gpi(n),i)
            gph => dataptr(gphi(n),i)
            rp  => dataptr(rhohalf(n),i)
            pp  => dataptr(pi(n),i)
            ph  => dataptr(phi(n),i)
            select case (dm)
            case (1)
               call hg_update_1d(proj_type, upn(:,1,1,1), ng_un, uon(:,1,1,1), ng_uo, &
                                 gpp(:,1,1,1), ng_gp, gph(:,1,1,1), ng_gh, rp(:,1,1,1), &
                                 ng_rh, pp(:,1,1,1), ng_p, ph(:,1,1,1), ng_h, lo, hi, dt)
            case (2)
               call hg_update_2d(proj_type, upn(:,:,1,:), ng_un, uon(:,:,1,:), ng_uo, &
                                 gpp(:,:,1,:), ng_gp, gph(:,:,1,:), ng_gh, rp(:,:,1,1), &
                                 ng_rh, pp(:,:,1,1), ng_p, ph(:,:,1,1), ng_h, lo, hi, dt)
            case (3)
               call hg_update_3d(proj_type, upn(:,:,:,:), ng_un, uon(:,:,:,:), ng_uo, &
                                 gpp(:,:,:,:), ng_gp, gph(:,:,:,:), ng_gh, rp(:,:,:,1), &
                                 ng_rh, pp(:,:,:,1), ng_p, ph(:,:,:,1), ng_h, lo, hi, dt)
            end select
         end do

         ! fill ghost cells for two adjacent grids at the same level
         ! this includes periodic domain boundary ghost cells
         call multifab_fill_boundary(pi(n))

      end do

      if (nlevs .eq. 1) then

         ! fill ghost cells for two adjacent grids at the same level
         ! this includes periodic domain boundary ghost cells
         call multifab_fill_boundary(unew(nlevs))
         call multifab_fill_boundary(gpi(nlevs))

         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc(unew(nlevs),1,1,dm,the_bc_level(nlevs))
         do i=1,dm
            call multifab_physbc(gpi(nlevs),i,foextrap_comp,1,the_bc_level(nlevs))
         end do

      else

         ! the loop over nlevs must count backwards to make sure the finer grids are 
         ! done first
         do n=nlevs,2,-1

            ! set level n-1 data to be the average of the level n data covering it
            call ml_cc_restriction(unew(n-1),unew(n),mla%mba%rr(n-1,:)) 
            call ml_cc_restriction(gpi(n-1),gpi(n),mla%mba%rr(n-1,:))

            ! fill level n ghost cells using interpolation from level n-1 data
            ! note that multifab_fill_boundary and multifab_physbc are called for
            ! both levels n-1 and n
            call multifab_fill_ghost_cells(unew(n),unew(n-1),nghost(unew(n)),mla%mba%rr(n-1,:), &
                                           the_bc_level(n-1),the_bc_level(n),1,1,dm, &
                                           fill_crse_input=.false.)
            do i=1,dm
               call multifab_fill_ghost_cells(gpi(n),gpi(n-1),1,mla%mba%rr(n-1,:), &
                                              the_bc_level(n-1),the_bc_level(n),i, &
                                              foextrap_comp,1,fill_crse_input=.false.)
            end do

         end do

      end if

      call destroy(bpt)

    end subroutine hg_update

    !   ****************************************************************************** !

    subroutine hg_update_1d(proj_type,unew,ng_un,uold,ng_uo,gpi,ng_gp,gphi,ng_gh, &
                            rhohalf,ng_rh,pi,ng_p,phi,ng_h,lo,hi,dt)

      use proj_parameters

      integer        , intent(in   ) :: ng_un,ng_uo,ng_gp,ng_gh,ng_rh,ng_p,ng_h
      integer        , intent(in   ) :: lo(:),hi(:)
      integer        , intent(in   ) :: proj_type
      real(kind=dp_t), intent(inout) ::    unew(lo(1)-ng_un:)
      real(kind=dp_t), intent(in   ) ::    uold(lo(1)-ng_uo:)
      real(kind=dp_t), intent(inout) ::     gpi(lo(1)-ng_gp:)
      real(kind=dp_t), intent(in   ) ::    gphi(lo(1)-ng_gh:)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-ng_rh:)
      real(kind=dp_t), intent(inout) ::      pi(lo(1)-ng_p :)
      real(kind=dp_t), intent(in   ) ::     phi(lo(1)-ng_h :)
      real(kind=dp_t), intent(in   ) :: dt

      !     Subtract off the density-weighted gradient.
      unew(lo(1):hi(1)) = unew(lo(1):hi(1)) - gphi(lo(1):hi(1))/rhohalf(lo(1):hi(1)) 

      if (proj_type .eq. pressure_iters_comp) &    ! unew held the projection of (ustar-uold)
           unew(lo(1):hi(1)) = uold(lo(1):hi(1)) + dt * unew(lo(1):hi(1))

      if ( (proj_type .eq. initial_projection_comp) .or. &
           (proj_type .eq. divu_iters_comp) ) then

         gpi = ZERO
         pi = ZERO

      else if (proj_type .eq. pressure_iters_comp) then

         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         gpi(lo(1):hi(1)  ) = gpi(lo(1):hi(1)  ) + gphi(lo(1):hi(1))
          pi(lo(1):hi(1)+1) =  pi(lo(1):hi(1)+1)  + phi(lo(1):hi(1)+1)

      else if (proj_type .eq. regular_timestep_comp) then

         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         gpi(lo(1):hi(1)  ) = (ONE/dt) * gphi(lo(1):hi(1))
          pi(lo(1):hi(1)+1) = (ONE/dt) *  phi(lo(1):hi(1)+1)

      end if

    end subroutine hg_update_1d

    !   ****************************************************************************** !

    subroutine hg_update_2d(proj_type,unew,ng_un,uold,ng_uo,gpi,ng_gp,gphi,ng_gh, &
                            rhohalf,ng_rh,pi,ng_p,phi,ng_h,lo,hi,dt)

      use proj_parameters

      integer        , intent(in   ) :: ng_un,ng_uo,ng_gp,ng_gh,ng_rh,ng_p,ng_h
      integer        , intent(in   ) :: lo(:),hi(:)
      integer        , intent(in   ) :: proj_type
      real(kind=dp_t), intent(inout) ::    unew(lo(1)-ng_un:,lo(2)-ng_un:,:)
      real(kind=dp_t), intent(in   ) ::    uold(lo(1)-ng_uo:,lo(2)-ng_uo:,:)
      real(kind=dp_t), intent(inout) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,:)
      real(kind=dp_t), intent(in   ) ::    gphi(lo(1)-ng_gh:,lo(2)-ng_gh:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-ng_rh:,lo(2)-ng_rh:)
      real(kind=dp_t), intent(inout) ::      pi(lo(1)-ng_p :,lo(2)-ng_p :)
      real(kind=dp_t), intent(in   ) ::     phi(lo(1)-ng_h :,lo(2)-ng_h :)
      real(kind=dp_t), intent(in   ) :: dt

      !     Subtract off the density-weighted gradient.
      unew(lo(1):hi(1),lo(2):hi(2),1) = unew(lo(1):hi(1),lo(2):hi(2),1) - &
          gphi(lo(1):hi(1),lo(2):hi(2),1)/rhohalf(lo(1):hi(1),lo(2):hi(2)) 
      unew(lo(1):hi(1),lo(2):hi(2),2) = unew(lo(1):hi(1),lo(2):hi(2),2) - &
          gphi(lo(1):hi(1),lo(2):hi(2),2)/rhohalf(lo(1):hi(1),lo(2):hi(2)) 

      if (proj_type .eq. pressure_iters_comp) &    ! unew held the projection of (ustar-uold)
           unew(lo(1):hi(1),lo(2):hi(2),:) = uold(lo(1):hi(1),lo(2):hi(2),:) + &
               dt * unew(lo(1):hi(1),lo(2):hi(2),:)

      if ( (proj_type .eq. initial_projection_comp) .or. &
           (proj_type .eq. divu_iters_comp) ) then

         gpi = ZERO
         pi = ZERO

      else if (proj_type .eq. pressure_iters_comp) then

         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         gpi(lo(1):hi(1)  ,lo(2):hi(2)  ,:) = gpi(lo(1):hi(1)  ,lo(2):hi(2)  ,:) + &
             gphi(lo(1):hi(1),lo(2):hi(2),:)
          pi(lo(1):hi(1)+1,lo(2):hi(2)+1)   =  pi(lo(1):hi(1)+1,lo(2):hi(2)+1)  +  &
             phi(lo(1):hi(1)+1,lo(2):hi(2)+1)

      else if (proj_type .eq. regular_timestep_comp) then

         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         gpi(lo(1):hi(1)  ,lo(2):hi(2)  ,:) = (ONE/dt) * gphi(lo(1):hi(1)  ,lo(2):hi(2)  ,:)
          pi(lo(1):hi(1)+1,lo(2):hi(2)+1)   = (ONE/dt) *  phi(lo(1):hi(1)+1,lo(2):hi(2)+1)

      end if

    end subroutine hg_update_2d

    !   ******************************************************************************* !

    subroutine hg_update_3d(proj_type,unew,ng_un,uold,ng_uo,gpi,ng_gp,gphi,ng_gh, &
                            rhohalf,ng_rh,pi,ng_p,phi,ng_h,lo,hi,dt)

      use proj_parameters

      integer        , intent(in   ) :: ng_un,ng_uo,ng_gp,ng_gh,ng_rh,ng_p,ng_h
      integer        , intent(in   ) :: lo(:),hi(:)
      integer        , intent(in   ) :: proj_type
      real(kind=dp_t), intent(inout) ::    unew(lo(1)-ng_un:,lo(2)-ng_un:,lo(3)-ng_un:,:)
      real(kind=dp_t), intent(in   ) ::    uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
      real(kind=dp_t), intent(inout) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
      real(kind=dp_t), intent(in   ) ::    gphi(lo(1)-ng_gh:,lo(2)-ng_gh:,lo(3)-ng_gh:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-ng_rh:,lo(2)-ng_rh:,lo(3)-ng_rh:)
      real(kind=dp_t), intent(inout) ::      pi(lo(1)-ng_p :,lo(2)-ng_p :,lo(3)-ng_p :)
      real(kind=dp_t), intent(in   ) ::     phi(lo(1)-ng_h :,lo(2)-ng_h :,lo(3)-ng_h :)
      real(kind=dp_t), intent(in   ) :: dt

      integer         :: i,j,k,m
      !
      ! Subtract off the density-weighted gradient
      !
      !$OMP PARALLEL PRIVATE(i,j,k,m)
      do m=1,3
         !$OMP DO
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  unew(i,j,k,m) = unew(i,j,k,m) - gphi(i,j,k,m)/rhohalf(i,j,k)
               end do
            end do
         end do
         !$OMP END DO NOWAIT
      end do
      !$OMP END PARALLEL

      if (proj_type .eq. pressure_iters_comp) then
         !
         ! unew held the projection of (ustar-uold)
         !
         !$OMP PARALLEL PRIVATE(i,j,k,m)
         do m=1,3
            !$OMP DO
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  do i=lo(1),hi(1)
                     unew(i,j,k,m) = uold(i,j,k,m) + dt*unew(i,j,k,m)
                  end do
               end do
            end do
            !$OMP END DO NOWAIT
         end do
         !$OMP END PARALLEL
      end if

      if ( (proj_type .eq. initial_projection_comp) .or. &
           (proj_type .eq. divu_iters_comp) ) then

         gpi = ZERO
         pi = ZERO

      else if (proj_type .eq. pressure_iters_comp) then
         !
         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         !
         !$OMP PARALLEL PRIVATE(i,j,k,m)
         do m=1,3
            !$OMP DO
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  do i=lo(1),hi(1)
                     gpi(i,j,k,m) = gpi(i,j,k,m) + gphi(i,j,k,m)
                  end do
               end do
            end do
            !$OMP END DO NOWAIT
         end do
         !$OMP END PARALLEL

         !$OMP PARALLEL DO PRIVATE(i,j,k)
            do k=lo(3),hi(3)+1
               do j=lo(2),hi(2)+1
                  do i=lo(1),hi(1)+1
                  pi(i,j,k) = pi(i,j,k) + phi(i,j,k)
               end do
            end do
         end do
         !$OMP END PARALLEL DO

      else if (proj_type .eq. regular_timestep_comp) then
         !
         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         !
         !$OMP PARALLEL PRIVATE(i,j,k,m)
         do m=1,3
            !$OMP DO
            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  do i=lo(1),hi(1)
                     gpi(i,j,k,m) = gphi(i,j,k,m) / dt
                  end do
               end do
            end do
            !$OMP END DO
         end do
         !$OMP END PARALLEL

         !$OMP PARALLEL DO PRIVATE(i,j,k)
         do k=lo(3),hi(3)+1
            do j=lo(2),hi(2)+1
               do i=lo(1),hi(1)+1
                  pi(i,j,k) = phi(i,j,k) / dt
               end do
            end do
         end do
         !$OMP END PARALLEL DO

      end if

    end subroutine hg_update_3d

  end subroutine hgproject

  ! ******************************************************************************** !

  subroutine hg_1d_solver(mla,rh,unew,rhohalf,phi,dx,the_bc_tower,divu_rhs)

    use bl_prof_module
    use bl_constants_module
    use ml_solve_module
    use nodal_divu_module
    use probin_module, only : nodal
    use probin_module, only : verbose, nodal
    use variables, only: press_comp

    type(ml_layout), intent(inout) :: mla
    type(multifab ), intent(inout) :: rh(:)
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: rhohalf(:)
    type(multifab ), intent(inout) :: phi(:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower

    type(multifab ), intent(in   ), optional :: divu_rhs(:)

    ! Local variables
    type(box     ) :: pd
    type(  layout) :: la
    type(mg_tower) :: mgt(mla%nlevel)

    type(multifab), allocatable :: coeffs(:)

    integer :: n, ns, dm, nlevs
    integer :: lo_inflow(1),hi_inflow(1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "hg_1d_solver")
    
    dm = mla%dim
    nlevs = mla%nlevel

    ns = 3

    do n = nlevs, 1, -1

       pd = layout_get_pd(mla%la(n))

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,press_comp), &
                           dh = dx(n,:), &
                           ns = ns, &
                           nodal = nodal)

    end do


    !! Fill coefficient array
    allocate(coeffs(nlevs))
    do n = nlevs,1,-1

       la = mla%la(n)
       pd = layout_get_pd(la)

       call multifab_build(coeffs(n), la, 1, 1)
       call setval(coeffs(n), 0.0_dp_t, 1, all=.true.)

       call mkcoeffs(rhohalf(n),coeffs(n))
       call multifab_fill_boundary(coeffs(n))

    end do

    ! Set the inflow array -- 1 if inflow, otherwise 0
    lo_inflow(1) = 0
    hi_inflow(1) = 0
    if (the_bc_tower%bc_tower_array(1)%phys_bc_level_array(0,1,1) == INLET) then
       lo_inflow(1) = 1
    end if
    if (the_bc_tower%bc_tower_array(1)%phys_bc_level_array(0,1,2) == INLET) then
       hi_inflow(1) = 1
    end if

    call divu(nlevs,mgt,unew,rh,mla%mba%rr,nodal,lo_inflow,hi_inflow)

    ! Do rh = rh - divu_rhs (this routine preserves rh=0 on
    !  nodes which have bc_dirichlet = true.
    if (present(divu_rhs)) &
       call subtract_divu_from_rh(nlevs,mgt,rh,divu_rhs)

!   call solve_1d(mla,mgt,rh,phi,mla%mba%rr)

    do n = nlevs,1,-1
       call multifab_fill_boundary(phi(n))
    end do

    do n = 1, nlevs
       call destroy(coeffs(n))
       call mg_tower_destroy(mgt(n))
    end do
    deallocate(coeffs)

    call destroy(bpt)

  end subroutine hg_1d_solver

  !   ********************************************************************************* !

  subroutine mkcoeffs(rho,coeffs)

    type(multifab) , intent(in   ) :: rho
    type(multifab) , intent(inout) :: coeffs

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    integer :: lo(rho%dim),hi(rho%dim)
    integer :: i,ng_r,ng_c,dm

    dm = get_dim(rho)

    ng_r = nghost(rho)
    ng_c = nghost(coeffs)

    do i = 1, nboxes(rho)
       if ( multifab_remote(rho, i) ) cycle
       lo = lwb(get_box(rho,i))
       hi = lwb(get_box(rho,i))
       rp => dataptr(rho   , i)
       cp => dataptr(coeffs, i)
       select case (dm)
       case (1)
          call mkcoeffs_1d(cp(:,1,1,1), ng_c, rp(:,1,1,1), lo, hi, ng_r)
       case (2)
          call mkcoeffs_2d(cp(:,:,1,1), ng_c, rp(:,:,1,1), lo, hi, ng_r)
       case (3)
          call mkcoeffs_3d(cp(:,:,:,1), ng_c, rp(:,:,:,1), lo, hi, ng_r)
       end select
    end do

  end subroutine mkcoeffs

  !   *********************************************************************************** !

  subroutine mkcoeffs_1d(coeffs,ng_c,rho,lo,hi,ng_r)

    use bl_constants_module

    integer                        :: ng_c,ng_r
    integer                        :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: coeffs(lo(1)-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(lo(1)-ng_r:)

    integer :: i

    do i = lo(1),hi(1)
       coeffs(i) = ONE / rho(i)
    end do

  end subroutine mkcoeffs_1d

  !   *********************************************************************************** !

  subroutine mkcoeffs_2d(coeffs,ng_c,rho,lo,hi,ng_r)

    use bl_constants_module

    integer                        :: ng_c,ng_r
    integer                        :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: coeffs(lo(1)-ng_c:,lo(2)-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(lo(1)-ng_r:,lo(2)-ng_r:)

    integer :: i,j

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          coeffs(i,j) = ONE / rho(i,j)
       end do
    end do

  end subroutine mkcoeffs_2d

  !   ********************************************************************************** !

  subroutine mkcoeffs_3d(coeffs,ng_c,rho,lo,hi,ng_r)

      use bl_constants_module

    integer                        :: ng_c,ng_r
    integer                        :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: coeffs(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)

    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             coeffs(i,j,k) = ONE / rho(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine mkcoeffs_3d

end module hgproject_module
