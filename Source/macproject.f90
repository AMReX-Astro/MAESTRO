! Project the edge-centered (MAC) velocities to satisfy the divergence
! constraint.
!
! We want to solve:  D [ (beta_0/rho) G phi ] = D ( beta_0 U ) - beta_0 S
!
! For the use_alt_energy_fix = T case, we instead solve:
!   D [ (beta_0**2/rho) G (phi/beta_0) ] = D ( beta_0 U ) - beta_0 S
! (note the extra beta_0)
!
! Here, phi is cell-centered.
! 
! after which we update the velocity as U = U - (1/rho) G phi
!
! basic outline:
!
!   -- multiply U by beta_0 (now the umac multifab holds: beta_0 U)
!
!   -- compute div ( beta_0 U)
!
!   -- store beta_0/rho (or beta_0**2/rho) in the beta multifab 
!      (mk_mac_coeffs followed by multiply either via mult_edge_by_1d_coeff 
!      or multifab_mult_mult)
!
!   -- solve the elliptic equation via multigrid
!
!   -- compute beta_0 U = beta_0 U - (beta_0/rho) G phi (mkumac subroutine)
!      (note it is beta_0 U that is updated since umac still holds beta_0 U)
!
!   -- divide out the beta_0 giving us the new U
!
module macproject_module

  use bl_types
  use bc_module
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use bl_constants_module

  implicit none

  private

  public :: macproject

contains 

  ! NOTE: this routine differs from that in varden because phi is passed in/out 
  !       rather than allocated here
  subroutine macproject(mla,umac,phi,rho,dx,the_bc_tower, &
                        divu_rhs,div_coeff_1d,div_coeff_1d_edge,div_coeff_cart_edge)

    use mac_hypre_module               , only : mac_hypre
    use mac_multigrid_module           , only : mac_multigrid
    use create_umac_grown_module       , only : create_umac_grown
    use multifab_physbc_edgevel_module , only : multifab_physbc_edgevel

    use geometry, only: spherical
    use probin_module, only: verbose, use_hypre, use_alt_energy_fix
    use variables, only: press_comp
    use ml_cc_restriction_module, only: ml_edge_restriction

    use mg_eps_module, only: eps_mac, eps_mac_max, mac_level_factor

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: umac(:,:)
    type(multifab ), intent(inout) :: phi(:)
    type(multifab ), intent(in   ) :: rho(:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower

    type(multifab ), intent(in   ), optional :: divu_rhs(:)
    real(dp_t)     , intent(in   ), optional :: div_coeff_1d(:,:)
    real(dp_t)     , intent(in   ), optional :: div_coeff_1d_edge(:,:)
    type(multifab ), intent(in   ), optional :: div_coeff_cart_edge(:,:)

    type(multifab)  :: rh(mla%nlevel),alpha(mla%nlevel),beta(mla%nlevel,mla%dim)
    type(bndry_reg) :: fine_flx(mla%nlevel)

!    real(dp_t)                   :: umac_norm(mla%nlevel)
    real(dp_t)                   :: umin,umax,vmin,vmax,wmin,wmax
    real(dp_t)                   :: rel_solver_eps
    real(dp_t)                   :: abs_solver_eps
    integer                      :: stencil_order,i,n,comp,dm,nlevs
    logical                      :: use_rhs, use_div_coeff_1d, use_div_coeff_cart_edge

    type(bl_prof_timer), save :: bpt

    call build(bpt, "macproject")

    dm = mla%dim
    nlevs = mla%nlevel

    use_rhs          = .false. ; if (present(divu_rhs)    ) use_rhs          = .true.
    use_div_coeff_1d = .false. ; if (present(div_coeff_1d)) use_div_coeff_1d = .true.
    use_div_coeff_cart_edge = .false. ; if (present(div_coeff_cart_edge)) use_div_coeff_cart_edge = .true.

    if (use_div_coeff_1d .and. use_div_coeff_cart_edge) then
       call bl_error('CANT HAVE 1D and 3D DIV_COEFF IN MACPROJECT')
    end if

    if (spherical .eq. 1 .and. use_div_coeff_1d) then
       call bl_error('CANT HAVE SPHERICAL .eq. 1 and 1D DIV_COEFF IN MACPROJECT')
    end if

    if (spherical .eq. 0 .and. use_div_coeff_cart_edge) then
       call bl_error('CANT HAVE SPHERICAL .eq. 0 and 3D DIV_COEFF IN MACPROJECT')
    end if

    stencil_order = 2

    do n = 1, nlevs
       call multifab_build(   rh(n), mla%la(n),  1, 0)
       call multifab_build(alpha(n), mla%la(n),  1, 1)
       do i = 1,dm
          call multifab_build_edge(beta(n,i),mla%la(n),1,1,i)
       end do
       call setval(alpha(n),ZERO,all=.true.)
    end do

    ! Print the norm of each component separately -- make sure to do
    !       this before we multiply the velocities by the div_coeff.
    if (verbose .eq. 1) then
       if (parallel_IOProcessor()) print *,''
       do n = 1,nlevs
          umax = multifab_max(umac(n,1))
          umin = multifab_min(umac(n,1))
          if (dm.eq.1) then
             if (parallel_IOProcessor()) then
                 write(6,1001) n, umax
                 write(6,1101) n, umin
             end if
          else if (dm.eq.2) then
             vmax = multifab_max(umac(n,2))
             vmin = multifab_min(umac(n,2))
             if (parallel_IOProcessor()) then
                 write(6,1002) n, umax, vmax
                 write(6,1102) n, umin, vmin
             end if
          else if (dm.eq.3) then
             vmax = multifab_max(umac(n,2))
             vmin = multifab_min(umac(n,2))
             wmax = multifab_max(umac(n,3))
             wmin = multifab_min(umac(n,3))
             if (parallel_IOProcessor()) then
                 write(6,1003) n, umax, vmax, wmax
                 write(6,1103) n, umin, vmin, wmin
             end if
          end if
       end do
    end if

    if (use_div_coeff_1d) then
       do n = 1,nlevs
          call mult_edge_by_1d_coeff(umac(n,:),div_coeff_1d(n,:),div_coeff_1d_edge(n,:),&
                                     .true.)
       end do
    else if (use_div_coeff_cart_edge) then
       do n=1,nlevs
          do comp=1,dm
             call multifab_mult_mult(umac(n,comp),div_coeff_cart_edge(n,comp))
          end do
       end do
    end if

    ! Compute umac_norm to be used inside the MG solver as part of a stopping criterion
!    umac_norm = -1.0_dp_t
!    do n = 1,nlevs
!       do i = 1,dm
!          umac_norm(n) = max(umac_norm(n),norm_inf(umac(n,i)))
!       end do
!    end do

    if (use_rhs) then
       call divumac(umac,rh,dx,mla%mba%rr,.true.,divu_rhs)
    else
       call divumac(umac,rh,dx,mla%mba%rr,.true.)
    end if

    ! first set beta = 1/rho
    call mk_mac_coeffs(mla,rho,beta)

    ! now make beta = beta_0/rho
    if (use_div_coeff_1d) then
       do n = 1,nlevs
          call mult_edge_by_1d_coeff(beta(n,:),div_coeff_1d(n,:), &
                                     div_coeff_1d_edge(n,:),.true.)
          if (use_alt_energy_fix) then
             call mult_edge_by_1d_coeff(beta(n,:),div_coeff_1d(n,:), &
                                        div_coeff_1d_edge(n,:),.true.)
          endif
       end do
    else if (use_div_coeff_cart_edge) then
       do n=1,nlevs
          do comp=1,dm
             call multifab_mult_mult(beta(n,comp),div_coeff_cart_edge(n,comp))
             if (use_alt_energy_fix) then
                call multifab_mult_mult(beta(n,comp),div_coeff_cart_edge(n,comp))                
             endif
          end do
       end do
    end if

    do n = 1,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    abs_solver_eps = -1.0_dp_t
!   if (present(umac_norm)) then
!      do n = 1,nlevs
!         abs_eps = max(abs_eps, umac_norm(n) / dx(n,1))
!      end do
!      abs_solver_eps = eps * abs_eps
!   end if

    rel_solver_eps = min(eps_mac*mac_level_factor**(nlevs-1), eps_mac_max)

    if (use_hypre) then
       call mac_hypre(mla,rh,phi,fine_flx,alpha,beta,dx,&
                      the_bc_tower,press_comp,stencil_order,mla%mba%rr,&
                      rel_solver_eps,abs_solver_eps)
    else
       call mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx,&
                          the_bc_tower,press_comp,stencil_order,&
                          rel_solver_eps,abs_solver_eps)
    endif

    call mkumac(mla,umac,phi,beta,fine_flx,dx,the_bc_tower)

    ! divide out the beta_0
    if (use_div_coeff_1d) then
       do n = 1,nlevs
          call mult_edge_by_1d_coeff(umac(n,:),div_coeff_1d(n,:),div_coeff_1d_edge(n,:), &
                                     .false.)
       end do
    else if (use_div_coeff_cart_edge) then
       do n=1,nlevs
          do comp=1,dm
             call multifab_div_div(umac(n,comp),div_coeff_cart_edge(n,comp))
          end do
       end do
    end if

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
      do i=1,dm
          call multifab_fill_boundary(umac(nlevs,i))
       enddo

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc_edgevel(umac(nlevs,:),the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
          do i=1,dm
             call ml_edge_restriction(umac(n-1,i),umac(n,i),mla%mba%rr(n-1,:),i)
          enddo
       end do

       do n=2,nlevs
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc_edgevel are called.
          call create_umac_grown(umac(n,:),umac(n-1,:), &
                                 the_bc_tower%bc_tower_array(n-1), &
                                 the_bc_tower%bc_tower_array(n), &
                                 n.eq.nlevs)
       end do

    end if
    
    ! Print the norm of each component separately -- make sure to do this after we divide
    !       the velocities by the div_coeff.                                 
    if (verbose .eq. 1) then
       do n = 1,nlevs
          umin = multifab_max(umac(n,1))
          umax = multifab_min(umac(n,1))
          if (dm.eq.1) then
             if (parallel_IOProcessor()) then
                 write(6,1001) n, umax
                 write(6,1101) n, umin
             end if
          else if (dm.eq.2) then
             vmin = multifab_max(umac(n,2))
             vmax = multifab_min(umac(n,2))
             if (parallel_IOProcessor()) then
                 write(6,1002) n, umax, vmax
                 write(6,1102) n, umin, vmin
             end if
          else if (dm.eq.3) then
             vmin = multifab_max(umac(n,2))
             vmax = multifab_min(umac(n,2))
             wmin = multifab_max(umac(n,3))
             wmax = multifab_min(umac(n,3))
             if (parallel_IOProcessor()) then
                 write(6,1003) n, umax, vmax, wmax
                 write(6,1103) n, umin, vmin, wmin
             end if
          end if
       end do
       if (parallel_IOProcessor()) print *,''
    end if

    do n = 1, nlevs
       call destroy(rh(n))
       call destroy(alpha(n))
       do i = 1,dm
          call destroy(beta(n,i))
       end do
    end do

    do n = 1,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

    call destroy(bpt)

1001 format('... level / max : MAC vels   ',i2,2x,e17.10)
1101 format('... level / min : MAC vels   ',i2,2x,e17.10)
1002 format('... level / max : MAC vels   ',i2,2x,e17.10,2x,e17.10)
1102 format('... level / min : MAC vels   ',i2,2x,e17.10,2x,e17.10)
1003 format('... level / max : MAC vels   ',i2,2x,e17.10,2x,e17.10,2x,e17.10)
1103 format('... level / min : MAC vels   ',i2,2x,e17.10,2x,e17.10,2x,e17.10)

  contains

    subroutine divumac(umac,rh,dx,ref_ratio,before,divu_rhs)

      use ml_cc_restriction_module, only: ml_cc_restriction, ml_edge_restriction
      use probin_module, only: verbose

      type(multifab) , intent(inout) :: umac(:,:)
      type(multifab) , intent(inout) :: rh(:)
      real(kind=dp_t), intent(in   ) :: dx(:,:)
      integer        , intent(in   ) :: ref_ratio(:,:)
      logical        , intent(in   ) :: before
      type(multifab ), intent(in   ), optional :: divu_rhs(:)

      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: vmp(:,:,:,:) 
      real(kind=dp_t), pointer :: wmp(:,:,:,:) 
      real(kind=dp_t), pointer :: rhp(:,:,:,:) 
      real(kind=dp_t)          :: rhmax
      integer :: i,lo(get_dim(rh(1))),hi(get_dim(rh(1)))
      integer :: ng_um,ng_rh,dm

      dm = get_dim(rh(1))

      ng_um = nghost(umac(1,1))
      ng_rh = nghost(rh(1))

      do n = nlevs,2,-1
         do i = 1,dm
            call ml_edge_restriction(umac(n-1,i),umac(n,i),ref_ratio(n-1,:),i)
         end do
      end do

      do n = 1,nlevs
         do i = 1, nfabs(rh(n))
            ump => dataptr(umac(n,1), i)
            rhp => dataptr(rh(n)    , i)
            lo =   lwb(get_box(rh(n), i))
            hi =   upb(get_box(rh(n), i))
            select case (dm)
            case (1)
               call divumac_1d(ump(:,1,1,1),ng_um,rhp(:,1,1,1),ng_rh,dx(n,:),lo,hi)
            case (2)
               vmp => dataptr(umac(n,2), i)
               call divumac_2d(ump(:,:,1,1),vmp(:,:,1,1),ng_um,rhp(:,:,1,1),ng_rh, &
                               dx(n,:),lo,hi)
            case (3)
               vmp => dataptr(umac(n,2), i)
               wmp => dataptr(umac(n,3), i)
               call divumac_3d(ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),ng_um, &
                               rhp(:,:,:,1),ng_rh,dx(n,:),lo,hi)
            end select
         end do
      end do

      !     NOTE: the sign convention is because the elliptic solver solves
      !            (alpha MINUS del dot beta grad) phi = RHS
      !            Here alpha is zero.

      ! Do rh = divu_rhs - rh
      if (present(divu_rhs)) then
         do n = 1, nlevs
            call multifab_sub_sub(rh(n),divu_rhs(n))
         end do
      end if
      ! ... or rh = -rh
      do n = 1, nlevs
         call multifab_mult_mult_s(rh(n),-ONE)
      end do

      rhmax = norm_inf(rh(nlevs))
      ! the loop over nlevs must count backwards to make sure the finer grids are done first
      do n = nlevs,2,-1
         ! set level n-1 data to be the average of the level n data covering it
         call ml_cc_restriction(rh(n-1),rh(n),ref_ratio(n-1,:))
         rhmax = max(rhmax,norm_inf(rh(n-1)))
      end do

      if (parallel_IOProcessor() .and. verbose .ge. 1) then
         if (before) then 
            write(6,1000) 
            write(6,1001) rhmax
         else
            write(6,1002) rhmax
         end if
      end if

1000  format(' ')
1001  format('... before mac_projection: max of [div (coeff * UMAC) - RHS)]',e15.8)
1002  format('...  after mac_projection: max of [div (coeff * UMAC) - RHS)]',e15.8)

    end subroutine divumac

    subroutine divumac_1d(umac,ng_um,rh,ng_rh,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_um,ng_rh
      real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_um:)
      real(kind=dp_t), intent(inout) ::   rh(lo(1)-ng_rh:)
      real(kind=dp_t), intent(in   ) ::   dx(:)

      integer :: i

      do i = lo(1),hi(1)
         rh(i) = (umac(i+1) - umac(i)) / dx(1) 
      end do

    end subroutine divumac_1d

    subroutine divumac_2d(umac,vmac,ng_um,rh,ng_rh,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_um,ng_rh
      real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_um:,lo(2)-ng_um:)
      real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_um:,lo(2)-ng_um:)
      real(kind=dp_t), intent(inout) ::   rh(lo(1)-ng_rh:,lo(2)-ng_rh:)
      real(kind=dp_t), intent(in   ) ::   dx(:)

      integer :: i,j

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rh(i,j) = (umac(i+1,j) - umac(i,j)) / dx(1) + &
                      (vmac(i,j+1) - vmac(i,j)) / dx(2)
         end do
      end do

    end subroutine divumac_2d

    subroutine divumac_3d(umac,vmac,wmac,ng_um,rh,ng_rh,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_um,ng_rh
      real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(in   ) :: wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(inout) ::   rh(lo(1)-ng_rh:,lo(2)-ng_rh:,lo(3)-ng_rh:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               rh(i,j,k) = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                           (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                           (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

    end subroutine divumac_3d

    subroutine mult_edge_by_1d_coeff(edge,div_coeff,div_coeff_edge,do_mult)

      type(multifab) , intent(inout) :: edge(:)
      real(dp_t)     , intent(in   ) :: div_coeff(0:)
      real(dp_t)     , intent(in   ) :: div_coeff_edge(0:)
      logical        , intent(in   ) :: do_mult

      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: vmp(:,:,:,:) 
      real(kind=dp_t), pointer :: wmp(:,:,:,:) 
      integer                  :: lo(get_dim(edge(1)))
      integer                  :: hi(get_dim(edge(1)))
      integer                  :: i,ng_um,dm

      dm = get_dim(edge(1))
      ng_um = nghost(edge(1))

      ! Multiply edge velocities by div coeff
      do i = 1, nfabs(edge(1))
         ump => dataptr(edge(1), i)
         lo =  lwb(get_box(edge(1), i))
         hi =  upb(get_box(edge(1), i))
         select case (dm)
         case (1)
            call mult_edge_by_1d_coeff_1d(ump(:,1,1,1), lo, hi, ng_um, &
                                          div_coeff_edge(lo(dm):), &
                                          do_mult)
         case (2)
            vmp => dataptr(edge(2), i)
            call mult_edge_by_1d_coeff_2d(ump(:,:,1,1), vmp(:,:,1,1), lo, hi, ng_um, &
                                          div_coeff(lo(dm):), div_coeff_edge(lo(dm):), &
                                          do_mult)
         case (3)
            vmp => dataptr(edge(2), i)
            wmp => dataptr(edge(3), i)
            call mult_edge_by_1d_coeff_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), lo, hi, ng_um, &
                                          div_coeff(lo(dm):), div_coeff_edge(lo(dm):), &
                                          do_mult)
         end select
      end do

    end subroutine mult_edge_by_1d_coeff

    subroutine mult_edge_by_1d_coeff_1d(uedge,lo,hi,ng_um,div_coeff_edge,do_mult)

      integer                        :: lo(:),hi(:)
      integer                        :: ng_um
      real(kind=dp_t), intent(inout) :: uedge(lo(1)-ng_um:)
      real(dp_t)     , intent(in   ) :: div_coeff_edge(lo(1):)
      logical        , intent(in   ) :: do_mult

      integer :: i

      if (do_mult) then
         do i = lo(1),hi(1)+1
            uedge(i) = uedge(i) * div_coeff_edge(i)
         end do
      else
         do i = lo(1),hi(1)+1
            uedge(i) = uedge(i) / div_coeff_edge(i)
         end do
      end if

    end subroutine mult_edge_by_1d_coeff_1d

    subroutine mult_edge_by_1d_coeff_2d(uedge,vedge,lo,hi,ng_um,div_coeff,div_coeff_edge,do_mult)

      integer                        :: lo(:),hi(:)
      integer                        :: ng_um
      real(kind=dp_t), intent(inout) :: uedge(lo(1)-ng_um:,lo(2)-ng_um:)
      real(kind=dp_t), intent(inout) :: vedge(lo(1)-ng_um:,lo(2)-ng_um:)
      real(dp_t)     , intent(in   ) :: div_coeff(lo(2):)
      real(dp_t)     , intent(in   ) :: div_coeff_edge(lo(2):)
      logical        , intent(in   ) :: do_mult

      integer :: j

      if (do_mult) then
         do j = lo(2),hi(2)
            uedge(lo(1):hi(1)+1,j) = uedge(lo(1):hi(1)+1,j) * div_coeff(j)
         end do
         do j = lo(2),hi(2)+1
            vedge(lo(1):hi(1),j) = vedge(lo(1):hi(1),j) * div_coeff_edge(j)
         end do
      else
         do j = lo(2),hi(2)
            uedge(lo(1):hi(1)+1,j) = uedge(lo(1):hi(1)+1,j) / div_coeff(j)
         end do
         do j = lo(2),hi(2)+1
            vedge(lo(1):hi(1),j) = vedge(lo(1):hi(1),j) / div_coeff_edge(j)
         end do
      end if

    end subroutine mult_edge_by_1d_coeff_2d

    subroutine mult_edge_by_1d_coeff_3d(uedge,vedge,wedge,lo,hi,ng_um,div_coeff,div_coeff_edge, &
                                        do_mult)

      integer                        :: lo(:),hi(:)
      integer                        :: ng_um
      real(kind=dp_t), intent(inout) :: uedge(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(inout) :: vedge(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(inout) :: wedge(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(dp_t)     , intent(in   ) :: div_coeff(lo(3):)
      real(dp_t)     , intent(in   ) :: div_coeff_edge(lo(3):)
      logical        , intent(in   ) :: do_mult

      integer :: k

      if (do_mult) then
         do k = lo(3),hi(3)
            uedge(lo(1):hi(1)+1,lo(2):hi(2),k) = uedge(lo(1):hi(1)+1,lo(2):hi(2),k) * div_coeff(k)
            vedge(lo(1):hi(1),lo(2):hi(2)+1,k) = vedge(lo(1):hi(1),lo(2):hi(2)+1,k) * div_coeff(k)
         end do
         do k = lo(3),hi(3)+1
            wedge(lo(1):hi(1),lo(2):hi(2),k) = wedge(lo(1):hi(1),lo(2):hi(2),k) * div_coeff_edge(k)
         end do
      else
         do k = lo(3),hi(3)
            uedge(lo(1):hi(1)+1,lo(2):hi(2),k) = uedge(lo(1):hi(1)+1,lo(2):hi(2),k) / div_coeff(k)
            vedge(lo(1):hi(1),lo(2):hi(2)+1,k) = vedge(lo(1):hi(1),lo(2):hi(2)+1,k) / div_coeff(k)
         end do
         do k = lo(3),hi(3)+1
            wedge(lo(1):hi(1),lo(2):hi(2),k) = wedge(lo(1):hi(1),lo(2):hi(2),k) / div_coeff_edge(k)
         end do
      end if

    end subroutine mult_edge_by_1d_coeff_3d

    subroutine mk_mac_coeffs(mla,rho,beta)

      use ml_cc_restriction_module, only: ml_edge_restriction

      type(ml_layout), intent(in   ) :: mla
      type(multifab ), intent(in   ) :: rho(:)
      type(multifab ), intent(inout) :: beta(:,:)

      real(kind=dp_t), pointer :: bxp(:,:,:,:) 
      real(kind=dp_t), pointer :: byp(:,:,:,:) 
      real(kind=dp_t), pointer :: bzp(:,:,:,:) 
      real(kind=dp_t), pointer :: rp(:,:,:,:) 
      integer :: i,ng_r,ng_b,lo(mla%dim),hi(mla%dim),dm

      dm = mla%dim
      ng_r = nghost(rho(1))
      ng_b = nghost(beta(1,1))

      do n = 1, nlevs
         do i = 1, nfabs(rho(n))
            rp => dataptr(rho(n) , i)
            bxp => dataptr(beta(n,1), i)
            lo = lwb(get_box(rho(n), i))
            hi = upb(get_box(rho(n), i))
            select case (dm)
            case (1)
               call mk_mac_coeffs_1d(bxp(:,1,1,1),ng_b, rp(:,1,1,1), &
                                     ng_r,lo,hi)
            case (2)
               byp => dataptr(beta(n,2), i)
               call mk_mac_coeffs_2d(bxp(:,:,1,1),byp(:,:,1,1),ng_b, rp(:,:,1,1), &
                                     ng_r,lo,hi)
            case (3)
               byp => dataptr(beta(n,2), i)
               bzp => dataptr(beta(n,3), i)
               call mk_mac_coeffs_3d(bxp(:,:,:,1),byp(:,:,:,1),bzp(:,:,:,1),&
                                     ng_b,rp(:,:,:,1),ng_r,lo,hi)
            end select
         end do
      end do

      ! Make sure that the fine edges average down onto the coarse edges.
      do n = nlevs,2,-1
         do i = 1,dm
            call ml_edge_restriction(beta(n-1,i),beta(n,i),mla%mba%rr(n-1,:),i)
         end do
      end do

    end subroutine mk_mac_coeffs

    subroutine mk_mac_coeffs_1d(betax,ng_b,rho,ng_r,lo,hi)

      integer :: ng_b,ng_r,lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: betax(lo(1)-ng_b:)
      real(kind=dp_t), intent(inout) ::   rho(lo(1)-ng_r:)

      integer :: i

      do i = lo(1),hi(1)+1
         betax(i) = TWO / (rho(i) + rho(i-1))
      end do

    end subroutine mk_mac_coeffs_1d

    subroutine mk_mac_coeffs_2d(betax,betay,ng_b,rho,ng_r,lo,hi)

      integer        , intent(in   ) :: ng_b,ng_r,lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: betax(lo(1)-ng_b:,lo(2)-ng_b:)
      real(kind=dp_t), intent(inout) :: betay(lo(1)-ng_b:,lo(2)-ng_b:)
      real(kind=dp_t), intent(in   ) ::   rho(lo(1)-ng_r:,lo(2)-ng_r:)

      integer :: i,j

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)+1
            betax(i,j) = TWO / (rho(i,j) + rho(i-1,j))
         end do
      end do

      do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)
            betay(i,j) = TWO / (rho(i,j) + rho(i,j-1))
         end do
      end do

    end subroutine mk_mac_coeffs_2d

    subroutine mk_mac_coeffs_3d(betax,betay,betaz,ng_b,rho,ng_r,lo,hi)

      integer        , intent(in   ) :: ng_b,ng_r,lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(inout) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(inout) :: betaz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(in   ) ::   rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)

      integer :: i,j,k

      !$OMP PARALLEL PRIVATE(i,j,k)
      !$OMP DO
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               betax(i,j,k) = TWO / (rho(i,j,k) + rho(i-1,j,k))
            end do
         end do
      end do
      !$OMP END DO NOWAIT
      !$OMP DO
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               betay(i,j,k) = TWO / (rho(i,j,k) + rho(i,j-1,k))
            end do
         end do
      end do
      !$OMP END DO NOWAIT
      !$OMP DO
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               betaz(i,j,k) = TWO / (rho(i,j,k) + rho(i,j,k-1))
            end do
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end subroutine mk_mac_coeffs_3d

    subroutine mkumac(mla,umac,phi,beta,fine_flx,dx,the_bc_tower)

      use variables, only: press_comp

      type(ml_layout),intent(in   ) :: mla
      type(multifab), intent(inout) :: umac(:,:)
      type(multifab), intent(in   ) ::  phi(:)
      type(multifab), intent(in   ) :: beta(:,:)
      type(bndry_reg),intent(in   ) :: fine_flx(:)
      real(dp_t)    , intent(in   ) :: dx(:,:)
      type(bc_tower), intent(in   ) :: the_bc_tower

      integer :: i,ng_um,ng_p,ng_b,lo(get_dim(phi(1))),hi(get_dim(phi(1))),dm

      type(bc_level)           :: bc
      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: vmp(:,:,:,:) 
      real(kind=dp_t), pointer :: wmp(:,:,:,:) 
      real(kind=dp_t), pointer :: php(:,:,:,:) 
      real(kind=dp_t), pointer :: bxp(:,:,:,:) 
      real(kind=dp_t), pointer :: byp(:,:,:,:) 
      real(kind=dp_t), pointer :: bzp(:,:,:,:) 
      real(kind=dp_t), pointer :: lxp(:,:,:,:) 
      real(kind=dp_t), pointer :: hxp(:,:,:,:) 
      real(kind=dp_t), pointer :: lyp(:,:,:,:) 
      real(kind=dp_t), pointer :: hyp(:,:,:,:) 
      real(kind=dp_t), pointer :: lzp(:,:,:,:) 
      real(kind=dp_t), pointer :: hzp(:,:,:,:) 

      dm = get_dim(phi(1))
      nlevs = size(phi)

      ng_um = nghost(umac(1,1))
      ng_p = nghost(phi(1))
      ng_b = nghost(beta(1,1))

      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, nfabs(phi(n))
            ump => dataptr(umac(n,1), i)
            php => dataptr( phi(n), i)
            bxp => dataptr(beta(n,1), i)
            lo  =  lwb(get_box(phi(n), i))
            hi  =  upb(get_box(phi(n), i))
            lxp => dataptr(fine_flx(n)%bmf(1,0), i)
            hxp => dataptr(fine_flx(n)%bmf(1,1), i)

            select case (dm)
            case (1)
               call mkumac_1d(ump(:,1,1,1), ng_um, & 
                              php(:,1,1,1), ng_p, &
                              bxp(:,1,1,1), ng_b, &
                              lxp(:,1,1,1),hxp(:,1,1,1), &
                              lo,hi,dx(n,:),bc%ell_bc_level_array(i,:,:,press_comp))
            case (2)
               vmp => dataptr(umac(n,2), i)
               byp => dataptr(beta(n,2), i)
               lyp => dataptr(fine_flx(n)%bmf(2,0), i)
               hyp => dataptr(fine_flx(n)%bmf(2,1), i)
               call mkumac_2d(ump(:,:,1,1),vmp(:,:,1,1),  ng_um, &
                              php(:,:,1,1),               ng_p, &
                              bxp(:,:,1,1), byp(:,:,1,1), ng_b,&
                              lxp(:,:,1,1),hxp(:,:,1,1),lyp(:,:,1,1),hyp(:,:,1,1), &
                              lo,hi,dx(n,:),bc%ell_bc_level_array(i,:,:,press_comp))

            case (3)
               vmp => dataptr(umac(n,2), i)
               wmp => dataptr(umac(n,3), i)
               byp => dataptr(beta(n,2), i)
               bzp => dataptr(beta(n,3), i)
               lyp => dataptr(fine_flx(n)%bmf(2,0), i)
               hyp => dataptr(fine_flx(n)%bmf(2,1), i)
               lzp => dataptr(fine_flx(n)%bmf(3,0), i)
               hzp => dataptr(fine_flx(n)%bmf(3,1), i)
               call mkumac_3d(ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1), ng_um,&
                              php(:,:,:,1),                             ng_p, &
                              bxp(:,:,:,1), byp(:,:,:,1), bzp(:,:,:,1), ng_b, &
                              lxp(:,:,:,1),hxp(:,:,:,1),lyp(:,:,:,1),hyp(:,:,:,1), &
                              lzp(:,:,:,1),hzp(:,:,:,1), &
                              lo,hi,dx(n,:),bc%ell_bc_level_array(i,:,:,press_comp))
            end select
         end do

         do i = 1,mla%dim
            call multifab_fill_boundary(umac(n,i))
         enddo

      end do
      
      do n = nlevs,2,-1
         do i = 1,mla%dim
            call ml_edge_restriction(umac(n-1,i),umac(n,i),mla%mba%rr(n-1,:),i)
         end do
      end do

    end subroutine mkumac

    subroutine mkumac_1d(umac,ng_um,phi,ng_p,betax,ng_b,lo_x_flx,hi_x_flx,lo,hi,dx,press_bc)

      integer        , intent(in   ) :: lo(:),hi(:)
      integer        , intent(in   ) :: ng_um,ng_p,ng_b
      real(kind=dp_t), intent(inout) ::  umac(lo(1)-ng_um:)
      real(kind=dp_t), intent(inout) ::   phi(lo(1)-ng_p:)
      real(kind=dp_t), intent(in   ) :: betax(lo(1)-ng_b:)
      real(kind=dp_t), intent(in   ) :: lo_x_flx(:)
      real(kind=dp_t), intent(in   ) :: hi_x_flx(:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      integer :: i

      ! At boundaries of grid
      umac(lo(1)  ) = umac(lo(1)  ) - lo_x_flx(1) * dx(1)
      umac(hi(1)+1) = umac(hi(1)+1) + hi_x_flx(1) * dx(1)

      ! In interior of grid
      do i = lo(1)+1, hi(1)
         umac(i) = umac(i) - betax(i) * (phi(i) - phi(i-1)) / dx(1)
      end do


    end subroutine mkumac_1d

    subroutine mkumac_2d(umac,vmac,ng_um,phi,ng_p,betax,betay,ng_b, &
                         lo_x_flx,hi_x_flx,lo_y_flx,hi_y_flx, &
                         lo,hi,dx,press_bc)

      integer        , intent(in   ) :: lo(:),hi(:)
      integer        , intent(in   ) :: ng_um,ng_p,ng_b
      real(kind=dp_t), intent(inout) ::  umac(lo(1)-ng_um:,lo(2)-ng_um:)
      real(kind=dp_t), intent(inout) ::  vmac(lo(1)-ng_um:,lo(2)-ng_um:)
      real(kind=dp_t), intent(inout) ::   phi(lo(1)-ng_p: ,lo(2)-ng_p:)
      real(kind=dp_t), intent(in   ) :: betax(lo(1)-ng_b: ,lo(2)-ng_b:)
      real(kind=dp_t), intent(in   ) :: betay(lo(1)-ng_b: ,lo(2)-ng_b:)
      real(kind=dp_t), intent(in   ) :: lo_x_flx(:,lo(2):), lo_y_flx(lo(1):,:)
      real(kind=dp_t), intent(in   ) :: hi_x_flx(:,lo(2):), hi_y_flx(lo(1):,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      real(kind=dp_t) :: gphix,gphiy
      integer :: i,j

      do j = lo(2),hi(2)

         umac(lo(1)  ,j) = umac(lo(1)  ,j) - lo_x_flx(1,j) * dx(1)
         umac(hi(1)+1,j) = umac(hi(1)+1,j) + hi_x_flx(1,j) * dx(1)

         do i = lo(1)+1,hi(1)
            gphix = (phi(i,j) - phi(i-1,j)) / dx(1)
            umac(i,j) = umac(i,j) - betax(i,j)*gphix
         end do

      end do

      do i = lo(1),hi(1)

         vmac(i,lo(2)  ) = vmac(i,lo(2)  ) - lo_y_flx(i,1) * dx(2)
         vmac(i,hi(2)+1) = vmac(i,hi(2)+1) + hi_y_flx(i,1) * dx(2)

         do j = lo(2)+1,hi(2)
            gphiy = (phi(i,j) - phi(i,j-1)) / dx(2)
            vmac(i,j) = vmac(i,j) - betay(i,j)*gphiy
         end do

      end do

    end subroutine mkumac_2d

    subroutine mkumac_3d(umac,vmac,wmac,   ng_um,&
                         phi,              ng_p, &
                         betax,betay,betaz,ng_b, &
                         lo_x_flx,hi_x_flx,lo_y_flx,hi_y_flx,lo_z_flx,hi_z_flx,&
                         lo,hi,dx,press_bc)

      integer        , intent(in   ) :: lo(:),hi(:)
      integer        , intent(in   ) :: ng_um,ng_p,ng_b
      real(kind=dp_t), intent(inout) :: umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(inout) :: vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(inout) :: wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(inout) ::  phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(in   ) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(in   ) :: betaz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(in   ) :: lo_x_flx(:,lo(2):,lo(3):), hi_x_flx(:,lo(2):,lo(3):)
      real(kind=dp_t), intent(in   ) :: lo_y_flx(lo(1):,:,lo(3):), hi_y_flx(lo(1):,:,lo(3):)
      real(kind=dp_t), intent(in   ) :: lo_z_flx(lo(1):,lo(2):,:), hi_z_flx(lo(1):,lo(2):,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      real(kind=dp_t) :: gphix,gphiy,gphiz
      integer :: i,j,k

      ! NOTE THAT THIS ASSUMES THAT LO, HI ARE THE GRID LIMITS, NOT TILE LIMITS!

      !$OMP PARALLEL PRIVATE(i,j,k,gphix,gphiy,gphiz)
      !$OMP DO
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            umac(lo(1)  ,j,k) = umac(lo(1)  ,j,k) - lo_x_flx(1,j,k) * dx(1)
            umac(hi(1)+1,j,k) = umac(hi(1)+1,j,k) + hi_x_flx(1,j,k) * dx(1)
            do i = lo(1)+1,hi(1)
               gphix = (phi(i,j,k) - phi(i-1,j,k)) / dx(1)
               umac(i,j,k) = umac(i,j,k) - betax(i,j,k)*gphix
            end do
         end do
      end do
      !$OMP END DO NOWAIT

      !$OMP DO
      do k = lo(3),hi(3)
         do i = lo(1),hi(1)
            vmac(i,lo(2)  ,k) = vmac(i,lo(2)  ,k) - lo_y_flx(i,1,k) * dx(2)
            vmac(i,hi(2)+1,k) = vmac(i,hi(2)+1,k) + hi_y_flx(i,1,k) * dx(2)
            do j = lo(2)+1,hi(2)
               gphiy = (phi(i,j,k) - phi(i,j-1,k)) / dx(2)
               vmac(i,j,k) = vmac(i,j,k) - betay(i,j,k)*gphiy
            end do
         end do
      end do
      !$OMP END DO NOWAIT

      !$OMP DO
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            wmac(i,j,lo(3)  ) = wmac(i,j,lo(3)  ) - lo_z_flx(i,j,1) * dx(3)
            wmac(i,j,hi(3)+1) = wmac(i,j,hi(3)+1) + hi_z_flx(i,j,1) * dx(3)
            do k = lo(3)+1,hi(3)
               gphiz = (phi(i,j,k) - phi(i,j,k-1)) / dx(3)
               wmac(i,j,k) = wmac(i,j,k) - betaz(i,j,k)*gphiz
            end do
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end subroutine mkumac_3d
  end subroutine macproject

end module macproject_module
