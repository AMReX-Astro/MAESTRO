module macproject_module

  use bl_types
  use bc_module
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use bl_constants_module
  use create_umac_grown_module
  use impose_phys_bcs_on_edges_module

  implicit none

  private

  public :: macproject

contains 

  ! NOTE: this routine differs from that in varden because phi is passed in/out 
  !       rather than allocated here
  subroutine macproject(mla,umac,phi,rho,dx,the_bc_tower, &
                        bc_comp,divu_rhs,div_coeff_1d,div_coeff_half_1d,div_coeff_3d)

    use mac_multigrid_module
    use geometry, only: dm, nlevs, spherical
    use probin_module, only: verbose, edge_nodal_flag

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: umac(:,:)
    type(multifab ), intent(inout) :: phi(:)
    type(multifab ), intent(in   ) :: rho(:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: bc_comp

    type(multifab ), intent(in   ), optional :: divu_rhs(:)
    real(dp_t)     , intent(in   ), optional :: div_coeff_1d(:,:)
    real(dp_t)     , intent(in   ), optional :: div_coeff_half_1d(:,:)
    type(multifab ), intent(in   ), optional :: div_coeff_3d(:)

    type(multifab)  :: rh(mla%nlevel),alpha(mla%nlevel),beta(mla%nlevel,dm)
    type(bndry_reg) :: fine_flx(2:mla%nlevel)

    real(dp_t)                   :: umac_norm(mla%nlevel)
    real(dp_t)                   :: unorm,vnorm,wnorm
    integer                      :: stencil_order,i,n
    logical                      :: use_rhs, use_div_coeff_1d, use_div_coeff_3d

    type(bl_prof_timer), save :: bpt

    call build(bpt, "macproject")

    use_rhs          = .false. ; if (present(divu_rhs)    ) use_rhs          = .true.
    use_div_coeff_1d = .false. ; if (present(div_coeff_1d)) use_div_coeff_1d = .true.
    use_div_coeff_3d = .false. ; if (present(div_coeff_3d)) use_div_coeff_3d = .true.

    if (use_div_coeff_1d .and. use_div_coeff_3d) then
       call bl_error('CANT HAVE 1D and 3D DIV_COEFF IN MACPROJECT')
    end if

    if (spherical .eq. 1 .and. use_div_coeff_1d) then
       call bl_error('CANT HAVE SPHERICAL .eq. 1 and 1D DIV_COEFF IN MACPROJECT')
    end if

    if (spherical .eq. 0 .and. use_div_coeff_3d) then
       call bl_error('CANT HAVE SPHERICAL .eq. 0 and 3D DIV_COEFF IN MACPROJECT')
    end if

    stencil_order = 2

    do n = 1, nlevs
       call multifab_build(   rh(n), mla%la(n),  1, 0)
       call multifab_build(alpha(n), mla%la(n),  1, 1)
       do i = 1,dm
          call multifab_build(beta(n,i),mla%la(n),1,1,nodal=edge_nodal_flag(i,:))
       end do
       call setval(alpha(n),ZERO,all=.true.)
    end do

    if (use_div_coeff_1d) then
       do n = 1,nlevs
          call mult_edge_by_1d_coeff(umac(n,:),div_coeff_1d(n,:),div_coeff_half_1d(n,:),&
                                     .true.)
       end do
    else if (use_div_coeff_3d) then
       do n = 1,nlevs
          call mult_edge_by_3d_coeff(umac(n,:),div_coeff_3d(n),ml_layout_get_pd(mla,n),&
                                     .true.)
       end do
    end if

    ! Compute umac_norm to be used inside the MG solver as part of a stopping criterion
    umac_norm = -1.0_dp_t
    do n = 1,nlevs
       do i = 1,dm
          umac_norm(n) = max(umac_norm(n),norm_inf(umac(n,i)))
       end do
    end do

    if (use_rhs) then
       call divumac(umac,rh,dx,mla%mba%rr,.true.,divu_rhs)
    else
       call divumac(umac,rh,dx,mla%mba%rr,.true.)
    end if

    ! Print the norm of each component separately
    if (verbose .eq. 1) then
       do n = 1,nlevs
          unorm = norm_inf(umac(n,1))
          vnorm = norm_inf(umac(n,2))
          if (dm.eq.3) wnorm = norm_inf(umac(n,3))
          if (parallel_IOProcessor()) then
            print *,'MAX OF UMAC AT LEVEL ',n,unorm
            print *,'MAX OF VMAC AT LEVEL ',n,vnorm
            if (dm.eq.3) print *,'MAX OF WMAC AT LEVEL ',n,wnorm
          end if
       end do
    end if

    call mk_mac_coeffs(mla,rho,beta,the_bc_tower)

    if (use_div_coeff_1d) then
       do n = 1,nlevs
          call mult_edge_by_1d_coeff(beta(n,:),div_coeff_1d(n,:),div_coeff_half_1d(n,:),&
                                     .true.)
       end do
    else if (use_div_coeff_3d) then
       do n = 1,nlevs
          call mult_edge_by_3d_coeff(beta(n,:),div_coeff_3d(n),ml_layout_get_pd(mla,n),&
                                     .true.)
       end do
    end if

    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    call mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx,&
                       the_bc_tower,bc_comp,stencil_order,mla%mba%rr,umac_norm)

    call mkumac(rh,umac,phi,beta,fine_flx,dx,the_bc_tower,bc_comp,mla%mba%rr)

    if (use_rhs) then
       call divumac(umac,rh,dx,mla%mba%rr,.false.,divu_rhs)
    else
       call divumac(umac,rh,dx,mla%mba%rr,.false.)
    end if

    ! Print the norm of each component separately
    if (verbose .eq. 1) then
       do n = 1,nlevs
          unorm = norm_inf(umac(n,1))
          vnorm = norm_inf(umac(n,2))
          if (dm.eq.3) wnorm = norm_inf(umac(n,3))
          if (parallel_IOProcessor()) then
            print *,'MAX OF UMAC AT LEVEL ',n,unorm
            print *,'MAX OF VMAC AT LEVEL ',n,vnorm
            if (dm.eq.3) print *,'MAX OF WMAC AT LEVEL ',n,wnorm
          end if
       end do
       if (parallel_IOProcessor()) print *,''
    end if

    if (use_div_coeff_1d) then
       do n = 1,nlevs
          call mult_edge_by_1d_coeff(umac(n,:),div_coeff_1d(n,:),div_coeff_half_1d(n,:), &
                                     .false.)
       end do
    else if (use_div_coeff_3d) then
       do n = 1,nlevs
          call mult_edge_by_3d_coeff(umac(n,:),div_coeff_3d(n),ml_layout_get_pd(mla,n), &
                                     .false.)
       end do
    end if

    if (nlevs .gt. 1) then
       do n=2,nlevs
          call create_umac_grown(n,umac(n,:),umac(n-1,:))
       end do
    else
       do n=1,nlevs
          do i=1,dm
             call multifab_fill_boundary(umac(n,i))
          end do
       end do
    end if

    ! This fills the same edges that create_umac_grown does but fills them from 
    !  physical boundary conditions rather than from coarser grids
    if (dm > 1) &
       call impose_phys_bcs_on_edges(rho,umac,the_bc_tower%bc_tower_array)
    
    do n = 1, nlevs
       call destroy(rh(n))
       call destroy(alpha(n))
       do i = 1,dm
          call destroy(beta(n,i))
       end do
    end do

    do n = 2,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

    call destroy(bpt)

  contains

    subroutine divumac(umac,rh,dx,ref_ratio,before,divu_rhs)

      use ml_restriction_module, only: ml_cc_restriction, ml_edge_restriction
      use probin_module, only: verbose
      use geometry, only: dm

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
      integer :: i,lo(dm),hi(dm)
      integer :: ng_um, ng_rh

      ng_um = umac(1,1)%ng
      ng_rh = rh(1)%ng

      do n = nlevs,2,-1
         do i = 1,dm
            call ml_edge_restriction(umac(n-1,i),umac(n,i),ref_ratio(n-1,:),i)
         end do
      end do

      do n = 1,nlevs
         do i = 1, rh(n)%nboxes
            if ( multifab_remote(rh(n), i) ) cycle
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

!$omp parallel do private(i,j,k)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               rh(i,j,k) = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                           (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                           (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)
            end do
         end do
      end do
!$omp end parallel do

    end subroutine divumac_3d

    subroutine mult_edge_by_1d_coeff(edge,div_coeff,div_coeff_half,do_mult)

      use geometry, only: dm

      type(multifab) , intent(inout) :: edge(:)
      real(dp_t)     , intent(in   ) :: div_coeff(0:)
      real(dp_t)     , intent(in   ) :: div_coeff_half(0:)
      logical        , intent(in   ) :: do_mult

      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: vmp(:,:,:,:) 
      real(kind=dp_t), pointer :: wmp(:,:,:,:) 
      integer                  :: lo(dm)
      integer                  :: i,ng_um

      ng_um = edge(1)%ng

      ! Multiply edge velocities by div coeff
      do i = 1, edge(1)%nboxes
         if ( multifab_remote(edge(1), i) ) cycle
         ump => dataptr(edge(1), i)
         lo =  lwb(get_box(edge(1), i))
         select case (dm)
         case (1)
            call mult_by_1d_coeff_1d(ump(:,1,1,1), ng_um, &
                                     div_coeff(lo(dm):), div_coeff_half(lo(dm):), do_mult)
         case (2)
            vmp => dataptr(edge(2), i)
            call mult_by_1d_coeff_2d(ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                     div_coeff(lo(dm):), div_coeff_half(lo(dm):), do_mult)
         case (3)
            vmp => dataptr(edge(2), i)
            wmp => dataptr(edge(3), i)
            call mult_by_1d_coeff_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                     div_coeff(lo(dm):), div_coeff_half(lo(dm):), do_mult)
         end select
      end do

    end subroutine mult_edge_by_1d_coeff

    subroutine mult_by_1d_coeff_1d(uedge,ng_um,div_coeff,div_coeff_half,do_mult)

      integer                        :: ng_um
      real(kind=dp_t), intent(inout) :: uedge(-ng_um:)
      real(dp_t)     , intent(in   ) :: div_coeff(0:)
      real(dp_t)     , intent(in   ) :: div_coeff_half(0:)
      logical        , intent(in   ) :: do_mult

      integer :: i,nx

      nx = size(uedge,dim=1)-2*ng_um-1

      if (do_mult) then
         do i = 0,nx
            uedge(i) = uedge(i) * div_coeff_half(i)
         end do
      else
         do i = 0,nx
            uedge(i) = uedge(i) / div_coeff_half(i)
         end do
      end if

    end subroutine mult_by_1d_coeff_1d

    subroutine mult_by_1d_coeff_2d(uedge,vedge,ng_um,div_coeff,div_coeff_half,do_mult)

      integer                        :: ng_Um
      real(kind=dp_t), intent(inout) :: uedge(-ng_um:,-ng_um:)
      real(kind=dp_t), intent(inout) :: vedge(-ng_um:,-ng_um:)
      real(dp_t)     , intent(in   ) :: div_coeff(0:)
      real(dp_t)     , intent(in   ) :: div_coeff_half(0:)
      logical        , intent(in   ) :: do_mult

      integer :: j,ny

      ny = size(uedge,dim=2)-2*ng_um

      if (do_mult) then
         do j = 0,ny-1
            uedge(:,j) = uedge(:,j) * div_coeff(j)
         end do
         do j = 0,ny
            vedge(:,j) = vedge(:,j) * div_coeff_half(j)
         end do
      else
         do j = 0,ny-1 
            uedge(:,j) = uedge(:,j) / div_coeff(j)
         end do
         do j = 0,ny
            vedge(:,j) = vedge(:,j) / div_coeff_half(j)
         end do
      end if

    end subroutine mult_by_1d_coeff_2d

    subroutine mult_by_1d_coeff_3d(uedge,vedge,wedge,ng_um,div_coeff,div_coeff_half,do_mult)

      integer                        :: ng_um
      real(kind=dp_t), intent(inout) :: uedge(-ng_um:,-ng_um:,-ng_um:)
      real(kind=dp_t), intent(inout) :: vedge(-ng_um:,-ng_um:,-ng_um:)
      real(kind=dp_t), intent(inout) :: wedge(-ng_um:,-ng_um:,-ng_um:)
      real(dp_t)     , intent(in   ) :: div_coeff(0:)
      real(dp_t)     , intent(in   ) :: div_coeff_half(0:)
      logical        , intent(in   ) :: do_mult

      integer :: k,nz

      nz = size(uedge,dim=3)-2*ng_um

      if (do_mult) then
!$omp parallel do private(k)
         do k = 0,nz-1 
            uedge(:,:,k) = uedge(:,:,k) * div_coeff(k)
            vedge(:,:,k) = vedge(:,:,k) * div_coeff(k)
         end do
!$omp end parallel do
!$omp parallel do private(k)
         do k = 0,nz
            wedge(:,:,k) = wedge(:,:,k) * div_coeff_half(k)
         end do
!$omp end parallel do
      else
!$omp parallel do private(k)
         do k = 0,nz-1 
            uedge(:,:,k) = uedge(:,:,k) / div_coeff(k)
            vedge(:,:,k) = vedge(:,:,k) / div_coeff(k)
         end do
!$omp end parallel do
!$omp parallel do private(k)
         do k = 0,nz
            wedge(:,:,k) = wedge(:,:,k) / div_coeff_half(k)
         end do
!$omp end parallel do
      end if

    end subroutine mult_by_1d_coeff_3d

    subroutine mult_edge_by_3d_coeff(edge,div_coeff,domain,do_mult)

      use geometry, only: dm

      type(multifab) , intent(inout) :: edge(:)
      type(multifab) , intent(in   ) :: div_coeff
      type(box)      , intent(in   ) :: domain
      logical        , intent(in   ) :: do_mult

      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: vmp(:,:,:,:) 
      real(kind=dp_t), pointer :: wmp(:,:,:,:) 
      real(kind=dp_t), pointer ::  dp(:,:,:,:) 
      integer :: i,lo(dm),hi(dm)
      integer :: domlo(dm),domhi(dm)

      domlo =  lwb(domain)
      domhi =  upb(domain)

      ! Multiply edge velocities by div coeff
      do i = 1, edge(1)%nboxes
         if ( multifab_remote(edge(1), i) ) cycle
         ump => dataptr(edge(1), i)
         vmp => dataptr(edge(2), i)
         wmp => dataptr(edge(3), i)
         dp => dataptr(div_coeff, i)
         lo =  lwb(get_box(edge(1), i))
         hi =  upb(get_box(edge(1), i))
         select case (dm)
         case (3)
            call mult_by_3d_coeff_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                     dp(:,:,:,1), lo, hi, domlo, domhi, do_mult)
         end select
      end do

    end subroutine mult_edge_by_3d_coeff

    subroutine mult_by_3d_coeff_3d(uedge,vedge,wedge,div_coeff,lo,hi,domlo,domhi,do_mult)

      integer        , intent(in   ) :: lo(:),hi(:),domlo(:),domhi(:)
      real(kind=dp_t), intent(inout) ::      uedge(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(kind=dp_t), intent(inout) ::      vedge(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(kind=dp_t), intent(inout) ::      wedge(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(dp_t)     , intent(in   ) :: div_coeff(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      logical        , intent(in   ) :: do_mult

      integer :: i,j,k

      if (do_mult) then

!$omp parallel private(i,j,k)
!$omp do
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1)+1,hi(1)
                  uedge(i,j,k) = uedge(i,j,k) * HALF * (div_coeff(i,j,k)+div_coeff(i-1,j,k))
               end do
               if (lo(1).eq.domlo(1)) then
                  uedge(lo(1),j,k) = uedge(lo(1),j,k) * div_coeff(lo(1),j,k)
               else
                  uedge(lo(1),j,k) = uedge(lo(1),j,k) * HALF * &
                       (div_coeff(lo(1),j,k)+div_coeff(lo(1)-1,j,k))
               end if
               if (hi(1).eq.domhi(1)) then
                  uedge(hi(1)+1,j,k) = uedge(hi(1)+1,j,k) * div_coeff(hi(1),j,k)
               else
                  uedge(hi(1)+1,j,k) = uedge(hi(1)+1,j,k) * HALF * &
                       (div_coeff(hi(1)+1,j,k)+div_coeff(hi(1),j,k))
               end if
            end do
         end do
!$omp end do nowait
!$omp do
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               do j = lo(2)+1,hi(2)
                  vedge(i,j,k) = vedge(i,j,k) * HALF * (div_coeff(i,j,k)+div_coeff(i,j-1,k))
               end do
               if (lo(2).eq.domlo(2)) then
                  vedge(i,lo(2),k) = vedge(i,lo(2),k) * div_coeff(i,lo(2),k)
               else
                  vedge(i,lo(2),k) = vedge(i,lo(2),k) * HALF * &
                       (div_coeff(i,lo(2),k)+div_coeff(i,lo(2)-1,k))
               end if
               if (hi(2).eq.domhi(2)) then
                  vedge(i,hi(2)+1,k) = vedge(i,hi(2)+1,k) * div_coeff(i,hi(2),k)
               else
                  vedge(i,hi(2)+1,k) = vedge(i,hi(2)+1,k) * HALF * &
                       (div_coeff(i,hi(2)+1,k)+div_coeff(i,hi(2),k))
               end if
            end do
         end do
!$omp end do nowait
!$omp do
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               do k = lo(3)+1,hi(3)
                  wedge(i,j,k) = wedge(i,j,k) * HALF * (div_coeff(i,j,k)+div_coeff(i,j,k-1))
               end do
               if (lo(3).eq.domlo(3)) then
                  wedge(i,j,lo(3)) = wedge(i,j,lo(3)) * div_coeff(i,j,lo(3))
               else
                  wedge(i,j,lo(3)) = wedge(i,j,lo(3)) * HALF * &
                       (div_coeff(i,j,lo(3))+div_coeff(i,j,lo(3)-1))
               end if
               if (hi(3).eq.domhi(3)) then
                  wedge(i,j,hi(3)+1) = wedge(i,j,hi(3)+1) * div_coeff(i,j,hi(3))
               else
                  wedge(i,j,hi(3)+1) = wedge(i,j,hi(3)+1) * HALF * &
                       (div_coeff(i,j,hi(3)+1)+div_coeff(i,j,hi(3)))
               end if
            end do
         end do
!$omp end do
!$omp end parallel

      else

!$omp parallel private(i,j,k)
!$omp do
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1)+1,hi(1)
                  uedge(i,j,k) = uedge(i,j,k) / ( HALF*(div_coeff(i,j,k)+div_coeff(i-1,j,k)))
               end do
               if (lo(1).eq.domlo(1)) then
                  uedge(lo(1),j,k) = uedge(lo(1),j,k) / div_coeff(lo(1),j,k)
               else
                  uedge(lo(1),j,k) = uedge(lo(1),j,k) / &
                       ( HALF * (div_coeff(lo(1),j,k)+div_coeff(lo(1)-1,j,k)))
               end if
               if (hi(1).eq.domhi(1)) then
                  uedge(hi(1)+1,j,k) = uedge(hi(1)+1,j,k) / div_coeff(hi(1),j,k)
               else
                  uedge(hi(1)+1,j,k) = uedge(hi(1)+1,j,k) / &
                       ( HALF * (div_coeff(hi(1)+1,j,k)+div_coeff(hi(1),j,k)))
               end if
            end do
         end do
!$omp end do nowait
!$omp do
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               do j = lo(2)+1,hi(2)
                  vedge(i,j,k) = vedge(i,j,k) / ( HALF*(div_coeff(i,j,k)+div_coeff(i,j-1,k)))
               end do
               if (lo(2).eq.domlo(2)) then
                  vedge(i,lo(2),k) = vedge(i,lo(2),k) / div_coeff(i,lo(2),k)
               else
                  vedge(i,lo(2),k) = vedge(i,lo(2),k) / &
                       ( HALF * (div_coeff(i,lo(2),k)+div_coeff(i,lo(2)-1,k)))
               end if
               if (hi(2).eq.domhi(2)) then
                  vedge(i,hi(2)+1,k) = vedge(i,hi(2)+1,k) / div_coeff(i,hi(2),k)
               else
                  vedge(i,hi(2)+1,k) = vedge(i,hi(2)+1,k) / &
                       ( HALF * (div_coeff(i,hi(2)+1,k)+div_coeff(i,hi(2),k)))
               end if
            end do
         end do
!$omp end do nowait
!$omp do
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               do k = lo(3)+1,hi(3)
                  wedge(i,j,k) = wedge(i,j,k) / ( HALF*(div_coeff(i,j,k)+div_coeff(i,j,k-1)))
               end do
               if (lo(3).eq.domlo(3)) then
                  wedge(i,j,lo(3)) = wedge(i,j,lo(3)) / div_coeff(i,j,lo(3))
               else
                  wedge(i,j,lo(3)) = wedge(i,j,lo(3)) / &
                       ( HALF * (div_coeff(i,j,lo(3))+div_coeff(i,j,lo(3)-1)))
               end if
               if (hi(3).eq.domhi(3)) then
                  wedge(i,j,hi(3)+1) = wedge(i,j,hi(3)+1) / div_coeff(i,j,hi(3))
               else
                  wedge(i,j,hi(3)+1) = wedge(i,j,hi(3)+1) / &
                       ( HALF * (div_coeff(i,j,hi(3)+1)+div_coeff(i,j,hi(3))))
               end if
            end do
         end do
!$omp end do
!$omp end parallel

      end if

    end subroutine mult_by_3d_coeff_3d

    subroutine mk_mac_coeffs(mla,rho,beta,the_bc_tower)

      use geometry, only: dm
      use ml_restriction_module, only: ml_edge_restriction

      type(ml_layout), intent(in   ) :: mla
      type(multifab ), intent(in   ) :: rho(:)
      type(multifab ), intent(inout) :: beta(:,:)
      type(bc_tower ), intent(in   ) :: the_bc_tower

      real(kind=dp_t), pointer :: bxp(:,:,:,:) 
      real(kind=dp_t), pointer :: byp(:,:,:,:) 
      real(kind=dp_t), pointer :: bzp(:,:,:,:) 
      real(kind=dp_t), pointer :: rp(:,:,:,:) 
      integer :: i,ng_r,ng_b,lo(dm),hi(dm)

      ng_r = rho(1)%ng
      ng_b = beta(1,1)%ng

      do n = 1, nlevs
         do i = 1, rho(n)%nboxes
            if ( multifab_remote(rho(n), i) ) cycle
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

      integer :: i,j

      do i = lo(1),hi(1)+1
         betax(i) = TWO / (rho(i) + rho(i-1))
      end do

    end subroutine mk_mac_coeffs_1d

    subroutine mk_mac_coeffs_2d(betax,betay,ng_b,rho,ng_r,lo,hi)

      integer :: ng_b,ng_r,lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: betax(lo(1)-ng_b:,lo(2)-ng_b:)
      real(kind=dp_t), intent(inout) :: betay(lo(1)-ng_b:,lo(2)-ng_b:)
      real(kind=dp_t), intent(inout) ::   rho(lo(1)-ng_r:,lo(2)-ng_r:)

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

      integer :: ng_b,ng_r,lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: betax(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(inout) :: betay(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(inout) :: betaz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(inout) ::   rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)

      integer :: i,j,k

!$omp parallel private(i,j,k)
!$omp do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               betax(i,j,k) = TWO / (rho(i,j,k) + rho(i-1,j,k))
            end do
         end do
      end do
!$omp end do nowait
!$omp do
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               betay(i,j,k) = TWO / (rho(i,j,k) + rho(i,j-1,k))
            end do
         end do
      end do
!$omp end do nowait
!$omp do
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               betaz(i,j,k) = TWO / (rho(i,j,k) + rho(i,j,k-1))
            end do
         end do
      end do
!$omp end do
!$omp end parallel

    end subroutine mk_mac_coeffs_3d

    subroutine mkumac(rh,umac,phi,beta,fine_flx,dx,the_bc_tower,press_comp,ref_ratio)

      use ml_restriction_module, only: ml_edge_restriction
      use geometry, only: dm, nlevs

      type(multifab), intent(inout) :: umac(:,:)
      type(multifab), intent(inout) ::   rh(:)
      type(multifab), intent(in   ) ::  phi(:)
      type(multifab), intent(in   ) :: beta(:,:)
      type(bndry_reg),intent(in   ) :: fine_flx(2:)
      real(dp_t)    , intent(in   ) :: dx(:,:)
      type(bc_tower), intent(in   ) :: the_bc_tower
      integer       , intent(in   ) :: press_comp
      integer       , intent(in   ) :: ref_ratio(:,:)

      integer :: i
      integer :: ng_um,ng_p,ng_b
      integer :: lo(dm),hi(dm)

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


      ng_um = umac(1,1)%ng
      ng_p = phi(1)%ng
      ng_b = beta(1,1)%ng

      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, rh(n)%nboxes
            if ( multifab_remote(rh(n), i) ) cycle
            ump => dataptr(umac(n,1), i)
            php => dataptr( phi(n), i)
            bxp => dataptr(beta(n,1), i)
            lo = lwb(get_box(phi(n), i))
            hi = upb(get_box(phi(n), i))
            select case (dm)
            case (1)
               call mkumac_1d(n,ump(:,1,1,1), ng_um, & 
                              php(:,1,1,1), ng_p, &
                              bxp(:,1,1,1), ng_b, &
                              dx(n,:),bc%ell_bc_level_array(i,:,:,press_comp))
               if (n > 1) then
                  lxp => dataptr(fine_flx(n)%bmf(1,0), i)
                  hxp => dataptr(fine_flx(n)%bmf(1,1), i)
                  call correct_mkumac_1d(ump(:,1,1,1),ng_um, &
                                         lxp(:,1,1,1),hxp(:,1,1,1),lo,hi,dx(n,:))
               end if
            case (2)
               vmp => dataptr(umac(n,2), i)
               byp => dataptr(beta(n,2), i)
               call mkumac_2d(n,ump(:,:,1,1),vmp(:,:,1,1), ng_um, & 
                              php(:,:,1,1), ng_p, &
                              bxp(:,:,1,1), byp(:,:,1,1), ng_b, &
                              dx(n,:),bc%ell_bc_level_array(i,:,:,press_comp))
               if (n > 1) then
                  lxp => dataptr(fine_flx(n)%bmf(1,0), i)
                  hxp => dataptr(fine_flx(n)%bmf(1,1), i)
                  lyp => dataptr(fine_flx(n)%bmf(2,0), i)
                  hyp => dataptr(fine_flx(n)%bmf(2,1), i)
                  call correct_mkumac_2d(ump(:,:,1,1),vmp(:,:,1,1),ng_um, &
                                         lxp(:,:,1,1),hxp(:,:,1,1),lyp(:,:,1,1),hyp(:,:,1,1), &
                                         lo,hi,dx(n,:))
               end if
            case (3)
               vmp => dataptr(umac(n,2), i)
               wmp => dataptr(umac(n,3), i)
               byp => dataptr(beta(n,2), i)
               bzp => dataptr(beta(n,3), i)
               call mkumac_3d(n,ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),ng_um,&
                              php(:,:,:,1), ng_p, &
                              bxp(:,:,:,1), byp(:,:,:,1), bzp(:,:,:,1), ng_b, &
                              dx(n,:),bc%ell_bc_level_array(i,:,:,press_comp))
               if (n > 1) then
                  lxp => dataptr(fine_flx(n)%bmf(1,0), i)
                  hxp => dataptr(fine_flx(n)%bmf(1,1), i)
                  lyp => dataptr(fine_flx(n)%bmf(2,0), i)
                  hyp => dataptr(fine_flx(n)%bmf(2,1), i)
                  lzp => dataptr(fine_flx(n)%bmf(3,0), i)
                  hzp => dataptr(fine_flx(n)%bmf(3,1), i)
                  call correct_mkumac_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                         lxp(:,:,:,1),hxp(:,:,:,1),lyp(:,:,:,1),hyp(:,:,:,1), &
                                         lzp(:,:,:,1),hzp(:,:,:,1),lo,hi,dx(n,:))
               end if
            end select
         end do
      end do

      do n = nlevs,2,-1
         do i = 1,dm
            call ml_edge_restriction(umac(n-1,i),umac(n,i),ref_ratio(n-1,:),i)
         end do
      end do

    end subroutine mkumac

    subroutine mkumac_1d(n,umac,ng_um,phi,ng_p,betax,ng_b,dx,press_bc)

      integer        , intent(in   ) :: n,ng_um,ng_p,ng_b
      real(kind=dp_t), intent(inout) :: umac(-ng_um:)
      real(kind=dp_t), intent(inout) ::  phi(-ng_p:)
      real(kind=dp_t), intent(in   ) :: betax(-ng_b:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      real(kind=dp_t) :: gphix
      integer :: i,nx,imin,imax

      nx = size(phi,dim=1) - 2

      if (n.eq.1) then
         imin = 0
         imax = nx
      else
         imin = 1
         imax = nx-1
      end if

      if (press_bc(1,1) == BC_NEU) then
         phi(-1) = phi(0)
      else if (press_bc(1,1) == BC_DIR) then
         phi(-1) = -TWO*phi(0) + THIRD * phi(1)
      end if

      if (press_bc(1,2) == BC_NEU) then
         phi(nx) = phi(nx-1)
      else if (press_bc(1,2) == BC_DIR) then
         phi(nx) = -TWO*phi(nx-1) + THIRD * phi(nx-2)
      end if

      do i = imin,imax
         gphix = (phi(i) - phi(i-1)) / dx(1)
         umac(i) = umac(i) - betax(i)*gphix
      end do

      ! Here we reset phi == 0 at BC_DIR to be used in later iteration if necessary
      if (press_bc(1,1) == BC_DIR) phi(-1) = ZERO
      if (press_bc(1,2) == BC_DIR) phi(nx) = ZERO

    end subroutine mkumac_1d

    subroutine correct_mkumac_1d(umac,ng_um,lo_x_flx,hi_x_flx,lo,hi,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_um
      real(kind=dp_t), intent(inout) :: umac(lo(1)-ng_um:)
      real(kind=dp_t), intent(in   ) :: lo_x_flx(:)
      real(kind=dp_t), intent(in   ) :: hi_x_flx(:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      umac(lo(1)  ) = umac(lo(1)  ) - lo_x_flx(1) * dx(1)
      umac(hi(1)+1) = umac(hi(1)+1) + hi_x_flx(1) * dx(1)

    end subroutine correct_mkumac_1d

    subroutine mkumac_2d(n,umac,vmac,ng_um,phi,ng_p,betax,betay,ng_b,dx,press_bc)

      integer        , intent(in   ) :: n,ng_um,ng_p,ng_b
      real(kind=dp_t), intent(inout) :: umac(-ng_um:,-ng_um:)
      real(kind=dp_t), intent(inout) :: vmac(-ng_um:,-ng_um:)
      real(kind=dp_t), intent(inout) ::  phi(-ng_p:,-ng_p:)
      real(kind=dp_t), intent(in   ) :: betax(-ng_b:,-ng_b:)
      real(kind=dp_t), intent(in   ) :: betay(-ng_b:,-ng_b:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      real(kind=dp_t) :: gphix,gphiy
      integer :: i,j,nx,ny
      integer :: imin,imax,jmin,jmax

      nx = size(phi,dim=1) - 2
      ny = size(phi,dim=2) - 2

      if (n.eq.1) then
         imin = 0
         imax = nx
         jmin = 0
         jmax = ny
      else
         imin = 1
         imax = nx-1
         jmin = 1
         jmax = ny-1
      end if

      if (press_bc(1,1) == BC_NEU) then
         do j = 0,ny-1
            phi(-1,j) = phi(0,j)
         end do
      else if (press_bc(1,1) == BC_DIR) then
         do j = 0,ny-1
            phi(-1,j) = -TWO*phi(0,j) + THIRD * phi(1,j)
         end do
      end if
      if (press_bc(1,2) == BC_NEU) then
         do j = 0,ny-1
            phi(nx,j) = phi(nx-1,j)
         end do
      else if (press_bc(1,2) == BC_DIR) then
         do j = 0,ny-1
            phi(nx,j) = -TWO*phi(nx-1,j) + THIRD * phi(nx-2,j)
         end do
      end if
      if (press_bc(2,1) == BC_NEU) then
         do i = 0,nx-1
            phi(i,-1) = phi(i,0)
         end do
      else if (press_bc(2,1) == BC_DIR) then
         do i = 0,nx-1
            phi(i,-1) = -TWO*phi(i,0) + THIRD * phi(i,1)
         end do
      end if
      if (press_bc(2,2) == BC_NEU) then
         do i = 0,nx-1
            phi(i,ny) = phi(i,ny-1)
         end do
      else if (press_bc(2,2) == BC_DIR) then
         do i = 0,nx-1
            phi(i,ny) = -TWO*phi(i,ny-1) + THIRD * phi(i,ny-2)
         end do
      end if

      do j = 0,ny-1
         do i = imin,imax
            gphix = (phi(i,j) - phi(i-1,j)) / dx(1)
            umac(i,j) = umac(i,j) - betax(i,j)*gphix
         end do
      end do

      do i = 0,nx-1
         do j = jmin,jmax
            gphiy = (phi(i,j) - phi(i,j-1)) / dx(2)
            vmac(i,j) = vmac(i,j) - betay(i,j)*gphiy
         end do
      end do

      ! Here we reset phi == 0 at BC_DIR to be used in later iteration if necessary
      if (press_bc(1,1) == BC_DIR) phi(-1,:) = ZERO
      if (press_bc(1,2) == BC_DIR) phi(nx,:) = ZERO
      if (press_bc(2,1) == BC_DIR) phi(:,-1) = ZERO
      if (press_bc(2,2) == BC_DIR) phi(:,ny) = ZERO

    end subroutine mkumac_2d

    subroutine correct_mkumac_2d(umac,vmac,ng_um, &
                                 lo_x_flx,hi_x_flx,lo_y_flx,hi_y_flx, &
                                 lo,hi,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_um
      real(kind=dp_t), intent(inout) :: umac(lo(1)-ng_um:,lo(2)-ng_um:)
      real(kind=dp_t), intent(inout) :: vmac(lo(1)-ng_um:,lo(2)-ng_um:)
      real(kind=dp_t), intent(in   ) :: lo_x_flx(:,lo(2):), lo_y_flx(lo(1):,:)
      real(kind=dp_t), intent(in   ) :: hi_x_flx(:,lo(2):), hi_y_flx(lo(1):,:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j

      do j = lo(2),hi(2)
         umac(lo(1)  ,j) = umac(lo(1)  ,j) - lo_x_flx(1,j) * dx(1)
         umac(hi(1)+1,j) = umac(hi(1)+1,j) + hi_x_flx(1,j) * dx(1)
      end do


      do i = lo(1),hi(1)
         vmac(i,lo(2)  ) = vmac(i,lo(2)  ) - lo_y_flx(i,1) * dx(2)
         vmac(i,hi(2)+1) = vmac(i,hi(2)+1) + hi_y_flx(i,1) * dx(2)
      end do

    end subroutine correct_mkumac_2d

    subroutine mkumac_3d(n,umac,vmac,wmac,ng_um,phi,ng_p, &
                         betax,betay,betaz,ng_b,dx,press_bc)

      integer        , intent(in   ) :: n,ng_um,ng_p,ng_b
      real(kind=dp_t), intent(inout) :: umac(-ng_um:,-ng_um:,-ng_um:)
      real(kind=dp_t), intent(inout) :: vmac(-ng_um:,-ng_um:,-ng_um:)
      real(kind=dp_t), intent(inout) :: wmac(-ng_um:,-ng_um:,-ng_um:)
      real(kind=dp_t), intent(inout) ::  phi(-ng_p:,-ng_p:,-ng_p:)
      real(kind=dp_t), intent(in   ) :: betax(-ng_b:,-ng_b:,-ng_b:)
      real(kind=dp_t), intent(in   ) :: betay(-ng_b:,-ng_b:,-ng_b:)
      real(kind=dp_t), intent(in   ) :: betaz(-ng_b:,-ng_b:,-ng_b:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      real(kind=dp_t) :: gphix,gphiy,gphiz
      integer :: i,j,k,nx,ny,nz
      integer :: imin,imax,jmin,jmax,kmin,kmax

      nx = size(phi,dim=1) - 2
      ny = size(phi,dim=2) - 2
      nz = size(phi,dim=3) - 2

      if (n.eq.1) then
         imin = 0
         imax = nx
         jmin = 0
         jmax = ny
         kmin = 0
         kmax = nz
      else
         imin = 1
         imax = nx-1
         jmin = 1
         jmax = ny-1
         kmin = 1
         kmax = nz-1
      end if

      if (press_bc(1,1) == BC_NEU) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(-1,j,k) = phi(0,j,k)
            end do
         end do
      else if (press_bc(1,1) == BC_DIR) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(-1,j,k) = -TWO*phi(0,j,k) + THIRD * phi(1,j,k)
            end do
         end do
      end if
      if (press_bc(1,2) == BC_NEU) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(nx,j,k) = phi(nx-1,j,k)
            end do
         end do
      else if (press_bc(1,2) == BC_DIR) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(nx,j,k) = -TWO*phi(nx-1,j,k) + THIRD * phi(nx-2,j,k)
            end do
         end do
      end if
      if (press_bc(2,1) == BC_NEU) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,-1,k) = phi(i,0,k)
            end do
         end do
      else if (press_bc(2,1) == BC_DIR) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,-1,k) = -TWO*phi(i,0,k) + THIRD * phi(i,1,k)
            end do
         end do
      end if
      if (press_bc(2,2) == BC_NEU) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,ny,k) = phi(i,ny-1,k)
            end do
         end do
      else if (press_bc(2,2) == BC_DIR) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,ny,k) = -TWO*phi(i,ny-1,k) + THIRD * phi(i,ny-2,k)
            end do
         end do
      end if
      if (press_bc(3,1) == BC_NEU) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,-1) = phi(i,j,0)
            end do
         end do
      else if (press_bc(3,1) == BC_DIR) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,-1) = -TWO*phi(i,j,0) + THIRD * phi(i,j,1)
            end do
         end do
      end if
      if (press_bc(3,2) == BC_NEU) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,nz) = phi(i,j,nz-1)
            end do
         end do
      else if (press_bc(3,2) == BC_DIR) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,nz) = -TWO*phi(i,j,nz-1) + THIRD * phi(i,j,nz-2)
            end do
         end do
      end if

!$omp parallel do private(i,j,k,gphix)
      do k = 0,nz-1 
         do j = 0,ny-1
            do i = imin,imax
               gphix = (phi(i,j,k) - phi(i-1,j,k)) / dx(1)
               umac(i,j,k) = umac(i,j,k) - betax(i,j,k)*gphix
            end do
         end do
      end do
!$omp end parallel do

!$omp parallel do private(i,j,k,gphiy)
      do k = 0,nz-1
         do j = jmin,jmax
            do i = 0,nx-1
               gphiy = (phi(i,j,k) - phi(i,j-1,k)) / dx(2)
               vmac(i,j,k) = vmac(i,j,k) - betay(i,j,k)*gphiy
            end do
         end do
      end do
!$omp end parallel do

!$omp parallel do private(i,j,k,gphiz)
      do k = kmin,kmax
         do j = 0,ny-1
            do i = 0,nx-1
               gphiz = (phi(i,j,k) - phi(i,j,k-1)) / dx(3)
               wmac(i,j,k) = wmac(i,j,k) - betaz(i,j,k)*gphiz
            end do
         end do
      end do
!$omp end parallel do

      ! Here we reset phi == 0 at BC_DIR to be used in later iteration if necessary
      if (press_bc(1,1) == BC_DIR) phi(-1,:,:) = ZERO
      if (press_bc(1,2) == BC_DIR) phi(nx,:,:) = ZERO
      if (press_bc(2,1) == BC_DIR) phi(:,-1,:) = ZERO
      if (press_bc(2,2) == BC_DIR) phi(:,ny,:) = ZERO
      if (press_bc(3,1) == BC_DIR) phi(:,:,-1) = ZERO
      if (press_bc(3,2) == BC_DIR) phi(:,:,nz) = ZERO

    end subroutine mkumac_3d

    subroutine correct_mkumac_3d(umac,vmac,wmac,ng_um, &
                                 lo_x_flx,hi_x_flx,lo_y_flx,hi_y_flx,lo_z_flx,hi_z_flx, &
                                 lo,hi,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_um
      real(kind=dp_t), intent(inout) :: umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(inout) :: vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(inout) :: wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
      real(kind=dp_t), intent(in   ) :: lo_x_flx(:,lo(2):,lo(3):), hi_x_flx(:,lo(2):,lo(3):)
      real(kind=dp_t), intent(in   ) :: lo_y_flx(lo(1):,:,lo(3):), hi_y_flx(lo(1):,:,lo(3):)
      real(kind=dp_t), intent(in   ) :: lo_z_flx(lo(1):,lo(2):,:), hi_z_flx(lo(1):,lo(2):,:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k

!$omp parallel do private(j,k)
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
         umac(lo(1)  ,j,k) = umac(lo(1)  ,j,k) - lo_x_flx(1,j,k) * dx(1)
         umac(hi(1)+1,j,k) = umac(hi(1)+1,j,k) + hi_x_flx(1,j,k) * dx(1)
      end do
      end do
!$omp end parallel do

!$omp parallel do private(i,k)
      do k = lo(3),hi(3)
      do i = lo(1),hi(1)
         vmac(i,lo(2)  ,k) = vmac(i,lo(2)  ,k) - lo_y_flx(i,1,k) * dx(2)
         vmac(i,hi(2)+1,k) = vmac(i,hi(2)+1,k) + hi_y_flx(i,1,k) * dx(2)
      end do
      end do
!$omp end parallel do

!$omp parallel do private(i,j)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         wmac(i,j,lo(3)  ) = wmac(i,j,lo(3)  ) - lo_z_flx(i,j,1) * dx(3)
         wmac(i,j,hi(3)+1) = wmac(i,j,hi(3)+1) + hi_z_flx(i,j,1) * dx(3)
      end do
      end do
!$omp end parallel do

    end subroutine correct_mkumac_3d

  end subroutine macproject

end module macproject_module
