!
! Compute eta_rho = Avg { rho' U dot e_r }  (see paper III, Eq. 30)
!
! We keep make three quantities here: 
!    etarho     is edge-centered
!    etarho_cc  is cell-centered
!
! For plane-parallel geometries, we compute etarho by averaging up 
! interface fluxes (etarho_flux) created in mkflux.  We compute etarho_cc
! from etarho.
!
! For spherical geometries, 
!      We construct a multifab containing {rho' (U dot e_r)} and 
!      use the average routine to put it in cell-centers 
!      on the base state to get etarho_cc.  We compute etarho from these 
!      cell-centered quantites by averaging to the center.  
!

module make_eta_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: make_etarho_planar, make_etarho_spherical

contains

  !---------------------------------------------------------------------------
  ! plane-parallel geometry routines
  !---------------------------------------------------------------------------

  subroutine make_etarho_planar(etarho_ec,etarho_cc,etarhoflux,mla)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, r_start_coord, r_end_coord, numdisjointchunks
    use restrict_base_module

    real(kind=dp_t)   , intent(  out) :: etarho_ec(:,0:)
    real(kind=dp_t)   , intent(  out) :: etarho_cc(:,0:)
    type(multifab)    , intent(in   ) :: etarhoflux(:)
    type(ml_layout)   , intent(in   ) :: mla

    ! local
    real(kind=dp_t), pointer :: efp(:,:,:,:)
    
    real(kind=dp_t) :: ncell(mla%nlevel)

    real(kind=dp_t) :: etarhosum_proc(0:nr_fine,mla%nlevel)
    real(kind=dp_t) ::      etarhosum(0:nr_fine,mla%nlevel)

    type(box) :: domain

    integer :: domlo(mla%dim),domhi(mla%dim),lo(mla%dim),hi(mla%dim)
    integer :: i,r,n,ng_e,dm,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_etarho")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_e = nghost(etarhoflux(1))

    etarhosum_proc  = ZERO
    etarhosum       = ZERO
    etarho_ec       = ZERO
    etarho_cc       = ZERO
    
    if (spherical .eq. 1) then
       call bl_error("ERROR: make_eta should not be called for spherical")
    end if
    
    do n=1,nlevs
       domain = layout_get_pd(mla%la(n))
       domlo  = lwb(domain)
       domhi  = upb(domain)

       if (dm .eq. 1) then
          ncell(n) = 1
       else if (dm .eq. 2) then
          ncell(n) = domhi(1)-domlo(1)+1
       else if (dm .eq. 3) then
          ncell(n) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
       end if

       do i=1,nfabs(etarhoflux(n))
          efp => dataptr(etarhoflux(n), i)
          lo =  lwb(get_box(etarhoflux(n), i))
          hi =  upb(get_box(etarhoflux(n), i))
          select case (dm)
          case (1)
             call sum_etarho_1d(n,lo,hi,efp(:,1,1,1),ng_e,etarhosum_proc(:,n))
          case (2)
             call sum_etarho_2d(n,lo,hi,efp(:,:,1,1),ng_e,etarhosum_proc(:,n))
          case (3)
             call sum_etarho_3d(n,lo,hi,efp(:,:,:,1),ng_e,etarhosum_proc(:,n))
          end select
       end do

       call parallel_reduce(etarhosum(:,n), etarhosum_proc(:,n), MPI_SUM)

       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)+1
             etarho_ec(n,r) = etarhosum(r,n) / dble(ncell(n))
          end do
       end do

    end do

    ! These calls shouldn't be needed since the planar algorithm doesn't use
    ! these outside of this function, but this is just to be safe in case
    ! things change in the future.
    call restrict_base(etarho_ec,.false.)
    call fill_ghost_base(etarho_ec,.false.)

    ! make the cell-centered etarho_cc by averaging etarho to centers
    do n=1,nlevs
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             etarho_cc(n,r) = HALF*(etarho_ec(n,r) + etarho_ec(n,r+1))
          enddo
       enddo
    enddo

    ! These calls shouldn't be needed since the planar algorithm only uses
    ! etarho_cc to make_psi, and then we fill ghost cells in make_psi, but
    ! this is just to be safe in case things change in the future
    call restrict_base(etarho_cc,.true.)
    call fill_ghost_base(etarho_cc,.true.)

    call destroy(bpt)

  end subroutine make_etarho_planar

!---------------------------------------------------------------------------

  subroutine sum_etarho_1d(n,lo,hi,etarhoflux,ng_e,etarhosum)

    use geometry, only: r_end_coord, numdisjointchunks

    integer         , intent(in   ) :: n, lo(:), hi(:), ng_e
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1)-ng_e:)
    real (kind=dp_t), intent(inout) :: etarhosum(0:)

    ! local
    integer :: i,k
    logical :: top_edge

    do i=lo(1),hi(1)
       etarhosum(i) = etarhoflux(i)
    end do

    ! we only add the contribution at the top edge if we are at the top of grid at a level
    ! this prevents double counting
    top_edge = .false.
    do k = 1,numdisjointchunks(n)
       if (hi(1) .eq. r_end_coord(n,k)) then
          top_edge = .true.
       end if
    end do

    if (top_edge) then
       i=hi(1)+1
       etarhosum(i) = etarhoflux(i)
    end if

  end subroutine sum_etarho_1d

!---------------------------------------------------------------------------

  subroutine sum_etarho_2d(n,lo,hi,etarhoflux,ng_e,etarhosum)

    use geometry, only: r_end_coord, numdisjointchunks

    integer         , intent(in   ) :: n, lo(:), hi(:), ng_e
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1)-ng_e:,lo(2)-ng_e:)
    real (kind=dp_t), intent(inout) :: etarhosum(0:)

    ! local
    integer :: i,j
    logical :: top_edge

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          etarhosum(j) = etarhosum(j) + etarhoflux(i,j)
       end do
    end do

    ! we only add the contribution at the top edge if we are at the top of grid at a level
    ! this prevents double counting
    top_edge = .false.
    do i=1,numdisjointchunks(n)
       if (hi(2) .eq. r_end_coord(n,i)) then
          top_edge = .true.
       end if
    end do
    if(top_edge) then
       j=hi(2)+1
       do i=lo(1),hi(1)
          etarhosum(j) = etarhosum(j) + etarhoflux(i,j)
       end do
    end if

  end subroutine sum_etarho_2d

!---------------------------------------------------------------------------

  subroutine sum_etarho_3d(n,lo,hi,etarhoflux,ng_e,etarhosum)

    use geometry, only: r_end_coord, numdisjointchunks

    integer         , intent(in   ) :: n,lo(:), hi(:), ng_e
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real (kind=dp_t), intent(inout) :: etarhosum(0:)

    ! local
    integer :: i,j,k
    logical :: top_edge

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etarhosum(k) = etarhosum(k) + etarhoflux(i,j,k)
          end do
       end do
    end do

    ! we only add the contribution at the top edge if we are at the top of the domain
    ! this prevents double counting
    top_edge = .false.
    do i=1,numdisjointchunks(n)
       if (hi(3) .eq. r_end_coord(n,i)) then
          top_edge = .true.
       end if
    end do
    if(top_edge) then
       k=hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etarhosum(k) = etarhosum(k) + etarhoflux(i,j,k)
          end do
       end do
    end if

  end subroutine sum_etarho_3d

  !---------------------------------------------------------------------------
  ! spherical routines
  !---------------------------------------------------------------------------

  subroutine make_etarho_spherical(sold,snew,umac,w0mac,rho0_old,rho0_new, &
                                   dx,normal,etarho_ec,etarho_cc,mla,the_bc_level)

    use bl_constants_module
    use geometry, only: spherical, nr_fine
    use variables
    use average_module
    use ml_restrict_fill_module

    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    type(multifab) , intent(in   ) :: sold(:), snew(:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:), rho0_new(:,0:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(  out) :: etarho_ec(:,0:)
    real(kind=dp_t), intent(  out) :: etarho_cc(:,0:)
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    type(multifab) :: eta_cart(mla%nlevel)
    
    real(kind=dp_t), pointer :: ep(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:), snp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)
    real(kind=dp_t), pointer :: wxp(:,:,:,:), wyp(:,:,:,:), wzp(:,:,:,:)
    real(kind=dp_t), pointer :: nop(:,:,:,:)

    integer :: n,i,lo(mla%dim),hi(mla%dim),ng_so,ng_sn,ng_um,ng_n,ng_e,ng_wm
    integer :: r,dm,nlevs

    dm = mla%dim
    nlevs = mla%nlevel

    if (spherical .eq. 0) then
       call bl_error("ERROR: make_eta_spherical should not be called for plane-parallel")
    end if

    ! construct a multifab containing  [ rho' (U dot e_r) ] 
    ! and another containing [ rho' ]
    ng_so = nghost(sold(1))
    ng_sn = nghost(snew(1))
    ng_um = nghost(umac(1,1))
    ng_n  = nghost(normal(1))
    ng_wm = nghost(w0mac(1,1))

    do n=1,nlevs

       call multifab_build( eta_cart(n), get_layout(sold(n)), 1, 1)

       ng_e = nghost(eta_cart(n))

       do i=1, nfabs(eta_cart(n))
          ep  => dataptr(eta_cart(n), i)
          sop => dataptr(sold(n), i)
          snp => dataptr(snew(n), i)
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)
          wmp => dataptr(umac(n,3), i)
          wxp => dataptr(w0mac(n,1), i)
          wyp => dataptr(w0mac(n,2), i)
          wzp => dataptr(w0mac(n,3), i)
          nop  => dataptr(normal(n), i)
          lo = lwb(get_box(eta_cart(n),i))
          hi = upb(get_box(eta_cart(n),i))
          call construct_eta_cart(sop(:,:,:,rho_comp), ng_so, &
                                  snp(:,:,:,rho_comp), ng_sn, &
                                  ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                  wxp(:,:,:,1), wyp(:,:,:,1), wzp(:,:,:,1), ng_wm, &
                                  nop(:,:,:,:), ng_n, ep(:,:,:,1), ng_e, &
                                  rho0_old(1,:), rho0_new(1,:), &
                                  dx(n,:), lo, hi)
       enddo

    enddo

    call ml_restrict_and_fill(nlevs, eta_cart, mla%mba%rr, the_bc_level, &
         icomp=1, bcomp=foextrap_comp, nc=1, ng=ng_e)

    ! compute etarho_cc as the average of eta_cart = [ rho' (U dot e_r) ]
    call average(mla,eta_cart,etarho_cc,dx,1)

    do n=1,nlevs
       call destroy(eta_cart(n))
    enddo

    ! put eta on base state edges
    ! note that in spherical the base state has no refinement
    ! the 0th value of etarho = 0, since U dot . e_r must be 
    ! zero at the center (since e_r is not defined there)
    etarho_ec(1,0) = ZERO
    do r=1,nr_fine-1
       etarho_ec(1,r) = HALF*(etarho_cc(1,r) + etarho_cc(1,r-1))
    enddo
    ! probably should do some better extrapolation here eventually
    etarho_ec(1,nr_fine) = etarho_cc(1,nr_fine-1)

  end subroutine make_etarho_spherical

  subroutine construct_eta_cart(rho_old, ng_so, rho_new, ng_sn, umac, vmac, wmac, ng_um, &
                                w0macx, w0macy, w0macz, ng_wm, &
                                normal, ng_n, eta_cart, ng_e, rho0_old, rho0_new, dx, lo, hi)

    use bl_constants_module
    use geometry, only: nr_fine
    use fill_3d_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_so, ng_sn, ng_um, ng_n, ng_e, ng_wm
    real(kind=dp_t), intent(in   ) ::       rho_old(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:)
    real(kind=dp_t), intent(in   ) ::       rho_new(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:)
    real(kind=dp_t), intent(in   ) ::          umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::          vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::          wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::        w0macx(lo(1)-ng_wm:,lo(2)-ng_wm:,lo(3)-ng_wm:)
    real(kind=dp_t), intent(in   ) ::        w0macy(lo(1)-ng_wm:,lo(2)-ng_wm:,lo(3)-ng_wm:)
    real(kind=dp_t), intent(in   ) ::        w0macz(lo(1)-ng_wm:,lo(2)-ng_wm:,lo(3)-ng_wm:)
    real(kind=dp_t), intent(in   ) ::        normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:)
    real(kind=dp_t), intent(inout) ::      eta_cart(lo(1)-ng_e :,lo(2)-ng_e :,lo(3)-ng_e :)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:), rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t) ::      rho0_nph(0:nr_fine-1)

    real(kind=dp_t), allocatable :: rho0_new_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: rho0_nph_cart(:,:,:,:)

    real(kind=dp_t) :: U_dot_er
    integer :: i,j,k,r

    allocate(rho0_new_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(rho0_nph_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    ! put the time-centered base state density on a Cartesian patch.
    do r = 0, nr_fine-1
       rho0_nph(r) = HALF*(rho0_old(r) + rho0_new(r))
    enddo

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,rho0_new,rho0_new_cart,lo,hi,dx,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,rho0_nph,rho0_nph_cart,lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k,U_dot_er)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             U_dot_er = HALF*(    umac(i,j,k) +   umac(i+1,j,k) &
                              + w0macx(i,j,k) + w0macx(i+1,j,k)) * normal(i,j,k,1) + &
                        HALF*(    vmac(i,j,k) +   vmac(i,j+1,k) &
                              + w0macy(i,j,k) + w0macy(i,j+1,k)) * normal(i,j,k,2) + &
                        HALF*(    wmac(i,j,k) +   wmac(i,j,k+1) &
                              + w0macz(i,j,k) + w0macz(i,j,k+1)) * normal(i,j,k,3)

             ! construct time-centered [ rho' (U dot e_r) ]
             eta_cart(i,j,k) = (HALF*(rho_old(i,j,k) + rho_new(i,j,k)) - &
                                rho0_nph_cart(i,j,k,1)) * U_dot_er

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(rho0_new_cart,rho0_nph_cart)
    
  end subroutine construct_eta_cart

end module make_eta_module
