! Compute eta_rho = Avg { rho' Utilde.e_r }  (see paper III, Eq. 30)
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
!      We construct a multifab containing ! { rho' Utilde.e_r } and 
!      use the average routine to put it in cell-centers 
!      on the base state to get etarho_cc.  We compute etarho from these 
!      cell-centered quantites by averaging to the center.  


module make_eta_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: make_etarho_planar, make_etarho_spherical, make_sprimebar_spherical

contains

  !---------------------------------------------------------------------------
  ! plane-parallel geometry routines
  !---------------------------------------------------------------------------

  subroutine make_etarho_planar(etarho_ec,etarho_cc,etarhoflux,mla)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, r_start_coord, r_end_coord, numdisjointchunks, &
         dr, dm, nlevs
    use restrict_base_module

    real(kind=dp_t)   , intent(  out) :: etarho_ec(:,0:)
    real(kind=dp_t)   , intent(  out) :: etarho_cc(:,0:)
    type(multifab)    , intent(in   ) :: etarhoflux(:)
    type(ml_layout)   , intent(in   ) :: mla

    ! local
    real(kind=dp_t), pointer :: efp(:,:,:,:)
    
    real(kind=dp_t) ::          ncell(nlevs,0:nr_fine)
    real(kind=dp_t) :: etarhosum_proc(nlevs,0:nr_fine)
    real(kind=dp_t) ::      etarhosum(nlevs,0:nr_fine)

    real(kind=dp_t) :: source_buffer(0:nr_fine)
    real(kind=dp_t) :: target_buffer(0:nr_fine)

    type(box) :: domain

    integer :: domlo(dm),domhi(dm)
    integer :: lo(dm),hi(dm)
    integer :: i,r,n,ng_e

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_etarho")

    ng_e = etarhoflux(1)%ng

    ncell           = ZERO
    etarhosum_proc  = ZERO
    etarhosum       = ZERO
    etarho_ec       = ZERO
    
    if (spherical .eq. 1) then
       call bl_error("ERROR: make_eta should not be called for spherical")
    end if
    
    do n=1,nlevs
       domain = layout_get_pd(mla%la(n))
       domlo = lwb(domain)
       domhi = upb(domain)

       if (dm .eq. 2) then
          ncell(n,:) = domhi(1)-domlo(1)+1
       else if(dm .eq. 3) then
          ncell(n,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
       end if

       do i=1,layout_nboxes(mla%la(n))
          if ( multifab_remote(etarhoflux(n), i) ) cycle
          efp => dataptr(etarhoflux(n), i)
          lo =  lwb(get_box(mla%la(n), i))
          hi =  upb(get_box(mla%la(n), i))
          select case (dm)
          case (2)
             call sum_etarho_2d(n,lo,hi,efp(:,:,1,1),ng_e,etarhosum_proc(n,:))
          case (3)
             call sum_etarho_3d(n,lo,hi,efp(:,:,:,1),ng_e,etarhosum_proc(n,:))
          end select
       end do

       ! gather etarhosum
       source_buffer = etarhosum_proc(n,:)
       call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
       etarhosum(n,:) = target_buffer

       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)+1
             etarho_ec(n,r) = etarhosum(n,r) / dble(ncell(n,r))
          end do
       end do

    end do

    call restrict_base(etarho_ec,.false.)

    ! make the cell-centered etarho_cc by averaging etarho to centers
    do n=1,nlevs
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             etarho_cc(n,r) = HALF*(etarho_ec(n,r) + etarho_ec(n,r+1))
          enddo
       enddo
    enddo

    call destroy(bpt)

  end subroutine make_etarho_planar

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

  subroutine make_etarho_spherical(sold,snew,umac,rho0_old,rho0_new, &
                                   dx,dt,normal,etarho_ec,etarho_cc,mla,the_bc_level)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, dm, nlevs
    use variables
    use average_module
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module

    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:), snew(:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:), rho0_new(:,0:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(  out) :: etarho_ec(:,0:)
    real(kind=dp_t), intent(  out) :: etarho_cc(:,0:)
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    type(multifab) :: eta_cart(mla%nlevel)
    
    real(kind=dp_t), pointer :: ep(:,:,:,:),rpp(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:), snp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    integer :: n,i,lo(dm),hi(dm),ng_so,ng_sn,ng_um,ng_n,ng_e
    integer :: r

    if (spherical .eq. 0) then
       call bl_error("ERROR: make_eta_spherical should not be called for plane-parallel")
    end if

    ! construct a multifab containing  [ rho' (Utilde . e_r) ] 
    ! and another containing [ rho' ]
    ng_so = sold(1)%ng
    ng_sn = snew(1)%ng
    ng_um = umac(1,1)%ng
    ng_n  = normal(1)%ng

    do n=1,nlevs

       call multifab_build(     eta_cart(n), sold(n)%la, 1, 1)

       ng_e = eta_cart(n)%ng

       do i=1,eta_cart(n)%nboxes
          if ( multifab_remote(eta_cart(n),i) ) cycle
          ep  => dataptr(eta_cart(n), i)
          sop => dataptr(sold(n), i)
          snp => dataptr(snew(n), i)
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)
          wmp => dataptr(umac(n,3), i)
          np  => dataptr(normal(n), i)
          lo = lwb(get_box(eta_cart(n),i))
          hi = upb(get_box(eta_cart(n),i))
          call construct_eta_cart(sop(:,:,:,rho_comp), ng_so, &
                                  snp(:,:,:,rho_comp), ng_sn, &
                                  ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                  np(:,:,:,:), ng_n, ep(:,:,:,1), ng_e, &
                                  rho0_old(1,:), rho0_new(1,:), &
                                  dx(n,:), lo, hi)
       enddo

    enddo

    ! fill eta_cart ghostcells
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(eta_cart(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(eta_cart(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(eta_cart(n-1)     ,eta_cart(n)     ,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(eta_cart(n),eta_cart(n-1), &
                                         ng_e,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         1,foextrap_comp,1)

       enddo

    end if
    
    ! compute etarho_cc as the average of eta_cart = [ rho' (Utilde . e_r) ]
    call average(mla,eta_cart,etarho_cc,dx,1)

    do n=1,nlevs
       call destroy(eta_cart(n))
    enddo

    ! put eta on base state edges -- here we are assuming that there
    ! is no refinement
       
    ! the 0th value of etarho = 0, since Utilde . e_r must be 
    ! zero at the center (since e_r is not defined there)
    etarho_ec(1,0) = ZERO
    do r=1,nr_fine-1
       etarho_ec(1,r) = HALF*(etarho_cc(1,r) + etarho_cc(1,r-1))
    enddo
    ! probably should do some better extrapolation here eventually
    etarho_ec(1,nr_fine) = etarho_cc(1,nr_fine-1)

  end subroutine make_etarho_spherical

  subroutine construct_eta_cart(rho_old, ng_so, rho_new, ng_sn, umac, vmac, wmac, ng_um, &
                                normal, ng_n, eta_cart, ng_e, rho0_old, rho0_new, dx, lo, hi)

    use bl_constants_module
    use geometry, only: nr_fine
    use fill_3d_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_so, ng_sn, ng_um, ng_n, ng_e
    real(kind=dp_t), intent(in   ) ::       rho_old(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:)
    real(kind=dp_t), intent(in   ) ::       rho_new(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:)
    real(kind=dp_t), intent(in   ) ::          umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::          vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::          wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::        normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:)
    real(kind=dp_t), intent(inout) ::      eta_cart(lo(1)-ng_e :,lo(2)-ng_e :,lo(3)-ng_e :)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:), rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t) ::      rho0_nph(0:nr_fine-1)
    real(kind=dp_t) :: rho0_new_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)
    real(kind=dp_t) :: rho0_nph_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)

    real(kind=dp_t) :: Utilde_dot_er
    integer :: i,j,k,r

    ! put the time-centered base state density on a Cartesian patch.
    do r = 0, nr_fine-1
       rho0_nph(r) = HALF*(rho0_old(r) + rho0_new(r))
    enddo

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,rho0_new,rho0_new_cart,lo,hi,dx,0,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,rho0_nph,rho0_nph_cart,lo,hi,dx,0,0)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             Utilde_dot_er = HALF*(umac(i,j,k) + umac(i+1,j,k)) * normal(i,j,k,1) + &
                             HALF*(vmac(i,j,k) + vmac(i,j+1,k)) * normal(i,j,k,2) + &
                             HALF*(wmac(i,j,k) + wmac(i,j,k+1)) * normal(i,j,k,3)

             ! construct time-centered [ rho' (Utilde . e_r) ]
             eta_cart(i,j,k) = (HALF*(rho_old(i,j,k) + rho_new(i,j,k)) - &
                                rho0_nph_cart(i,j,k,1)) * Utilde_dot_er

          enddo
       enddo
    enddo

  end subroutine construct_eta_cart

  subroutine make_sprimebar_spherical(s,comp,s0,dx,sprimebar,mla,the_bc_level)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, dm, nlevs
    use variables
    use average_module
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module

    integer        , intent(in   ) :: comp
    type(multifab) , intent(in   ) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(  out) :: sprimebar(:,0:)
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    type(multifab) :: sprime(mla%nlevel)
    
    real(kind=dp_t), pointer :: spp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    integer :: n,i,lo(dm),hi(dm),ng_sp,ng_s

    if (spherical .eq. 0) then
       call bl_error("ERROR: make_sprimebar_spherical shouldn't be called in plane-parallel")
    end if

    do n=1,nlevs

       call multifab_build(sprime(n),s(n)%la,1,1)

       ng_sp = sprime(n)%ng
       ng_s = s(n)%ng

       do i=1,sprime(n)%nboxes
          if ( multifab_remote(sprime(n),i) ) cycle
          spp  => dataptr(sprime(n), i)
          sp => dataptr(s(n), i)
          lo = lwb(get_box(sprime(n),i))
          hi = upb(get_box(sprime(n),i))
          call construct_sprime(sp(:,:,:,comp),ng_s,spp(:,:,:,1),ng_sp,s0(1,:),dx(n,:),lo,hi)
       enddo

    enddo

    ! fill sprime ghostcells
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(sprime(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(sprime(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(sprime(n-1),sprime(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(sprime(n),sprime(n-1),ng_sp,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n),1,foextrap_comp,1)
       enddo

    end if

    call average(mla,sprime,sprimebar,dx,1)

    do n=1,nlevs
       call destroy(sprime(n))
    enddo

  end subroutine make_sprimebar_spherical

  subroutine construct_sprime(s,ng_s,sprime,ng_sp,s0,dx,lo,hi)

    use bl_constants_module
    use geometry, only: nr_fine
    use fill_3d_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_sp
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(inout) :: sprime(lo(1)-ng_sp:,lo(2)-ng_sp:,lo(3)-ng_sp:)
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t) :: s0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)

    integer :: i,j,k

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0,s0_cart,lo,hi,dx,0,0)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             sprime(i,j,k) = s(i,j,k) - s0_cart(i,j,k,1)
 
          enddo
       enddo
    enddo

  end subroutine construct_sprime

end module make_eta_module
