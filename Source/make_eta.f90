! Compute eta_rho = Avg { rho' Utilde.e_r }  (see paper III, Eq. 30)
!
! For plane-parallel geometries, we compute eta_rho by averaging up 
! interface fluxes (etarho_flux) created in mkflux.
!
! For spherical geometries, we construct a multifab containing 
! { rho' Utilde.e_r } and use the average routine to put it in cell-centers 
! on the base state.
!
! We keep make two quantities here: etarho is edge-centered and etarho_cc
! is cell-centered.  Only spherical (in make_w0) needs the cell-centered
! quantity.

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
  subroutine make_etarho_planar(nlevs,etarho,etarho_cc,etarhoflux,mla)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, r_start_coord, r_end_coord, numdisjointchunks
    use restrict_base_module

    integer           , intent(in   ) :: nlevs
    real(kind=dp_t)   , intent(inout) :: etarho(:,0:)
    real(kind=dp_t)   , intent(inout) :: etarho_cc(:,0:)
    type(multifab)    , intent(inout) :: etarhoflux(:)
    type(ml_layout)   , intent(inout) :: mla

    ! local
    real(kind=dp_t), pointer :: efp(:,:,:,:)
    
    real(kind=dp_t), allocatable :: ncell(:,:)
    real(kind=dp_t), allocatable :: etarhosum_proc(:,:)
    real(kind=dp_t), allocatable :: etarhosum(:,:)
    real(kind=dp_t), allocatable :: etarhopert_proc(:,:)
    real(kind=dp_t), allocatable :: etarhopert(:,:)
    real(kind=dp_t), allocatable :: source_buffer(:)
    real(kind=dp_t), allocatable :: target_buffer(:)
    logical                      :: fine_grids_span_domain_width

    type(box) :: domain

    integer :: domlo(mla%dim),domhi(mla%dim)
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: i,r,rpert,n,dm,rr,ng_e

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_etarho")

    dm = mla%dim
    ng_e = etarhoflux(1)%ng

    ! ncell is a function of r only for spherical
    allocate(ncell          (nlevs,0:nr_fine)) 
    allocate(etarhosum_proc (nlevs,0:nr_fine))
    allocate(etarhosum      (nlevs,0:nr_fine))
    allocate(etarhopert_proc(nlevs,0:nr_fine))
    allocate(etarhopert     (nlevs,0:nr_fine))

    allocate(source_buffer(0:nr_fine))
    allocate(target_buffer(0:nr_fine))

    ncell        = ZERO
    etarhosum_proc  = ZERO
    etarhosum       = ZERO
    etarhopert_proc = ZERO
    etarhopert      = ZERO
    etarho          = ZERO

    if (spherical .eq. 0) then

       fine_grids_span_domain_width = .true.

       if (fine_grids_span_domain_width) then
          
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
                   call sum_etarho_2d(n,lo,hi,domhi,efp(:,:,1,1),ng_e,etarhosum_proc(n,:))
                case (3)
                   call sum_etarho_3d(n,lo,hi,domhi,efp(:,:,:,1),ng_e,etarhosum_proc(n,:))
                end select
             end do
             
             ! gather etarhosum
             source_buffer = etarhosum_proc(n,:)
             call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
             etarhosum(n,:) = target_buffer
             
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n,i),r_end_coord(n,i)+1
                   etarho(n,r) = etarhosum(n,r) / dble(ncell(n,r))
                end do
             end do
             
          end do
          
       else
    
          domain = layout_get_pd(mla%la(1))
          domlo = lwb(domain)
          domhi = upb(domain)
          
          if (dm .eq. 2) then
             ncell(1,:) = domhi(1)-domlo(1)+1
          else if(dm .eq. 3) then
             ncell(1,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
          end if
          
          ! the first step is to compute etarho assuming the coarsest level 
          ! is the only level in existence
          do i=1,layout_nboxes(mla%la(1))
             if ( multifab_remote(etarhoflux(1), i) ) cycle
             efp => dataptr(etarhoflux(1), i)
             lo =  lwb(get_box(mla%la(1), i))
             hi =  upb(get_box(mla%la(1), i))
             select case (dm)
             case (2)
                call sum_etarho_2d(1,lo,hi,domhi,efp(:,:,1,1),ng_e,etarhosum_proc(1,:))
             case (3)
                call sum_etarho_3d(1,lo,hi,domhi,efp(:,:,:,1),ng_e,etarhosum_proc(1,:))
             end select
          end do
          
          ! gather etarhosum
          source_buffer = etarhosum_proc(1,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          etarhosum(1,:) = target_buffer
          
          do r=0,r_end_coord(1,1)
             etarho(1,r) = etarhosum(1,r) / dble(ncell(1,r))
          end do
          
          ! now we compute etarho at the finer levels
          do n=2,nlevs
             
             rr = mla%mba%rr(n-1,dm)
             
             if (mla%mba%rr(n-1,1) .ne. mla%mba%rr(n-1,dm) .or. &
                  mla%mba%rr(n-1,2) .ne. mla%mba%rr(n-1,dm)) then
                print*,"ERROR: In make_etarho, refinement ratio in each direction must match"
                stop
             endif
             
             domain = layout_get_pd(mla%la(n))
             domlo  = lwb(domain)
             domhi  = upb(domain)
             
             if (dm .eq. 2) then
                ncell(n,:) = domhi(1)-domlo(1)+1
             else if (dm .eq. 3) then
                ncell(n,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
             end if
             
             ! on faces that exist at the next coarser level, etarho is the same since
             ! the fluxes have been restricted.  We copy these values of etarho directly and
             ! copy scaled values of etarhosum directly.
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n-1,i),r_end_coord(n-1,i)+1
                   etarho   (n,r*rr) = etarho   (n-1,r)
                   etarhosum(n,r*rr) = etarhosum(n-1,r)*rr**(dm-1)
                end do
             end do
             
             ! on faces that do not exist at the next coarser level, we use linear
             ! interpolation to get etarhosum at these faces.
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n-1,i),r_end_coord(n-1,i)
                   do rpert=1,rr-1
                      etarhosum(n,r*rr+rpert) = &
                           dble(rpert)/dble(rr)*(etarhosum(n,r*rr)) + &
                           dble(rr-rpert)/dble(rr)*(etarhosum(n,(r+1)*rr))
                   end do
                end do
             end do
             
             ! compute etarhopert_proc on faces that do not exist at the coarser level
             do i=1,layout_nboxes(mla%la(n))
                if ( multifab_remote(etarhoflux(n), i) ) cycle
                efp  => dataptr(etarhoflux(n), i)
                lo =  lwb(get_box(mla%la(n), i))
                hi =  upb(get_box(mla%la(n), i))
                select case (dm)
                case (2)
                   call compute_etarhopert_2d(lo,hi,efp(:,:,1,1),ng_e,etarhosum_proc(1,:),rr)
                case (3)
                   call compute_etarhopert_3d(lo,hi,efp(:,:,:,1),ng_e,etarhosum_proc(1,:),rr)
                end select
             end do
             
             ! gather etarhopert for rhoh
             source_buffer = etarhopert_proc(n,:)
             call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
             etarhopert(n,:) = target_buffer
             
             ! update etasum on faces that do not exist at the coarser level
             ! then recompute eta on these faces
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n-1,i),r_end_coord(n-1,i)
                   do rpert=1,rr-1
                      etarhosum(n,r*rr+rpert) = &
                           etarhosum(n,r*rr+rpert) + etarhopert(n,r*rr+rpert)
                      etarho(n,r*rr+rpert) = &
                           etarhosum(n,r*rr+rpert)/dble(ncell(n,r*rr+rpert))
                   end do
                end do
             end do

          end do ! end loop over levels
          
       end if

       call restrict_base(nlevs,etarho,.false.)
       
    else

       call bl_error("ERROR: make_eta should not be called for spherical")

    end if

    ! make the cell-centered etarho_cc by averaging etarho to centers
    do n=1,nlevs
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             etarho_cc(n,r) = HALF*(etarho(n,r) + etarho(n,r+1))
          enddo
       enddo
    enddo

    deallocate(ncell)
    deallocate(etarhosum_proc,etarhosum)
    deallocate(etarhopert_proc,etarhopert)
    deallocate(source_buffer,target_buffer)

    call destroy(bpt)

  end subroutine make_etarho_planar

  subroutine sum_etarho_2d(n,lo,hi,domhi,etarhoflux,ng_e,etarhosum)

    use geometry, only: r_end_coord, numdisjointchunks

    integer         , intent(in   ) :: n, lo(:), hi(:), domhi(:), ng_e
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

  subroutine sum_etarho_3d(n,lo,hi,domhi,etarhoflux,ng_e,etarhosum)

    use geometry, only: r_end_coord, numdisjointchunks

    integer         , intent(in   ) :: n,lo(:), hi(:), domhi(:), ng_e
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

  subroutine compute_etarhopert_2d(lo,hi,etarhoflux,ng_e,etarhopert,rr)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), rr, ng_e
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1)-ng_e:,lo(2)-ng_e:)
    real (kind=dp_t), intent(inout) :: etarhopert(0:)
    
    ! local
    integer         :: i,j,ipert,jpert
    real(kind=dp_t) :: loavg,hiavg,crseavg

    do j=lo(2),hi(2)-rr,rr
       do jpert=1,rr-1

          do i=lo(1),hi(1)-rr,rr

             loavg = ZERO
             hiavg = ZERO
             do ipert=0,rr-1
                loavg = loavg + etarhoflux(i+ipert,j)
                hiavg = hiavg + etarhoflux(i+ipert,j+rr)
             end do
             loavg = loavg / dble(rr)
             hiavg = hiavg / dble(rr)
             crseavg = dble(jpert)/dble(rr)*loavg + dble(rr-jpert)/dble(rr)*hiavg
             do ipert=0,rr-1
                etarhopert(j+jpert) = etarhoflux(i+ipert,j+jpert) - crseavg
             end do
             
          end do

       end do
    end do

  end subroutine compute_etarhopert_2d

  subroutine compute_etarhopert_3d(lo,hi,etarhoflux,ng_e,etarhopert,rr)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), rr, ng_e
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real (kind=dp_t), intent(inout) :: etarhopert(0:)
    
    ! local
    integer         :: i,j,k,ipert,jpert,kpert
    real(kind=dp_t) :: loavg,hiavg,crseavg

    do k=lo(3),hi(3)-rr,rr
       do kpert=1,rr-1

          do j=lo(2),hi(2)-rr,rr
             do i=lo(1),hi(1)-rr,rr

                loavg = ZERO
                hiavg = ZERO
                do ipert=0,rr-1
                   do jpert=0,rr-1
                      loavg = loavg + etarhoflux(i+ipert,j+jpert,k)
                      hiavg = hiavg + etarhoflux(i+ipert,j+jpert,k+rr)
                   end do
                end do
                loavg = loavg / dble(rr**2)
                hiavg = hiavg / dble(rr**2)
                crseavg = dble(kpert)/dble(rr)*loavg + dble(rr-kpert)/dble(rr)*hiavg
                do ipert=0,rr-1
                   do jpert=0,rr-1
                      etarhopert(k+kpert) = &
                           etarhoflux(i+ipert,j+jpert,k+kpert) - crseavg
                   end do
                end do

             end do
          end do

       end do
    end do

  end subroutine compute_etarhopert_3d


  !---------------------------------------------------------------------------
  ! spherical routines
  !---------------------------------------------------------------------------
  subroutine make_etarho_spherical(nlevs,sold,snew,umac,rho0_old,rho0_new, &
                                   dx,normal,etarho,etarho_cc,mla,the_bc_level)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, nr
    use variables
    use average_module
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:), snew(:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:), rho0_new(:,0:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(  out) :: etarho(:,0:)
    real(kind=dp_t), intent(  out) :: etarho_cc(:,0:)
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    type(multifab) :: eta_cart(mla%nlevel)
    
    real(kind=dp_t), pointer :: ep(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:), snp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    integer :: n,i,lo(sold(1)%dim),hi(sold(1)%dim),ng_so,ng_sn,ng_um,ng_n,ng_e
    integer :: r


    ! construct a multifab containing  [ rho' (Utilde . e_r) ]
    ng_so = sold(1)%ng
    ng_sn = snew(1)%ng
    ng_um = umac(1,1)%ng    ! here we are assuming all components have the 
                            ! same # of ghostcells
    ng_n = normal(1)%ng

    do n=1,nlevs

       call multifab_build(eta_cart(n), sold(n)%la, 1, 1)

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

          ng_e = eta_cart(n)%ng 

          call construct_eta_cart(n, sop(:,:,:,rho_comp), ng_so, &
                                  snp(:,:,:,rho_comp), ng_sn, &
                                  ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                  np(:,:,:,:), ng_n, ep(:,:,:,1), ng_e, &
                                  rho0_old(n,:), rho0_new(n,:), &
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
          call ml_cc_restriction(eta_cart(n-1)    ,eta_cart(n)    ,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(eta_cart(n),eta_cart(n-1), &
                                         ng_e,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         1,foextrap_comp,1)
       enddo

    end if
    
    
    ! average 
    call average(mla,eta_cart,etarho_cc,dx,1)


    do n=1,nlevs
       call destroy(eta_cart(n))
    enddo

    ! put eta on base state edges -- here we are assuming that there
    ! is no refinement
    do n=1,nlevs
       
       ! the 0th value of etarho = 0, since Utilde . e_r must be 
       ! zero at the center (since e_r is not defined there)
       etarho(n,0) = ZERO
       do r=1, nr(n)-1
          etarho(n,r) = HALF*(etarho_cc(n,r) + etarho_cc(n,r-1))
       enddo

       ! probably should do some better extrapolation here eventually
       etarho(n,nr(n)) = etarho_cc(n,nr(n)-1)

    enddo

  end subroutine make_etarho_spherical

  subroutine construct_eta_cart(n, rho_old, ng_so, rho_new, ng_sn, &
                                umac, vmac, wmac, ng_um, &
                                normal, ng_n, eta_cart, ng_e, &
                                rho0_old, rho0_new, dx, lo, hi)

    use bl_constants_module
    use geometry, only: nr_fine
    use fill_3d_module

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_so, ng_sn, ng_um, ng_n, ng_e
    real(kind=dp_t), intent(in   ) ::  rho_old(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:)
    real(kind=dp_t), intent(in   ) ::  rho_new(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:)
    real(kind=dp_t), intent(in   ) ::     umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)    
    real(kind=dp_t), intent(in   ) ::     vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)    
    real(kind=dp_t), intent(in   ) ::     wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)    
    real(kind=dp_t), intent(in   ) ::   normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:)
    real(kind=dp_t), intent(inout) :: eta_cart(lo(1)-ng_e :,lo(2)-ng_e :,lo(3)-ng_e :)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:), rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t), allocatable :: rho0_nph(:)
    real(kind=dp_t), allocatable :: rho0_cart(:,:,:,:)

    real(kind=dp_t) :: Utilde_dot_er
    integer :: i,j,k,r

    ! put the time-centered base state density on a Cartesian patch.
    allocate(rho0_nph(0:nr_fine-1))
    do r = 0, nr_fine-1
       rho0_nph(r) = HALF*(rho0_old(r) + rho0_new(r))
    enddo

    allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,rho0_nph,rho0_cart,lo,hi,dx,0,0)


    ! construct time-centered [ rho' (Utilde . e_r) ]
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             Utilde_dot_er = HALF*(umac(i,j,k) + umac(i+1,j,k)) * normal(i,j,k,1) + &
                             HALF*(vmac(i,j,k) + vmac(i,j+1,k)) * normal(i,j,k,2) + &
                             HALF*(wmac(i,j,k) + wmac(i,j,k+1)) * normal(i,j,k,3)

             eta_cart(i,j,k) = (HALF*(rho_old(i,j,k) + rho_new(i,j,k)) - &
                                rho0_cart(i,j,k,1)) * Utilde_dot_er

          enddo
       enddo
    enddo

    deallocate (rho0_cart,rho0_nph)

  end subroutine construct_eta_cart

end module make_eta_module
