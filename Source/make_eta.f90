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
  subroutine make_etarho_planar(nlevs,etarho,etarhoflux,mla)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, r_start_coord, r_end_coord
    use restrict_base_module

    integer           , intent(in   ) :: nlevs
    real(kind=dp_t)   , intent(inout) :: etarho(:,0:)
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
    integer :: i,r,rpert,n,dm,rr

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_etarho")

    dm = mla%dim

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
                   call sum_etarho_coarsest_2d(lo,hi,domhi,efp(:,:,1,1),etarhosum_proc(n,:))
                case (3)
                   call sum_etarho_coarsest_3d(lo,hi,domhi,efp(:,:,:,1),etarhosum_proc(n,:))
                end select
             end do
             
             ! gather etarhosum
             source_buffer = etarhosum_proc(n,:)
             call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
             etarhosum(n,:) = target_buffer
             
             do r=r_start_coord(n),r_end_coord(n)+1
                etarho(n,r) = etarhosum(n,r) / dble(ncell(n,r))
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
                call sum_etarho_coarsest_2d(lo,hi,domhi,efp(:,:,1,1),etarhosum_proc(1,:))
             case (3)
                call sum_etarho_coarsest_3d(lo,hi,domhi,efp(:,:,:,1),etarhosum_proc(1,:))
             end select
          end do
          
          ! gather etarhosum
          source_buffer = etarhosum_proc(1,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          etarhosum(1,:) = target_buffer
          
          do r=0,r_end_coord(1)
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
             do r=r_start_coord(n-1),r_end_coord(n-1)+1
                etarho   (n,r*rr) = etarho   (n-1,r)
                etarhosum(n,r*rr) = etarhosum(n-1,r)*rr**(dm-1)
             end do
             
             ! on faces that do not exist at the next coarser level, we use linear
             ! interpolation to get etarhosum at these faces.
             do r=r_start_coord(n-1),r_end_coord(n-1)
                do rpert=1,rr-1
                   etarhosum(n,r*rr+rpert) = &
                        dble(rpert)/dble(rr)*(etarhosum(n,r*rr)) + &
                        dble(rr-rpert)/dble(rr)*(etarhosum(n,(r+1)*rr))
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
                   call compute_etarhopert_2d(lo,hi,efp(:,:,1,1),etarhosum_proc(1,:),rr)
                case (3)
                   call compute_etarhopert_3d(lo,hi,efp(:,:,:,1),etarhosum_proc(1,:),rr)
                end select
             end do
             
             ! gather etarhopert for rhoh
             source_buffer = etarhopert_proc(n,:)
             call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
             etarhopert(n,:) = target_buffer
             
             ! update etasum on faces that do not exist at the coarser level
             ! then recompute eta on these faces
             do r=r_start_coord(n-1),r_end_coord(n-1)
                do rpert=1,rr-1
                   etarhosum(n,r*rr+rpert) = &
                        etarhosum(n,r*rr+rpert) + etarhopert(n,r*rr+rpert)
                   etarho(n,r*rr+rpert) = &
                        etarhosum(n,r*rr+rpert)/dble(ncell(n,r*rr+rpert))
                end do
             end do
             
          end do ! end loop over levels
          
       end if

       call restrict_base(nlevs,etarho,.false.)
       
    else

       call bl_error("ERROR: make_eta should not be called for spherical")

    end if

    deallocate(ncell)
    deallocate(etarhosum_proc,etarhosum)
    deallocate(etarhopert_proc,etarhopert)
    deallocate(source_buffer,target_buffer)

    call destroy(bpt)

  end subroutine make_etarho

  subroutine sum_etarho_coarsest_2d(lo,hi,domhi,etarhoflux,etarhosum)

    integer         , intent(in   ) :: lo(:), hi(:), domhi(:)
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1):,lo(2):)
    real (kind=dp_t), intent(inout) :: etarhosum(0:)

    ! local
    integer :: i,j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          etarhosum(j) = etarhosum(j) + etarhoflux(i,j)
       end do
    end do

    ! we only add the contribution at the top edge if we are at the top of the domain
    ! this prevents double counting
    if(hi(2) .eq. domhi(2)) then
       j=hi(2)+1
       do i=lo(1),hi(1)
          etarhosum(j) = etarhosum(j) + etarhoflux(i,j)
       end do
    end if

  end subroutine sum_etarho_coarsest_2d

  subroutine sum_etarho_coarsest_3d(lo,hi,domhi,etarhoflux,etarhosum)

    integer         , intent(in   ) :: lo(:), hi(:), domhi(:)
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1):,lo(2):,lo(3):)
    real (kind=dp_t), intent(inout) :: etarhosum(0:)

    ! local
    integer :: i,j,k

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etarhosum(k) = etarhosum(k) + etarhoflux(i,j,k)
          end do
       end do
    end do

    ! we only add the contribution at the top edge if we are at the top of the domain
    ! this prevents double counting
    if(hi(3) .eq. domhi(3)) then
       k=hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etarhosum(k) = etarhosum(k) + etarhoflux(i,j,k)
          end do
       end do
    end if

  end subroutine sum_etarho_coarsest_3d

  subroutine compute_etarhopert_2d(lo,hi,etarhoflux,etarhopert,rr)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), rr
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1):,lo(2):)
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

  subroutine compute_etarhopert_3d(lo,hi,etarhoflux,etarhopert,rr)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), rr
    real (kind=dp_t), intent(in   ) :: etarhoflux(lo(1):,lo(2):,lo(3):)
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
                                   dx,normal,etarho,mla)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, r_end_coord, nr
    use variables
    use average_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:), snew(:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:), rho0_new(:,0:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(  out) :: etarho(:,0:)
    type(ml_layout), intent(in   ) :: mla

    type(multifab) :: eta_cart(mla%nlevel)
    
    real(kind=dp_t), pointer :: ep(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:), snp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    integer :: n,i,lo(sold(1)%dim),hi(sold(1)%dim),ng_s
    integer :: r

    real(kind=dp_t), allocatable :: etarho_cc(:,:)

    allocate (etarho_cc(nlevs,0:nr_fine-1))


    ! construct a multifab containing  [ rho' (Utilde . e_r) ]
    ng_s = sold(1)%ng

    do n=1,nlevs

       call multifab_build(eta_cart(n), sold(n)%la, 1, 0)

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
          hi = lwb(get_box(eta_cart(n),i))

          call construct_eta_cart(n, sop(:,:,:,rho_comp), snp(:,:,:,rho_comp), &
                                  ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                  np(:,:,:,:), ep(:,:,:,1), &
                                  rho0_old(n,:), rho0_new(n,:), &
                                  dx(n,:), lo, hi, ng_s)
          
       enddo

    enddo
    
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

  subroutine construct_eta_cart(n, rho_old, rho_new, &
                                umac, vmac, wmac, &
                                normal, eta_cart, &
                                rho0_old, rho0_new, dx, lo, hi, ng)

    use bl_constants_module
    use geometry, only: nr_fine
    use fill_3d_module

    integer        , intent(in   ) :: n,lo(:),hi(:),ng
    real(kind=dp_t), intent(in   ) ::  rho_old(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind=dp_t), intent(in   ) ::  rho_new(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind=dp_t), intent(in   ) ::     umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)    
    real(kind=dp_t), intent(in   ) ::     vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)    
    real(kind=dp_t), intent(in   ) ::     wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)    
    real(kind=dp_t), intent(in   ) ::   normal(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real(kind=dp_t), intent(inout) :: eta_cart(lo(1):,lo(2):,lo(3):)
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
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,rho0_nph,rho0_cart,lo,hi,dx,0)


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

    deallocate (rho0_cart)

  end subroutine construct_eta_cart

end module make_eta_module
