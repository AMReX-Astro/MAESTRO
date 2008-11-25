! given a multifab of data (phi), average down to a base state quantity,
! phibar.  If we are in plane-parallel, the averaging is at constant
! height.  If we are spherical, then the averaging is done at constant
! radius.  

module average_module
  
  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private
  public :: average, average_one_level

contains

  subroutine average(mla,phi,phibar,dx,incomp)

    use geometry, only: nr_fine, r_start_coord, r_end_coord, spherical, numdisjointchunks, &
         dm, nlevs
    use bl_prof_module
    use bl_constants_module
    use restrict_base_module
    use probin_module, only: drdxfac

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: incomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    real(kind=dp_t), pointer     :: pp(:,:,:,:)
    logical, pointer             :: mp(:,:,:,:)

    type(box)                    :: domain
    integer                      :: domlo(dm),domhi(dm),lo(dm),hi(dm)
    integer                      :: i,j,r,n,ng,rr,nr_crse
    real(kind=dp_t)              :: w_lo, w_hi, del_w, wsix, theta

    real(kind=dp_t) ::   ncell_grid(nlevs,0:nr_fine-1)
    real(kind=dp_t) ::   ncell_proc(nlevs,0:nr_fine-1)
    real(kind=dp_t) ::        ncell(nlevs,0:nr_fine-1)
    real(kind=dp_t) ::  phisum_proc(nlevs,0:nr_fine-1)
    real(kind=dp_t) ::       phisum(nlevs,0:nr_fine-1)
    real(kind=dp_t) :: phipert_proc(nlevs,0:nr_fine-1)
    real(kind=dp_t) ::      phipert(nlevs,0:nr_fine-1)

    real(kind=dp_t) :: source_buffer(0:nr_fine-1)
    real(kind=dp_t) :: target_buffer(0:nr_fine-1)

    real(kind=dp_t), allocatable :: ncell_crse(:)
    real(kind=dp_t), allocatable :: phibar_crse(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average")

    ng = phi(1)%ng

    phibar = ZERO

    ncell        = ZERO
    ncell_proc   = ZERO
    phisum       = ZERO       
    phisum_proc  = ZERO
    phipert      = ZERO
    phipert_proc = ZERO

    rr = 2

    if (spherical .eq. 0) then
       
       do n=1,nlevs

          domain = layout_get_pd(phi(n)%la)
          domlo  = lwb(domain)
          domhi  = upb(domain)

          if (dm .eq. 2) then
             ncell(n,:) = domhi(1)-domlo(1)+1
          else if (dm .eq. 3) then
             ncell(n,:) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
          end if

          do i=1,phi(n)%nboxes
             if ( multifab_remote(phi(n), i) ) cycle
             pp => dataptr(phi(n), i)
             lo =  lwb(get_box(phi(n), i))
             hi =  upb(get_box(phi(n), i))
             select case (dm)
             case (2)
                call sum_phi_2d(pp(:,:,1,:),phisum_proc(n,:),lo,hi,ng,incomp)
             case (3)
                call sum_phi_3d(pp(:,:,:,:),phisum_proc(n,:),lo,hi,ng,incomp)
             end select
          end do

          ! gather phisum
          source_buffer = phisum_proc(n,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          phisum(n,:) = target_buffer
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)
                phibar(n,r) = phisum(n,r) / dble(ncell(n,r))
             end do
          end do

       end do

       call fill_ghost_base(phibar,.true.)
       call restrict_base(phibar,.true.)

    else if(spherical .eq. 1) then

       ! The spherical case is tricky because the base state only exists at 
       ! one level with dr = dx(nlevs) / drdxfac.

       ! Therefore, for each level, we compute the contribution to each radial bin
       ! Then we sum the contributions to a coarse radial bin,
       ! then interpolate to fill the fine radial bin.

       ! phisum(nlevs,:,:) will be the volume weighted sum over all levels.
       ! ncell(nlevs,:) will be the volume weighted number of cells over 
       ! all levels.

       ! We make sure to use mla%mask to not double count cells, i.e.,
       ! we only sum up cells that are not covered by finer cells.
       ! we use the convention that a cell volume of 1 corresponds to 
       ! dx(n=1)**3

       ! First we compute ncell(nlevs,:) and phisum(nlevs,:,:) as if the 
       ! finest level were the only level in existence.
       ! Then, we add contributions from each coarser cell that is not 
       ! covered by a finer cell.

       do n=nlevs,1,-1

          do i=1,phi(n)%nboxes
             if ( multifab_remote(phi(n), i) ) cycle
             pp => dataptr(phi(n), i)
             lo =  lwb(get_box(phi(n), i))
             hi =  upb(get_box(phi(n), i))
             ncell_grid(n,:) = ZERO
             if (n .eq. nlevs) then
                call sum_phi_3d_sphr(n,pp(:,:,:,:),phisum_proc(n,:), &
                                     lo,hi,ng,dx(n,:),ncell_grid(n,:),incomp)
             else
                mp => dataptr(mla%mask(n), i)
                call sum_phi_3d_sphr(n,pp(:,:,:,:),phisum_proc(n,:), &
                                     lo,hi,ng,dx(n,:),ncell_grid(n,:),incomp, &
                                     mp(:,:,:,1))
             end if

             ncell_proc(n,:) = ncell_proc(n,:) + ncell_grid(n,:)
          end do

          call parallel_reduce(ncell(n,:), ncell_proc(n,:), MPI_SUM)

          source_buffer = phisum_proc(n,:)
          call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
          phisum(n,:) = target_buffer

          if (n .ne. nlevs) then
             ncell(nlevs,:) = ncell(nlevs,:) + ncell(n,:)
             do r=0,nr_fine-1
                phisum(nlevs,r) = phisum(nlevs,r) + phisum(n,r)
             end do
          end if

       end do

       if (drdxfac .ne. 1) then

          if ( mod(nr_fine,drdxfac) .eq. 0 ) then
             nr_crse = nr_fine / drdxfac
          else
             nr_crse = nr_fine / drdxfac + 1
          end if

          allocate( ncell_crse(0:nr_crse-1))
          allocate(phibar_crse(-2:nr_crse+1))

          do r = 0, nr_crse-1

             phibar_crse(r) = 0.d0
              ncell_crse(r) = 0.d0
   
             ! Sum fine data onto the crse grid
             do j = drdxfac*r, min(drdxfac*r+(drdxfac-1),nr_fine-1)
                phibar_crse(r) = phibar_crse(r) + phisum(nlevs,j)
                 ncell_crse(r) =  ncell_crse(r) +  ncell(nlevs,j)
             end do

             ! Now compute the average
             if (ncell_crse(r) .gt. ZERO) then
                phibar_crse(r) = phibar_crse(r) / ncell_crse(r)
             else if (r .eq. nr_crse-3) then
                phibar_crse(r) = phibar_crse(nr_crse-4)
             else if (r .eq. nr_crse-2) then
                phibar_crse(r) = phibar_crse(nr_crse-3)
             else if (r .eq. nr_crse-1) then
                phibar_crse(r) = phibar_crse(nr_crse-2)
             else 
                call bl_error("ERROR: ncell_crse = 0 in average")
             end if
          end do

          ! Reflect (even) across origin
          phibar_crse(-1) = phibar_crse(0)
          phibar_crse(-2) = phibar_crse(1)

          ! Extend at high r
          phibar_crse(nr_crse  ) = phibar_crse(nr_crse-1)
          phibar_crse(nr_crse+1) = phibar_crse(nr_crse-1)

          ! Put the average back onto the fine grid
          do r = 0, nr_crse-1
   
             w_lo = ( 7.d0 * (phibar_crse(r  ) + phibar_crse(r-1)) &
                     -1.d0 * (phibar_crse(r+1) + phibar_crse(r-2)) ) / 12.d0
             w_hi = ( 7.d0 * (phibar_crse(r  ) + phibar_crse(r+1)) &
                     -1.d0 * (phibar_crse(r-1) + phibar_crse(r+2)) ) / 12.d0

             del_w = w_hi - w_lo

             wsix = 6.d0 * ( phibar_crse(r) - 0.5d0 * (w_lo + w_hi) )

             do j = 0, min(drdxfac-1,nr_fine-drdxfac*r-1)
                ! piecewise constant
                ! phibar(nlevs,drdxfac*r+j) = phibar_crse(r)
   
                ! parabolic interpolation
                theta = dble(j) / dble(drdxfac)
                phibar(1,drdxfac*r+j) = w_lo + theta * del_w + &
                                       theta * (1.d0 - theta) * wsix

             end do

          end do

       else 

          ! if drdxfac = 1 then divide the total phisum by the number of cells to get phibar
          do r=0,nr_fine-1
             if (ncell(nlevs,r) .gt. ZERO) then
                phibar(1,r) = phisum(nlevs,r) / ncell(nlevs,r)
             else if (r .eq. nr_fine-1 .and. ncell(nlevs,r) .eq. ZERO) then
                phibar(1,r) = phibar(nlevs,r-1)
             else
                call bl_error("ERROR: ncell_crse = 0 in average")
             end if
          end do

       end if
       
       ! temporary hack for the case where the outermost radial bin average 
       ! to zero because there is no contribution from any Cartesian cell 
       ! that lies in this bin.  This needs to be addressed - perhaps in the 
       ! definition of nr_fine in varden.f90 for spherical problems.
       if (ncell(nlevs,nr_fine-1) .eq. ZERO) then
          phibar(1,nr_fine-1) = phibar(1,nr_fine-2)
       end if

    endif

    call destroy(bpt)

  end subroutine average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_2d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i,j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          phisum(j) = phisum(j) + phi(i,j,incomp)
       end do
    end do

  end subroutine sum_phi_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_3d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i,j,k

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phisum(k) = phisum(k) + phi(i,j,k,incomp)
          end do
       end do
    end do

  end subroutine sum_phi_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_3d_sphr(n,phi,phisum,lo,hi,ng,dx,ncell,incomp,mask)

    use geometry, only: dr, center
    use ml_layout_module
    use bl_constants_module

    integer         , intent(in   )           :: n, lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   )           :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout)           :: phisum(0:)
    real (kind=dp_t), intent(in   )           :: dx(:)
    real (kind=dp_t), intent(inout)           :: ncell(0:)
    logical         , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    ! local
    real (kind=dp_t) :: x, y, z, radius
    integer          :: i, j, k, idx
    real (kind=dp_t) :: cell_weight
    logical          :: cell_valid

    cell_weight = ONE / EIGHT**(n-1)

    do k=lo(3),hi(3)
       z = (dble(k) + HALF)*dx(3) - center(3)

       do j=lo(2),hi(2)
          y = (dble(j) + HALF)*dx(2) - center(2)

          do i=lo(1),hi(1)
             x = (dble(i) + HALF)*dx(1) - center(1)

             cell_valid = .true.
             if ( present(mask) ) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             end if

             if (cell_valid) then
                radius = sqrt(x**2 + y**2 + z**2)
                idx  = int(radius / dr(1))
                
                phisum(idx) = phisum(idx) + cell_weight*phi(i,j,k,incomp)
                ncell(idx) = ncell(idx) + cell_weight
             end if

          end do
       end do
    end do

  end subroutine sum_phi_3d_sphr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_one_level(n,phi,phibar,incomp)

    use geometry, only: nr_fine, nr, spherical, dm
    use bl_prof_module
    use bl_constants_module

    integer        , intent(in   ) :: n,incomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar(:,0:)

    ! local
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    type(box)                :: domain
    integer                  :: domlo(dm),domhi(dm)
    integer                  :: lo(dm),hi(dm)
    integer                  :: i,r,ng,ncell

    real(kind=dp_t) ::   phisum_proc(0:nr_fine-1)
    real(kind=dp_t) ::        phisum(0:nr_fine-1)
    real(kind=dp_t) :: source_buffer(0:nr_fine-1)
    real(kind=dp_t) :: target_buffer(0:nr_fine-1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average_one_level")

    ng = phi(1)%ng

    phibar = ZERO

    ncell        = ZERO
    phisum       = ZERO       
    phisum_proc  = ZERO

    if (spherical .eq. 0) then
       
       domain = layout_get_pd(phi(n)%la)
       domlo  = lwb(domain)
       domhi  = upb(domain)

       if (dm .eq. 2) then
          ncell = domhi(1)-domlo(1)+1
       else if (dm .eq. 3) then
          ncell = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
       end if

       do i=1,phi(n)%nboxes
          if ( multifab_remote(phi(n), i) ) cycle
          pp => dataptr(phi(n), i)
          lo =  lwb(get_box(phi(n), i))
          hi =  upb(get_box(phi(n), i))
          select case (dm)
          case (2)
             call sum_phi_2d(pp(:,:,1,:),phisum_proc,lo,hi,ng,incomp)
          case (3)
             call sum_phi_3d(pp(:,:,:,:),phisum_proc,lo,hi,ng,incomp)
          end select
       end do

       ! gather phisum
       source_buffer = phisum_proc
       call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
       phisum = target_buffer
       do r=0,nr(n)-1
          phibar(n,r) = phisum(r) / dble(ncell)
       end do

    else if(spherical .eq. 1) then

       call bl_error("average_one_level not written for multilevel spherical")

    end if

    call destroy(bpt)

  end subroutine average_one_level

end module average_module
