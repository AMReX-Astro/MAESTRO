
module horizontal_average_module
  
  use bl_types        , only: dp_t
  use box_module      , only: lwb, upb, box
  use multifab_module , only: multifab, multifab_remote, dataptr, &
                              get_layout, get_box, nghost, nboxes, get_dim
  use layout_module   , only: get_pd
  use ml_layout_module, only: ml_layout
  use parallel        , only: parallel_reduce, MPI_SUM
  use bl_error_module , only: bl_error

  implicit none

  private

  public :: make_horizontalaverage, average_irreg, average_one_level

contains

  subroutine make_horizontalaverage(plot_file_name,mla,s,u,rho0,div_coeff,w0,dx,time)

    use geometry, only: nr_fine, nr_irreg, r_start_coord, r_end_coord, spherical, &
                        numdisjointchunks, dr
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use bl_constants_module, only: ZERO, HALF
    use restrict_base_module, only: restrict_base, fill_ghost_base
    use probin_module, only: max_levs, drdxfac, prob_lo
    use geometry, only: nr
    use bl_IO_module, only: unit_new
    use parallel     , only: parallel_IOProcessorNode, parallel_IOProcessor

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: u(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: time
    character(len=*),intent(in   ) :: plot_file_name

    ! local
    real(kind=dp_t), pointer     :: pu(:,:,:,:)
    real(kind=dp_t), pointer     :: ps(:,:,:,:)
    logical, pointer             :: mp(:,:,:,:)

    type(box)                    :: domain
    integer                      :: domlo(mla%dim),domhi(mla%dim),lo(mla%dim),hi(mla%dim)
    integer                      :: max_rcoord(mla%nlevel),rcoord(mla%nlevel)
    integer                      :: stencil_coord(mla%nlevel)
    integer                      :: dm,nlevs,i,j,r,n,ng_u,ng_s,min_all,min_lev
    real(kind=dp_t)              :: radius
    integer                      :: un

    character(len=256)           :: filename
    real(kind=dp_t)              :: rloc

    integer, allocatable ::  ncell_proc(:,:)
    integer, allocatable ::       ncell(:,:)

    integer, allocatable :: which_lev(:)

    real(kind=dp_t), allocatable :: kesum_proc(:,:), pesum_proc(:,:)
    real(kind=dp_t), allocatable ::      kesum(:,:), pesum(:,:)
    real(kind=dp_t), allocatable :: kebar(:,:), pebar(:,:)

    real(kind=dp_t), allocatable :: radii(:,:)

    logical, save :: firstCall_io = .true.

    type(bl_prof_timer), save :: bpt

    logical :: limit

    call build(bpt, "horizontal_average")

    dm = mla%dim
    nlevs = mla%nlevel

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! NOTE: The indices for ncell_proc, etc., are switched in this subroutine
    !       to have the number of levels second because due to memory ordering, 
    !       it will run much faster for larger problems
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (spherical .eq. 1) then

    else

       allocate(ncell_proc (0:nr_fine-1,nlevs))
       allocate(ncell      (0:nr_fine-1,nlevs))
       allocate(kesum_proc(0:nr_fine-1,nlevs))
       allocate(kesum     (0:nr_fine-1,nlevs))
       allocate(pesum_proc(0:nr_fine-1,nlevs))
       allocate(pesum     (0:nr_fine-1,nlevs))
       allocate(kebar     (nlevs,0:nr_fine-1))
       allocate(pebar     (nlevs,0:nr_fine-1))

    end if

    ng_s = nghost(s(1))
    ng_u = nghost(u(1))

    kebar        = ZERO
    pebar        = ZERO
    ncell        = 0
    ncell_proc   = 0
    kesum        = ZERO     
    pesum        = ZERO 
    kesum_proc   = ZERO
    pesum_proc   = ZERO

    if (spherical .eq. 0) then
       
       ! The plane-parallel case is straightforward.
       ! Simply average all the cells values at a particular height.

       do n=1,nlevs

          domain = get_pd(get_layout(u(n)))
          domlo  = lwb(domain)
          domhi  = upb(domain)

          if (dm .eq. 1) then
             ncell(:,n) = 1
          else if (dm .eq. 2) then
             ncell(:,n) = domhi(1)-domlo(1)+1
          else if (dm .eq. 3) then
             ncell(:,n) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
          end if

          do i=1, nboxes(u(n))
             if ( multifab_remote(u(n), i) ) cycle
             pu => dataptr(u(n), i)
             ps => dataptr(s(n), i)
             lo =  lwb(get_box(u(n), i))
             hi =  upb(get_box(u(n), i))
             select case (dm)
            ! case (1)
            !    call sum_phi_1d(pp(:,1,1,:),phisum_proc(:,n),lo,hi,ng,incomp)
             case (2)
                call sum_energy_2d(ps(:,:,1,:),pu(:,:,1,:),rho0(n,:),div_coeff(n,:), &
                                   w0(n,:),dx(n,:),kesum_proc(:,n),pesum_proc(:,n), &
                                   lo,hi,ng_u,ng_s)
            ! case (3)
            !    call sum_phi_3d(pp(:,:,:,:),phisum_proc(:,n),lo,hi,ng,incomp)
             end select
          end do

          call parallel_reduce(kesum(:,n), kesum_proc(:,n), MPI_SUM, &
                               proc = parallel_IOProcessorNode())
          call parallel_reduce(pesum(:,n), pesum_proc(:,n), MPI_SUM, &
                               proc = parallel_IOProcessorNode())

          ! compute ebar by normalizing esum
          if (parallel_IOProcessor()) then
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n,i),r_end_coord(n,i)
                   kebar(n,r) = kesum(r,n) / dble(ncell(r,n))
                   pebar(n,r) = pesum(r,n) / dble(ncell(r,n))
                end do
             end do
          end if

       end do

       call restrict_base(kebar,.true.)
       call restrict_base(pebar,.true.)

    else if(spherical .eq. 1) then

   end if

   !---------------------------------------------------------------
   ! writing the data
   !---------------------------------------------------------------

 999 format("# job name: ",a)
1000 format(1x,3(g24.10,1x))
1001 format("# time: ",(g24.10,1x))
1002 format("#",3(a24,1x))

   if (parallel_IOProcessor()) then
      un = unit_new()
      filename = "horizontal_averages/" // plot_file_name
      open(unit=un, file=filename,status="new")

      write (un, *) " "
      write (un, 1001) time
      write (un, 1002) "location","KE","PE"

      do r=0,nr(1)-1
         rloc=prob_lo(2)+(dble(r)+HALF)*dr(1)
         write (un, 1000) rloc, kebar(1,r), pebar(1,r)
      end do

      close(un)
   end if

   call destroy(bpt)

 contains

   subroutine cubic_interp(x,x0,x1,x2,x3,y,y0,y1,y2,y3)

     real(kind=dp_t), intent(in   ) :: x,x0,x1,x2,x3,y0,y1,y2,y3
     real(kind=dp_t), intent(  out) :: y

     y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
          + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1) &
          + ( ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0) &
          -((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1))/(x3-x1) ) / (x3-x0) &
          *(x-x0)*(x-x1)*(x-x2)

     if (y .gt. max(y0,y1,y2,y3)) y = max(y0,y1,y2,y3)
     if (y .lt. min(y0,y1,y2,y3)) y = min(y0,y1,y2,y3)

   end subroutine cubic_interp

   subroutine quad_interp(x,x0,x1,x2,y,y0,y1,y2,limit)

     real(kind=dp_t), intent(in   ) :: x,x0,x1,x2,y0,y1,y2
     real(kind=dp_t), intent(  out) :: y
     logical,         intent(in   ) :: limit

     y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
          + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1)


     if (limit) then
        if (y .gt. max(y0,y1,y2)) y = max(y0,y1,y2)
        if (y .lt. min(y0,y1,y2)) y = min(y0,y1,y2)
     end if

   end subroutine quad_interp

   subroutine lin_interp(x,x0,x1,y,y0,y1)

     real(kind=dp_t), intent(in   ) :: x,x0,x1,y0,y1
     real(kind=dp_t), intent(  out) :: y

     y = y0 + (y1-y0)/(x1-x0)*(x-x0)

   end subroutine lin_interp

 end subroutine make_horizontalaverage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_irreg(mla,phi,phibar_irreg,dx,incomp)

    use geometry, only: nr_fine, nr_irreg, r_start_coord, r_end_coord, spherical, &
         numdisjointchunks, dr
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use bl_constants_module, only: ZERO
    use restrict_base_module, only: restrict_base, fill_ghost_base

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: incomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar_irreg(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    real(kind=dp_t), pointer     :: pp(:,:,:,:)
    logical, pointer             :: mp(:,:,:,:)

    integer                      :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer                      :: i,r,n,ng

    integer, allocatable ::  ncell_proc(:,:)
    integer, allocatable ::       ncell(:,:)

    real(kind=dp_t), allocatable :: phisum_proc(:,:)
    real(kind=dp_t), allocatable ::      phisum(:,:)
    real(kind=dp_t), allocatable :: radii(:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average_irreg")

    dm = mla%dim
    nlevs = mla%nlevel

    if (spherical .eq. 0) then
       call bl_error("average_irreg only written for spherical")
    end if
    
    allocate(ncell_proc ( 0:nr_irreg  ,nlevs))
    allocate(ncell      (-1:nr_irreg  ,nlevs))
    allocate(phisum_proc( 0:nr_irreg  ,nlevs))
    allocate(phisum     (-1:nr_irreg  ,nlevs))
    allocate(radii      (-1:nr_irreg+1,nlevs))
    !
    ! radii contains every possible distance that a cell-center at the finest
    ! level can map into
    !
    !$OMP PARALLEL PRIVATE(r,n)
    do n=1,nlevs
       !$OMP DO
       do r=0,nr_irreg
          radii(r,n) = sqrt(0.75d0+2.d0*r)*dx(n,1)
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

    radii(nr_irreg+1,:) = 1.d99

    ng = nghost(phi(1))

    phibar_irreg = ZERO
    ncell        = 0
    ncell_proc   = 0
    phisum       = ZERO       
    phisum_proc  = ZERO
    !
    ! For spherical, we construct a 1D array, phisum, that has space
    ! allocated for every possible radius that a cell-center at the finest
    ! level can map into.  The radius locations have been precomputed and stored
    ! in radii(:).
    
    ! For cells at the non-finest level, map a weighted contribution into the nearest
    ! bin in phisum.
    !
    do n=nlevs,1,-1

       do i=1, nboxes(phi(n))
          if ( multifab_remote(phi(n), i) ) cycle
          pp => dataptr(phi(n), i)
          lo =  lwb(get_box(phi(n), i))
          hi =  upb(get_box(phi(n), i))

          if (n .eq. nlevs) then
             call sum_phi_3d_sphr(radii(0:,n),nr_irreg,pp(:,:,:,:),phisum_proc(:,n), &
                                  lo,hi,ng,dx(n,:),ncell_proc(:,n),incomp)
          else
             ! we include the mask so we don't double count; i.e., we only consider
             ! cells that we can "see" when constructing the sum
             mp => dataptr(mla%mask(n), i)
             call sum_phi_3d_sphr(radii(0:,n),nr_irreg,pp(:,:,:,:),phisum_proc(:,n), &
                                  lo,hi,ng,dx(n,:),ncell_proc(:,n),incomp, &
                                  mp(:,:,:,1))
          end if
       end do

       call parallel_reduce(ncell (0:,n), ncell_proc (:,n), MPI_SUM)
       call parallel_reduce(phisum(0:,n), phisum_proc(:,n), MPI_SUM)
    end do

    ! compute phibar_irreg
    do n=1,nlevs
       do r=0,nr_irreg
          if (ncell(r,n) .ne. 0.d0) then
             phibar_irreg(n,r) = phisum(r,n) / dble(ncell(r,n))
          end if
       end do
    end do

   call destroy(bpt)

 end subroutine average_irreg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_1d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i

    do i=lo(1),hi(1)
       phisum(i) = phi(i,incomp)
    end do

  end subroutine sum_phi_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_energy_2d(s,u,rhob,div_coeff,w0,dx,kesum,pesum,lo,hi,ng_u,ng_s)

    use variables, only: rho_comp
    use bl_constants_module, only: HALF, ZERO
    use probin_module, only: prob_lo, grav_const, rho0, nn, K
    use geometry, only: nr, dr

    integer         , intent(in   ) :: lo(:), hi(:), ng_u, ng_s
    real (kind=dp_t), intent(in   ) :: u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rhob(0:), div_coeff(0:), w0(0:)
    real (kind=dp_t), intent(inout) :: kesum(0:), pesum(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !***** local variables
    real (kind=dp_t)  :: velsqr, rho1, y
    real (kind=dp_t)  :: N2(lo(2):hi(2)), rho_over_beta(lo(2):hi(2))
    real (kind=dp_t)  :: entropy(lo(2):hi(2))
    integer :: i,j

    integer           :: r
    real (kind=dp_t)  :: rloc,rmax

    ! calculate N2

    do j = lo(2), hi(2)
       rho_over_beta(j) = rhob(j) / div_coeff(j)
       entropy(j) = log(rho_over_beta(j))
    enddo

    do j = lo(2), hi(2)-1
       N2(j) = grav_const*(entropy(j+1)-entropy(j))/dx(2)
    enddo
    N2(hi(2)) = N2(hi(2)-1)

    do j = lo(2),hi(2)
       y = prob_lo(2) + (dble(j) + HALF) * dx(2)
       do i = lo(1),hi(1)
          rho1   = s(i,j,rho_comp) - rhob(j)
          velsqr = u(i,j,1)**2 + u(i,j,2)**2

          kesum(j) = kesum(j) + HALF * rhob(j) * velsqr
          pesum(j) = pesum(j) + HALF * grav_const**2 * rhob(j) / N2(j) &
                     * (rho1/rhob(j))**2
       end do
    end do

  end subroutine sum_energy_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_3d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phisum(k) = phisum(k) + phi(i,j,k,incomp)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine sum_phi_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_3d_sphr(radii,nr_irreg,phi,phisum,lo,hi,ng,dx,ncell,incomp,mask)

    use geometry, only: dr, center
    use probin_module, only: prob_lo
    use bl_constants_module, only: HALF

    integer         , intent(in   )           :: lo(:), hi(:), ng, incomp, nr_irreg
    real (kind=dp_t), intent(in   )           :: radii(0:)
    real (kind=dp_t), intent(in   )           :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout)           :: phisum(0:)
    real (kind=dp_t), intent(in   )           :: dx(:)
    integer         , intent(inout)           :: ncell(0:)
    logical         , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    ! local
    real (kind=dp_t) :: x, y, z, radius
    integer          :: i, j, k, index
    logical          :: cell_valid

    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k) + HALF)*dx(3) - center(3)
       
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j) + HALF)*dx(2) - center(2)
          
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i) + HALF)*dx(1) - center(1)
             
             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if ( present(mask) ) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             end if
             
             if (cell_valid) then

                ! compute distance to center
                radius = sqrt(x**2 + y**2 + z**2)
                
                ! figure out which radii index this point maps into
                index = ((radius / dx(1))**2 - 0.75d0) / 2.d0
                
                ! due to roundoff error, need to ensure that we are in the proper radial bin
                if (index .lt. nr_irreg) then
                   if (abs(radius-radii(index)) .gt. abs(radius-radii(index+1))) then
                      index = index+1
                   end if
                end if
                
                phisum(index) = phisum(index) + phi(i,j,k,incomp)
                ncell(index)  = ncell(index) + 1

             end if
             
          end do
       end do
    end do

  end subroutine sum_phi_3d_sphr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_one_level(n,phi,phibar,incomp)

    use geometry, only: nr_fine, nr, spherical
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use bl_constants_module, only: ZERO

    integer        , intent(in   ) :: n,incomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar(:,0:)

    ! local
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    type(box)                :: domain
    integer                  :: domlo(get_dim(phi(1))),domhi(get_dim(phi(1)))
    integer                  :: lo(get_dim(phi(1))),hi(get_dim(phi(1)))
    integer                  :: i,r,ng,ncell,dm

    real(kind=dp_t) ::   phisum_proc(0:nr_fine-1)
    real(kind=dp_t) ::        phisum(0:nr_fine-1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average_one_level")

    dm = get_dim(phi(1))
    ng = nghost(phi(1))

    phibar = ZERO

    ncell        = 0
    phisum       = ZERO       
    phisum_proc  = ZERO

    if (spherical .eq. 0) then
       
       domain = get_pd(get_layout(phi(n)))
       domlo  = lwb(domain)
       domhi  = upb(domain)

       if (dm .eq. 1) then
          ncell = 1
       else if (dm .eq. 2) then
          ncell = domhi(1)-domlo(1)+1
       else if (dm .eq. 3) then
          ncell = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
       end if

       do i=1, nboxes(phi(n))
          if ( multifab_remote(phi(n), i) ) cycle
          pp => dataptr(phi(n), i)
          lo =  lwb(get_box(phi(n), i))
          hi =  upb(get_box(phi(n), i))
          select case (dm)
          case (1)
             call sum_phi_1d(pp(:,1,1,:),phisum_proc,lo,hi,ng,incomp)
          case (2)
!             call sum_phi_2d(pp(:,:,1,:),phisum_proc,lo,hi,ng,incomp)
          case (3)
             call sum_phi_3d(pp(:,:,:,:),phisum_proc,lo,hi,ng,incomp)
          end select
       end do

       call parallel_reduce(phisum, phisum_proc, MPI_SUM)

       do r=0,nr(n)-1
          phibar(n,r) = phisum(r) / dble(ncell)
       end do

    else if(spherical .eq. 1) then

       call bl_error("average_one_level not written for multilevel spherical")

    end if

    call destroy(bpt)

  end subroutine average_one_level

end module horizontal_average_module
