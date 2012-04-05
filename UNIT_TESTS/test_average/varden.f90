subroutine varden()

  use BoxLib
  use f2kcli
  use list_box_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use base_state_module
  use ml_restriction_module
  use bc_module
  use define_bc_module
  use bl_space, ONLY: MAX_SPACEDIM
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use variables
  use geometry
  use network
  use average_module
  use eos_module
  use fill_3d_module
  use probin_module
  use runtime_init_module
  use bl_constants_module
  use initialize_module
  use multifab_physbc_module
  use multifab_fill_ghost_module

  implicit none

  real(dp_t) :: dist,lenx,leny,lenz,max_dist

  integer :: i,n,r
  integer :: dm, nlevs

  type(ml_layout) :: mla

  type(multifab), allocatable :: phi(:)

  real(kind=dp_t), pointer :: nrp(:,:,:,:)
  real(kind=dp_t), pointer :: pp(:,:,:,:)
  real(kind=dp_t), pointer :: dx(:,:)

  integer, allocatable :: lo(:),hi(:)

  type(box) :: domain
  integer   :: domhi(MAX_SPACEDIM)

  type(layout)      :: la
  type(ml_boxarray) :: mba

  real(dp_t), allocatable :: phi_exact(:,:)
  real(dp_t), allocatable :: phi_avg(:,:)

  type(bc_tower) ::  the_bc_tower

  call runtime_init()
  call init_spherical()
  call init_center()

  call init_variables()

  call network_init()
  call eos_init(use_eos_coulomb=use_eos_coulomb)

  call read_a_hgproj_grid(mba, test_set)

  call ml_layout_build(mla,mba,pmask)

  ! check for proper nesting
  if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
     call bl_error('fixed_grids not properly nested')
  end if

  ! initialize nlevs
  nlevs = mla%nlevel

  if (nlevs .ne. max_levs) then
     call bl_error('varden.f90: nlevs .ne. max_levs not supported yet')
  end if

  nlevs_radial = merge(1, nlevs, spherical .eq. 1)
  if ( parallel_IOProcessor() ) &
       print *, 'nlevs = ', nlevs

  ! initialize dm
  dm = mla%dim

  ! initialize boundary conditions
  call initialize_bc(the_bc_tower,nlevs,pmask)
  do n = 1,nlevs
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! allocate states
  allocate(phi(nlevs))

  ! build states
  do n = 1,nlevs
     call multifab_build(phi(n),mla%la(n),1,1)
     call setval(phi(n),0.0_dp_t,all=.true.)
  end do

  ! initialize_dx
  call initialize_dx(dx,mba,nlevs)

  ! now that we have dx we can initialize nr_fine and dr_fine
  if (spherical .eq. 1) then
     
     dr_fine = dx(nlevs,1) / dble(drdxfac)
     
     if (.not. octant) then
        lenx = HALF * (prob_hi(1) - prob_lo(1))
        leny = HALF * (prob_hi(2) - prob_lo(2))
        lenz = HALF * (prob_hi(3) - prob_lo(3))
     else
        lenx = prob_hi(1) - prob_lo(1)
        leny = prob_hi(2) - prob_lo(2)
        lenz = prob_hi(3) - prob_lo(3)
     end if
     
     max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
     nr_fine = int(max_dist / dr_fine) + 1

     ! compute nr_irreg
     domain = layout_get_pd(phi(nlevs)%la)
     domhi  = upb(domain)+1
     if (.not. octant) then
        nr_irreg = (3*(domhi(1)/2-0.5d0)**2-0.75d0)/2.d0
     else
        nr_irreg = (3*(domhi(1)-0.5d0)**2-0.75d0)/2.d0
     endif
     
  else
     
     nr_fine = extent(mla%mba%pd(nlevs),dm)
     dr_fine = (prob_hi(dm)-prob_lo(dm)) / dble(nr_fine)
     
  end if

  ! create numdisjointchunks, r_start_coord, r_end_coord
  call init_multilevel(phi)

  ! now that we have nr_fine and dr_fine we can create nr, dr, r_cc_loc, r_edge_loc
  call init_radial(nlevs,mba)

  allocate( phi_exact(nlevs,0:nr_fine-1))
  allocate( phi_avg  (nlevs,0:nr_fine-1))

  phi_exact(:,:) = ZERO
  phi_avg  (:,:) = ZERO

  do r=0,nr(1)-1
     dist = (dble(r)+HALF)*dr(1)
     phi_exact(1,r) = exp(-dist**2/0.1d0)
  end do

  la = mla%la(1)

  allocate(lo(dm))
  allocate(hi(dm))

  do n=1,nlevs
     do i = 1, phi(n)%nboxes
        if ( multifab_remote(phi(n),i) ) cycle
        pp => dataptr(phi(n),i)
        lo =  lwb(get_box(phi(n),i))
        hi =  upb(get_box(phi(n),i))
        select case (dm)
        case (2)
        case (3)
           if (spherical .eq. 1) then
              call initgaussian_3d_sphr(phi_exact(1,:), pp(:,:,:,:), phi(n)%ng, &
                                        lo, hi, dx(n,:))
           else
           end if
        end select
     end do
  enddo

  if (dump_phi_plotfile) then
     call fabio_ml_multifab_write_d(phi,mla%mba%rr(:,1),"phi_pltfile")
  endif


  if (nlevs .eq. 1) then

     ! fill ghost cells for two adjacent grids at the same level
     ! this includes periodic domain boundary ghost cells
     call multifab_fill_boundary(phi(nlevs))

     ! fill non-periodic domain boundary ghost cells
     call multifab_physbc(phi(nlevs),1,foextrap_comp,1,the_bc_tower%bc_tower_array(nlevs))

  else

     ! the loop over nlevs must count backwards to make sure the finer grids are done first
     do n=nlevs,2,-1

        ! set level n-1 data to be the average of the level n data covering it
        call ml_cc_restriction(phi(n-1),phi(n),mla%mba%rr(n-1,:))

        ! fill level n ghost cells using interpolation from level n-1 data
        ! note that multifab_fill_boundary and multifab_physbc are called for
        ! both levels n-1 and n
        call multifab_fill_ghost_cells(phi(n),phi(n-1),phi(n)%ng,mla%mba%rr(n-1,:), &
                                       the_bc_tower%bc_tower_array(n-1), &
                                       the_bc_tower%bc_tower_array(n), &
                                       1,foextrap_comp,1,fill_crse_input=.false.)

     enddo

  end if

  ! now that we are initialized, try averaging the state to 1-d
  ! and compare to the base state
  if ( parallel_IOProcessor() ) &
       print *, 'averaging...'

  call average(mla,phi,phi_avg,dx,1)

  if ( parallel_IOProcessor() ) &
       print *, 'done'

  ! compute the error against the base state
  if ( parallel_IOProcessor() ) then
     open (unit=10, file="phi.error")
     write (10,*) "r_cc_loc, phi_exact, phi_avg, phi_exact-phi_avg, (phi_exact-phi_avg)/phi_exact"
     do n=1,nlevs_radial
        do r=0,nr(n)-1
           write (10,1000) r_cc_loc(n,r), phi_exact(n,r), phi_avg(n,r), &
                phi_exact(n,r)-phi_avg(n,r), &
                (phi_exact(n,r)-phi_avg(n,r))/phi_exact(n,r)
        enddo
     enddo
     close (10)
  endif

1000 format(1x,6(g24.16))

  do n = 1,nlevs
     call destroy(phi(n))
  end do

  call destroy(mla)
  call destroy(mba)

  deallocate(phi)
  deallocate(phi_exact,phi_avg,lo,hi)

  call bc_tower_destroy(the_bc_tower)

  call runtime_close()

  deallocate(dr,r_cc_loc,r_edge_loc,r_start_coord,r_end_coord,nr,numdisjointchunks,dx)

  contains

    subroutine initgaussian_3d_sphr(phi_exact,phi,ng,lo,hi,dx)

      use probin_module, only: prob_lo, perturb_model

      real (kind = dp_t), intent(in   ) :: phi_exact(:)
      integer           , intent(in   ) :: lo(:),hi(:),ng
      real (kind = dp_t), intent(inout) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) :: dx(:)

      !     Local variables
      integer         :: i,j,k
      real(kind=dp_t) :: x,y,z,dist

      ! initialize (rho h) using the EOS
!      do k = lo(3), hi(3)
!         z = (dble(k)+0.5d0)*dx(3) - center(3)
!         do j = lo(2), hi(2)
!            y = (dble(j)+0.5d0)*dx(2) - center(2)
!            do i = lo(1), hi(1)
!               x = (dble(i)+0.5d0)*dx(1) - center(1)
!
!               dist = sqrt(x**2 + y**2 + z**2)
!               phi(i,j,k,1) = exp(-dist**2/0.1d0)
!
!            enddo
!         enddo
!      enddo

      call put_1d_array_on_cart_3d_sphr(.false.,.false.,phi_exact(:), &
                                        phi(:,:,:,1:),lo,hi,dx,ng)

    end subroutine initgaussian_3d_sphr

end subroutine varden
