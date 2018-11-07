module make_pi_cc_module

  use bl_types

  implicit none

contains

  !---------------------------------------------------------------------------
  ! make_pi_cc
  !---------------------------------------------------------------------------
  subroutine make_pi_cc(mla,pi,pi_cc,comp,the_bc_level,beta0)

    ! average the nodal pi to the cell-centers and normalize
    ! it such that it integrates to 0.  Store the result in pi_cc

    use multifab_module
    use define_bc_module
    use ml_layout_module
    use bc_module
    use ml_restrict_fill_module
    use variables , only: foextrap_comp

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: pi(:)
    type(multifab) , intent(inout) :: pi_cc(:)
    integer        , intent(in   ) :: comp
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(in   ) :: beta0(:)

    real(kind=dp_t), pointer :: ppn(:,:,:,:)
    real(kind=dp_t), pointer :: ppc(:,:,:,:)
    real(kind=dp_t), pointer ::  bp(:,:,:,:)
    logical,         pointer ::  mp(:,:,:,:)

    integer :: i,n,ng_pn,ng_pc,ng_b,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t) :: ncell_proc(mla%nlevel), ncell(mla%nlevel)
    real(kind=dp_t) :: pisum_proc(mla%nlevel), pisum(mla%nlevel)

    real(kind=dp_t) :: weight,avg

    dm = mla%dim
    nlevs = mla%nlevel

    ncell      = 0.d0
    pisum      = 0.d0
    ncell_proc = 0.d0
    pisum_proc = 0.d0
    ng_pn      = nghost(pi(1))
    ng_pc      = nghost(pi_cc(1))
    ng_b       = nghost(beta0(1))

    do n=1,nlevs
       weight = 2.d0**(dm*(n-1))
       do i=1,nfabs(pi_cc(n))
          ppn => dataptr(pi(n), i)
          ppc => dataptr(pi_cc(n), i)
          bp  => dataptr(beta0(n), i)
          lo  =  lwb(get_box(pi_cc(n), i))
          hi  =  upb(get_box(pi_cc(n), i))
          select case (dm)
          case (1)
             if (n .eq. nlevs) then
                call make_pi_cc_1d(weight,ppn(:,1,1,1),ng_pn,ppc(:,1,1,comp),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n),bp(:,1,1,1),ng_b)
             else
                mp => dataptr(mla%mask(n), i)
                call make_pi_cc_1d(weight,ppn(:,1,1,1),ng_pn,ppc(:,1,1,comp),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n),bp(:,1,1,1),ng_b,mp(:,1,1,1))
             end if
          case (2)
             if (n .eq. nlevs) then
                call make_pi_cc_2d(weight,ppn(:,:,1,1),ng_pn,ppc(:,:,1,comp),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n),bp(:,:,1,1),ng_b)
             else
                mp => dataptr(mla%mask(n), i)
                call make_pi_cc_2d(weight,ppn(:,:,1,1),ng_pn,ppc(:,:,1,comp),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n),bp(:,:,1,1),ng_b,mp(:,:,1,1))
             end if
          case (3)
             if (n .eq. nlevs) then
                call make_pi_cc_3d(weight,ppn(:,:,:,1),ng_pn,ppc(:,:,:,comp),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n),bp(:,:,:,1),ng_b)
             else
                mp => dataptr(mla%mask(n), i)
                call make_pi_cc_3d(weight,ppn(:,:,:,1),ng_pn,ppc(:,:,:,comp),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n),bp(:,:,:,1),ng_b,mp(:,:,:,1))
             end if
          end select
       end do
    end do

    call parallel_reduce(ncell, ncell_proc, MPI_SUM)
    call parallel_reduce(pisum, pisum_proc, MPI_SUM)

    ! now ncell will contain the total number of cells over all levels
    ! now picum will contain the sum of (volume weighted) pi over all levels
    do n=2,nlevs
       ncell(1) = ncell(1) + ncell(n)
       pisum(1) = pisum(1) + pisum(n)
    end do

    ! divide the sum by the number of cells
    avg = pisum(1)/ncell(1)

    ! if there are no outlet boundary conditions, normalize pi_cc so the
    ! sum over the domain is zero
    if (.not.(any(the_bc_level(1)%phys_bc_level_array(:,:,:) .eq. OUTLET))) then
       do n=1,nlevs
          call multifab_sub_sub_s_c(pi_cc(n),comp,avg,1,ng_pc)
       end do
    end if

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,pi_cc,mla%mba%rr,the_bc_level, &
                                icomp=1, &
                                bcomp=foextrap_comp, &
                                nc=1, &
                                ng=pi_cc(1)%ng, &
                                same_boundary=.false.)

  end subroutine make_pi_cc

  subroutine make_pi_cc_1d(weight,pi,ng_pn,pi_cc,ng_pc,lo,hi,ncell,pisum,beta0,ng_b,mask)

    use probin_module, only: use_alt_energy_fix

    real (kind=dp_t), intent(in   )           :: weight
    integer         , intent(in   )           :: lo(:), hi(:), ng_pn, ng_pc, ng_b
    real (kind=dp_t), intent(in   )           ::    pi(lo(1)-ng_pn:)
    real (kind=dp_t), intent(inout)           :: pi_cc(lo(1)-ng_pc:)
    real (kind=dp_t), intent(in   )           :: beta0(lo(1)-ng_b:)
    real (kind=dp_t), intent(inout)           :: ncell,pisum
    logical         , intent(in   ), optional ::  mask(lo(1):hi(1))

    ! local
    integer :: i

    logical :: cell_valid

    do i=lo(1),hi(1)

       pi_cc(i) = (pi(i) + pi(i+1)) / 2.d0

       if (use_alt_energy_fix) then
          pi_cc(i) = pi_cc(i)*beta0(i)
       endif

       ! make sure the cell isn't covered by finer cells
       cell_valid = .true.
       if (present(mask)) then
          cell_valid = mask(i)
       endif
       
       if (cell_valid) then
          pisum = pisum + weight*pi_cc(i)
          ncell = ncell + weight
       end if
       
    end do

  end subroutine make_pi_cc_1d

  subroutine make_pi_cc_2d(weight,pi,ng_pn,pi_cc,ng_pc,lo,hi,ncell,pisum,beta0,ng_b,mask)

    use probin_module, only: use_alt_energy_fix

    real (kind=dp_t), intent(in   )           :: weight
    integer         , intent(in   )           :: lo(:), hi(:), ng_pn, ng_pc, ng_b
    real (kind=dp_t), intent(in   )           ::    pi(lo(1)-ng_pn:,lo(2)-ng_pn:)
    real (kind=dp_t), intent(inout)           :: pi_cc(lo(1)-ng_pc:,lo(2)-ng_pc:)
    real (kind=dp_t), intent(in   )           :: beta0(lo(1)-ng_b: ,lo(2)-ng_b:)
    real (kind=dp_t), intent(inout)           :: ncell,pisum
    logical         , intent(in   ), optional ::  mask(lo(1):hi(1) ,lo(2):hi(2))

    ! local
    integer :: i,j

    logical :: cell_valid

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          pi_cc(i,j) = (pi(i,j) + pi(i+1,j) + pi(i,j+1) + pi(i+1,j+1)) / 4.d0

          if (use_alt_energy_fix) then
             pi_cc(i,j) = pi_cc(i,j)*beta0(i,j)
          endif

          ! make sure the cell isn't covered by finer cells
          cell_valid = .true.
          if (present(mask)) then
             cell_valid = mask(i,j)
          endif

          if (cell_valid) then
             pisum = pisum + weight*pi_cc(i,j)
             ncell = ncell + weight
          end if
             
       end do
    end do

  end subroutine make_pi_cc_2d

  subroutine make_pi_cc_3d(weight,pi,ng_pn,pi_cc,ng_pc,lo,hi,ncell,pisum,beta0,ng_b,mask)

    use probin_module, only: use_alt_energy_fix

    real(kind=dp_t), intent(in   )           :: weight
    integer        , intent(in   )           :: lo(:), hi(:), ng_pn, ng_pc, ng_b
    real(kind=dp_t), intent(in   )           ::    pi(lo(1)-ng_pn:,lo(2)-ng_pn:,lo(3)-ng_pn:)
    real(kind=dp_t), intent(inout)           :: pi_cc(lo(1)-ng_pc:,lo(2)-ng_pc:,lo(3)-ng_pc:)
    real(kind=dp_t), intent(inout)           :: beta0(lo(1)-ng_b: ,lo(2)-ng_b: ,lo(3)-ng_b:)
    real(kind=dp_t), intent(inout)           :: ncell,pisum
    logical        , intent(in   ), optional ::  mask(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    ! local
    integer :: i,j,k

    logical :: cell_valid

    !$OMP PARALLEL DO PRIVATE(i,j,k,cell_valid) reduction(+:pisum,ncell)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             
             pi_cc(i,j,k) = (pi(i,j,k) + pi(i+1,j,k) + pi(i,j+1,k) + pi(i,j,k+1) &
                  + pi(i+1,j+1,k) + pi(i+1,j,k+1) + pi(i,j+1,k+1) + pi(i+1,j+1,k+1)) / 8.d0

             if (use_alt_energy_fix) then
                pi_cc(i,j,k) = pi_cc(i,j,k)*beta0(i,j,k)
             endif
             
             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if (present(mask)) then
                cell_valid = mask(i,j,k)
             endif

             if (cell_valid) then
                pisum = pisum + weight*pi_cc(i,j,k)
                ncell = ncell + weight
             end if
             
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine make_pi_cc_3d

end module make_pi_cc_module
