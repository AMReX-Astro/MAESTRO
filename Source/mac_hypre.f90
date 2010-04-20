module mac_hypre_module
 
  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use bl_constants_module
  use impose_phys_bcs_on_edges_module
 
  implicit none
 
  private
 
  public :: mac_hypre
 
contains
iii
  subroutine mac_hypre(mla,rh,phi,fine_flx,alpha,beta,dx,the_bc_tower,bc_comp, &
                       stencil_order,ref_ratio,umac_norm)
    use mg_module
!   use coeffs_module
    use probin_module, only : verbose
    use geometry, only: dm, nlevs

    include 'HYPREf.h'

    type(ml_layout), intent(in   )        :: mla
    integer        , intent(in   )        :: stencil_order
    integer        , intent(in   )        :: ref_ratio(:,:)
    real(dp_t)     , intent(in)           :: dx(:,:)
    type(bc_tower) , intent(in)           :: the_bc_tower
    integer        , intent(in   )        :: bc_comp
    type(multifab) , intent(in   )        :: alpha(:), beta(:)
    type(multifab) , intent(inout)        ::    rh(:),  phi(:)
    type(bndry_reg), intent(inout)        :: fine_flx(2:)
    real(dp_t)     , intent(in), optional :: umac_norm(:)

    integer         :: lo(dm),hi(dm)
    integer         :: pdlo(dm),pdhi(dm)
    integer, parameter :: ns = 7

    integer         :: i, j, k, n, nb

    integer(kind=8) :: grid
    integer(kind=8) :: A,b,x
    integer(kind=8) :: stencil
    integer(kind=8) :: solver,precond

    integer         :: ierr
    integer         :: stencil_indices(ns)
    integer         :: offsets(dm,ns)
    type(box)       :: pd
    double precision :: tol
    double precision, allocatable :: values(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    if (nlevs > 1) &
       call bl_error('mac_hypre: not set up for nlevs > 1')

    do i = 1, ns
       stencil_indices(i) = i-1
    enddo

    ! Center
    offsets(1,1) = 0
    offsets(2,1) = 0
    offsets(3,1) = 0

    !Left
    offsets(1,2) = -1
    offsets(2,2) = 0
    offsets(3,2) = 0

    !Right
    offsets(1,3) = 1
    offsets(2,3) = 0
    offsets(3,3) = 0

    !Down(j)
    offsets(1,4) = 0
    offsets(2,4) = -1
    offsets(3,4) = 0

    !Up(j)
    offsets(1,5) = 0
    offsets(2,5) = 1
    offsets(3,5) = 0

    !Down(k)
    offsets(1,6) = 0
    offsets(2,6) = 0
    offsets(3,6) = -1

    !Up(k)
    offsets(1,7) = 0
    offsets(2,7) = 0
    offsets(3,7) = 1

    call HYPRE_StructGridCreate(MPI_COMM_WORLD,dm,grid,ierr) 

    do n = 1,nlevs
      do i = 1, rh(n)%nboxes
         if ( multifab_remote(rh(n), i) ) cycle
         lo =  lwb(get_box(rh(n), i))
         hi =  upb(get_box(rh(n), i))
         call HYPRE_StructGridSetExtents(grid,lo,hi,ierr)
      end do
    end do

    call HYPRE_StructGridAssemble(grid,ierr)

    call HYPRE_StructStencilCreate(dm,ns,stencil,ierr)

    do i = 1, ns
       call HYPRE_StructStencilSetElement(stencil,i-1,offsets(1,i),ierr);
    end do

    call HYPRE_StructMatrixCreate(MPI_COMM_WORLD,grid,stencil,A,ierr)
    call HYPRE_StructMatrixInitialize(A,ierr)

    do n = 1,nlevs

      pd = layout_get_pd(mla%la(n))
      pdlo =  lwb(pd)
      pdhi =  upb(pd)
      do nb = 1, rh(n)%nboxes
         if ( multifab_remote(rh(n), nb) ) cycle
         lo =  lwb(get_box(rh(n), nb))
         hi =  upb(get_box(rh(n), nb))

         allocate(values(1:ns,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
         values(:,:,:,:) = 0.d0

         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)

            if (i.ne.pdlo(1)) &
              values(2,i,j,k) = -1.0

            if (i.ne.pdhi(1)) &
              values(3,i,j,k) = -1.0

            if (j.ne.pdlo(2)) &
              values(4,i,j,k) = -1.0

            if (j.ne.pdhi(2)) &
              values(5,i,j,k) = -1.0

            if (k.ne.pdlo(3)) &
              values(6,i,j,k) = -1.0

            if (k.ne.pdhi(3)) &
              values(7,i,j,k) = -1.0

            values(1,i,j,k) = -(values(2,i,j,k)+values(3,i,j,k)+values(4,i,j,k)+ &
                                values(5,i,j,k)+values(6,i,j,k)+values(7,i,j,k))

         end do
         end do
         end do

         call HYPRE_StructMatrixSetBoxValues(A,lo,hi,ns,stencil_indices,values,ierr)
         deallocate(values)

      end do
    end do

    call HYPRE_StructMatrixAssemble(A,ierr)

    ! Create an empty vector object
    call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, b, ierr)
    call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, x, ierr)
 
    ! Indicate that the vector coefficients are ready to be set
    call HYPRE_StructVectorInitialize(b,ierr)
    call HYPRE_StructVectorInitialize(x,ierr)

    do n = 1,nlevs
      do i = 1, rh(n)%nboxes
         if ( multifab_remote(rh(n), i) ) cycle
         lo =  lwb(get_box(rh(n), i))
         hi =  upb(get_box(rh(n), i))

         ! Set RHS
         allocate(values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
         values(1,:,:,:) = 0.0
         if (i .eq. 1) then
            values(1,lo(1)+3,lo(2)+3,lo(3)+3) =  1.0
            values(1,lo(1)+4,lo(2)+4,lo(3)+4) = -1.0
         end if
         call HYPRE_StructVectorSetBoxValues(b, lo, hi, values, ierr)

         ! Set initial guess
         values(1,:,:,:) = 0.0
         call HYPRE_StructVectorSetBoxValues(x, lo, hi, values, ierr)

         deallocate(values)

      end do
    end do

    ! This is a collective call finalizing the vector assembly.
    !  The vectors are now ``ready to be used''
    call HYPRE_StructVectorAssemble(b,ierr)
    call HYPRE_StructVectorAssemble(x,ierr)

!   Create an empty PCG Struct solver
    call HYPRE_StructPCGCreate(MPI_COMM_WORLD, solver, ierr)
!   Set PCG parameters
    tol = 1.d-12
    call HYPRE_StructPCGSetTol(solver, tol, ierr)
    call HYPRE_StructPCGSetPrintLevel(solver, 2, ierr)
    call HYPRE_StructPCGSetMaxIter(solver, 50, ierr)
 
!   Create the Struct PFMG solver for use as a preconditioner
    call HYPRE_StructPFMGCreate(MPI_COMM_WORLD, precond, ierr)
 
!   Set PFMG parameters
    call HYPRE_StructPFMGSetMaxIter(precond, 1, ierr)
    call HYPRE_StructPFMGSetTol(precond, 0.0, ierr)
    call HYPRE_StructPFMGSetZeroGuess(precond, ierr)
    call HYPRE_StructPFMGSetNumPreRelax(precond, 2, ierr)
    call HYPRE_StructPFMGSetNumPostRelax(precond, 2, ierr)
 
!   Non-Galerkin coarse grid (more efficient for this problem)
    call HYPRE_StructPFMGSetRAPType(precond, 1, ierr)
 
!   R/B Gauss-Seidel
    call HYPRE_StructPFMGSetRelaxType(precond, 2, ierr)
 
!   Skip relaxation on some levels (more efficient for this problem)
    call HYPRE_StructPFMGSetSkipRelax(precond, 1, ierr)

!   Set preconditioner (PFMG = 1) and solve
    call HYPRE_StructPCGSetPrecond(solver, 1, precond, ierr)
    call HYPRE_StructPCGSetup(solver, A, b, x, ierr)
    call HYPRE_StructPCGSolve(solver, A, b, x, ierr)

!   Free memory
!   call HYPRE_StructPFMGDestroy(precond, ierr)
!   call HYPRE_StructPCGDestroy(solver);
!   call HYPRE_StructGridDestroy(grid);
!   call HYPRE_StructStencilDestroy(stencil);
!   call HYPRE_StructMatrixDestroy(A);
!   call HYPRE_StructVectorDestroy(b);
!   call HYPRE_StructVectorDestroy(x);

    call build(bpt, "mac_hypre")

    call destroy(bpt)

  end subroutine mac_hypre
end module mac_hypre_module
