! t1 is a test example based on
!
!   http://www.radford.edu/~thompson/vodef90web/vodef90source/Double/Prologue_Examples/example1.f90
!
! Some notes from example1.f90:
!
!   The problem is from chemical kinetics, and consists of the following
!   three rate equations:
!
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
!
!   on the interval from t = 0.0d0 to t = 4.d10, with initial
!   conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
!
! The original version of this test distributed with BDF solved the same
! problem, except:
!   * we'll evolve two solutions / initial conditions at the same time
!   * the Jacobian will be computed using the first solution only
!
! This version has been adapted for developing and testing OpenACC acceleration
! of VBDF for use in Maestro's nuclear reaction integration.  We solve the
! problem as in example1.f90 but do so for a vector of hydro cell data.  This
! mimics what we would do in Maestro, where loops over hydro cells are common.


!module feval
!  !use bdf
!  use bl_types
!  implicit none
!  integer, parameter :: neq = 3
!  integer, parameter :: npt = 1
!contains
!  subroutine f_rhs_vec(neq, npt, y, t, ydot, upar)
!    !$acc routine seq
!    integer,  intent(in   ) :: neq, npt
!    real(dp_t), intent(in   ) :: y(neq,npt), t
!    real(dp_t), intent(  out) :: ydot(neq,npt)
!    real(dp_t), intent(inout), optional :: upar(:,:)
!
!    !For the purposes of t1, npt=1
!    ydot(1,1) = -.04d0*y(1,1) + 1.d4*y(2,1)*y(3,1)
!    ydot(3,1) = 3.e7*y(2,1)*y(2,1)
!    ydot(2,1) = -ydot(1,1) - ydot(3,1)
!  end subroutine f
!
!  subroutine jac_vec(neq, npt, y, t, pd, upar)
!    !$acc routine seq
!    integer,  intent(in   ) :: neq, npt
!    real(dp_t), intent(in   ) :: y(neq,npt), t
!    real(dp_t), intent(  out) :: pd(neq,neq,npt)
!    real(dp_t), intent(inout), optional :: upar(:,:)
!
!    !For the purposes of t1, npt=1
!    pd(1,1,1) = -.04d0
!    pd(1,2,1) = 1.d4*y(3,1)
!    pd(1,3,1) = 1.d4*y(2,1)
!    pd(2,1,1) = .04d0
!    pd(2,3,1) = -pd(1,3,1)
!    pd(3,2,1) = 6.e7*y(2,1)
!    pd(2,2,1) = -pd(1,2,1) - pd(3,2,1)
!  end subroutine J
!end module feval

program test
   use bdf
   use bl_types
   implicit none

   integer, parameter :: NEQ = 3
   integer, parameter :: NPT = 1
   integer, parameter :: NCELLS = 64 !32**3 = 32768 is a common grid size in
                                     !Maestro problems
   integer, parameter :: MAX_ORDER = 3
   type(bdf_ts)  :: ts(NCELLS)
   real(dp_t) :: rtol(NEQ), atol(NEQ), dt
   real(dp_t) :: y0(NEQ,NPT), t0, y1(NEQ,NPT), t1, state(NCELLS,NEQ)
   real(dp_t), allocatable :: upar(:,:)

   integer :: i, ierr, navg

   !Build ncells of state data and timestepper objects,
   !create a device version of the ts array
   allocate(upar(1, NPT))
   !$acc enter data create(ts)
   do i = 1, NCELLS
      state(i,:) = [ 1.d0, 0.d0, 0.d0 ]
      rtol = 1.d-4
      atol = [ 1.d-8, 1.d-14, 1.d-6 ]
      call bdf_ts_build(ts(i), NEQ, NPT, rtol, atol, MAX_ORDER, upar)
      ts(i)%temp_data = ts(i)%A(0,0)

      !Now that we've built the ts(i) object on the host, we need to update the
      !device.  For now, this means updating each member that is not allocatable
      !or a pointer.  Were we to do something like
      !   acc update device(ts(i))
      !my understanding is that this would overwrite the device's pointer
      !address with the host's.  By updating each individual member, we instead
      !copy host data into the corresponding device memory.
      !$acc update device(       &
      !$acc    ts(i)%neq,        &
      !$acc    ts(i)%npt,        &
      !$acc    ts(i)%max_order,  &
      !$acc    ts(i)%max_steps,  &
      !$acc    ts(i)%max_iters,  &
      !$acc    ts(i)%verbose,    &
      !$acc    ts(i)%dt_min,     &
      !$acc    ts(i)%eta_min,    &
      !$acc    ts(i)%eta_max,    &
      !$acc    ts(i)%eta_thresh, &
      !$acc    ts(i)%max_j_age,  &
      !$acc    ts(i)%max_p_age,  &
      !$acc    ts(i)%debug,      &
      !$acc    ts(i)%dump_unit,  &
      !$acc    ts(i)%t,          &
      !$acc    ts(i)%t1,         &
      !$acc    ts(i)%dt,         &
      !$acc    ts(i)%dt_nwt,     &
      !$acc    ts(i)%k,          &
      !$acc    ts(i)%n,          &
      !$acc    ts(i)%j_age,      &
      !$acc    ts(i)%p_age,      &
      !$acc    ts(i)%k_age,      &
      !$acc    ts(i)%tq,         &
      !$acc    ts(i)%tq2save,    &
      !$acc    ts(i)%temp_data,  &
      !$acc    ts(i)%refactor,   &
      !$acc    ts(i)%nfe,        &
      !$acc    ts(i)%nje,        &
      !$acc    ts(i)%nlu,        &
      !$acc    ts(i)%nit,        &
      !$acc    ts(i)%nse,        &
      !$acc    ts(i)%ncse,       &
      !$acc    ts(i)%ncit,       &
      !$acc    ts(i)%ncdtmin)

      !For the dynamic data, we have to copy it in.  This will create each of
      !these allocatables on the device, attach them to the appropriate user
      !type, and copy the host data onto the device.
      !$acc enter data copyin(   &
      !$acc    ts(i)%rtol,       &
      !$acc    ts(i)%atol,       &
      !$acc    ts(i)%J,          &
      !$acc    ts(i)%P,          &
      !$acc    ts(i)%z(:,:,0:),  &
      !$acc    ts(i)%z0(:,:,0:), &
      !$acc    ts(i)%h(0:),      &
      !$acc    ts(i)%l(0:),      &
      !$acc    ts(i)%shift(0:),  &
      !$acc    ts(i)%upar,       &
      !$acc    ts(i)%y,          &
      !$acc    ts(i)%yd,         &
      !$acc    ts(i)%rhs,        &
      !$acc    ts(i)%e,          &
      !$acc    ts(i)%e1,         &
      !$acc    ts(i)%ewt,        &
      !$acc    ts(i)%b,          &
      !$acc    ts(i)%ipvt,       &
      !$acc    ts(i)%A(0:,0:))
   enddo
   !$acc enter data copyin(state(:,:))

   t0 = 0.d0
   t1 = 0.4d0
   dt = 1.d-8
   navg=0
   print *, 'state in: ', state(1,:)
   print *, 'A in:     ', ts(1)%A(0,0)

   !Have the GPU loop over state data, with the intention of having each
   !CUDA core execute the acc seq routine bdf_advance on a cell of hydro data

   !$acc parallel loop gang vector reduction(+:navg) private(y0,y1,ierr) &
   !$acc    present(ts, state)
   do i = 1, NCELLS
      !print *, 't, y1(1), y1(2), y1(3), ierr, message'
      !print *, t1, y1(:,1), ierr, errors(ierr)

      !ts(i)%temp_data = ts(i)%rtol(1)
      ts(i)%temp_data = -3.5
      y0(:,NPT) = state(i,:)

      call bdf_advance(ts(i), NEQ, NPT, y0, t0, y1, t1, dt, &
         .true., .false., ierr, .true.)

      state(i,:) = y1(:,NPT)

      navg = navg + ts(i)%n
      !print *, 'td: ', ts%temp_data
      !if (ierr /= BDF_ERR_SUCCESS) exit
   end do


   !Clean up ncells of state data and timestepper objects
   !$acc exit data copyout(state(:,:))
   !TODO: You will have to explicitly update all data members you want on host
   !$acc update host(ts(1)%temp_data)
   print *, 'state out 1: ', state(1,:)
   print *, 'state out 2: ', state(2,:)
   print *, 'temp data 1: ', ts(1)%temp_data
   print *, 'temp data 2: ', ts(2)%temp_data
   do i=1, NCELLS
      !$acc exit data delete(     &
      !$acc    ts(i)%rtol(:),     &
      !$acc    ts(i)%atol(:),     &
      !$acc    ts(i)%J(:,:,:),    &
      !$acc    ts(i)%P(:,:,:),    &
      !$acc    ts(i)%z(:,:,:),    &
      !$acc    ts(i)%z0(:,:,:),   &
      !$acc    ts(i)%h(:),        &
      !$acc    ts(i)%l(:),        &
      !$acc    ts(i)%shift(:),    &
      !$acc    ts(i)%upar(:,:),   &
      !$acc    ts(i)%y(:,:),      &
      !$acc    ts(i)%yd(:,:),     &
      !$acc    ts(i)%rhs(:,:),    &
      !$acc    ts(i)%e(:,:),      &
      !$acc    ts(i)%e1(:,:),     &
      !$acc    ts(i)%ewt(:,:),    &
      !$acc    ts(i)%b(:,:),      &
      !$acc    ts(i)%ipvt(:,:),   &
      !$acc    ts(i)%A(:,:))
      call bdf_ts_destroy(ts(i))
   end do
   deallocate(upar)

   !WARNING! Do *not* do copyout on ts(:), it'll break
   !$acc exit data delete(ts(:))

   !TODO: Either rewrite to get this info for each cell, or delete
   !print *, ''
   !print *, 'max stats for last interval'
   !navg = navg / NCELLS
   print *, 'total number of steps taken, NCELLS: ', navg, NCELLS
   !print *, 'number of function evals   ', ts%nfe
   !print *, 'number of jacobian evals   ', ts%nje
   !print *, 'number of lu decomps       ', ts%nlu
   !print *, 'number of solver iterations', ts%nit
   !print *, 'number of solver errors    ', ts%nse

end program test
