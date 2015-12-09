! vec.bdf.f90 is the fully vectorized version of Matt Emmett's BDF code.
! This means all points (e.g. hydro cells) have their own time, timestep, etc,
! whereas the original bdf had limited vectorization in which all points
! shared a single Jacobian and stepped at the same rate.  This version
! also includes adaptations to integrate bdf into BoxLib/Maestro.
!
!
! BDF (backward differentiation formula) time-stepping routines.
!
! See
!
!   1. VODE: A variable-coefficient ODE solver; Brown, Byrne, and
!      Hindmarsh; SIAM J. Sci. Stat. Comput., vol. 10, no. 5, pp.
!      1035-1051, 1989.
!
!   2. An alternative implementation of variable step-size multistep
!      formulas for stiff ODES; Jackson and Sacks-Davis; ACM
!      Trans. Math. Soft., vol. 6, no. 3, pp. 295-318, 1980.
!
!   3. A polyalgorithm for the numerical solution of ODEs; Byrne and
!      Hindmarsh; ACM Trans. Math. Soft., vol. 1, no. 1, pp. 71-96,
!      1975.
!

! TODO:
!  1) In fully vecorizing BDF into VBDF, many loops over npts have been
!  introduced.  We need to make sure the memory access pattern is cache
!  efficient.
module bdf

  use bl_types
  use bl_error_module
  use parallel

  implicit none

  real(dp_t), private, parameter :: one  = 1.0_dp_t
  real(dp_t), private, parameter :: two  = 2.0_dp_t
  real(dp_t), private, parameter :: half = 0.5_dp_t

  integer, parameter :: bdf_max_iters = 666666666

  integer, parameter :: BDF_ERR_SUCCESS  = 0
  integer, parameter :: BDF_ERR_SOLVER   = 1
  integer, parameter :: BDF_ERR_MAXSTEPS = 2
  integer, parameter :: BDF_ERR_DTMIN    = 3

  character(len=64), parameter :: errors(0:3) = [ &
       'Success.                                                ', &
       'Newton solver failed to converge several times in a row.', &
       'Too many steps were taken.                              ', &
       'Minimum time-step reached several times in a row.       ' ]

  !
  ! bdf time-stepper
  !
  type :: bdf_ts

     integer  :: neq                      ! number of equations (degrees of freedom) per point
     integer  :: npt                      ! number of points
     integer  :: max_order                ! maximum order (1 to 6)
     integer  :: max_steps                ! maximum allowable number of steps
     integer  :: max_iters                ! maximum allowable number of newton iterations
     integer  :: verbose                  ! verbosity level
     real(dp_t) :: dt_min                   ! minimum allowable step-size
     real(dp_t) :: eta_min                  ! minimum allowable step-size shrink factor
     real(dp_t) :: eta_max                  ! maximum allowable step-size growth factor
     real(dp_t) :: eta_thresh               ! step-size growth threshold
     integer  :: max_j_age                ! maximum age of jacobian
     integer  :: max_p_age                ! maximum age of newton iteration matrix

     logical  :: debug
     integer  :: dump_unit

     real(dp_t), pointer :: rtol(:)         ! realtive tolerances
     real(dp_t), pointer :: atol(:)         ! absolute tolerances

     ! state
     real(dp_t), pointer :: t(:)          ! current time
     real(dp_t) :: t1                     ! final time
     real(dp_t), pointer :: dt(:)         ! current time step
     real(dp_t), pointer :: dt_nwt(:)     ! dt used when building newton iteration matrix
     integer, pointer :: k(:)             ! current order
     integer  :: n                        ! current step
     !TODO: As of now, j_age is a de-facto scalar.  
     !   Should decide if vectorized age boosts performance
     integer, pointer :: j_age(:)         ! age of jacobian
     integer, pointer :: p_age(:)         ! age of newton iteration matrix
     integer,   pointer  :: k_age(:)      ! number of steps taken at current order
     real(dp_t), pointer :: tq(:,:)    ! error coefficients (test quality)
     real(dp_t), pointer :: tq2save(:)
     logical, pointer  :: refactor(:)

     real(dp_t), pointer :: J(:,:,:)        ! jacobian matrix
     real(dp_t), pointer :: P(:,:,:)        ! newton iteration matrix
     real(dp_t), pointer :: z(:,:,:)        ! nordsieck histroy array, indexed as (dof, p, n)
     real(dp_t), pointer :: z0(:,:,:)       ! nordsieck predictor array
     real(dp_t), pointer :: h(:,:)          ! time steps, h = [ h_n, h_{n-1}, ..., h_{n-k} ]
     real(dp_t), pointer :: l(:,:)          ! predictor/corrector update coefficients (0:max_order,npt)
     real(dp_t), pointer :: upar(:,:)       ! array of user parameters (passed to
                                            ! user's Jacobian and f)
     real(dp_t), pointer :: y(:,:)          ! current y
     real(dp_t), pointer :: yd(:,:)         ! current \dot{y}
     real(dp_t), pointer :: rhs(:,:)        ! solver rhs (NOT the f in xdot = f)
     real(dp_t), pointer :: e(:,:)          ! accumulated correction
     real(dp_t), pointer :: e1(:,:)         ! accumulated correction, previous step
     real(dp_t), pointer :: ewt(:,:)        ! cached error weights
     real(dp_t), pointer :: b(:,:)          ! solver work space
     integer,  pointer :: ipvt(:,:)         ! pivots (neq,npts)
     integer,  pointer :: A(:,:)            ! pascal matrix

     ! counters
     integer, pointer :: nfe(:)             ! number of function evaluations for each pt (npt)
     integer, pointer :: nje(:)             ! number of jacobian evaluations
     integer, pointer :: nlu(:)             ! number of factorizations
     integer, pointer :: nit(:)             ! number of non-linear solver iterations
     integer, pointer :: nse(:)             ! number of non-linear solver errors
     integer, pointer :: ncse(:)            ! number of consecutive non-linear solver errors
     integer, pointer :: ncit(:)            ! number of current non-linear solver iterations
     integer, pointer :: ncdtmin(:)         ! number of consecutive times we tried to shrink beyound the minimum time step

  end type bdf_ts

  private :: &
       rescale_timestep, decrease_order, increase_order, &
       alpha0, alphahat0, xi_j, xi_star_inv, ewts, norm, eye_r, eye_i, factorial
  !public subroutines: bdf_advance, bdf_update, bdf_predict, bdf_solve, bdf_check
  !                    bdf_correct, bdf_dump, bdf_adjust, bdf_reset, print_y
  !                    bdf_ts_build, bdf_ts_destroy, bdf_wrap

contains

  !
  ! Wrapper of the vectorized BDF (VBDF) integrator that mirrors the interface of DVODE.
  ! It translates DVODE input into the equivalent VBDF input and wraps
  ! DVODE-style interfaces with VBDF-style interfaces.
  !
  ! This will be the quickest way to replace DVODE with VBDF, but there will be
  ! no performance benefit.  This is intended for debugging and comparing VBDF
  ! with DVODE.
  !
  ! See the DVODE source code's extensive comments for an explanation of this
  ! interface.
  !
  
  subroutine bdf_wrap(f, neq, y, t, tout, itol, rtol, atol, itask, &
       istate, iopt, rwork, lrw, iwork, liw, jac, mf,    &
       rpar, ipar)
     integer,         intent(in   ) :: neq, itol, itask, iopt, &
                                       lrw, liw, mf
     integer,         intent(inout) :: istate
     integer,         intent(in   ) :: iwork(liw), ipar(:)
     real(kind=dp_t), intent(in   ) :: tout, rtol(:), atol(:), &
                                       rwork(lrw)
     real(kind=dp_t), intent(inout) :: y(neq), t, rpar(:)
     external f, Jac 
     integer, parameter :: NPT = 1         !For DVODE-style calls there's no concept of npt>1
     integer, parameter :: MAX_ORDER = 3   !This is arbitrary, should investigate other values
     logical, parameter :: RESET = .true.  !.true. means we want to initialize the bdf_ts object
     logical, parameter :: REUSE = .false. !.false. means don't reuse the Jacobian
     integer, parameter :: MF_ANALYTIC_JAC = 21
     real(kind=dp_t), parameter :: DT0 = 1.0d-9 !Initial dt to be used in getting from 
                                                !t to tout.  Also arbitrary,
                                                !multiple values should be
                                                !explored.
     type(bdf_ts)    :: ts
     logical         :: first_call
     integer         :: ierr
     real(kind=dp_t) :: y0(neq,NPT), y1(neq,NPT)
     real(kind=dp_t),allocatable :: upar(:,:)

     ! Check user input
     if(mf .ne. MF_ANALYTIC_JAC) then
        call bl_error("ERROR in BDF integrator: mf != MF_ANALYTIC_JAC not yet supported")
     endif

     ! Build the bdf_ts time-stepper object
     allocate(upar(size(rpar),NPT))
     upar(:,NPT) = rpar(:)
     call bdf_ts_build(ts, neq, NPT, rtol, atol, MAX_ORDER, upar)
     ts%dt_min     = 1.0e-20_dp_t 
     
     ! Translate DVODE args into args for bdf_advance
     y0(:,NPT) = y(:)
     call bdf_advance(ts, f_wrap, Jac_wrap, neq, NPT, y0, t, y1, tout, &
                      DT0, RESET, REUSE, ierr, initial_call=.true.)
     t = tout !BDF is designed to always end at tout, 
              !set t to tout to mimic the output behavior of DVODE
     
     istate = ierr
     if (istate .eq. 0) then
        y(:) = y1(:,NPT)
        rpar(:) = upar(:,NPT)
     else
        call bl_error("bdf_advance returned error!: ", errors(istate))
     endif

     ! Cleanup
     call bdf_ts_destroy(ts)

     contains
       ! Wraps the DVODE-style f in a BDF-style interface
       ! ASSUMPTION: All t(:) are the same
       subroutine f_wrap(neq, y, t, yd, upar)
          integer,  intent(in   ) :: neq
          real(kind=dp_t), intent(in   ) :: y(neq), t
          real(kind=dp_t), intent(  out) :: yd(neq)
          real(kind=dp_t), intent(inout), optional :: upar(:)

          real(kind=dp_t), allocatable :: rpar(:)
          integer :: ipar(2) !Dummy array to match DVODE interface

          ipar = -1

          if (present(upar)) then
             allocate(rpar(size(upar)))
          endif

          rpar(:) = upar(:)
          call f(neq, t, y, yd, rpar, ipar)
          upar(:) = rpar(:)
       end subroutine f_wrap

       ! Wraps the DVODE-style Jacobian in a BDF-style interface
       ! ASSUMPTION: All t(:) are the same
       subroutine Jac_wrap(neq, y, t, J, upar)
          integer,  intent(in   ) :: neq
          real(kind=dp_t), intent(in   ) :: y(neq), t
          real(kind=dp_t), intent(  out) :: J(neq, neq)
          real(kind=dp_t), intent(inout), optional :: upar(:)

          real(kind=dp_t), allocatable :: rpar(:)
          integer :: ipar(2), ml, mu

          ml = -1
          mu = -1
          ipar = -1

          if (present(upar)) then
             allocate(rpar(size(upar)))
          endif
          
          rpar(:) = upar(:)
          call Jac(neq, t, y, ml, mu, J, neq, rpar, ipar)
          upar(:) = rpar(:)
       end subroutine Jac_wrap
  end subroutine bdf_wrap

  !
  ! Advance system from t0 to t1.
  !
  subroutine bdf_advance(ts, f, Jac, neq, npt, y0, t0, y1, t1, dt0, reset, reuse, ierr, initial_call)
     type(bdf_ts), intent(inout) :: ts
     integer,      intent(in   ) :: neq, npt
     real(dp_t),   intent(in   ) :: y0(neq,npt), t0, t1, dt0
     real(dp_t),   intent(  out) :: y1(neq,npt)
     logical,      intent(in   ) :: reset, reuse
     integer,      intent(  out) :: ierr
     logical,      intent(in   ), optional :: initial_call
     interface
        subroutine f(neq, y, t, yd, upar)
           import dp_t
           integer,    intent(in   ) :: neq
           real(dp_t), intent(in   ) :: y(neq), t
           real(dp_t), intent(  out) :: yd(neq)
           real(dp_t), intent(inout), optional :: upar(:)
        end subroutine f
        subroutine Jac(neq, y, t, J, upar)
           import dp_t
           integer,    intent(in   ) :: neq
           real(dp_t), intent(in   ) :: y(neq), t
           real(dp_t), intent(  out) :: J(neq, neq)
           real(dp_t), intent(inout), optional :: upar(:)
        end subroutine Jac
     end interface

     integer  :: k, p, m
     logical  :: retry, linitial

     linitial = .false.; if (present(initial_call)) linitial = initial_call

     if (reset) call bdf_reset(ts, f, y0, dt0, reuse)

     ierr = BDF_ERR_SUCCESS

     ts%t1 = t1; ts%t = t0; ts%ncse = 0; ts%ncdtmin = 0;
     do k = 1, bdf_max_iters + 1
        if (ts%n > ts%max_steps .or. k > bdf_max_iters) then
           ierr = BDF_ERR_MAXSTEPS; return
        end if

        if (k == 1) &
             call bdf_dump(ts)
        call bdf_update(ts)                ! update various coeffs (l, tq) based on time-step history
        call bdf_predict(ts)               ! predict nordsieck array using pascal matrix
        if(linitial .and. k == 1) then
           !This is the initial solve, so use the user's initial value, 
           !not the predicted value.
           do p = 1, ts%npt
              do m = 1, ts%neq
                 !Overwrite the predicted z0 with the user's y0
                 ts%z0(m,p,0) = ts%y(m,p)
              end do
           end do
        endif
        call bdf_solve(ts, f, Jac)         ! solve for y_n based on predicted y and yd
        call bdf_check(ts, retry, ierr)    ! check for solver errors and test error estimate

        if (ierr /= BDF_ERR_SUCCESS) return
        if (retry) cycle

        call bdf_correct(ts)               ! new solution looks good, correct history and advance

        call bdf_dump(ts)
        !TODO: As of now, points that reach t1 are just skipped in the 
        !      various calculation routines.
        !      Should we attempt to prune these points out and only
        !      iterate on points that haven't reached t1?
        if (minval(ts%t) >= t1) exit
        
        call bdf_adjust(ts)                ! adjust step-size/order
     end do

     !TODO: Handle how to display dt, k now that it's vector
     if (ts%verbose > 0) &
          print '("BDF: n:",i6,", fe:",i6,", je: ",i3,", lu: ",i3,", &
          &it: ",i3,", se: ",i3,", min(dt): ",e15.8,", min(k): ",i2)', &
          ts%n, maxval(ts%nfe), maxval(ts%nje), maxval(ts%nlu), maxval(ts%nit), maxval(ts%nse), minval(ts%dt), minval(ts%k)

     y1 = ts%z(:,:,0)
  end subroutine bdf_advance

  !
  ! Compute Nordsieck update coefficients l and error coefficients tq.
  !
  ! Regarding the l coefficients, see section 5, and in particular
  ! eqn. 5.2, of Jackson and Sacks-Davis (1980).
  !
  ! Regarding the error coefficients tq, these have been adapted from
  ! cvode.  The tq array is indexed as:
  !
  !  tq(-1) coeff. for order k-1 error est.
  !  tq(0)  coeff. for order k error est.
  !  tq(1)  coeff. for order k+1 error est.
  !  tq(2)  coeff. for order k+1 error est. (used for e_{n-1})
  !
  ! Note:
  !
  !   1. The input vector t = [ t_n, t_{n-1}, ... t_{n-k} ] where we
  !      are advancing from step n-1 to step n.
  !
  !   2. The step size h_n = t_n - t_{n-1}.
  !
  subroutine bdf_update(ts)
    type(bdf_ts), intent(inout) :: ts

    integer  :: j, p
    real(dp_t) :: a0, a0hat, a1, a2, a3, a4, a5, a6, xistar_inv, xi_inv, c

    ts%l  = 0
    ts%tq = 0

    ! compute l vector
    do p = 1, ts%npt
       ts%l(0,p) = 1
       ts%l(1,p) = xi_j(ts%h(:,p), 1)
       if (ts%k(p) > 1) then
          do j = 2, ts%k(p)-1
             ts%l(:,p) = ts%l(:,p) + eoshift(ts%l(:,p), -1) / xi_j(ts%h(:,p), j)
          end do
          ts%l(:,p) = ts%l(:,p) + eoshift(ts%l(:,p), -1) * xi_star_inv(ts%k(p), ts%h(:,p))
       end if
    enddo

    ! compute error coefficients (adapted from cvode)
    do p = 1, ts%npt
       a0hat = alphahat0(ts%k(p), ts%h(:,p))
       a0    = alpha0(ts%k(p))

       xi_inv     = one
       xistar_inv = one
       if (ts%k(p) > 1) then
          xi_inv     = one / xi_j(ts%h(:,p), ts%k(p))
          xistar_inv = xi_star_inv(ts%k(p), ts%h(:,p))
       end if

       a1 = one - a0hat + a0
       a2 = one + ts%k(p) * a1
       ts%tq(0,p) = abs(a1 / (a0 * a2))
       ts%tq(2,p) = abs(a2 * xistar_inv / (ts%l(ts%k(p),p) * xi_inv))
       if (ts%k(p) > 1) then
          c  = xistar_inv / ts%l(ts%k(p),p)
          a3 = a0 + one / ts%k(p)
          a4 = a0hat + xi_inv
          ts%tq(-1,p) = abs(c * (one - a4 + a3) / a3)
       else
          ts%tq(-1,p) = one
       end if

       xi_inv = ts%h(0,p) / sum(ts%h(0:ts%k(p),p))
       a5 = a0 - one / (ts%k(p)+1)
       a6 = a0hat - xi_inv
       ts%tq(1,p) = abs((one - a6 + a5) / a2 / (xi_inv * (ts%k(p)+2) * a5))
    enddo

    call ewts(ts)
  end subroutine bdf_update

  !
  ! Predict (apply Pascal matrix).
  !
  subroutine bdf_predict(ts)
    type(bdf_ts), intent(inout) :: ts
    integer :: i, j, m, p, pp
    do p = 1, ts%npt
      do i = 0, ts%k(p)
          ts%z0(:,p,i) = 0
          do j = i, ts%k(p)
             do m = 1, ts%neq
                ts%z0(m,p,i) = ts%z0(m,p,i) + ts%A(i,j) * ts%z(m,p,j)
             end do
          end do
       end do
    end do
  end subroutine bdf_predict

  !
  ! Solve "y_n - dt f(y_n,t) = y - dt yd" for y_n where y and yd are
  ! predictors from the Nordsieck form.
  !
  ! Newton iteration is:
  !   solve:   P x = -c G(y(k)) for x
  !   update:  y(k+1) = y(k) + x
  ! where
  !   G(y) = y - dt * f(y,t) - rhs
  !
  subroutine bdf_solve(ts, f, Jac)
    type(bdf_ts), intent(inout) :: ts
    interface
       subroutine f(neq, y, t, yd, upar)
         import dp_t
         integer,  intent(in   ) :: neq
         real(dp_t), intent(in   ) :: y(neq), t
         real(dp_t), intent(  out) :: yd(neq)
         real(dp_t), intent(inout), optional :: upar(:)
       end subroutine f
       subroutine Jac(neq, y, t, J, upar)
         import dp_t
         integer,  intent(in   ) :: neq
         real(dp_t), intent(in   ) :: y(neq), t
         real(dp_t), intent(  out) :: J(neq, neq)
         real(dp_t), intent(inout), optional :: upar(:)
       end subroutine Jac
    end interface

    integer  :: k, m, n, p, info
    real(dp_t) :: c, dt_adj, dt_rat, inv_l1 
    logical  :: rebuild, iterating

    do p = 1, ts%npt
       do m = 1, ts%neq
          inv_l1 = 1.0_dp_t / ts%l(1,p)
          ts%e(m,p)   = 0
          ts%rhs(m,p) = ts%z0(m,p,0) - ts%z0(m,p,1) * inv_l1
          ts%y(m,p)   = ts%z0(m,p,0)
       end do

       dt_adj = ts%dt(p) / ts%l(1,p)
       dt_rat = dt_adj / ts%dt_nwt(p)

       if (ts%p_age(p) > ts%max_p_age) ts%refactor(p) = .true.
       if (dt_rat < 0.7d0 .or. dt_rat > 1.429d0) ts%refactor(p) = .true.

       do k = 1, ts%max_iters
          iterating = .true.
          ! build iteration matrix and factor
          if (ts%refactor(p)) then
             rebuild = .true.
             if (ts%ncse(p) == 0 .and. ts%j_age(p) < ts%max_j_age) rebuild = .false.
             if (ts%ncse(p) > 0  .and. (dt_rat < 0.2d0 .or. dt_rat > 5.d0)) rebuild = .false.

             if (rebuild) then
                call Jac(ts%neq, ts%y(:,p), ts%t(p), ts%J(:,:,p), ts%upar(:,p))
                ts%nje(p)   = ts%nje(p) + 1
                ts%j_age(p) = 0
             end if

             call eye_r(ts%P(:,:,p))
            
             do m = 1, ts%neq
                do n = 1, ts%neq
                   ts%P(n,m,p) = ts%P(n,m,p) - dt_adj * ts%J(n,m,p)
                end do
             end do
             call dgefa(ts%P(:,:,p), ts%neq, ts%neq, ts%ipvt(:,p), info)
             ! lapack      call dgetrf(neq, neq, ts%P, neq, ts%ipvt, info)

             ts%nlu(p) = ts%nlu(p) + 1 !The number of times dgefa was called above, LU decomps
             
             ts%dt_nwt(p) = dt_adj
             ts%p_age(p)  = 0
             ts%refactor(p)  = .false.
          end if

          c = 2 * ts%dt_nwt(p) / (dt_adj + ts%dt_nwt(p))

          call f(ts%neq, ts%y(:,p), ts%t(p), ts%yd(:,p), ts%upar(:,p))
          ts%nfe(p) = ts%nfe(p) + 1

          ! solve using factorized iteration matrix
          do m = 1, ts%neq
             ts%b(m,p) = c * (ts%rhs(m,p) - ts%y(m,p) + dt_adj * ts%yd(m,p))
          end do
          call dgesl(ts%P(:,:,p), ts%neq, ts%neq, ts%ipvt(:,p), ts%b(:,p), 0)
          ! lapack   call dgetrs ('N', neq, 1, ts%P, neq, ts%ipvt, ts%b, neq, info)
          ts%nit = ts%nit + 1

          do m = 1, ts%neq
             ts%e(m,p) = ts%e(m,p) + ts%b(m,p)
             ts%y(m,p) = ts%z0(m,p,0) + ts%e(m,p)
          end do
          if (norm(ts%b(:,p), ts%ewt(:,p)) < one) iterating = .false.
          if (.not. iterating) exit
          if (ts%t(p) >= ts%t1) exit
       end do

       ts%ncit(p) = k; ts%p_age(p) = ts%p_age(p) + 1; ts%j_age(p) = ts%j_age(p) + 1
    end do
  end subroutine bdf_solve

  !
  ! Check error estimates.
  !
  subroutine bdf_check(ts, retry, err)
    type(bdf_ts), intent(inout) :: ts
    logical,      intent(out)   :: retry
    integer,      intent(out)   :: err

    real(dp_t) :: error, eta
    integer    :: p

    retry = .false.; err = BDF_ERR_SUCCESS

    ! if solver failed many times, bail
    if (maxval(ts%ncit) >= ts%max_iters .and. maxval(ts%ncse) > 7) then
       err = BDF_ERR_SOLVER
       return
    end if

    ! if solver failed to converge, shrink dt and try again
    ! TODO: As of now, we retry if ANY point failed
    do p = 1, ts%npt
      if (ts%ncit(p) >= ts%max_iters) then
         ts%refactor(p) = .true.; ts%nse(p) = ts%nse(p) + 1; ts%ncse(p) = ts%ncse(p) + 1
         call rescale_timestep(ts, 0.25d0, p)
         retry = .true.
      endif
    enddo
    if (retry) return
    ts%ncse = 0

    ! if local error for a point is too large, shrink dt and try again
    do p = 1, ts%npt
       error = ts%tq(0,p) * norm(ts%e(:,p), ts%ewt(:,p))
       if (error > one) then
          eta = one / ( (6.d0 * error) ** (one / ts%k(p)) + 1.d-6 )
          call rescale_timestep(ts, eta, p)
          retry = .true.
          !TODO: Discuss if we want to use original epsilon implementation
          !if (ts%dt(p) < ts%dt_min + epsilon(ts%dt_min)) ts%ncdtmin = ts%ncdtmin + 1
          if (ts%dt(p) < ts%dt_min) ts%ncdtmin(p) = ts%ncdtmin(p) + 1
          if (ts%ncdtmin(p) > 7) err = BDF_ERR_DTMIN
       end if
    end do
    if (retry) return
    ts%ncdtmin = 0
  end subroutine bdf_check

  !
  ! Correct (apply l coeffs) and advance step.
  !
  subroutine bdf_correct(ts)
    type(bdf_ts), intent(inout) :: ts
    integer :: i, m, p

    do p = 1, ts%npt
       if (ts%t(p) >= ts%t1) cycle !No need to advance
       do i = 0, ts%k(p)
          do m = 1, ts%neq
             ts%z(m,p,i) = ts%z0(m,p,i) + ts%e(m,p) * ts%l(i,p)
          end do
       end do
       ts%h(:,p)   = eoshift(ts%h(:,p), -1)
       ts%t(p)     = ts%t(p) + ts%dt(p)
       ts%h(0,p)   = ts%dt(p)
       ts%k_age(p) = ts%k_age(p) + 1
    end do

    ts%n     = ts%n + 1
  end subroutine bdf_correct


  !
  ! Dump (for debugging)...
  !
  subroutine bdf_dump(ts)
    type(bdf_ts), intent(inout) :: ts
    integer :: i, m, p

    if (.not. ts%debug) return
    write(ts%dump_unit,*) ts%t
    write(ts%dump_unit,*) ts%z(:,:,0)
  end subroutine bdf_dump

  !
  ! Adjust step-size/order to maximize step-size.
  !
  subroutine bdf_adjust(ts)
    type(bdf_ts), intent(inout) :: ts

    real(dp_t) :: c, error(ts%npt), eta(-1:1,ts%npt), rescale, etamax(ts%npt), etaminmax, delta(ts%npt)
    integer  :: p

    rescale = 0

    do p = 1, ts%npt
       ! compute eta(k-1), eta(k), eta(k+1)
       eta(:,p) = 0
       error(p)  = ts%tq(0,p) * norm(ts%e(:,p), ts%ewt(:,p))
       eta(0,p) = one / ( (6.d0 * error(p)) ** (one / ts%k(p)) + 1.d-6 )
       if (ts%k_age(p) > ts%k(p)) then
          if (ts%k(p) > 1) then
             error(p)     = ts%tq(-1,p) * norm(ts%z(:,p,ts%k(p)), ts%ewt(:,p))
             eta(-1,p) = one / ( (6.d0 * error(p)) ** (one / ts%k(p)) + 1.d-6 )
          end if
          if (ts%k(p) < ts%max_order) then
             c = (ts%tq(2,p) / ts%tq2save(p)) * (ts%h(0,p) / ts%h(2,p)) ** (ts%k(p)+1)
             error  = ts%tq(1,p) * norm(ts%e(:,p) - c * ts%e1(:,p), ts%ewt(:,p))
             eta(1,p) = one / ( (10.d0 * error(p)) ** (one / (ts%k(p)+2)) + 1.d-6 )
          end if
          ts%k_age(p) = 0
       end if

       ! choose which eta will maximize the time step
       etamax(p) = 0
       if (eta(-1,p) > etamax(p)) then
          etamax(p) = eta(-1,p)
          delta(p)  = -1
       end if
       if (eta(1,p) > etamax(p)) then
          etamax(p) = eta(1,p)
          delta(p)  = 1
       end if
       if (eta(0,p) > etamax(p)) then
          etamax(p) = eta(0,p)
          delta(p)  = 0
       end if

       ! optimize order
       rescale = 0
       if (etamax(p) > ts%eta_thresh) then
          if (delta(p) == -1) then
             call decrease_order(ts,p)
          else if (delta(p) == 1) then
             call increase_order(ts,p)
          end if
          rescale = etamax(p)
       end if

       ! optimize timestep
       if (ts%t(p) + ts%dt(p) > ts%t1) then
          rescale = (ts%t1 - ts%t(p)) / ts%dt(p)
          call rescale_timestep(ts, rescale, p, .true.)
       else if (rescale /= 0) then
          call rescale_timestep(ts, rescale, p)
       end if
   
       ! save for next step (needed to compute eta(1))
       ts%e1(:,p) = ts%e(:,p)
       ts%tq2save(p) = ts%tq(2,p)

    end do
    
  end subroutine bdf_adjust

  !
  ! Reset counters, set order to one, init Nordsieck history array.
  !
  subroutine bdf_reset(ts, f, y0, dt, reuse)
    type(bdf_ts), intent(inout) :: ts
    real(dp_t),     intent(in   ) :: y0(ts%neq, ts%npt), dt
    logical,      intent(in   ) :: reuse
    interface
       subroutine f(neq, y, t, yd, upar)
         import dp_t
         integer,  intent(in   ) :: neq
         real(dp_t), intent(in   ) :: y(neq), t
         real(dp_t), intent(  out) :: yd(neq)
         real(dp_t), intent(inout), optional :: upar(:)
       end subroutine f
    end interface
    integer :: p

    ts%nfe = 0
    ts%nje = 0
    ts%nlu = 0
    ts%nit = 0
    ts%nse = 0

    ts%y  = y0
    ts%dt = dt
    ts%n  = 1
    ts%k  = 1

    ts%h        = dt
    ts%dt_nwt   = dt
    ts%refactor = .true.

    do p = 1, ts%npt
      call f(ts%neq, ts%y(:,p), ts%t(p), ts%yd(:,p), ts%upar(:,p))
      ts%nfe(p) = ts%nfe(p) + 1
    enddo

    ts%z(:,:,0) = ts%y
    ts%z(:,:,1) = dt * ts%yd

    ts%k_age = 0
    if (.not. reuse) then
       ts%j_age = ts%max_j_age + 1
       ts%p_age = ts%max_p_age + 1
    else
       ts%j_age = 0
       ts%p_age = 0
    end if

  end subroutine bdf_reset

  !
  ! Rescale time-step.
  !
  ! This consists of:
  !   1. bound eta to honor eta_min, eta_max, and dt_min
  !   2. scale dt and adjust time array t accordingly
  !   3. rescale Nordsieck history array
  !
  subroutine rescale_timestep(ts, eta_in, p_in, force_in)
    type(bdf_ts), intent(inout)           :: ts
    real(dp_t),   intent(in   )           :: eta_in
    integer,      intent(in   )           :: p_in
    logical,      intent(in   ), optional :: force_in

    real(dp_t) :: eta
    integer  :: i
    logical  :: force

    force = .false.; if (present(force_in)) force = force_in

    if (force) then
       eta = eta_in
    else
       eta = max(eta_in, ts%dt_min / ts%dt(p_in), ts%eta_min)
       eta = min(eta, ts%eta_max)

       if (ts%t(p_in) + eta*ts%dt(p_in) > ts%t1) then
          eta = (ts%t1 - ts%t(p_in)) / ts%dt(p_in)
       end if
    end if

    ts%dt(p_in)  = eta * ts%dt(p_in)
    ts%h(0,p_in) = ts%dt(p_in)

    do i = 1, ts%k(p_in)
       ts%z(:,p_in,i) = eta**i * ts%z(:,p_in,i)
    end do
  end subroutine rescale_timestep

  !
  ! Decrease order.
  !
  subroutine decrease_order(ts, p)
    type(bdf_ts), intent(inout) :: ts
    integer,      intent(in   ) :: p
    integer  :: j
    real(dp_t) :: c(0:6)

    if (ts%k(p) > 2) then
       c = 0
       c(2) = 1
       do j = 1, ts%k(p)-2
          c = eoshift(c, -1) + c * xi_j(ts%h(:,p), j)
       end do

       do j = 2, ts%k(p)-1
          ts%z(:,p,j) = ts%z(:,p,j) - c(j) * ts%z(:,p,ts%k(p))
       end do
    end if

    ts%z(:,p,ts%k(p)) = 0
    ts%k(p) = ts%k(p) - 1
  end subroutine decrease_order

  !
  ! Increase order.
  !
  subroutine increase_order(ts,p)
    type(bdf_ts), intent(inout) :: ts
    integer,      intent(in   ) :: p
    integer  :: j
    real(dp_t) :: c(0:6)

    c = 0
    c(2) = 1
    do j = 1, ts%k(p)-2
       c = eoshift(c, -1) + c * xi_j(ts%h(:,p), j)
    end do

    ts%z(:,p,ts%k(p)+1) = 0
    do j = 2, ts%k(p)+1
       ts%z(:,p,j) = ts%z(:,p,j) + c(j) * ts%e(:,p)
    end do

    ts%k(p) = ts%k(p) + 1
  end subroutine increase_order

  !
  ! Return $\alpha_0$.
  !
  function alpha0(k) result(a0)
    integer,  intent(in) :: k
    real(dp_t) :: a0
    integer  :: j
    a0 = -1
    do j = 2, k
       a0 = a0 - one / j
    end do
  end function alpha0

  !
  ! Return $\hat{\alpha}_{n,0}$.
  !
  function alphahat0(k, h) result(a0)
    integer,  intent(in) :: k
    real(dp_t), intent(in) :: h(0:k)
    real(dp_t) :: a0
    integer  :: j
    a0 = -1
    do j = 2, k
       a0 = a0 - h(0) / sum(h(0:j-1))
    end do
  end function alphahat0

  !
  ! Return 1 / $\xi^*_k$.
  !
  ! Note that a lot of simplifications can be made to the formula for
  ! $\xi^*_k$ that appears in Jackson and Sacks-Davis.
  !
  function xi_star_inv(k, h) result(xii)
    integer,  intent(in) :: k
    real(dp_t), intent(in) :: h(0:)
    real(dp_t) :: xii, hs
    integer  :: j
    hs = 0.0_dp_t
    xii = -alpha0(k)
    do j = 0, k-2
       hs  = hs + h(j)
       xii = xii - h(0) / hs
    end do
  end function xi_star_inv

  !
  ! Return $\xi_j$.
  !
  function xi_j(h, j) result(xi)
    integer,  intent(in) :: j
    real(dp_t), intent(in) :: h(0:)
    real(dp_t) :: xi
    xi = sum(h(0:j-1)) / h(0)
  end function xi_j

  !
  ! Pre-compute error weights.
  !
  subroutine ewts(ts)
    type(bdf_ts), intent(inout) :: ts
    integer :: m, p
    do p = 1, ts%npt
       do m = 1, ts%neq
          ts%ewt(m,p) = one / (ts%rtol(m) * abs(ts%y(m,p)) + ts%atol(m))
       end do
    end do
  end subroutine ewts

  subroutine print_y(ts)
    type(bdf_ts), intent(in) :: ts
    integer :: p
    do p = 1, ts%npt
       print *, ts%y(:,p)
    end do
  end subroutine print_y

  !
  ! Compute weighted norm of y.
  !
  function norm(y, ewt) result(r)
    real(dp_t), intent(in) :: y(1:), ewt(1:)
    real(dp_t) :: r
    integer :: m, n
    n = size(y)
    r = 0.0_dp_t
    do m = 1, n
       r = r + (y(m)*ewt(m))**2
    end do
    r = sqrt(r/n)
  end function norm

  !
  ! Copy over all point-based values
  ! TODO: We may not need this
  !
  subroutine bdf_ts_ptcopy(ts_src, p_src, ts_dst, p_dst)
     type(bdf_ts), intent(inout) :: ts_src, ts_dst
     integer,      intent(in   ) :: p_src,  p_dst

     if(ts_dst%max_order /= ts_src%max_order) then
        call bl_error('unequal max_orders!')
     endif

     ts_dst%t(p_dst)       = ts_src%t(p_src)
     ts_dst%dt(p_dst)      = ts_src%dt(p_src)
     ts_dst%dt_nwt(p_dst)  = ts_src%dt_nwt(p_src)
     ts_dst%k(p_dst)       = ts_src%k(p_src)
     ts_dst%j_age(p_dst)   = ts_src%j_age(p_src)
     ts_dst%p_age(p_dst)   = ts_src%p_age(p_src)
     ts_dst%k_age(p_dst)   = ts_src%k_age(p_src)
     ts_dst%tq(-1:2,p_dst) = ts_src%tq(-1:2,p_src)
     !ts_dst%tq(:,p_dst) = ts_src%tq(:,p_src)
     ts_dst%tq2save(p_dst) = ts_src%tq2save(p_src)
     ts_dst%J(:,:,p_dst)   = ts_src%J(:,:,p_src)
     ts_dst%P(:,:,p_dst)   = ts_src%P(:,:,p_src)
     !ts_dst%z(:,p_dst,0:ts_dst%max_order) = &
     !   ts_src%z(:,p_src,0:ts_src%max_order)
     ts_dst%z(:,p_dst,:) = ts_src%z(:,p_src,:)
     !ts_dst%z0(:,p_dst,0:ts_dst%max_order) = &
     !   ts_src%z0(:,p_src,0:ts_src%max_order)
     ts_dst%z0(:,p_dst,:) = ts_src%z0(:,p_src,:)
     !ts_dst%h(0:ts_dst%max_order,p_dst)     = ts_src%h(0:ts_src%max_order,p_src)
     ts_dst%h(:,p_dst)     = ts_src%h(:,p_src)
     !ts_dst%l(0:ts_dst%max_order,p_dst)     = ts_src%l(0:ts_src%max_order,p_src)
     ts_dst%l(:,p_dst)     = ts_src%l(:,p_src)
     ts_dst%upar(:,p_dst)  = ts_src%upar(:,p_src)
     ts_dst%y(:,p_dst)     = ts_src%y(:,p_src)
     ts_dst%yd(:,p_dst)    = ts_src%yd(:,p_src)
     ts_dst%rhs(:,p_dst)   = ts_src%rhs(:,p_src)
     ts_dst%e(:,p_dst)     = ts_src%e(:,p_src)
     ts_dst%e1(:,p_dst)    = ts_src%e1(:,p_src)
     ts_dst%ewt(:,p_dst)   = ts_src%ewt(:,p_src)
     ts_dst%b(:,p_dst)     = ts_src%b(:,p_src)
     ts_dst%ipvt(:,p_dst)  = ts_src%ipvt(:,p_src)
  end subroutine bdf_ts_ptcopy

  !
  ! Copy over global values.
  ! TODO: We may not need this
  !
  subroutine bdf_ts_globalcopy(ts_src, ts_dst)
     type(bdf_ts), intent(inout) :: ts_src, ts_dst

     ts_dst%neq        = ts_src%neq         
     ts_dst%npt        = ts_src%npt       
     ts_dst%max_order  = ts_src%max_order 
     ts_dst%max_steps  = ts_src%max_steps 
     ts_dst%max_iters  = ts_src%max_iters 
     ts_dst%verbose    = ts_src%verbose   
     ts_dst%dt_min     = ts_src%dt_min    
     ts_dst%eta_min    = ts_src%eta_min   
     ts_dst%eta_max    = ts_src%eta_max   
     ts_dst%eta_thresh = ts_src%eta_thresh
     ts_dst%max_j_age  = ts_src%max_j_age 
     ts_dst%max_p_age  = ts_src%max_p_age 
     ts_dst%debug      = ts_src%debug
     ts_dst%dump_unit  = ts_src%dump_unit
     ts_dst%rtol       = ts_src%rtol      
     ts_dst%atol       = ts_src%atol      
     ts_dst%t1         = ts_src%t1        
     ts_dst%n          = ts_src%n         
     ts_dst%refactor   = ts_src%refactor
     ts_dst%A          = ts_src%A         
     ts_dst%nfe        = ts_src%nfe       
     ts_dst%nje        = ts_src%nje       
     ts_dst%nlu        = ts_src%nlu       
     ts_dst%nit        = ts_src%nit       
     ts_dst%nse        = ts_src%nse       
     ts_dst%ncse       = ts_src%ncse      
     ts_dst%ncit       = ts_src%ncit      
     ts_dst%ncdtmin    = ts_src%ncdtmin   
  end subroutine bdf_ts_globalcopy

  !
  ! Build/destroy BDF time-stepper.
  !
  subroutine bdf_ts_build(ts, neq, npt, rtol, atol, max_order, upar)
    type(bdf_ts), intent(inout) :: ts
    integer,      intent(in   ) :: max_order, neq, npt
    real(dp_t),     intent(in   ) :: rtol(neq), atol(neq)
    real(dp_t),     intent(in   ), optional :: upar(:,:)

    integer :: k, U(max_order+1, max_order+1), Uk(max_order+1, max_order+1)

    allocate(ts%rtol(neq))
    allocate(ts%atol(neq))
    allocate(ts%t(npt))
    allocate(ts%dt(npt))
    allocate(ts%dt_nwt(npt))
    allocate(ts%k(npt))
    allocate(ts%j_age(npt))
    allocate(ts%p_age(npt))
    allocate(ts%k_age(npt))
    allocate(ts%tq(-1:2, npt))
    allocate(ts%tq2save(npt))
    allocate(ts%refactor(npt))
    allocate(ts%z(neq, npt, 0:max_order))
    allocate(ts%z0(neq, npt, 0:max_order))
    allocate(ts%l(0:max_order, npt))
    allocate(ts%h(0:max_order, npt))
    allocate(ts%P(neq, neq, npt))
    allocate(ts%J(neq, neq, npt))
    allocate(ts%y(neq, npt))
    allocate(ts%yd(neq, npt))
    allocate(ts%rhs(neq, npt))
    allocate(ts%e(neq, npt))
    allocate(ts%e1(neq, npt))
    allocate(ts%ewt(neq, npt))
    allocate(ts%b(neq, npt))
    allocate(ts%ipvt(neq,npt))
    allocate(ts%A(0:max_order, 0:max_order))
    allocate(ts%nfe(npt))
    allocate(ts%nje(npt))
    allocate(ts%nlu(npt))
    allocate(ts%nit(npt))
    allocate(ts%nse(npt))
    allocate(ts%ncse(npt))
    allocate(ts%ncit(npt))
    allocate(ts%ncdtmin(npt))

    if(present(upar)) then
      allocate(ts%upar(size(upar,1),npt))
      ts%upar = upar
    else
      nullify(ts%upar)
    endif

    ts%neq        = neq
    ts%npt        = npt
    ts%max_order  = max_order
    ts%max_steps  = 1000000
    ts%max_iters  = 10
    ts%verbose    = 0
    !ts%dt_min     = 1.e-6*epsilon(ts%dt_min) !TODO: Keep this as default?
    ts%dt_min     = epsilon(ts%dt_min) !TODO: Keep this as default?
    ts%eta_min    = 0.2_dp_t
    ts%eta_max    = 10.0_dp_t
    ts%eta_thresh = 1.50_dp_t
    ts%max_j_age  = 50
    ts%max_p_age  = 20

    ts%k = -1

    ts%rtol = rtol
    ts%atol = atol

    ts%J  = 0
    ts%P  = 0
    ts%yd = 0

    ts%j_age = 666666666
    ts%p_age = 666666666
    ts%k_age = 666666666

    ts%debug = .false.

    ! build pascal matrix A using A = exp(U)
    U = 0
    do k = 1, max_order
       U(k,k+1) = k
    end do
    Uk = U
    call eye_i(ts%A)
    do k = 1, max_order+1
       ts%A  = ts%A + Uk / factorial(k)
       Uk = matmul(U, Uk)
    end do
  end subroutine bdf_ts_build

  subroutine bdf_ts_destroy(ts)
    type(bdf_ts), intent(inout) :: ts
    integer :: p
    deallocate(ts%h,ts%l,ts%ewt,ts%rtol,ts%atol)
    deallocate(ts%t,ts%dt,ts%dt_nwt,ts%k,ts%j_age)
    deallocate(ts%p_age,ts%k_age,ts%tq,ts%tq2save)
    deallocate(ts%refactor,ts%y,ts%yd,ts%z,ts%z0,ts%A)
    deallocate(ts%nfe,ts%nje,ts%nlu,ts%nit,ts%nse,ts%ncse)
    deallocate(ts%ncit,ts%ncdtmin)
    deallocate(ts%P,ts%J,ts%rhs,ts%e,ts%e1,ts%b,ts%ipvt)
    if(associated(ts%upar)) then
      deallocate(ts%upar)
    endif
  end subroutine bdf_ts_destroy

  !
  ! Various misc. helper functions
  !
  subroutine eye_r(A)
    real(dp_t), intent(inout) :: A(:,:)
    integer :: i
    A = 0
    do i = 1, size(A, 1)
       A(i,i) = 1.0_dp_t
    end do
  end subroutine eye_r
  subroutine eye_i(A)
    integer, intent(inout) :: A(:,:)
    integer :: i
    A = 0
    do i = 1, size(A, 1)
       A(i,i) = 1
    end do
  end subroutine eye_i
  recursive function factorial(n) result(r)
    integer, intent(in) :: n
    integer :: r
    if (n == 1) then
       r = 1
    else
       r = n * factorial(n-1)
    end if
  end function factorial

end module bdf
