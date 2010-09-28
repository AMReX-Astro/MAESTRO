! compute beta_0, the coefficient in our constraint equation,
! div{beta_0 U} = beta_0 S

module make_div_coeff_module

  use bl_types

  implicit none

  private

  public :: make_div_coeff

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_div_coeff(div_coeff,rho0,p0,gamma1bar,grav_center)

    use bl_constants_module
    use geometry, only: nr_fine, dr, anelastic_cutoff_coord, r_start_coord, &
         r_end_coord, nr, numdisjointchunks, nlevs_radial
    use restrict_base_module
    use probin_module, only: beta_type
!FIXME
    use make_grav_module

    real(kind=dp_t), intent(  out) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:), p0(:,0:), gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_center(:,0:)

    ! local
    integer :: r, n, i, refrat, j
    real(kind=dp_t) :: integral
    real(kind=dp_t) :: beta0_edge(nlevs_radial,0:nr_fine)
    real(kind=dp_t) :: lambda, mu, nu, kappa
    real(kind=dp_t) :: denom, coeff1, coeff2, coeff3
    real(kind=dp_t) :: del,dpls,dmin,slim,sflag
    real(kind=dp_t) :: offset

    ! FIXME debugging
    real(kind=dp_t) :: dbdr, rho_cc, lambda2, lambda3
    real(kind=dp_t), allocatable :: rho_edge(:,:),grav_edge(:,:)
!    real(kind=dp_t) :: rho_lin(:), p0_lin(:), gamma_lin(:), grav_lin(:)


    div_coeff = 0.d0

    if (beta_type .eq. 1) then

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Compute beta0 on the edges and average to the center      
       !
       ! Multilevel Outline:
       !
       ! First, compute beta0 on edges and centers at level 1 only
       ! Obtain the starting value from rho0 at the bottom of the domain.
       ! do n=2,nlevs_radial
       !   Compute beta0 on edges and centers at level n
       !   Obtain the starting value of beta0_edge_lo from the coarser grid
       !   if n>1, compare the difference between beta0 at the top of level n to the
       !           corresponding point on level n-1
       !   do i=n-1,1,-1
       !     Offset the centered beta on level i above this point so the total integral 
       !      is consistent
       !     Redo the anelastic cutoff part
       !   end do
       ! end do
       ! call restrict_base and fill_ghost_base
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       

! FIXME debugging
       open(unit=11,file='HSE_712_beta.dat',form = "formatted", access = "sequential",action="write")
       write(11,*)'# This has linear gravity with the switch from'
       write(11,*)'#    beta_edge(0) = rho0(0)  to' 
       write(11,*)'#    beta_cell-centered(0) = rho0(0)'
       write(11,*)'# And one-sided slopes in the first and last cells'  
       open(unit=12,file='HSE_712_lin.dat',form = "formatted", access = "sequential",action="write")
       write(12,*)'# This file contains the linearized p, g, rho, gamma used'
       write(12,*)'# to integrate beta'
       write(12,*)'# This has linear gravity with'
       write(12,*)'#    beta_cell-centered(0) = rho0(0)'
       write(12,*)'# instead of  beta_edge(0) = rho0(0)' 
       write(12,*)'# And one-sided slopes in the first and last cells'  
       open(unit=13,file='HSE_712_deriv.dat',form = "formatted", access = "sequential",action="write")
       write(13,*)'# This file contains expressions which should be'
       write(13,*)'# analyticaly the same in isentropic regions.'
       write(13,*)'# Linearized p, g, rho, gamma were used to integrate beta'
       write(13,*)'# This has linear gravity with'
       write(13,*)'#    beta_cell-centered(0) = rho0(0)'
       write(13,*)'# instead of  beta_edge(0) = rho0(0)' 
       write(13,*)'# And one-sided slopes in the first and last cells'  
       write(13,*)'#'
       write(13,*)'# r_index   r   1/beta dbeta/dr   1/gamma*p0 dp0dr   ', &
            '1/rho0 drho0dr   rho0*g/(gamma* p0)'

       allocate(grav_edge(nlevs_radial,0:nr_fine))
       call make_grav_edge(grav_edge,rho0)
       allocate(rho_edge(nlevs_radial,0:nr_fine))
       ! worry about multi-levle later
       rho_edge(1,0) = rho0(1,0)
       do r = 1, nr_fine-1
          rho_edge(1,r) = HALF*(rho0(1,r+1)+rho0(1,r))
       end do
!!!!!!!!!!!

       do n=1,nlevs_radial

          do j=1,numdisjointchunks(n)

             ! Compute beta0 on edges and centers at level n

             if (n .eq. 1) then

                ! original choice
                beta0_edge(1,0) = rho0(1,0) 

                ! try instead requiring that 
                ! beta_cell-centered(0) = rho0(0) 
!                integral = -abs(grav_center(n,0))*rho0(n,0)*dr(n)/(p0(n,0)*gamma1bar(n,0))

!                 ! try one-sided difference for slopes rather than slope = 0
!                 r = 0
!                 lambda = (rho0(n,r+1) - rho0(n,r))/dr(n)
!                 mu =     (gamma1bar(n,r+1) - gamma1bar(n,r))/dr(n)
!                 nu =     (p0(n,r+1) - p0(n,r))/dr(n)
!                 kappa =  (grav_center(n,r+1) - grav_center(n,r))/dr(n)

!                 if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
!                      ( mu*p0(n,r) - nu*gamma1bar(n,r)) .eq. ZERO .or. &
!                      ((gamma1bar(n,r) + HALF*mu*dr(n))/ &
!                      (gamma1bar(n,r) - HALF*mu*dr(n))) .le. ZERO .or. &
!                      ((p0(n,r) + HALF*nu*dr(n))/ &
!                      (p0(n,r) - HALF*nu*dr(n))) .le. ZERO) then

!                    integral = -abs(grav_center(n,r))*rho0(n,r)*dr(n)/(p0(n,r)*gamma1bar(n,r))
!                    write(*,*)r, ': using constant integrand'
!                 else 

!                    denom = mu*p0(n,r) - nu*gamma1bar(n,r) 
!                    coeff1 = (lambda*gamma1bar(n,r) - mu*rho0(n,r)) * &
!                             (kappa *gamma1bar(n,r) + mu*abs(grav_center(n,r))) / &
!                             (mu*mu*denom)
!                    coeff2 = (lambda*p0(n,r) - nu*rho0(n,r))* &
!                             (-kappa*p0(n,r) - nu*abs(grav_center(n,r))) / &
!                             (nu*nu*denom)
!                    coeff3 = kappa*lambda / (mu*nu)

!                    integral =  &
!                         coeff1*log( (gamma1bar(n,r) + HALF*mu*dr(n))/ &
!                                     (gamma1bar(n,r) - HALF*mu*dr(n)) ) + &
!                         coeff2*log( (p0(n,r) + HALF*nu*dr(n))/ &
!                                     (p0(n,r) - HALF*nu*dr(n)) ) + &
!                         coeff3*dr(n)
!                 endif

!                 beta0_edge(1,0) = TWO*rho0(1,0) / (ONE + exp(integral)) 

             else
                ! Obtain the starting value of beta0_edge_lo from the coarser grid
                beta0_edge(n,r_start_coord(n,j)) = beta0_edge(n-1,r_start_coord(n,j)/2)
             end if

             do r=r_start_coord(n,j),r_end_coord(n,j)

                if (r .eq. 0 .or. r .eq. nr(n)-1) then

                   lambda = ZERO
                   mu = ZERO
                   nu = ZERO
                   kappa = ZERO
                ! use one-sided differences                   
!                 if (r .eq. 0 ) then

!                    lambda = (rho0(n,r+1) - rho0(n,r))/dr(n)
!                    mu =     (gamma1bar(n,r+1) - gamma1bar(n,r))/dr(n)
!                    nu =     (p0(n,r+1) - p0(n,r))/dr(n)
!                    kappa =  (grav_center(n,r+1) - grav_center(n,r))/dr(n)

!                 else if ( r .eq. nr(n)-1) then

!                    lambda = (rho0(n,r) - rho0(n,r-1))/dr(n)
!                    mu =     (gamma1bar(n,r) - gamma1bar(n,r-1))/dr(n)
!                    nu =     (p0(n,r) - p0(n,r-1))/dr(n)
!                    kappa =  (grav_center(n,r) - grav_center(n,r-1))/dr(n)  

                else

                   ! decide between second order and fourth order
                   if( .true.) then

                      !second order slopes

                      ! edge based
!                       del    = HALF* (rho_edge(n,r+1) - rho_edge(n,r-1))/dr(n)
!                       dpls   = TWO * (rho_edge(n,r+1) - rho_edge(n,r  ))/dr(n)
!                       dmin   = TWO * (rho_edge(n,r  ) - rho_edge(n,r-1))/dr(n)
!                       slim   = min(abs(dpls), abs(dmin))
!                       slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
!                       sflag  = sign(ONE,del)
!                       lambda = sflag*min(slim,abs(del))

!                       del   = HALF* (gamma1bar(n,r+1) - gamma1bar(n,r-1))/dr(n)
!                       dpls  = TWO * (gamma1bar(n,r+1) - gamma1bar(n,r  ))/dr(n)
!                       dmin  = TWO * (gamma1bar(n,r  ) - gamma1bar(n,r-1))/dr(n)
!                       slim  = min(abs(dpls), abs(dmin))
!                       slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
!                       sflag = sign(ONE,del)
!                       mu    = sflag*min(slim,abs(del))

!                       del   = HALF* (p0(n,r+1) - p0(n,r-1))/dr(n)
!                       dpls  = TWO * (p0(n,r+1) - p0(n,r  ))/dr(n)
!                       dmin  = TWO * (p0(n,r  ) - p0(n,r-1))/dr(n)
!                       slim  = min(abs(dpls), abs(dmin))
!                       slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
!                       sflag = sign(ONE,del)
!                       nu    = sflag*min(slim,abs(del))

!                       del   = HALF* (grav_edge(n,r+1) - grav_edge(n,r-1))/dr(n)
!                       dpls  = TWO * (grav_edge(n,r+1) - grav_edge(n,r  ))/dr(n)
!                       dmin  = TWO * (grav_edge(n,r  ) - grav_edge(n,r-1))/dr(n)
!                       slim  = min(abs(dpls), abs(dmin))
!                       slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
!                       sflag = sign(ONE,del)
!                       kappa = sflag*min(slim,abs(del))

                      !cell-centered
                      del    = HALF* (rho0(n,r+1) - rho0(n,r-1))/dr(n)
                      dpls   = TWO * (rho0(n,r+1) - rho0(n,r  ))/dr(n)
                      dmin   = TWO * (rho0(n,r  ) - rho0(n,r-1))/dr(n)
                      slim   = min(abs(dpls), abs(dmin))
                      slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag  = sign(ONE,del)
                      lambda = sflag*min(slim,abs(del))

                      del   = HALF* (gamma1bar(n,r+1) - gamma1bar(n,r-1))/dr(n)
                      dpls  = TWO * (gamma1bar(n,r+1) - gamma1bar(n,r  ))/dr(n)
                      dmin  = TWO * (gamma1bar(n,r  ) - gamma1bar(n,r-1))/dr(n)
                      slim  = min(abs(dpls), abs(dmin))
                      slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(ONE,del)
                      mu    = sflag*min(slim,abs(del))

                      del   = HALF* (p0(n,r+1) - p0(n,r-1))/dr(n)
                      dpls  = TWO * (p0(n,r+1) - p0(n,r  ))/dr(n)
                      dmin  = TWO * (p0(n,r  ) - p0(n,r-1))/dr(n)
                      slim  = min(abs(dpls), abs(dmin))
                      slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(ONE,del)
                      nu    = sflag*min(slim,abs(del))

                      del   = HALF* (grav_center(n,r+1) - grav_center(n,r-1))/dr(n)
                      dpls  = TWO * (grav_center(n,r+1) - grav_center(n,r  ))/dr(n)
                      dmin  = TWO * (grav_center(n,r  ) - grav_center(n,r-1))/dr(n)
                      slim  = min(abs(dpls), abs(dmin))
                      slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(ONE,del)
                      kappa = sflag*min(slim,abs(del))

                   else
                      
!                       do i = is-2,ie+2 
!                          dxscr(i,cen) = half*(s(i+1,comp)-s(i-1,comp))
!                          dmin = two*(s(i  ,comp)-s(i-1,comp))
!                          dpls = two*(s(i+1,comp)-s(i  ,comp))
!                          dxscr(i,lim)= min(abs(dmin),abs(dpls))
!                          dxscr(i,lim) = merge(dxscr(i,lim),zero,dpls*dmin.gt.ZERO)
!                          dxscr(i,flag) = sign(one,dxscr(i,cen))
!                          dxscr(i,fromm)= dxscr(i,flag)*min(dxscr(i,lim), &
!                               abs(dxscr(i,cen)))
!                       enddo

!                       do i = is-1,ie+1 
!                          ds = two * two3rd * dxscr(i,cen) - &
!                               sixth * (dxscr(i+1,fromm) + dxscr(i-1,fromm)) 
!                          slx(i,comp) = dxscr(i,flag)*min(abs(ds),dxscr(i,lim))
!                       enddo
                    end if

                end if

                ! FOR CONST GRAV
                ! catch the cases where the analytic formula would have us 
                !  dividing by zero or taking log(-#)
!                 if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
!                      (nu*gamma1bar(n,r) - mu*p0(n,r)) .eq. ZERO .or. &
!                      ((gamma1bar(n,r) + HALF*mu*dr(n))/ &
!                      (gamma1bar(n,r) - HALF*mu*dr(n))) .le. ZERO .or. &
!                      ((p0(n,r) + HALF*nu*dr(n))/ &
!                      (p0(n,r) - HALF*nu*dr(n))) .le. ZERO) then

!                    integral = -abs(grav_center(n,r))*rho0(n,r)*dr(n)/(p0(n,r)*gamma1bar(n,r))

!                 else 

!                    denom = nu*gamma1bar(n,r) - mu*p0(n,r)
!                    coeff1 = lambda*gamma1bar(n,r)/mu - rho0(n,r)
!                    coeff2 = lambda*p0(n,r)/nu - rho0(n,r)

!                    integral = (-abs(grav_center(n,r))/denom)* &
!                         (coeff1*log( (gamma1bar(n,r) + HALF*mu*dr(n))/ &
!                         (gamma1bar(n,r) - HALF*mu*dr(n))) - &
!                         coeff2*log( (p0(n,r) + HALF*nu*dr(n))/ &
!                         (p0(n,r) - HALF*nu*dr(n))) )

!                endif

                ! FOR CELL_CENTERED RHO, G
                ! catch the cases where the analytic formula would have us 
                !  dividing by zero or taking log(-#)
                if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
                     ( mu*p0(n,r) - nu*gamma1bar(n,r)) .eq. ZERO .or. &
                     ((gamma1bar(n,r) + HALF*mu*dr(n))/ &
                     (gamma1bar(n,r) - HALF*mu*dr(n))) .le. ZERO .or. &
                     ((p0(n,r) + HALF*nu*dr(n))/ &
                     (p0(n,r) - HALF*nu*dr(n))) .le. ZERO) then

                   integral = -abs(grav_center(n,r))*rho0(n,r)*dr(n)/(p0(n,r)*gamma1bar(n,r))
                   write(*,*)r, ': using constant integrand'
                else 

                   denom = mu*p0(n,r) - nu*gamma1bar(n,r) 
                   coeff1 = (lambda*gamma1bar(n,r) - mu*rho0(n,r)) * &
                            (kappa *gamma1bar(n,r) + mu*abs(grav_center(n,r))) / &
                            (mu*mu*denom)
                   coeff2 = (lambda*p0(n,r) - nu*rho0(n,r))* &
                            (-kappa*p0(n,r) - nu*abs(grav_center(n,r))) / &
                            (nu*nu*denom)
                   coeff3 = kappa*lambda / (mu*nu)

                   integral =  &
                        coeff1*log( (gamma1bar(n,r) + HALF*mu*dr(n))/ &
                                    (gamma1bar(n,r) - HALF*mu*dr(n)) ) + &
                        coeff2*log( (p0(n,r) + HALF*nu*dr(n))/ &
                                    (p0(n,r) - HALF*nu*dr(n)) ) + &
                        coeff3*dr(n)
                endif

                ! FOR EDGE BASED RHO, G
! NEED TO CHANGE TO INTEGRAL FROM 0 TO dr RATHER THAN -dr/2 TO dr/2 
                ! catch the cases where the analytic formula would have us 
                !  dividing by zero or taking log(-#)
!                 if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
!                      ( mu*p0(n,r) - nu*gamma1bar(n,r)) .eq. ZERO .or. &
!                      ((gamma1bar(n,r) + HALF*mu*dr(n))/ &
!                      (gamma1bar(n,r) - HALF*mu*dr(n))) .le. ZERO .or. &
!                      ((p0(n,r) + HALF*nu*dr(n))/ &
!                      (p0(n,r) - HALF*nu*dr(n))) .le. ZERO) then

!                    integral = -abs(grav_edge(n,r))*rho_edge(n,r)*dr(n) / &
!                         (p0(n,r)*gamma1bar(n,r))
!                    write(*,*)r, ': using constant integrand'
!                 else 

!                    denom = mu*p0(n,r) - nu*gamma1bar(n,r) 
!                    coeff1 = (lambda*gamma1bar(n,r) - mu* &
!                               (rho_edge(n,r)+lambda*dr(n)*HALF)) * &
!                             (kappa *gamma1bar(n,r) + mu* &
!                               (abs(grav_edge(n,r))+kappa*dr(n)*HALF) ) / &
!                             (mu*mu*denom)
!                    coeff2 = (lambda*p0(n,r) - nu* &
!                               (rho_edge(n,r)+lambda*dr(n)*HALF)) * &
!                             (-kappa*p0(n,r) - nu* &
!                               (abs(grav_edge(n,r))+kappa*dr(n)*HALF) ) / &
!                             (nu*nu*denom)
!                    coeff3 = kappa*lambda / (mu*nu)

!                    integral =  &
!                         coeff1*log( (gamma1bar(n,r) + HALF*mu*dr(n))/ &
!                                     (gamma1bar(n,r) - HALF*mu*dr(n)) ) + &
!                         coeff2*log( (p0(n,r) + HALF*nu*dr(n))/ &
!                                     (p0(n,r) - HALF*nu*dr(n)) ) + &
!                         coeff3*dr(n)
!                 endif


                beta0_edge(n,r+1) = beta0_edge(n,r) * exp(integral)
                div_coeff(n,r) = HALF*(beta0_edge(n,r) + beta0_edge(n,r+1))

!!!! FIXME
                ! edges
!                 write(12,1000) r, (r+1)*dr(n), &
!                      rho_edge(n,r) + HALF*lambda*dr(n), p0(n,r) + HALF*nu*dr(n), &
!                      gamma1bar(n,r) + HALF*mu*dr(n), &
!                      grav_edge(n,r) + HALF*kappa*dr(n)
!                 write(12,1000) r, r*dr(n), &
!                      rho_edge(n,r) - HALF*lambda*dr(n), p0(n,r) - HALF*nu*dr(n), &
!                      gamma1bar(n,r) - HALF*mu*dr(n), &
!                      grav_edge(n,r) - HALF*kappa*dr(n)

!                 if ( r.eq.0 ) then
!                    dbdr = 0
!                 else
!                    dbdr  = (div_coeff(n,r  ) - div_coeff(n,r-1))/dr(n)
!                 endif
! ! This would actually be a better idea
! !                dbdr  = (beta0_edge(n,r+1) - beta0_edge(n,r  ))/dr(n)

!                 rho_cc = (rho_edge(n,r+1) + rho_edge(n,r  ))/TWO
!                 lambda2 = (rho_edge(n,r+1) - rho_edge(n,r  ))/dr(n)

!                 del    = HALF* (rho0(n,r+1) - rho0(n,r-1))/dr(n)
!                 dpls   = TWO * (rho0(n,r+1) - rho0(n,r  ))/dr(n)
!                 dmin   = TWO * (rho0(n,r  ) - rho0(n,r-1))/dr(n)
!                 slim   = min(abs(dpls), abs(dmin))
!                 slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
!                 sflag  = sign(ONE,del)
!                 lambda3 = sflag*min(slim,abs(del))

!                 write(13,1000) r, (r+HALF)*dr(n), r*dr(n), &
!                      dbdr/div_coeff(n,r), &
!                      (beta0_edge(n,r+1)-beta0_edge(n,r))/ &
!                         (dr(n)*div_coeff(n,r)),&
!                      nu/(gamma1bar(n,r)*p0(n,r)), &
!                      lambda/rho_edge(n,r), lambda2/rho_cc, lambda3/rho0(n,r), &
!                      -(rho_edge(n,r)+rho_edge(n,r))* &
!                          (abs(grav_edge(n,r))+abs(grav_edge(n,r)))/ &
!                          (FOUR*gamma1bar(n,r)*p0(n,r)), &
!                      dr(n)

                ! cell-centered
                write(12,1000) r, (r+1)*dr(n), &
                     rho0(n,r) + HALF*lambda*dr(n), p0(n,r) + HALF*nu*dr(n), &
                     gamma1bar(n,r) + HALF*mu*dr(n), &
                     grav_center(n,r) + HALF*kappa*dr(n)
                write(12,1000) r, r*dr(n), &
                     rho0(n,r) - HALF*lambda*dr(n), p0(n,r) - HALF*nu*dr(n), &
                     gamma1bar(n,r) - HALF*mu*dr(n), &
                     grav_center(n,r) - HALF*kappa*dr(n)

                if ( r.eq.0 ) then
                   dbdr = 0
                else
                   dbdr  = (div_coeff(n,r  ) - div_coeff(n,r-1))/dr(n)
                endif
! This would actually be a better idea
!                dbdr  = (beta0_edge(n,r+1) - beta0_edge(n,r  ))/dr(n)
                write(13,1000) r, (r+HALF)*dr(n), dbdr/div_coeff(n,r), &
                     nu/(gamma1bar(n,r)*p0(n,r)), lambda/rho0(n,r), &
                     -rho0(n,r)*abs(grav_center(n,r))/(gamma1bar(n,r)*p0(n,r))
!!!!!!!!!!!!

! try different averaging
! this did not seem to be a good idea
! compute ln(beta_edge)
!                beta0_edge(n,r+1) = beta0_edge(n,r) - integral
! compute beta at cell-centers
!                div_coeff(n,r) = exp(HALF*(beta0_edge(n,r) + beta0_edge(n,r+1)))

             end do

! fixme debugging
! didn't think at all about the multi-level case
!              r = r_end_coord(n,j)
!              div_coeff(n,r) = HALF*(beta0_edge(n,r) + beta0_edge(n,r+1))
!              r = r_start_coord(n,j)
!              div_coeff(n,r) = HALF*(beta0_edge(n,r) + beta0_edge(n,r+1))

!              do r=r_start_coord(n,j)+1,r_end_coord(n,j)-1
!                 div_coeff(n,r) = ( SEVEN*(beta0_edge(n,r) + beta0_edge(n,r+1)) &
!                                  - beta0_edge(n,r-1) & 
!                                  - beta0_edge(n,r+2))/TWELVE
!              enddo
!!!!!!!!!!!!!

             do r = anelastic_cutoff_coord(n),r_end_coord(n,j)
                div_coeff(n,r) = div_coeff(n,r-1) * (rho0(n,r)/rho0(n,r-1))
             end do

! FIXME debugging!
             write(13,*)
             write(13,*)

             ! edge
!              do r=r_start_coord(n,j),r_end_coord(n,j)

!                 if ( r.eq.0 ) then
!                    dbdr = (div_coeff(n,r+1) - div_coeff(n,r))/dr(n)
!                 else if ( r.eq.nr(n)-1 ) then
!                    dmin  =(div_coeff(n,r  ) - div_coeff(n,r-1))/dr(n) 
!                 else
!                    del   = HALF* (div_coeff(n,r+1) - div_coeff(n,r-1))/dr(n)
!                    dpls  = TWO * (div_coeff(n,r+1) - div_coeff(n,r  ))/dr(n)
!                    dmin  = TWO * (div_coeff(n,r  ) - div_coeff(n,r-1))/dr(n)
!                    slim  = min(abs(dpls), abs(dmin))
!                    slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
!                    sflag = sign(ONE,del)
!                    dbdr    = sflag*min(slim,abs(del))                      
!                 endif
!                 write(13,1000) r, (r+HALF)*dr(n), dbdr/div_coeff(n,r)

!                 write(11,1000)r,(r+HALF)*dr(n), grav_edge(n,r),p0(n,r), &
!                      gamma1bar(n,r),rho_edge(n,r),div_coeff(n,r), &
!                      (rho_edge(n,r)-div_coeff(n,r))/rho_edge(n,r), &
!                      -rho_edge(n,r)*abs(grav_edge(n,r))/(gamma1bar(n,r)*p0(n,r))
!              end do

             ! cell-centered
             do r=r_start_coord(n,j),r_end_coord(n,j)

                if ( r.eq.0 ) then
                   dbdr = (div_coeff(n,r+1) - div_coeff(n,r))/dr(n)
                else if ( r.eq.nr(n)-1 ) then
                   dmin  =(div_coeff(n,r  ) - div_coeff(n,r-1))/dr(n) 
                else
                   del   = HALF* (div_coeff(n,r+1) - div_coeff(n,r-1))/dr(n)
                   dpls  = TWO * (div_coeff(n,r+1) - div_coeff(n,r  ))/dr(n)
                   dmin  = TWO * (div_coeff(n,r  ) - div_coeff(n,r-1))/dr(n)
                   slim  = min(abs(dpls), abs(dmin))
                   slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                   sflag = sign(ONE,del)
                   dbdr    = sflag*min(slim,abs(del))                      
                endif
                write(13,1000) r, (r+HALF)*dr(n), dbdr/div_coeff(n,r)

                write(11,1000)r,(r+HALF)*dr(n), grav_center(n,r),p0(n,r), &
                     gamma1bar(n,r),rho0(n,r),div_coeff(n,r), &
                     (rho0(n,r)-div_coeff(n,r))/rho0(n,r), &
                     -rho0(n,r)*abs(grav_center(n,r))/(gamma1bar(n,r)*p0(n,r))
             end do

             stop
1000 format(i8,32(e30.20,1x))
!!!!!!!!

             if (n .gt. 1) then

                ! Compare the difference between beta0 at the top of level n to the 
                ! corresponding point on level n-1
                offset = beta0_edge(n,r_end_coord(n,j)+1) &
                     - beta0_edge(n-1,(r_end_coord(n,j)+1)/2)

                do i=n-1,1,-1

                   refrat = 2**(n-i)

                   ! Offset the centered beta on level i above this point so the total 
                   ! integral is consistent
                   do r=r_end_coord(n,j)/refrat+1,nr(i)
                      div_coeff(i,r) = div_coeff(i,r) + offset
                   end do

                   ! Redo the anelastic cutoff part
                   do r=anelastic_cutoff_coord(i),nr(i)
                      if (rho0(i,r-1) /= ZERO) then
                         div_coeff(i,r) = div_coeff(i,r-1) * (rho0(i,r)/rho0(i,r-1))
                      endif
                   end do

                   ! This next piece of coded is needed for the case when the anelastic 
                   ! cutoff coordinate lives on level n.  We first average div_coeff from 
                   ! level i+1 to level i in the region between the anelastic cutoff and 
                   ! the top of grid n.  Then recompute div_coeff at level i above the top 
                   ! of grid n.
                   if (r_end_coord(n,j) .ge. anelastic_cutoff_coord(n)) then

                      do r=anelastic_cutoff_coord(i),(r_end_coord(n,j)+1)/refrat-1
                         div_coeff(i,r) = HALF*(div_coeff(i+1,2*r)+div_coeff(i+1,2*r+1))
                      end do

                      do r=(r_end_coord(n,j)+1)/refrat,nr(i)
                         div_coeff(i,r) = div_coeff(i,r-1) * (rho0(i,r)/rho0(i,r-1))
                      end do

                   end if

                end do ! end loop over i=n-1,1,-1

             end if ! end if (n .gt. 1)

          end do ! end loop over disjoint chunks

       end do ! end loop over levels

       ! zero the div_coeff where there is no corresponding full state array
       do n=2,nlevs_radial
          do j=1,numdisjointchunks(n)
             if (j .eq. numdisjointchunks(n)) then
                do r=r_end_coord(n,j)+1,nr(n)-1
                   div_coeff(n,r) = ZERO
                end do
             else
                do r=r_end_coord(n,j)+1,r_start_coord(n,j+1)-1
                   div_coeff(n,r) = ZERO
                end do
             end if
          end do
       end do

    else if (beta_type .eq. 2) then

       ! beta_0 = rho_0
       do n=1,nlevs_radial
          do j=1,numdisjointchunks(n)
             do r=r_start_coord(n,j),r_end_coord(n,j)
                div_coeff(n,r) = rho0(n,r)
             end do
          end do
       end do

    else if (beta_type .eq. 3) then

       ! beta_0 = 1.d0
       do n=1,nlevs_radial
          do j=1,numdisjointchunks(n)
             do r=r_start_coord(n,j),r_end_coord(n,j)
                div_coeff = 1.d0
             end do
          end do
       end do

    end if

    call restrict_base(div_coeff,.true.)
    call fill_ghost_base(div_coeff,.true.)

  end subroutine make_div_coeff

end module make_div_coeff_module
