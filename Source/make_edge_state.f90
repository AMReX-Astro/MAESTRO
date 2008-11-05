! make_edge_state constructs the edge state of a variable, using a 
! second-order Taylor expansion in space (through dx/2) and time 
! (though dt/2). 
!
! If is_vel = .true., then we are computing the edge states for the
! velocity.
!
! If velpred = 1 (then is_vel should also be 1), and we are computing 
! the edge states for the MAC.  In this case, we only need the normal 
! edge states (i.e. the x-edge state for u, the y-edge state for v, ...)
!
! If velpred = 0, then we are computing all edge states for each 
! variable.  This is what is done for the final updates of the state 
! variables and velocity.

module make_edge_state_module

  use bl_types

  implicit none

  private

  public :: make_edge_state_1d
  
contains

   subroutine make_edge_state_1d(s,sedgex,w0,force,dx,dt)

     use geometry, only: r_start_coord, r_end_coord, nr_fine, nr, numdisjointchunks, nlevs
     use probin_module, only: slope_order
     use bl_constants_module
     use variables, only: rel_eps
     
     real(kind=dp_t), intent(in   ) ::      s(:,0:)
     real(kind=dp_t), intent(inout) :: sedgex(:,0:)
     real(kind=dp_t), intent(in   ) ::   w0(:,0:)
     real(kind=dp_t), intent(in   ) ::  force(:,0:)
     real(kind=dp_t), intent(in   ) :: dx(:),dt
     
     real(kind=dp_t) :: dmin,dpls,ds,del,slim,sflag
     real(kind=dp_t) :: ubardth, dth, savg, u
     
     integer :: r,lo,hi,n,i

     integer        , parameter :: cen = 1, lim = 2, flag = 3, fromm = 4
     real(kind=dp_t), parameter :: fourthirds = 4.0_dp_t / 3.0_dp_t
        
     real(kind=dp_t) :: slopex(nlevs, 0:nr_fine-1)
     real(kind=dp_t) ::    s_l(nlevs,-1:nr_fine+1)
     real(kind=dp_t) ::    s_r(nlevs,-1:nr_fine+1)
     real(kind=dp_t) ::  dxscr(nlevs, 0:nr_fine-1,4)

     dth = HALF*dt

     ! compute slopes
     do n=1,nlevs

        do i=1,numdisjointchunks(n)
           
           lo = r_start_coord(n,i)
           hi = r_end_coord(n,i)

           if (slope_order .eq. 0) then

              slopex(n,:) = ZERO

           else if (slope_order .eq. 2) then

              if (n .eq. 1) then ! second order slopes for coarsest level

                 ! do standard limiting on interior cells
                 do r = lo+1,hi-1
                    del = half*(s(n,r+1) - s(n,r-1))
                    dpls = two*(s(n,r+1) - s(n,r  ))
                    dmin = two*(s(n,r  ) - s(n,r-1))
                    slim = min(abs(dpls), abs(dmin))
                    slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                    sflag = sign(one,del)
                    slopex(n,r)= sflag*min(slim,abs(del))
                 enddo

                 ! set slopes next to domain boundaries to zero
                 slopex(n,lo) = ZERO
                 slopex(n,hi) = ZERO

              else ! second order slopes for non-coarsest levels

                 do r = lo,hi
                    if (r .eq. 0 .or. r .eq. nr(n)-1) then
                       ! set slopes next to domain boundaries to zero
                       slopex(n,r) = ZERO
                    else
                       ! do standard limiting on interior cells
                       del = half*(s(n,r+1) - s(n,r-1))
                       dpls = two*(s(n,r+1) - s(n,r  ))
                       dmin = two*(s(n,r  ) - s(n,r-1))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                       sflag = sign(one,del)
                       slopex(n,r)= sflag*min(slim,abs(del))
                    end if
                 enddo

              end if

           else if (slope_order .eq. 4) then

              if (n .eq. 1) then ! fourth order slopes for coarsest level

                 do r = lo+1,hi-1
                    dxscr(n,r,cen) = half*(s(n,r+1)-s(n,r-1))
                    dpls = two*(s(n,r+1)-s(n,r  ))
                    dmin = two*(s(n,r  )-s(n,r-1))
                    dxscr(n,r,lim)= min(abs(dmin),abs(dpls))
                    dxscr(n,r,lim) = merge(dxscr(n,r,lim),zero,dpls*dmin.gt.ZERO)
                    dxscr(n,r,flag) = sign(one,dxscr(n,r,cen))
                    dxscr(n,r,fromm)= dxscr(n,r,flag)*min(dxscr(n,r,lim),abs(dxscr(n,r,cen)))
                 enddo

                 dxscr(n,lo,fromm) = ZERO
                 dxscr(n,hi,fromm) = ZERO

                 ! fourth order limited slopes on interior
                 do r = lo+1,hi-1
                    ds = fourthirds * dxscr(n,r,cen) - &
                         sixth * (dxscr(n,r+1,fromm) + dxscr(n,r-1,fromm))
                    slopex(n,r) = dxscr(n,r,flag)*min(abs(ds),dxscr(n,r,lim))
                 enddo

                 ! set slopes adjacent to domain boundaries to zero
                 slopex(n,lo) = ZERO
                 slopex(n,hi) = ZERO

              else ! fourth order slopes for non-coarsest levels

                 do r=lo,hi

                    if (r .eq. 0 .or. r .eq. nr(n)-1) then
                       dxscr(n,r,fromm) = ZERO
                    else
                       dxscr(n,r,cen) = half*(s(n,r+1)-s(n,r-1))
                       dpls = two*(s(n,r+1)-s(n,r  ))
                       dmin = two*(s(n,r  )-s(n,r-1))
                       dxscr(n,r,lim)= min(abs(dmin),abs(dpls))
                       dxscr(n,r,lim) = merge(dxscr(n,r,lim),zero,dpls*dmin.gt.ZERO)
                       dxscr(n,r,flag) = sign(one,dxscr(n,r,cen))
                       dxscr(n,r,fromm)= &
                            dxscr(n,r,flag)*min(dxscr(n,r,lim),abs(dxscr(n,r,cen)))
                    end if

                 end do

                 do r=lo,hi

                    if (r .eq. 0 .or. r .eq. nr(n)-1) then
                       ! set slopes adjacent to domain boundaries to zero
                       slopex(n,r) = ZERO
                    else if (r .eq. r_start_coord(n,i) .or. r .eq. r_end_coord(n,i)) then
                       ! drop order to second-order limited differences
                       del = half*(s(n,r+1) - s(n,r-1))
                       dpls = two*(s(n,r+1) - s(n,r  ))
                       dmin = two*(s(n,r  ) - s(n,r-1))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                       sflag = sign(one,del)
                       slopex(n,r)= sflag*min(slim,abs(del))
                    else
                       ! fourth-order limited slopes on interior
                       ds = fourthirds * dxscr(n,r,cen) - &
                            sixth * (dxscr(n,r+1,fromm) + dxscr(n,r-1,fromm))
                       slopex(n,r) = dxscr(n,r,flag)*min(abs(ds),dxscr(n,r,lim))
                    end if

                 end do

              end if ! which level

           end if ! slope order

        end do ! disjoint chunks

     end do ! end compute slopes

     ! compute s_l and s_r
     do n=1,nlevs

        do i=1,numdisjointchunks(n)
           
           lo = r_start_coord(n,i)
           hi = r_end_coord(n,i)
           
           do r = lo,hi
              
              u = HALF * (w0(n,r) + w0(n,r+1))
              ubardth = dth*u/dx(n)

              
              s_l(n,r+1)= s(n,r) + (HALF-ubardth)*slopex(n,r) + dth * force(n,r)
              s_r(n,r  )= s(n,r) - (HALF+ubardth)*slopex(n,r) + dth * force(n,r)
              
           end do

        end do
        
     end do ! end compute s_l and s_r

     ! compute edge states from s_l and s_r
     do n=1,nlevs

        do i=1,numdisjointchunks(n)
           
           lo = r_start_coord(n,i)
           hi = r_end_coord(n,i)

           ! if we are not at the finest level
           ! copy in the s_r and s_l states from the next finer level at the c-f interface
           if (n .ne. nlevs) then
              s_r(n,r_start_coord(n+1,i)/2) = s_r(n+1,r_start_coord(n+1,i))
              s_l(n,(r_end_coord(n+1,i)+1)/2) = s_l(n+1,r_end_coord(n+1,i)+1)
           end if

           ! if we are not at the coarsest level
           ! copy in the s_l and s_r states from the next coarser level at the c-f interface
           if (n .ne. 1) then
              s_l(n,lo) = s_l(n-1,lo/2)
              s_r(n,hi+1) = s_r(n-1,(hi+1)/2)
           end if

           do r=lo,hi+1
              if (r .eq. 0) then
                 sedgex(n,r) = s_r(n,r)
              else if (r .eq. nr(n)) then
                 sedgex(n,r) = s_l(n,r)
              else
                 sedgex(n,r)=merge(s_l(n,r),s_r(n,r),w0(n,r).gt.ZERO)
                 savg = HALF*(s_r(n,r) + s_l(n,r))
                 sedgex(n,r)=merge(savg,sedgex(n,r),abs(w0(n,r)) .lt. rel_eps)
              end if
           end do

        end do
        
     end do ! end compute edge state from s_l and s_r
     
   end subroutine make_edge_state_1d
   
 end module make_edge_state_module
