! The sponge acts to damp the velocities at the edge of the star.

module sponge_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module

  implicit none

  real(dp_t), save :: r_sp, r_md, r_tp
  real(dp_t), save :: r_sp_outer, r_tp_outer

  ! the sponge_start_density should be the density below which the
  ! sponge first turns on.  Different problems may compute this in
  ! different ways (i.e. not using sponge_center_density and
  ! sponge_start_factor), so we provide this public module variable to
  ! ensure that the rest of the code always knows at what density the
  ! sponge begins.
  real(dp_t), save, public :: sponge_start_density
  
  private

  public :: init_sponge, make_sponge

contains

  subroutine init_sponge(rho0,dx,prob_lo_r)

    ! The sponge has a HALF * ( 1 - cos( (r - r_sp)/L)) profile, where
    ! the width, L, is r_tp - r_sp.
    !
    ! The center of the sponge, r_md, is set to the radius where r =
    ! sponge_center_density
    !
    ! The start of the sponge, r_sp, (moving outward from the center)
    ! is the radius where r = sponge_start_factor * sponge_center_density
    ! 
    ! The top of the sponge is then 2 * r_md - r_tp

    use geometry, only: dr, r_end_coord, spherical, polar
    use bl_constants_module
    use probin_module, only: verbose, sponge_start_factor, sponge_center_density

    real(kind=dp_t), intent(in   ) :: rho0(0:),prob_lo_r
    real(kind=dp_t), intent(in   ) :: dx(:)

    real (kind = dp_t) :: rloc
    real (kind = dp_t) :: r_top
    integer            :: r

    r_top = prob_lo_r + dble(r_end_coord(1,1)+1) * dr(1)
    r_sp = r_top

    sponge_start_density = sponge_start_factor*sponge_center_density

    ! set r_sp
    do r=0,r_end_coord(1,1)
       rloc = prob_lo_r + (dble(r)+HALF) * dr(1)
       if (rho0(r) < sponge_start_density) then
          r_sp = rloc
          exit
       endif
    enddo

    ! set r_md
    r_md = r_top
    do r=0,r_end_coord(1,1)
       rloc = prob_lo_r + (dble(r)+HALF) * dr(1)
       if (rho0(r) < sponge_center_density) then
          r_md = rloc
          exit
       endif
    enddo

    ! set r_tp
    r_tp = TWO * r_md - r_sp

    ! outer sponge parameters used for spherical problems
    if (spherical .eq. 1 .or. polar .eq. 1) then
       r_sp_outer = r_tp
       r_tp_outer = r_sp_outer + 4.d0 * dx(2)
    end if

    if ( parallel_IOProcessor() .and. verbose .ge. 1) write(6,1000) r_sp, r_tp
    if (spherical .eq. 1 .or. polar .eq. 1) then
       if ( parallel_IOProcessor() .and. verbose .ge. 1) write(6,1001) r_sp_outer, r_tp_outer
    end if
    if ( parallel_IOProcessor() .and. verbose .ge. 1) print*,""

1000 format('inner sponge: r_sp      , r_tp      : ',e20.12,2x,e20.12)
1001 format('outer sponge: r_sp_outer, r_tp_outer: ',e20.12,2x,e20.12)

  end subroutine init_sponge

  subroutine make_sponge(sponge,dx,dt,mla)

    use bl_constants_module
    use ml_cc_restriction_module, only: ml_cc_restriction

    type(multifab) , intent(inout) :: sponge(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(ml_layout), intent(in   ) :: mla

    ! Local variables
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer :: i,n,ng_sp
    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs

    dm = mla%dim
    nlevs = mla%nlevel

    ng_sp = nghost(sponge(1))

    do n=1,nlevs
       do i = 1, nfabs(sponge(n))
          sp => dataptr(sponge(n), i)
          lo =  lwb(get_box(sponge(n), i))
          hi =  upb(get_box(sponge(n), i))
          select case (dm)
          case (2)
             call mk_sponge_2d(sp(:,:,1,1),ng_sp,lo,hi,dx(n,:),dt)
          case (3)
             call mk_sponge_3d(sp(:,:,:,1),ng_sp,lo,hi,dx(n,:),dt)
          end select
       end do
    end do

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       call ml_cc_restriction(sponge(n-1),sponge(n),mla%mba%rr(n-1,:))
    end do

  end subroutine make_sponge

  subroutine mk_sponge_1d(sponge,ng_sp,lo,hi,dx,dt)

    use bl_constants_module
    use probin_module, only: prob_lo, sponge_kappa

    integer        , intent(in   ) ::  lo(:),hi(:), ng_sp
    real(kind=dp_t), intent(inout) :: sponge(lo(1)-ng_sp:)
    real(kind=dp_t), intent(in   ) ::     dx(:),dt

    integer         :: i
    real(kind=dp_t) :: y,smdamp

    sponge = ONE

    do i = lo(1),hi(1)
       y = prob_lo(1) + (dble(i)+HALF)*dx(1)

       if (y >= r_sp) then
          if (y < r_tp) then
             smdamp = HALF*(ONE - cos(M_PI*(y - r_sp)/(r_tp - r_sp)))
          else
             smdamp = ONE
          endif
          sponge(i) = ONE / (ONE + dt * smdamp* sponge_kappa)
       endif

    end do

  end subroutine mk_sponge_1d

  subroutine mk_sponge_2d(sponge,ng_sp,lo,hi,dx,dt)

    use bl_constants_module
    use probin_module, only: prob_lo, sponge_kappa, sponge_mode
    use geometry, only: polar, center

    integer        , intent(in   ) ::  lo(:),hi(:), ng_sp
    real(kind=dp_t), intent(inout) :: sponge(lo(1)-ng_sp:,lo(2)-ng_sp:)
    real(kind=dp_t), intent(in   ) ::     dx(:),dt

    integer         :: i,j
    real(kind=dp_t) :: x,y,r,smdamp

    sponge = ONE
    
    if (polar .eq. 0) then
        do j = lo(2),hi(2)
        y = prob_lo(2) + (dble(j)+HALF)*dx(2)
            
        if (y >= r_sp) then
            if (y < r_tp) then
                smdamp = HALF*(ONE - cos(M_PI*(y - r_sp)/(r_tp - r_sp)))
            else
                smdamp = ONE
            endif
            sponge(:,j) = ONE / (ONE + dt * smdamp* sponge_kappa)
        endif

        end do
    
    else
        if (sponge_mode .ne. 2) then
        ! the classical sponge, damping at a certain radial density
            !$OMP PARALLEL DO PRIVATE(i,j,x,y,r,smdamp)
                do j = lo(2),hi(2)
                    y = prob_lo(2) + (dble(j)+HALF)*dx(2)
                    
                    do i = lo(1),hi(1)
                        x = prob_lo(1) + (dble(i)+HALF)*dx(1)

                        r = sqrt( (x-center(1))**2 + (y-center(2))**2 )
                        
                        ! Inner sponge: damps velocities at edge of star
                        if (r >= r_sp) then
                            if (r < r_tp) then
                                smdamp = HALF*(ONE - cos(M_PI*(r - r_sp)/(r_tp - r_sp)))
                            else
                                smdamp = ONE
                            endif
                            sponge(i,j) = ONE / (ONE + dt * smdamp * sponge_kappa)
                        endif

                        ! Outer sponge: damps velocities in the corners of the domain
                        if (r >= r_sp_outer) then
                            if (r < r_tp_outer) then
                                smdamp = HALF * &
                                    (ONE - cos(M_PI*(r - r_sp_outer)/(r_tp_outer - r_sp_outer)))
                            else
                                smdamp = ONE
                            endif
                            sponge(i,j) = sponge(i,j) / &
                                (ONE + dt * smdamp * 10.d0 * sponge_kappa)
                        endif

                    end do
                end do
            !$OMP END PARALLEL DO
       else if (sponge_mode .eq. 2) then  
        ! the classical sponge, damping all velocity components equally, but in a "box" form
            !$OMP PARALLEL DO PRIVATE(i,j,x,y,r,smdamp)
                do j = lo(2),hi(2)
                    y = prob_lo(2) + (dble(j)+HALF)*dx(2)
                    y = abs(y-center(2))
                    
                    do i = lo(1),hi(1)
                        x = prob_lo(1) + (dble(i)+HALF)*dx(1)
                        x = abs(x-center(1))
                        ! Inner sponge: damps velocities at edge of star
                        if (x >= r_sp .and. x>y) then
                            if (x < r_tp) then
                                smdamp = HALF*(ONE - cos(M_PI*(x - r_sp)/(r_tp - r_sp)))
                            else
                                smdamp = ONE
                            endif
                            sponge(i,j) = ONE / (ONE + dt * smdamp * sponge_kappa)
                        
                        else if (y>= r_sp) then
                            if (y < r_tp) then
                                smdamp = HALF*(ONE - cos(M_PI*(y - r_sp)/(r_tp - r_sp)))
                            else
                                smdamp = ONE
                            endif
                            sponge(i,j) = ONE / (ONE + dt * smdamp * sponge_kappa)
                        endif

                        ! Outer sponge: damps velocities in the corners of the domain
                        if (x >= r_sp_outer .and. x>y) then
                            if (x < r_tp_outer) then
                                smdamp = HALF * &
                                    (ONE - cos(M_PI*(x - r_sp_outer)/(r_tp_outer - r_sp_outer)))
                            else
                                smdamp = ONE
                            endif
                            sponge(i,j) = sponge(i,j) / &
                                (ONE + dt * smdamp * 10.d0 * sponge_kappa)
                        else if (y>=r_sp_outer) then
                            if (y < r_tp_outer) then
                                smdamp = HALF * &
                                    (ONE - cos(M_PI*(y - r_sp_outer)/(r_tp_outer - r_sp_outer)))
                            else
                                smdamp = ONE
                            endif
                            sponge(i,j) = sponge(i,j) / &
                                (ONE + dt * smdamp * 10.d0 * sponge_kappa)
                        endif

                    end do
                end do
            !$OMP END PARALLEL DO    

       else
            call bl_error('Error: sponge mode not defined')
       endif
    end if

  end subroutine mk_sponge_2d

  subroutine mk_sponge_3d(sponge,ng_sp,lo,hi,dx,dt)

    use geometry, only: spherical, center
    use bl_constants_module
    use probin_module, only: prob_lo, sponge_kappa

    integer        , intent(in   ) :: lo(:),hi(:), ng_sp
    real(kind=dp_t), intent(inout) :: sponge(lo(1)-ng_sp:,lo(2)-ng_sp:,lo(3)-ng_sp:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    integer         :: i,j,k
    real(kind=dp_t) :: x,y,z,r,smdamp

    sponge = ONE

    if (spherical .eq. 0) then

       !$OMP PARALLEL DO PRIVATE(k,z,smdamp)
       do k = lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+HALF)*dx(3)

          if (z >= r_sp) then
             if (z < r_tp) then
                smdamp = HALF*(ONE - cos(M_PI*(z - r_sp)/(r_tp - r_sp)))
             else
                smdamp = ONE
             endif
             sponge(:,:,k) = ONE / (ONE + dt * smdamp* sponge_kappa)
          endif

       end do
       !$OMP END PARALLEL DO

    else

       !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,r,smdamp)
       do k = lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+HALF)*dx(3)

          do j = lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+HALF)*dx(2)
             
             do i = lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+HALF)*dx(1)

                r = sqrt( (x-center(1))**2 + (y-center(2))**2 + (z-center(3))**2 )
                
                ! Inner sponge: damps velocities at edge of star
                if (r >= r_sp) then
                   if (r < r_tp) then
                      smdamp = HALF*(ONE - cos(M_PI*(r - r_sp)/(r_tp - r_sp)))
                   else
                      smdamp = ONE
                   endif
                   sponge(i,j,k) = ONE / (ONE + dt * smdamp * sponge_kappa)
                endif

                ! Outer sponge: damps velocities in the corners of the domain
                if (r >= r_sp_outer) then
                   if (r < r_tp_outer) then
                      smdamp = HALF * &
                           (ONE - cos(M_PI*(r - r_sp_outer)/(r_tp_outer - r_sp_outer)))
                   else
                      smdamp = ONE
                   endif
                   sponge(i,j,k) = sponge(i,j,k) / &
                        (ONE + dt * smdamp * 10.d0 * sponge_kappa)
                endif

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    end if

  end subroutine mk_sponge_3d

end module sponge_module
