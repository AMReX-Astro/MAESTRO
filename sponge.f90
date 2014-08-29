! a module for storing the geometric information so we don't have to pass it

module sponge_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  real(dp_t), save :: topsponge_lo_r, topsponge_hi_r, botsponge_lo_r, botsponge_hi_r

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

    use geometry, only: dr, r_end_coord
    use bl_constants_module
    use probin_module, only: anelastic_cutoff, prob_hi, sponge_start_factor

    real(kind=dp_t), intent(in   ) :: rho0(0:),prob_lo_r
    real(kind=dp_t), intent(in   ) :: dx(:)

    real (kind = dp_t) :: rloc
    real (kind = dp_t) :: r_top
    integer            :: r

    r_top = prob_lo_r + dble(r_end_coord(1,1)+1) * dr(1)
    topsponge_lo_r = r_top

    ! we do things differently than what is in Source/sponge.f90 
    sponge_start_density = sponge_start_factor*anelastic_cutoff

    do r=0,r_end_coord(1,1)
       rloc = prob_lo_r + (dble(r)+HALF) * dr(1)
       if (rho0(r) < sponge_start_density) then
          topsponge_lo_r = rloc
          exit
       endif
    enddo

    do r=0,r_end_coord(1,1)
       rloc = prob_lo_r + (dble(r)+HALF) * dr(1)
       if (rho0(r) < anelastic_cutoff) then
          topsponge_hi_r = rloc
          exit
       endif
    enddo

    do r=0,r_end_coord(1,1)
       rloc = prob_lo_r + (dble(r)+HALF) * dr(1)
       if (rho0(r) < 6.d7) then
          botsponge_lo_r = rloc
          exit
       endif
    enddo

    do r=0,r_end_coord(1,1)
       rloc = prob_lo_r + (dble(r)+HALF) * dr(1)
       if (rho0(r) < 5.d7) then
          botsponge_hi_r = rloc
          exit
       endif
    enddo


    if ( parallel_IOProcessor() ) write(6,1000) topsponge_lo_r, topsponge_hi_r

1000 format('inner sponge: topsponge_lo_r      , topsponge_hi_r      : ',e20.12,2x,e20.12)

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
    integer :: lo(get_dim(sponge(1))),hi(get_dim(sponge(1))),dm,nlevs

    dm = get_dim(sponge(1))
    nlevs = mla%nlevel

    ng_sp = sponge(1)%ng

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

  subroutine mk_sponge_2d(sponge,ng_sp,lo,hi,dx,dt)

    use bl_constants_module
    use probin_module, only: prob_lo, xrb_use_bottom_sponge, sponge_min

    integer        , intent(in   ) :: lo(:),hi(:),ng_sp
    real(kind=dp_t), intent(inout) :: sponge(lo(1)-ng_sp:,lo(2)-ng_sp:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    integer         :: j
    real(kind=dp_t) :: y

    if (xrb_use_bottom_sponge) then

       do j = lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+HALF)*dx(2)
          if(y .le. botsponge_lo_r) then
             sponge(:,j) = sponge_min
          else if(y .le. botsponge_hi_r) then
             sponge(:,j) = -HALF*(ONE-sponge_min) &
                  * cos(M_PI*(y-botsponge_lo_r)/(botsponge_hi_r-botsponge_lo_r)) &
                  + HALF*(ONE+sponge_min)
          else if(y .le. topsponge_lo_r) then
             sponge(:,j) = ONE
          else if (y .le. topsponge_hi_r) then
             sponge(:,j) = HALF*(ONE-sponge_min) &
                  * cos(M_PI*(y-topsponge_lo_r)/(topsponge_hi_r-topsponge_lo_r)) &
                  + HALF*(ONE+sponge_min)
          else
             sponge(:,j) = sponge_min
          end if
       end do

    else

       do j = lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+HALF)*dx(2)
          if(y .le. topsponge_lo_r) then
             sponge(:,j) = ONE
          else if (y .le. topsponge_hi_r) then
             sponge(:,j) = HALF*(ONE-sponge_min) &
                  * cos(M_PI*(y-topsponge_lo_r)/(topsponge_hi_r-topsponge_lo_r)) &
                  + HALF*(ONE+sponge_min)
          else
             sponge(:,j) = sponge_min
          end if
       end do

    end if

  end subroutine mk_sponge_2d

  subroutine mk_sponge_3d(sponge,ng_sp,lo,hi,dx,dt)

    use bl_constants_module
    use probin_module, only: prob_lo, xrb_use_bottom_sponge, sponge_min

    integer        , intent(in   ) :: lo(:),hi(:),ng_sp
    real(kind=dp_t), intent(inout) :: sponge(lo(1)-ng_sp:,lo(2)-ng_sp:,lo(3)-ng_sp:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    integer         :: k
    real(kind=dp_t) :: z

    if (xrb_use_bottom_sponge) then

       do k = lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+HALF)*dx(3)
          if(z .le. botsponge_lo_r) then
             sponge(:,:,k) = sponge_min
          else if(z .le. botsponge_hi_r) then
             sponge(:,:,k) = -HALF*(ONE-sponge_min) &
                  * cos(M_PI*(z-botsponge_lo_r)/(botsponge_hi_r-botsponge_lo_r)) &
                  + HALF*(ONE+sponge_min)
          else if(z .le. topsponge_lo_r) then
             sponge(:,:,k) = ONE
          else if (z .le. topsponge_hi_r) then
             sponge(:,:,k) = HALF*(ONE-sponge_min) &
                  * cos(M_PI*(z-topsponge_lo_r)/(topsponge_hi_r-topsponge_lo_r)) &
                  + HALF*(ONE+sponge_min)
          else
             sponge(:,:,k) = sponge_min
          end if
       end do

    else

       do k = lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+HALF)*dx(3)
          if(z .le. topsponge_lo_r) then
             sponge(:,:,k) = ONE
          else if (z .le. topsponge_hi_r) then
             sponge(:,:,k) = HALF*(ONE-sponge_min) &
                  * cos(M_PI*(z-topsponge_lo_r)/(topsponge_hi_r-topsponge_lo_r)) &
                  + HALF*(ONE+sponge_min)
          else
             sponge(:,:,k) = sponge_min
          end if
       end do

    end if

   end subroutine mk_sponge_3d

end module sponge_module
