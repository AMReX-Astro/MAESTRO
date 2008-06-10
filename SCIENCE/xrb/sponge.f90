! a module for storing the geometric information so we don't have to pass it

module sponge_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  real(dp_t), save :: topsponge_lo_r, topsponge_hi_r, botsponge_lo_r, botsponge_hi_r

  private

  public :: init_sponge, make_sponge

contains

  subroutine init_sponge(rho0,prob_hi,dx,prob_lo_r)

    use geometry, only: dr, r_end_coord
    use bl_constants_module
    use probin_module, only: anelastic_cutoff

    real(kind=dp_t), intent(in   ) :: rho0(0:),prob_lo_r
    real(kind=dp_t), intent(in   ) :: prob_hi(:),dx(:)

    real (kind = dp_t) :: rloc
    real (kind = dp_t) :: r_top
    integer            :: r

    r_top = prob_lo_r + dble(r_end_coord(1)+1) * dr(1)
    topsponge_lo_r = r_top

    do r=0,r_end_coord(1)
       rloc = prob_lo_r + (dble(r)+HALF) * dr(1)
       if (rho0(r) < 25.d0*anelastic_cutoff) then
          topsponge_lo_r = rloc
          exit
       endif
    enddo

    do r=0,r_end_coord(1)
       rloc = prob_lo_r + (dble(r)+HALF) * dr(1)
       if (rho0(r) < anelastic_cutoff) then
          topsponge_hi_r = rloc
          exit
       endif
    enddo

    do r=0,r_end_coord(1)
       rloc = prob_lo_r + (dble(r)+HALF) * dr(1)
       if (rho0(r) < 6.d7) then
          botsponge_lo_r = rloc
          exit
       endif
    enddo

    do r=0,r_end_coord(1)
       rloc = prob_lo_r + (dble(r)+HALF) * dr(1)
       if (rho0(r) < 5.d7) then
          botsponge_hi_r = rloc
          exit
       endif
    enddo


    if ( parallel_IOProcessor() ) write(6,1000) topsponge_lo_r, topsponge_hi_r

1000 format('inner sponge: topsponge_lo_r      , topsponge_hi_r      : ',e20.12,2x,e20.12)

  end subroutine init_sponge

  subroutine make_sponge(nlevs,sponge,dx,dt,mla)

    use bl_constants_module
    use ml_restriction_module, only: ml_cc_restriction

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: sponge(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(ml_layout), intent(in   ) :: mla

    ! Local variables
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer :: i,dm,n
    integer :: lo(sponge(1)%dim),hi(sponge(1)%dim)

    dm = sponge(1)%dim

    do n=1,nlevs

       do i = 1, sponge(n)%nboxes
          if ( multifab_remote(sponge(n), i) ) cycle
          sp => dataptr(sponge(n), i)
          lo =  lwb(get_box(sponge(n), i))
          hi =  upb(get_box(sponge(n), i))
          select case (dm)
          case (2)
             call mk_sponge_2d(sp(:,:,1,1),lo,hi,dx(n,:),dt)
          case (3)
             call mk_sponge_3d(sp(:,:,:,1),lo,hi,dx(n,:),dt)
          end select
       end do

    end do

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       call ml_cc_restriction(sponge(n-1),sponge(n),mla%mba%rr(n-1,:))
    end do

  end subroutine make_sponge

  subroutine mk_sponge_2d(sponge,lo,hi,dx,dt)

    use bl_constants_module
    use probin_module, only: prob_lo_y, use_xrb_bottom_sponge

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sponge(lo(1):,lo(2):)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    integer         :: j
    real(kind=dp_t) :: y,spongemin

    spongemin = 0.01d0

    if (use_xrb_bottom_sponge) then

       do j = lo(2),hi(2)
          y = prob_lo_y + (dble(j)+HALF)*dx(2)
          if(y .le. botsponge_lo_r) then
             sponge(:,j) = spongemin
          else if(y .le. botsponge_hi_r) then
             sponge(:,j) = -HALF*(ONE-spongemin) &
                  * cos(M_PI*(y-botsponge_lo_r)/(botsponge_hi_r-botsponge_lo_r)) &
                  + HALF*(ONE+spongemin)
          else if(y .le. topsponge_lo_r) then
             sponge(:,j) = ONE
          else if (y .le. topsponge_hi_r) then
             sponge(:,j) = HALF*(ONE-spongemin) &
                  * cos(M_PI*(y-topsponge_lo_r)/(topsponge_hi_r-topsponge_lo_r)) &
                  + HALF*(ONE+spongemin)
          else
             sponge(:,j) = spongemin
          end if
       end do

    else

       do j = lo(2),hi(2)
          y = prob_lo_y + (dble(j)+HALF)*dx(2)
          if(y .le. topsponge_lo_r) then
             sponge(:,j) = ONE
          else if (y .le. topsponge_hi_r) then
             sponge(:,j) = HALF*(ONE-spongemin) &
                  * cos(M_PI*(y-topsponge_lo_r)/(topsponge_hi_r-topsponge_lo_r)) &
                  + HALF*(ONE+spongemin)
          else
             sponge(:,j) = spongemin
          end if
       end do

    end if

  end subroutine mk_sponge_2d

  subroutine mk_sponge_3d(sponge,lo,hi,dx,dt)

    use bl_constants_module
    use probin_module, only: prob_lo_z, use_xrb_bottom_sponge

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sponge(lo(1):,lo(2):,lo(3):)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    integer         :: k
    real(kind=dp_t) :: z,spongemin

    spongemin = 0.01d0

    if (use_xrb_bottom_sponge) then

       do k = lo(3),hi(3)
          z = prob_lo_z + (dble(k)+HALF)*dx(3)
          if(z .le. botsponge_lo_r) then
             sponge(:,:,k) = spongemin
          else if(z .le. botsponge_hi_r) then
             sponge(:,:,k) = -HALF*(ONE-spongemin) &
                  * cos(M_PI*(z-botsponge_lo_r)/(botsponge_hi_r-botsponge_lo_r)) &
                  + HALF*(ONE+spongemin)
          else if(z .le. topsponge_lo_r) then
             sponge(:,:,k) = ONE
          else if (z .le. topsponge_hi_r) then
             sponge(:,:,k) = HALF*(ONE-spongemin) &
                  * cos(M_PI*(z-topsponge_lo_r)/(topsponge_hi_r-topsponge_lo_r)) &
                  + HALF*(ONE+spongemin)
          else
             sponge(:,:,k) = spongemin
          end if
       end do

    else

       do k = lo(3),hi(3)
          z = prob_lo_z + (dble(k)+HALF)*dx(3)
          if(z .le. topsponge_lo_r) then
             sponge(:,:,k) = ONE
          else if (z .le. topsponge_hi_r) then
             sponge(:,:,k) = HALF*(ONE-spongemin) &
                  * cos(M_PI*(z-topsponge_lo_r)/(topsponge_hi_r-topsponge_lo_r)) &
                  + HALF*(ONE+spongemin)
          else
             sponge(:,:,k) = spongemin
          end if
       end do

    end if

   end subroutine mk_sponge_3d

end module sponge_module
