! a module for storing the geometric information so we don't have to pass it
 
module sponge_module
 
  use bl_types
  use bl_constants_module
  use multifab_module
  use geometry
  use variables
 
  implicit none
 
  real(dp_t), save :: r_sp, r_md, r_tp
  real(dp_t), save :: alpha
 
contains

  subroutine init_sponge (s0,anelastic_cutoff)

    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    real(kind=dp_t), intent(in   ) :: anelastic_cutoff

    real (kind = dp_t) :: x, y, z, r
    real (kind = dp_t) :: r_top
    integer            :: j,nr

    nr = size(s0,dim=1)
    r_top = dble(nr) * dr

    r_sp = r_top
    do j = 0, nr-1
       r = (dble(j)+HALF) * dr
       if (s0(j,rho_comp) < 10.d0*anelastic_cutoff) then
          r_sp = r
          exit
       endif
    enddo

      r_md = r_top
      do j = 0,nr-1
         r = (dble(j)+HALF) * dr
         if (s0(j,rho_comp) < anelastic_cutoff) then
            r_md = r
            exit
         endif
      enddo

      r_tp = 2.d0 * r_md - r_sp

!     alpha = 100.d0
      alpha = 1000.d0

      print *, 'sponge : r_sp, r_md, r_tp, r_top ', r_sp, r_md, r_tp, r_top

  end subroutine init_sponge

  subroutine make_sponge (sponge,dx,dt)

      type(multifab) , intent(inout) :: sponge
      real(kind=dp_t), intent(in   ) :: dx(:),dt

      ! Local variables
      real(kind=dp_t), pointer::  sp(:,:,:,:)
      integer :: i,lo(sponge%dim),hi(sponge%dim),dm

      dm = sponge%dim

      do i = 1, sponge%nboxes
         if ( multifab_remote(sponge, i) ) cycle
          sp => dataptr(sponge, i)
          lo =  lwb(get_box(sponge, i))
          hi =  upb(get_box(sponge, i))
         select case (dm)
            case (2)
              call mk_sponge_2d(sp(:,:,1,1),lo,hi,dx,dt)
            case (3)
              call mk_sponge_3d(sp(:,:,:,1),lo,hi,dx,dt)
         end select
      end do

  end subroutine make_sponge

  subroutine mk_sponge_2d(sponge,lo,hi,dx,dt)

      integer        , intent(in   ) ::  lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: sponge(lo(1):,lo(2):)
      real(kind=dp_t), intent(in   ) ::     dx(:),dt

      integer         :: i,j
      real(kind=dp_t) :: y,smdamp

      sponge = ONE

      do j = lo(2),hi(2)
        y = (dble(j)+HALF)*dx(2)

          if (y >= r_sp) then
             if (y < r_tp) then
               smdamp = HALF*(ONE - cos(M_PI*(y - r_sp)/(r_tp - r_sp)))
             else
               smdamp = ONE
             endif
             sponge(:,j) = ONE / (ONE + dt * smdamp* alpha)
          endif

      end do

  end subroutine mk_sponge_2d

  subroutine mk_sponge_3d(sponge,lo,hi,dx,dt)

      integer        , intent(in   ) ::  lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: sponge(lo(1):,lo(2):,lo(3):)
      real(kind=dp_t), intent(in   ) ::     dx(:),dt

      integer         :: i,j,k
      real(kind=dp_t) :: x,y,z,r,smdamp

      sponge = ONE

      if (spherical .eq. 0) then
        do k = lo(3),hi(3)
          z = (dble(k)+HALF)*dx(3)
          if (z >= r_sp) then
            if (z < r_tp) then
              smdamp = HALF*(ONE - cos(M_PI*(z - r_sp)/(r_tp - r_sp)))
            else
              smdamp = ONE
            endif
            sponge(:,:,k) = ONE / (ONE + dt * smdamp* alpha)
          end if
        end do

      else

        do k = lo(3),hi(3)
          z = (dble(k)+HALF)*dx(3)
          do j = lo(2),hi(2)
            y = (dble(j)+HALF)*dx(2)
            do i = lo(1),hi(1)
              x = (dble(i)+HALF)*dx(1)
  
              r = sqrt( (x-center(1))**2 + (y-center(2))**2 + (z-center(3))**2 )
              if (r >= r_sp) then
                 if (r < r_tp) then
                   smdamp = HALF*(ONE - cos(M_PI*(r - r_sp)/(r_tp - r_sp)))
                 else
                   smdamp = ONE
                 endif
                 sponge(i,j,k) = ONE / (ONE + dt * smdamp* alpha)
              endif

            end do
          end do
        end do

      end if

  end subroutine mk_sponge_3d

end module sponge_module
