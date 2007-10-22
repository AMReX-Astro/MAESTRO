module slope_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_module

  implicit none

contains

      subroutine slopex_2d(s,slx,lo,ng,nvar,bc,slope_order)

      integer        , intent(in   ) :: lo(2),ng,nvar
      real(kind=dp_t), intent(in   ) ::   s(lo(1)-ng:, lo(2)-ng:,:)
      real(kind=dp_t), intent(  out) :: slx(lo(1)- 1:, lo(2)- 1:,:) 
      integer, intent(in) ::  bc(:,:,:)
      integer, intent(in) ::  slope_order

!     Local variables
      integer :: hi(2)
      integer :: is,js,ie,je
      integer :: i,j,iv
      integer :: cen,lim,flag,fromm
      real(kind=dp_t) del,slim,sflag
      real(kind=dp_t) dpls,dmin,ds

      real(kind=dp_t), allocatable :: dxscr(:,:)

      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flag = 3 )
      parameter( fromm = 4 )

      hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

      is = lo(1)
      js = lo(2)
      ie = hi(1)
      je = hi(2)

      allocate(dxscr(is-2:ie+2,4))

!     HERE DOING 1ST ORDER
      if (slope_order .eq. 0) then
        slx = zero

!     HERE DOING 2ND ORDER
      else if (slope_order .eq. 2) then

        do iv=1,nvar 
          do j = js-1,je+1 
            do i = is-1,ie+1 
              del = half*(s(i+1,j,iv) - s(i-1,j,iv))
              dpls = two*(s(i+1,j,iv) - s(i  ,j,iv))
              dmin = two*(s(i  ,j,iv) - s(i-1,j,iv))
              slim = min(abs(dpls), abs(dmin))
              slim = merge(slim, zero, dpls*dmin.gt.ZERO)
              sflag = sign(one,del)
              slx(i,j,iv)= sflag*min(slim,abs(del))
            enddo

            if (bc(1,1,iv) .eq. EXT_DIR  .or. bc(1,1,iv) .eq. HOEXTRAP) then

              slx(is-1,j,iv) = zero
              del = (s(is+1,j,iv)+three*s(is,j,iv)- &
                     four*s(is-1,j,iv) ) * third
              dpls = two*(s(is+1,j,iv) - s(is  ,j,iv))
              dmin = two*(s(is  ,j,iv) - s(is-1,j,iv))
              slim = min(abs(dpls), abs(dmin))
              slim = merge(slim, zero, dpls*dmin.gt.ZERO)
              sflag = sign(one,del)
              slx(is,j,iv)= sflag*min(slim,abs(del))

            endif

            if (bc(1,2,iv) .eq. EXT_DIR  .or. bc(1,2,iv) .eq. HOEXTRAP) then

              slx(ie+1,j,iv) = zero
              del = -(s(ie-1,j,iv)+three*s(ie,j,iv)- &
                      four*s(ie+1,j,iv) ) * third
              dpls = two*(s(ie  ,j,iv) - s(ie-1,j,iv))
              dmin = two*(s(ie+1,j,iv) - s(ie  ,j,iv))
              slim = min(abs(dpls), abs(dmin))
              slim = merge(slim, zero, dpls*dmin.gt.ZERO)
              sflag = sign(one,del)
              slx(ie,j,iv)= sflag*min(slim,abs(del))

            endif
          enddo
          enddo

      else 

!     HERE DOING 4TH ORDER
      do iv=1,nvar 
        do j = js-1,je+1 

          do i = is-2,ie+2 
            dxscr(i,cen) = half*(s(i+1,j,iv)-s(i-1,j,iv))
            dmin = two*(s(i  ,j,iv)-s(i-1,j,iv))
            dpls = two*(s(i+1,j,iv)-s(i  ,j,iv))
            dxscr(i,lim)= min(abs(dmin),abs(dpls))
            dxscr(i,lim) = merge(dxscr(i,lim),zero,dpls*dmin.gt.ZERO)
            dxscr(i,flag) = sign(one,dxscr(i,cen))
            dxscr(i,fromm)= dxscr(i,flag)*min(dxscr(i,lim), &
                            abs(dxscr(i,cen)))
          enddo

          do i = is-1,ie+1 
            ds = two * two3rd * dxscr(i,cen) - &
                 sixth * (dxscr(i+1,fromm) + dxscr(i-1,fromm)) 
            slx(i,j,iv) = dxscr(i,flag)*min(abs(ds),dxscr(i,lim))
          enddo

          if (bc(1,1,iv) .eq. EXT_DIR  .or. bc(1,1,iv) .eq. HOEXTRAP) then

            slx(is-1,j,iv) = zero

            del = -sixteen/fifteen*s(is-1,j,iv) + half*s(is,j,iv) + &
                            two3rd*s(is+1,j,iv) - tenth*s(is+2,j,iv)
            dmin = two*(s(is  ,j,iv)-s(is-1,j,iv))
            dpls = two*(s(is+1,j,iv)-s(is  ,j,iv))
            slim = min(abs(dpls), abs(dmin))
            slim = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag = sign(one,del)
            slx(is,j,iv)= sflag*min(slim,abs(del))

!           Recalculate the slope at is+1 using the revised dxscr(is,fromm)
            dxscr(is,fromm) = slx(is,j,iv)
            ds = two * two3rd * dxscr(is+1,cen) - &
                 sixth * (dxscr(is+2,fromm) + dxscr(is,fromm))
            slx(is+1,j,iv) = dxscr(is+1,flag)*min(abs(ds),dxscr(is+1,lim))

          endif

          if (bc(1,2,iv) .eq. EXT_DIR  .or. bc(1,2,iv) .eq. HOEXTRAP) then

            slx(ie+1,j,iv) = zero

            del = -( -sixteen/fifteen*s(ie+1,j,iv) + half*s(ie,j,iv) +  &
                               two3rd*s(ie-1,j,iv) - tenth*s(ie-2,j,iv) )
            dmin = two*(s(ie  ,j,iv)-s(ie-1,j,iv))
            dpls = two*(s(ie+1,j,iv)-s(ie  ,j,iv))
            slim = min(abs(dpls), abs(dmin))
            slim = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag = sign(one,del)
            slx(ie,j,iv)= sflag*min(slim,abs(del))

!           Recalculate the slope at ie-1 using the revised dxscr(ie,fromm)
            dxscr(ie,fromm) = slx(ie,j,iv)
            ds = two * two3rd * dxscr(ie-1,cen) - &
                 sixth * (dxscr(ie-2,fromm) + dxscr(ie,fromm))
            slx(ie-1,j,iv) = dxscr(ie-1,flag)*min(abs(ds),dxscr(ie-1,lim))

          endif
        enddo
      enddo

      endif

      end subroutine slopex_2d

      subroutine slopey_2d(s,sly,lo,ng,nvar,bc,slope_order)


      integer, intent(in) :: lo(:),ng,nvar
      integer, intent(in) :: bc(:,:,:)
      integer, intent(in) :: slope_order

      real(kind=dp_t), intent( in) ::     s(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(out) ::   sly(lo(1)- 1:,lo(2)- 1:,:)
      real(kind=dp_t), allocatable :: dyscr(:,:)

      real(kind=dp_t) :: dpls,dmin,ds
      real(kind=dp_t) :: del,slim,sflag
      integer :: hi(2)
      integer :: is,js,ie,je,i,j,iv

      integer cen,lim,flag,fromm
      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flag = 3 )
      parameter( fromm = 4 )

      hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

      is = lo(1)
      js = lo(2)
      ie = hi(1)
      je = hi(2)

      allocate(dyscr(lo(2)-2:hi(2)+2,4))

!     HERE DOING 1ST ORDER
      if (slope_order .eq. 0) then
        sly = zero

!     HERE DOING 2ND ORDER
      else if (slope_order .eq. 2) then

        do iv=1,nvar 
          do j = js-1,je+1 
            do i = is-1,ie+1 

              del  = half*(s(i,j+1,iv) - s(i,j-1,iv))
              dpls = two *(s(i,j+1,iv) - s(i,j  ,iv))
              dmin = two *(s(i,j  ,iv) - s(i,j-1,iv))
              slim = min(abs(dpls),abs(dmin))
              slim = merge(slim, zero, dpls*dmin.gt.ZERO)
              sflag = sign(one,del)
              sly(i,j,iv)= sflag*min(slim,abs(del))

            enddo
          enddo

          if (bc(2,1,iv) .eq. EXT_DIR .or. bc(2,1,iv) .eq. HOEXTRAP) then

            do i = is-1,ie+1 
              sly(i,js-1,iv) = zero
              del = (s(i,js+1,iv)+three*s(i,js,iv)- &
                     four*s(i,js-1,iv)) * third
              dpls = two*(s(i,js+1,iv) - s(i,js  ,iv))
              dmin = two*(s(i,js  ,iv) - s(i,js-1,iv))
              slim = min(abs(dpls), abs(dmin))
              slim = merge(slim, zero, dpls*dmin.gt.ZERO)
              sflag = sign(one,del)
              sly(i,js,iv)= sflag*min(slim,abs(del))
            enddo

          endif

          if (bc(2,2,iv) .eq. EXT_DIR .or. bc(2,2,iv) .eq. HOEXTRAP) then

            do i = is-1, ie+1 
              sly(i,je+1,iv) = zero
              del = -(s(i,je-1,iv)+three*s(i,je,iv)- &
                      four*s(i,je+1,iv)) * third
              dpls = two*(s(i,je+1,iv) - s(i,je ,iv))
              dmin = two*(s(i,je  ,iv) - s(i,je-1,iv))
              slim = min(abs(dpls), abs(dmin))
              slim = merge(slim, zero, dpls*dmin.gt.ZERO)
              sflag = sign(one,del)
              sly(i,je,iv)= sflag*min(slim,abs(del))
            enddo

          endif
        enddo

      else 

!     HERE DOING 4TH ORDER

      do iv=1,nvar 
        do i = is-1,ie+1 
          do j = js-2,je+2 
            dyscr(j,cen) = half*(s(i,j+1,iv)-s(i,j-1,iv))
            dmin = two*(s(i,j  ,iv)-s(i,j-1,iv))
            dpls = two*(s(i,j+1,iv)-s(i,j  ,iv))
            dyscr(j,lim)  = min(abs(dmin),abs(dpls))
            dyscr(j,lim)  = merge(dyscr(j,lim),zero,dpls*dmin.gt.ZERO)
            dyscr(j,flag) = sign(one,dyscr(j,cen))
            dyscr(j,fromm)= dyscr(j,flag)*min(dyscr(j,lim),abs(dyscr(j,cen)))
          enddo

          do j = js-1,je+1 
            ds = two * two3rd * dyscr(j,cen) -  &
                 sixth * (dyscr(j+1,fromm) + dyscr(j-1,fromm))
            sly(i,j,iv) = dyscr(j,flag)*min(abs(ds),dyscr(j,lim))
          enddo

          if (bc(2,1,iv) .eq. EXT_DIR .or. bc(2,1,iv) .eq. HOEXTRAP) then

            sly(i,js-1,iv) = zero
            del = -sixteen/fifteen*s(i,js-1,iv) +  half*s(i,js ,iv) +  &
                            two3rd*s(i,js+1,iv) - tenth*s(i,js+2,iv)
            dmin = two*(s(i,js  ,iv)-s(i,js-1,iv))
            dpls = two*(s(i,js+1,iv)-s(i,js  ,iv))
            slim = min(abs(dpls), abs(dmin))
            slim = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag = sign(one,del)
            sly(i,js,iv)= sflag*min(slim,abs(del))

!           Recalculate the slope at js+1 using the revised dyscr(js,fromm)
            dyscr(js,fromm) = sly(i,js,iv)
            ds = two * two3rd * dyscr(js+1,cen) - &
                 sixth * (dyscr(js+2,fromm) + dyscr(js,fromm))
            sly(i,js+1,iv) = dyscr(js+1,flag)*min(abs(ds),dyscr(js+1,lim))

          endif

          if (bc(2,2,iv) .eq. EXT_DIR .or. bc(2,2,iv) .eq. HOEXTRAP) then

            sly(i,je+1,iv) = zero
            del = -( -sixteen/fifteen*s(i,je+1,iv) +  half*s(i,je  ,iv) + &
                               two3rd*s(i,je-1,iv) - tenth*s(i,je-2,iv) )
            dmin = two*(s(i,je ,iv)-s(i,je-1,iv))
            dpls = two*(s(i,je+1,iv)-s(i,je ,iv))
            slim = min(abs(dpls), abs(dmin))
            slim = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag = sign(one,del)
            sly(i,je,iv)= sflag*min(slim,abs(del))

!           Recalculate the slope at js+1 using the revised dyscr(js,fromm)
            dyscr(je,fromm) = sly(i,je,iv)
            ds = two * two3rd * dyscr(je-1,cen) -  &
                 sixth * (dyscr(je-2,fromm) + dyscr(je,fromm))
            sly(i,je-1,iv) = dyscr(je-1,flag)*min(abs(ds),dyscr(je-1,lim))

          endif

        enddo
      enddo

      endif

      end subroutine slopey_2d

      subroutine slopez_3d(s,slz,lo,ng,nvar,bc,slope_order)

      integer, intent(in) :: lo(:),ng,nvar
      integer, intent(in) :: bc(:,:,:)
      integer, intent(in) :: slope_order

      real(kind=dp_t), intent( in) ::     s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real(kind=dp_t), intent(out) ::   slz(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real(kind=dp_t), allocatable :: dzscr(:,:)

      real(kind=dp_t) :: dpls,dmin,ds
      real(kind=dp_t) :: del,slim,sflag
      integer :: hi(3)
      integer :: is,js,ks,ie,je,ke,i,j,k,iv

      integer cen,lim,flag,fromm
      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flag = 3 )
      parameter( fromm = 4 )

      hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
      hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

      is = lo(1)
      js = lo(2)
      ks = lo(3)
      ie = hi(1)
      je = hi(2)
      ke = hi(3)

      allocate(dzscr(lo(3)-2:hi(3)+2,4))

!     HERE DOING 1ST ORDER
      if (slope_order .eq. 0) then
        slz = zero

!     HERE DOING 2ND ORDER
      else if (slope_order .eq. 2) then

        do iv=1,nvar 
          do k = ks-1,ke+1 
          do j = js-1,je+1 
            do i = is-1,ie+1 

              del  = half*(s(i,j,k+1,iv) - s(i,j,k-1,iv))
              dpls = two *(s(i,j,k+1,iv) - s(i,j,k  ,iv))
              dmin = two *(s(i,j,k  ,iv) - s(i,j,k-1,iv))
              slim = min(abs(dpls),abs(dmin))
              slim = merge(slim, zero, dpls*dmin.gt.ZERO)
              sflag = sign(one,del)
              slz(i,j,k,iv)= sflag*min(slim,abs(del))

            enddo
          enddo
          enddo

          if (bc(3,1,iv) .eq. EXT_DIR .or. bc(3,1,iv) .eq. HOEXTRAP) then
            do j = js-1,je+1 
            do i = is-1,ie+1 
              slz(i,j,ks-1,iv) = zero
              del = (s(i,j,ks+1,iv)+three*s(i,j,ks,iv)- &
                     four*s(i,j,ks-1,iv)) * third
              dpls = two*(s(i,j,ks+1,iv) - s(i,j,ks  ,iv))
              dmin = two*(s(i,j,ks  ,iv) - s(i,j,ks-1,iv))
              slim = min(abs(dpls), abs(dmin))
              slim = merge(slim, zero, dpls*dmin.gt.ZERO)
              sflag = sign(one,del)
              slz(i,j,ks,iv)= sflag*min(slim,abs(del))
            enddo
            enddo

          endif

          if (bc(3,2,iv) .eq. EXT_DIR .or. bc(3,2,iv) .eq. HOEXTRAP) then

            do j = js-1, je+1 
            do i = is-1, ie+1 
              slz(i,j,ke+1,iv) = zero
              del = -(s(i,j,ke-1,iv)+three*s(i,j,ke,iv)- &
                      four*s(i,j,ke+1,iv)) * third
              dpls = two*(s(i,j,ke+1,iv) - s(i,j,ke ,iv))
              dmin = two*(s(i,j,ke  ,iv) - s(i,j,ke-1,iv))
              slim = min(abs(dpls), abs(dmin))
              slim = merge(slim, zero, dpls*dmin.gt.ZERO)
              sflag = sign(one,del)
              slz(i,j,ke,iv)= sflag*min(slim,abs(del))
            enddo
            enddo

          endif
        enddo

      else 

!     HERE DOING 4TH ORDER

      do iv=1,nvar 
        do j = js-1,je+1 
        do i = is-1,ie+1 
          do k = ks-2,ke+2 
            dzscr(k,cen) = half*(s(i,j,k+1,iv)-s(i,j,k-1,iv))
            dmin = two*(s(i,j,k  ,iv)-s(i,j,k-1,iv))
            dpls = two*(s(i,j,k+1,iv)-s(i,j,k  ,iv))
            dzscr(k,lim)  = min(abs(dmin),abs(dpls))
            dzscr(k,lim)  = merge(dzscr(k,lim),zero,dpls*dmin.gt.ZERO)
            dzscr(k,flag) = sign(one,dzscr(k,cen))
            dzscr(k,fromm)= dzscr(k,flag)*min(dzscr(k,lim),abs(dzscr(k,cen)))
          enddo

          do k = ks-1,ke+1 
            ds = two * two3rd * dzscr(k,cen) -  &
                 sixth * (dzscr(k+1,fromm) + dzscr(k-1,fromm))
            slz(i,j,k,iv) = dzscr(k,flag)*min(abs(ds),dzscr(k,lim))
          enddo

          if (bc(3,1,iv) .eq. EXT_DIR .or. bc(3,1,iv) .eq. HOEXTRAP) then

            slz(i,j,ks-1,iv) = zero
            del = -sixteen/fifteen*s(i,j,ks-1,iv) +  half*s(i,j,ks ,iv) +  &
                            two3rd*s(i,j,ks+1,iv) - tenth*s(i,j,ks+2,iv)
            dmin = two*(s(i,j,ks  ,iv)-s(i,j,ks-1,iv))
            dpls = two*(s(i,j,ks+1,iv)-s(i,j,ks  ,iv))
            slim = min(abs(dpls), abs(dmin))
            slim = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag = sign(one,del)
            slz(i,j,ks,iv)= sflag*min(slim,abs(del))

!           Recalculate the slope at js+1 using the revised dzscr(js,fromm)
            dzscr(ks,fromm) = slz(i,j,ks,iv)
            ds = two * two3rd * dzscr(ks+1,cen) - &
                 sixth * (dzscr(ks+2,fromm) + dzscr(ks,fromm))
            slz(i,j,ks+1,iv) = dzscr(ks+1,flag)*min(abs(ds),dzscr(ks+1,lim))

          endif

          if (bc(3,2,iv) .eq. EXT_DIR .or. bc(3,2,iv) .eq. HOEXTRAP) then

            slz(i,j,ke+1,iv) = zero
            del = -( -sixteen/fifteen*s(i,j,ke+1,iv) +  half*s(i,j,ke  ,iv) + &
                               two3rd*s(i,j,ke-1,iv) - tenth*s(i,j,ke-2,iv) )
            dmin = two*(s(i,j,ke ,iv)-s(i,j,ke-1,iv))
            dpls = two*(s(i,j,ke+1,iv)-s(i,j,ke ,iv))
            slim = min(abs(dpls), abs(dmin))
            slim = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag = sign(one,del)
            slz(i,j,ke,iv)= sflag*min(slim,abs(del))

!           Recalculate the slope at ks+1 using the revised dzscr(ks,fromm)
            dzscr(ke,fromm) = slz(i,j,ke,iv)
            ds = two * two3rd * dzscr(ke-1,cen) -  &
                 sixth * (dzscr(ke-2,fromm) + dzscr(ke,fromm))
            slz(i,j,ke-1,iv) = dzscr(ke-1,flag)*min(abs(ds),dzscr(ke-1,lim))

          endif
        enddo
        enddo
      enddo

      endif

      end subroutine slopez_3d

end module slope_module
