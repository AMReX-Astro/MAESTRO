module bds_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: bds, bds_velpred

contains

  subroutine bds(s,sedge,umac,force,dx,dt,is_vel,the_bc_level, &
                 start_scomp,start_bccomp,num_comp,is_conservative,mla)

    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: force(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    logical        , intent(in   ) :: is_vel
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    logical        , intent(in   ) :: is_conservative
    type(ml_layout), intent(in   ) :: mla

    ! Local
    type(bl_prof_timer), save :: bpt

    ! this will hold slx, sly, slxy, etc.
    type(multifab) :: slope(mla%nlevel)

    real(kind=dp_t), pointer ::    sop(:,:,:,:)
    real(kind=dp_t), pointer ::   sepx(:,:,:,:)
    real(kind=dp_t), pointer ::   sepy(:,:,:,:)
    real(kind=dp_t), pointer ::   sepz(:,:,:,:)
    real(kind=dp_t), pointer :: slopep(:,:,:,:)
    real(kind=dp_t), pointer ::  umacp(:,:,:,:)
    real(kind=dp_t), pointer ::  vmacp(:,:,:,:)
    real(kind=dp_t), pointer ::  wmacp(:,:,:,:)
    real(kind=dp_t), pointer ::     fp(:,:,:,:)

    integer :: dm,ng_s,ng_c,ng_u,ng_se,ng_f
    integer :: n,i,comp,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    call build(bpt, "bds")

    nlevs = mla%nlevel
    dm = mla%dim

    if (dm .eq. 2) then
       ! 3 components and 1 ghost cell
       ! component 1 = slx
       ! component 2 = sly
       ! component 3 = slxy
       do n=1,nlevs
          call multifab_build(slope(n),mla%la(n),3,1)
       end do
    else if (dm .eq. 3) then
       ! 7 components and 1 ghost cell
       ! component 1 = slx
       ! component 2 = sly
       ! component 3 = slz
       ! component 4 = slxy
       ! component 5 = slxz
       ! component 6 = slyz
       ! component 7 = slxyz
       do n=1,nlevs
          call multifab_build(slope(n),mla%la(n),7,1)
       end do
    end if

    ng_s = s(1)%ng
    ng_c = slope(1)%ng
    ng_u = umac(1,1)%ng
    ng_f = force(1)%ng
    ng_se = sedge(1,1)%ng

    do comp=start_scomp,start_scomp+num_comp-1
       do n=1,nlevs
          do i = 1, nfabs(s(n))

             sop    => dataptr(s(n), i)
             sepx   => dataptr(sedge(n,1), i)
             sepy   => dataptr(sedge(n,2), i)
             slopep => dataptr(slope(n), i)
             umacp  => dataptr(umac(n,1), i)
             vmacp  => dataptr(umac(n,2), i)
             fp     => dataptr(force(n) ,i)
             lo =  lwb(get_box(s(n), i))
             hi =  upb(get_box(s(n), i))
             select case (dm)
             case (2)
                call bdsslope_2d(lo, hi, sop(:,:,1,comp), ng_s, &
                                 slopep(:,:,1,:), ng_c, dx(n,:)) 

                call bdsconc_2d(lo, hi, sop(:,:,1,comp), ng_s, &
                                slopep(:,:,1,:), ng_c, &
                                umacp(:,:,1,1), vmacp(:,:,1,1), ng_u, &
                                fp(:,:,1,comp), ng_f, &
                                sepx(:,:,1,comp), sepy(:,:,1,comp), ng_se, &
                                dx(n,:), dt, is_conservative)
             case (3)
                wmacp  => dataptr(umac(n,3), i)
                sepz   => dataptr(sedge(n,3), i)
                call bdsslope_3d(lo, hi, sop(:,:,:,comp), ng_s, &
                                 slopep(:,:,:,:), ng_c, dx(n,:))

                call bdsconc_3d(lo, hi, sop(:,:,:,comp), ng_s, &
                                slopep(:,:,:,:), ng_c, &
                                umacp(:,:,:,1), vmacp(:,:,:,1), wmacp(:,:,:,1), ng_u, &
                                fp(:,:,:,comp), ng_f, &
                                sepx(:,:,:,comp), sepy(:,:,:,comp), sepz(:,:,:,comp), ng_se, &
                                dx(n,:), dt, is_conservative)
             end select
          end do ! loop over boxes
       end do ! loop over levels
    end do ! loop over components

    do n=1,nlevs
       call multifab_destroy(slope(n))
    end do

    call destroy(bpt)

  end subroutine bds

  subroutine bdsslope_2d(lo,hi,s,ng_s,slope,ng_c,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_c
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :)
    real(kind=dp_t), intent(inout) ::  slope(lo(1)-ng_c :,lo(2)-ng_c :,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local variables
    real(kind=dp_t), allocatable :: sint(:,:)

    real(kind=dp_t) :: diff(4)
    real(kind=dp_t) :: smin(4)
    real(kind=dp_t) :: smax(4)
    real(kind=dp_t) :: sc(4)

    real(kind=dp_t) :: hx,hy,eps
    real(kind=dp_t) :: sumloc,redfac,redmax,div,kdp,sumdif,sgndif
    integer         :: i,j,ll,mm

    ! nodal with one ghost cell
    allocate(sint(lo(1)-ng_c:hi(1)+ng_c+1,lo(2)-ng_c:hi(2)+ng_c+1))

    hx = dx(1)
    hy = dx(2)

    eps = 1.d-10

    ! bicubic interpolation to corner points
    ! (i,j,k) refers to lower corner of cell
    do j = lo(2)-ng_c,hi(2)+ng_c+1
       do i = lo(1)-ng_c,hi(1)+ng_c+1
          sint(i,j) = (s(i-2,j-2) + s(i-2,j+1) + s(i+1,j-2) + s(i+1,j+1) &
               - 7.d0*(s(i-2,j-1) + s(i-2,j  ) + s(i-1,j-2) + s(i  ,j-2) + & 
                       s(i-1,j+1) + s(i  ,j+1) + s(i+1,j-1) + s(i+1,j  )) &
              + 49.d0*(s(i-1,j-1) + s(i  ,j-1) + s(i-1,j  ) + s(i  ,j  )) ) / 144.d0
       enddo
    enddo

    do j = lo(2)-ng_c,hi(2)+ng_c
       do i = lo(1)-ng_c,hi(1)+ng_c

          ! compute initial estimates of slopes from unlimited corner points

          ! sx
          slope(i,j,1) = 0.5d0*(sint(i+1,j+1) + sint(i+1,j  ) - &
                                sint(i  ,j+1) - sint(i  ,j  ) ) / hx

          ! sy
          slope(i,j,2) = 0.5d0*(sint(i+1,j+1) - sint(i+1,j  ) + &
                                sint(i  ,j+1) - sint(i  ,j  ) ) / hy

          ! sxy
          slope(i,j,3) = ( sint(i+1,j+1) - sint(i+1,j  ) &
                          -sint(i  ,j+1) + sint(i  ,j  ) ) / (hx*hy)

          ! ++ / sint(i+1,j+1)
          sc(4) = s(i,j) + 0.5d0*(hx*slope(i,j,1) + hy*slope(i,j,2))  &
               + 0.25d0*hx*hy*slope(i,j,3)

          ! +- / sint(i+1,j  )
          sc(3) = s(i,j) + 0.5d0*(hx*slope(i,j,1) - hy*slope(i,j,2))  &
               - 0.25d0*hx*hy*slope(i,j,3)

          ! -+ / sint(i  ,j+1)
          sc(2) = s(i,j) - 0.5d0*(hx*slope(i,j,1) - hy*slope(i,j,2)) &
               - 0.25d0*hx*hy*slope(i,j,3)

          ! -- / sint(i  ,j  )
          sc(1) = s(i,j) - 0.5d0*(hx*slope(i,j,1) + hy*slope(i,j,2)) &
               + 0.25d0*hx*hy*slope(i,j,3)

          ! enforce max/min bounds
          smin(4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
          smax(4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))

          smin(3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
          smax(3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))

          smin(2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
          smax(2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))

          smin(1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
          smax(1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))

          do mm=1,4
             sc(mm) = max(min(sc(mm), smax(mm)), smin(mm))
          enddo

          ! iterative loop
          do ll = 1,3 
             sumloc = 0.25d0*(sc(4) + sc(3) + sc(2) + sc(1))
             sumdif = (sumloc - s(i,j))*4.d0
             sgndif = sign(1.d0,sumdif)

             do mm=1,4
                diff(mm) = (sc(mm) - s(i,j))*sgndif
             enddo

             kdp = 0

             do mm=1,4
                if (diff(mm) .gt. eps) then
                   kdp = kdp+1
                end if
             end do

             do mm = 1,4 
                if (kdp.lt.1) then 
                   div = 1.d0
                else
                   div = dble(kdp)
                end if

                if (diff(mm).gt.eps) then
                   redfac = sumdif*sgndif/div
                   kdp = kdp-1
                else
                   redfac = 0.d0
                end if

                if (sgndif .gt. 0.d0) then
                   redmax = sc(mm) - smin(mm)
                else
                   redmax = smax(mm) - sc(mm)
                end if

                redfac = min(redfac,redmax)
                sumdif = sumdif - redfac*sgndif
                sc(mm) = sc(mm) - redfac*sgndif
             enddo
          enddo

          ! final slopes

          ! sx
          slope(i,j,1) = 0.5d0*( sc(4) + sc(3) &
                                -sc(1) - sc(2))/hx

          ! sy
          slope(i,j,2) = 0.5d0*( sc(4) + sc(2) &
                                -sc(1) - sc(3))/hy

          ! sxy
          slope(i,j,3) = ( sc(1) + sc(4) &
                          -sc(2) - sc(3) ) / (hx*hy)

       enddo
    enddo

    deallocate(sint)

  end subroutine bdsslope_2d

  subroutine bdsslope_3d(lo,hi,s,ng_s,slope,ng_c,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_c
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t), intent(inout) :: slope(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local variables
    real(kind=dp_t), allocatable :: sint(:,:,:)

    real(kind=dp_t) :: diff(8)
    real(kind=dp_t) :: smin(8)
    real(kind=dp_t) :: smax(8)
    real(kind=dp_t) :: sc(8)

    real(kind=dp_t) :: c1,c2,c3,c4
    real(kind=dp_t) :: hx,hy,hz,eps
    real(kind=dp_t) :: sumloc,redfac,redmax,div,kdp,sumdif,sgndif
    integer         :: i,j,k,ll,mm

    ! nodal with one ghost cell
    allocate(sint(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+2))

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    eps = 1.d-10

    c1 = (343.d0/1728.d0)
    c2 = (49.d0 /1728.d0)
    c3 = (7.d0  /1728.d0)
    c4 = (1.d0  /1728.d0)

    ! tricubic interpolation to corner points
    ! (i,j,k) refers to lower corner of cell
    do k = lo(3)-1,hi(3)+2
       do j = lo(2)-1,hi(2)+2
          do i = lo(1)-1,hi(1)+2
             sint(i,j,k) = c1*( s(i  ,j  ,k  ) + s(i-1,j  ,k  ) + s(i  ,j-1,k  ) &
                               +s(i  ,j  ,k-1) + s(i-1,j-1,k  ) + s(i-1,j  ,k-1) &
                               +s(i  ,j-1,k-1) + s(i-1,j-1,k-1) ) &
                          -c2*( s(i-1,j  ,k+1) + s(i  ,j  ,k+1) + s(i-1,j-1,k+1) &
                               +s(i  ,j-1,k+1) + s(i-1,j+1,k  ) + s(i  ,j+1,k  ) &
                               +s(i-2,j  ,k  ) + s(i+1,j  ,k  ) + s(i-2,j-1,k  ) &
                               +s(i+1,j-1,k  ) + s(i-1,j-2,k  ) + s(i  ,j-2,k  ) &
                               +s(i-1,j+1,k-1) + s(i  ,j+1,k-1) + s(i-2,j  ,k-1) &
                               +s(i+1,j  ,k-1) + s(i-2,j-1,k-1) + s(i+1,j-1,k-1) &
                               +s(i-1,j-2,k-1) + s(i  ,j-2,k-1) + s(i-1,j  ,k-2) &
                               +s(i  ,j  ,k-2) + s(i-1,j-1,k-2) + s(i  ,j-1,k-2) ) &
                          +c3*( s(i-1,j+1,k+1) + s(i  ,j+1,k+1) + s(i-2,j  ,k+1) &
                               +s(i+1,j  ,k+1) + s(i-2,j-1,k+1) + s(i+1,j-1,k+1) &
                               +s(i-1,j-2,k+1) + s(i  ,j-2,k+1) + s(i-2,j+1,k  ) &
                               +s(i+1,j+1,k  ) + s(i-2,j-2,k  ) + s(i+1,j-2,k  ) &
                               +s(i-2,j+1,k-1) + s(i+1,j+1,k-1) + s(i-2,j-2,k-1) &
                               +s(i+1,j-2,k-1) + s(i-1,j+1,k-2) + s(i  ,j+1,k-2) &
                               +s(i-2,j  ,k-2) + s(i+1,j  ,k-2) + s(i-2,j-1,k-2) &
                               +s(i+1,j-1,k-2) + s(i-1,j-2,k-2) + s(i  ,j-2,k-2) ) &
                          -c4*( s(i-2,j+1,k+1) + s(i+1,j+1,k+1) + s(i-2,j-2,k+1) &
                               +s(i+1,j-2,k+1) + s(i-2,j+1,k-2) + s(i+1,j+1,k-2) &
                               +s(i-2,j-2,k-2) + s(i+1,j-2,k-2) )
          enddo
       enddo
    enddo

    do k = lo(3)-1,hi(3)+1
       do j = lo(2)-1,hi(2)+1
          do i = lo(1)-1,hi(1)+1 

             ! compute initial estimates of slopes from unlimited corner points

             ! sx
             slope(i,j,k,1) = 0.25d0*( ( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  ) &
                                        +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1)) &
                                      -( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  ) &
                                        +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1)) ) / hx

             ! sy
             slope(i,j,k,2) = 0.25d0*( ( sint(i  ,j+1,k  ) + sint(i+1,j+1,k  ) &
                                        +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                                      -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                                        +sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1)) ) / hy

             ! sz
             slope(i,j,k,3) = 0.25d0*( ( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1) &
                                        +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                                      -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                                        +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  )) ) / hz

             ! sxy
             slope(i,j,k,4) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i  ,j  ,k+1) &
                                       +sint(i+1,j+1,k  ) + sint(i+1,j+1,k+1)) &
                                     -( sint(i+1,j  ,k  ) + sint(i+1,j  ,k+1) &
                                       +sint(i  ,j+1,k  ) + sint(i  ,j+1,k+1)) ) / (hx*hy)

             ! sxz
             slope(i,j,k,5) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  ) &
                                       +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1)) &
                                     -( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  ) &
                                       +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1)) ) / (hx*hz)

             ! syz
             slope(i,j,k,6) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                                       +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                                     -( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1) &
                                       +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  )) ) / (hy*hz)

             ! sxyz
             slope(i,j,k,7) = (-sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) + sint(i  ,j+1,k  ) &
                               +sint(i  ,j  ,k+1) - sint(i+1,j+1,k  ) - sint(i+1,j  ,k+1) &
                               -sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) ) / (hx*hy*hz)

             ! +++ / sint(i+1,j+1,k+1)
             sc(8) = s(i,j,k) &
                  +0.5d0*( hx*slope(i,j,k,1)+hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                  +0.25d0*( hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                  +0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! ++- / sint(i+1,j+1,k  )
             sc(7) = s(i,j,k) &
                  +0.5d0*( hx*slope(i,j,k,1)+hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                  +0.25d0*( hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                  -0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! +-+ / sint(i+1,j  ,k+1)
             sc(6) = s(i,j,k) &
                  +0.5d0*( hx*slope(i,j,k,1)-hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                  +0.25d0*(-hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                  -0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! +-- / sint(i+1,j  ,k  )
             sc(5) = s(i,j,k) &
                  +0.5d0*( hx*slope(i,j,k,1)-hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                  +0.25d0*(-hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                  +0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! -++ / sint(i  ,j+1,k+1)
             sc(4) = s(i,j,k) &
                  +0.5d0*(-hx*slope(i,j,k,1)+hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                  +0.25d0*(-hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                  -0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! -+- / sint(i  ,j+1,k  )
             sc(3) = s(i,j,k) &
                  +0.5d0*(-hx*slope(i,j,k,1)+hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                  +0.25d0*(-hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                  +0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! --+ / sint(i  ,j  ,k+1)
             sc(2) = s(i,j,k) &
                  +0.5d0*(-hx*slope(i,j,k,1)-hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                  +0.25d0*( hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                  +0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! ---/ sint(i  ,j  ,k  )
             sc(1) = s(i,j,k) &
                  +0.5d0*(-hx*slope(i,j,k,1)-hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                  +0.25d0*( hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                  -0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! enforce max/min bounds
             smin(8) = min(s(i,j,k),s(i+1,j,k),s(i,j+1,k),s(i,j,k+1), &
                           s(i+1,j+1,k),s(i+1,j,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
             smax(8) = max(s(i,j,k),s(i+1,j,k),s(i,j+1,k),s(i,j,k+1), &
                           s(i+1,j+1,k),s(i+1,j,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))

             smin(7) = min(s(i,j,k-1),s(i+1,j,k-1),s(i,j+1,k-1),s(i,j,k), &
                           s(i+1,j+1,k-1),s(i+1,j,k),s(i,j+1,k),s(i+1,j+1,k))
             smax(7) = max(s(i,j,k-1),s(i+1,j,k-1),s(i,j+1,k-1),s(i,j,k), &
                           s(i+1,j+1,k-1),s(i+1,j,k),s(i,j+1,k),s(i+1,j+1,k))

             smin(6) = min(s(i,j-1,k),s(i+1,j-1,k),s(i,j,k),s(i,j-1,k+1), &
                           s(i+1,j,k),s(i+1,j-1,k+1),s(i,j,k+1),s(i+1,j,k+1))
             smax(6) = max(s(i,j-1,k),s(i+1,j-1,k),s(i,j,k),s(i,j-1,k+1), &
                           s(i+1,j,k),s(i+1,j-1,k+1),s(i,j,k+1),s(i+1,j,k+1))

             smin(5) = min(s(i,j-1,k-1),s(i+1,j-1,k-1),s(i,j,k-1),s(i,j-1,k), &
                           s(i+1,j,k-1),s(i+1,j-1,k),s(i,j,k),s(i+1,j,k))
             smax(5) = max(s(i,j-1,k-1),s(i+1,j-1,k-1),s(i,j,k-1),s(i,j-1,k), &
                           s(i+1,j,k-1),s(i+1,j-1,k),s(i,j,k),s(i+1,j,k))

             smin(4) = min(s(i-1,j,k),s(i,j,k),s(i-1,j+1,k),s(i-1,j,k+1), &
                           s(i,j+1,k),s(i,j,k+1),s(i-1,j+1,k+1),s(i,j+1,k+1))
             smax(4) = max(s(i-1,j,k),s(i,j,k),s(i-1,j+1,k),s(i-1,j,k+1), &
                           s(i,j+1,k),s(i,j,k+1),s(i-1,j+1,k+1),s(i,j+1,k+1))

             smin(3) = min(s(i-1,j,k-1),s(i,j,k-1),s(i-1,j+1,k-1),s(i-1,j,k), &
                           s(i,j+1,k-1),s(i,j,k),s(i-1,j+1,k),s(i,j+1,k))
             smax(3) = max(s(i-1,j,k-1),s(i,j,k-1),s(i-1,j+1,k-1),s(i-1,j,k), &
                           s(i,j+1,k-1),s(i,j,k),s(i-1,j+1,k),s(i,j+1,k))

             smin(2) = min(s(i-1,j-1,k),s(i,j-1,k),s(i-1,j,k),s(i-1,j-1,k+1), &
                           s(i,j,k),s(i,j-1,k+1),s(i-1,j,k+1),s(i,j,k+1))
             smax(2) = max(s(i-1,j-1,k),s(i,j-1,k),s(i-1,j,k),s(i-1,j-1,k+1), &
                           s(i,j,k),s(i,j-1,k+1),s(i-1,j,k+1),s(i,j,k+1))

             smin(1) = min(s(i-1,j-1,k-1),s(i,j-1,k-1),s(i-1,j,k-1),s(i-1,j-1,k), &
                           s(i,j,k-1),s(i,j-1,k),s(i-1,j,k),s(i,j,k))
             smax(1) = max(s(i-1,j-1,k-1),s(i,j-1,k-1),s(i-1,j,k-1),s(i-1,j-1,k), &
                           s(i,j,k-1),s(i,j-1,k),s(i-1,j,k),s(i,j,k))

             do mm=1,8
                sc(mm) = max(min(sc(mm), smax(mm)), smin(mm))
             enddo
             
             ! iterative loop
             do ll = 1,3 
                sumloc = 0.125d0*(sc(1)+sc(2)+sc(3)+sc(4)+sc(5)+sc(6)+sc(7)+sc(8))
                sumdif = (sumloc - s(i,j,k))*8.d0
                sgndif = sign(1.d0,sumdif)

                do mm=1,8
                   diff(mm) = (sc(mm) - s(i,j,k))*sgndif
                enddo

                kdp = 0

                do mm=1,8
                   if (diff(mm) .gt. eps) then
                      kdp = kdp+1
                   end if
                end do

                do mm = 1,8
                   if (kdp.lt.1) then 
                      div = 1.d0
                   else
                      div = dble(kdp)
                   end if

                   if (diff(mm).gt.eps) then
                      redfac = sumdif*sgndif/div
                      kdp = kdp-1
                   else
                      redfac = 0.d0
                   end if

                   if (sgndif .gt. 0.d0) then
                      redmax = sc(mm) - smin(mm)
                   else
                      redmax = smax(mm) - sc(mm)
                   end if

                   redfac = min(redfac,redmax)
                   sumdif = sumdif - redfac*sgndif
                   sc(mm) = sc(mm) - redfac*sgndif
                enddo
             enddo

             ! final slopes

             ! sx
             slope(i,j,k,1) = 0.25d0*( ( sc(5) + sc(7) &
                                        +sc(6) + sc(8)) &
                                      -( sc(1) + sc(3) &
                                        +sc(2) + sc(4)) ) / hx

             ! sy
             slope(i,j,k,2) = 0.25d0*( ( sc(3) + sc(7) &
                                        +sc(4) + sc(8)) &
                                      -( sc(1) + sc(5) &
                                        +sc(2) + sc(6)) ) / hy

             ! sz
             slope(i,j,k,3) = 0.25d0*( ( sc(2) + sc(6) &
                                        +sc(4) + sc(8)) &
                                      -( sc(1) + sc(5) &
                                        +sc(3) + sc(7)) ) / hz

             ! sxy
             slope(i,j,k,4) = 0.5d0*( ( sc(1) + sc(2) &
                                       +sc(7) + sc(8)) &
                                     -( sc(5) + sc(6) &
                                       +sc(3) + sc(4)) ) / (hx*hy)

             ! sxz
             slope(i,j,k,5) = 0.5d0*( ( sc(1) + sc(3) &
                                       +sc(6) + sc(8)) &
                                     -( sc(5) + sc(7) &
                                       +sc(2) + sc(4)) ) / (hx*hz)

             ! syz
             slope(i,j,k,6) = 0.5d0*( ( sc(1) + sc(5) &
                                       +sc(4) + sc(8)) &
                                     -( sc(2) + sc(6) &
                                       +sc(3) + sc(7)) ) / (hy*hz)

             ! sxyz
             slope(i,j,k,7) = (-sc(1) + sc(5) + sc(3) &
                               +sc(2) - sc(7) - sc(6) &
                               -sc(4) + sc(8) ) / (hx*hy*hz)

          enddo
       enddo
    enddo


  end subroutine bdsslope_3d

  subroutine bdsconc_2d(lo,hi,s,ng_s,slope,ng_c,umac,vmac,ng_u,force,ng_f, &
                        sedgex,sedgey,ng_se,dx,dt,is_conservative)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_c,ng_u,ng_f,ng_se
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :)
    real(kind=dp_t), intent(in   ) ::  slope(lo(1)-ng_c :,lo(2)-ng_c :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-ng_u :,lo(2)-ng_u :)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-ng_u :,lo(2)-ng_u :)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:,lo(2)-ng_se:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_se:,lo(2)-ng_se:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_conservative

    ! local variables
    integer i,j,ioff,joff,ll

    real(kind=dp_t), allocatable ::   ux(:,:)
    real(kind=dp_t), allocatable ::   vy(:,:)
    real(kind=dp_t), allocatable :: divu(:,:)


    real(kind=dp_t) :: isign,jsign,hx,hy
    real(kind=dp_t) :: del(2),p1(2),p2(2),p3(2)
    real(kind=dp_t) :: val1,val2,val3
    real(kind=dp_t) :: u,v,gamma
    real(kind=dp_t) :: dt2,dt3,half

    allocate(  ux(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(  vy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(divu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))

    hx = dx(1)
    hy = dx(2)

    dt2 = dt/2.d0
    dt3 = dt/3.d0

    half = 0.5d0

    ! compute cell-centered ux and vy
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          ux(i,j) = (umac(i+1,j) - umac(i,j)) / hx
          vy(i,j) = (vmac(i,j+1) - vmac(i,j)) / hy
          divu(i,j) = ux(i,j) + vy(i,j)
       end do
    end do

    ! compute sedgex on x-faces
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute sedgex without transverse corrections
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (umac(i,j) .gt. 0) then
             isign = 1.d0
             ioff = -1
          else
             isign = -1.d0
             ioff = 0
          endif

          ! centroid of rectangular volume
          del(1) = isign*0.5d0*hx - 0.5d0*umac(i,j)*dt
          del(2) = 0.d0
          call eval_2d(s(i+ioff,j),slope(i+ioff,j,:),del,sedgex(i,j))

          ! source term
          if (is_conservative) then
             sedgex(i,j) = sedgex(i,j)*(1.d0 - dt2*ux(i+ioff,j)) + dt2*force(i+ioff,j)
          else
             sedgex(i,j) = sedgex(i,j)*(1.d0 + dt2*vy(i+ioff,j)) + dt2*force(i+ioff,j)
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute \Gamma^{y+}
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (vmac(i+ioff,j+1) .gt. 0) then
             jsign = 1.d0
             joff = 0
          else
             jsign = -1.d0
             joff = 1
          endif

          u = 0.d0
          if (umac(i,j)*umac(i,j+joff) .gt. 0) then
             u = umac(i,j+joff)
          endif

          p1(1) = isign*0.5d0*hx
          p1(2) = jsign*0.5d0*hy

          p2(1) = isign*0.5d0*hx - umac(i,j)*dt
          p2(2) = jsign*0.5d0*hy

          p3(1) = isign*0.5d0*hx - u*dt
          p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1)*dt

          do ll=1,2
             del(ll) = (p2(ll)+p3(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val1)

          do ll=1,2
             del(ll) = (p1(ll)+p3(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val2)

          do ll=1,2
             del(ll) = (p1(ll)+p2(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val3)

          ! average these centroid values to get the average value
          gamma = (val1+val2+val3)/3.d0

          ! source term
          if (is_conservative) then
             gamma = gamma*(1.d0 - dt3*divu(i+ioff,j+joff))
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! correct sedgex with \Gamma^{y+}
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          gamma = gamma * vmac(i+ioff,j+1)
          sedgex(i,j) = sedgex(i,j) - dt*gamma/(2.d0*hy)
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute \Gamma^{y-}
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          if (vmac(i+ioff,j) .gt. 0) then
             jsign = 1.d0
             joff = -1
          else
             jsign = -1.d0
             joff = 0
          endif

          u = 0.d0
          if (umac(i,j)*umac(i,j+joff) .gt. 0) then
             u = umac(i,j+joff)
          endif

          p1(1) = isign*0.5d0*hx
          p1(2) = jsign*0.5d0*hy

          p2(1) = isign*0.5d0*hx - umac(i,j)*dt
          p2(2) = jsign*0.5d0*hy

          p3(1) = isign*0.5d0*hx - u*dt
          p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j)*dt

          do ll=1,2
             del(ll) = (p2(ll)+p3(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val1)

          do ll=1,2
             del(ll) = (p1(ll)+p3(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val2)

          do ll=1,2
             del(ll) = (p1(ll)+p2(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val3)

          ! average these centroid values to get the average value
          gamma = (val1+val2+val3)/3.d0

          ! source term
          if (is_conservative) then
             gamma = gamma*(1.d0 - dt3*divu(i+ioff,j+joff))
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! correct sedgex with \Gamma^{y-}
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          gamma = gamma * vmac(i+ioff,j)
          sedgex(i,j) = sedgex(i,j) + dt*gamma/(2.d0*hy)

       end do
    end do

    ! compute sedgey on y-faces    
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute sedgey without transverse corrections
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! centroid of rectangular volume
          if (vmac(i,j) .gt. 0) then
             jsign = 1.d0
             joff = -1
          else
             jsign = -1.d0
             joff = 0
          endif

          del(1) = 0.d0
          del(2) = jsign*0.5d0*hy - 0.5d0*vmac(i,j)*dt
          call eval_2d(s(i,j+joff),slope(i,j+joff,:),del,sedgey(i,j))

          ! source term
          if (is_conservative) then
             sedgey(i,j) = sedgey(i,j)*(1.d0 - dt2*vy(i,j+joff)) + dt2*force(i,j+joff)
          else
             sedgey(i,j) = sedgey(i,j)*(1.d0 + dt2*ux(i,j+joff)) + dt2*force(i,j+joff)
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute \Gamma^{x+} without corner corrections
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (umac(i+1,j+joff) .gt. 0) then
             isign = 1.d0
             ioff = 0
          else
             isign = -1.d0
             ioff = 1
          endif

          v = 0.d0
          if (vmac(i,j)*vmac(i+ioff,j) .gt. 0) then
             v = vmac(i+ioff,j)
          endif

          p1(1) = isign*0.5d0*hx
          p1(2) = jsign*0.5d0*hy

          p2(1) = isign*0.5d0*hx
          p2(2) = jsign*0.5d0*hy - vmac(i,j)*dt

          p3(1) = isign*0.5d0*hx - umac(i+1,j+joff)*dt
          p3(2) = jsign*0.5d0*hy - v*dt

          do ll=1,2
             del(ll) = (p2(ll)+p3(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val1)

          do ll=1,2
             del(ll) = (p1(ll)+p3(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val2)

          do ll=1,2
             del(ll) = (p1(ll)+p2(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val3)

          ! average these centroid values to get the average value
          gamma = (val1+val2+val3)/3.d0

          ! source term
          if (is_conservative) then
             gamma = gamma*(1.d0 - dt3*divu(i+ioff,j+joff))
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! correct sedgey with \Gamma^{x+}
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
          gamma = gamma * umac(i+1,j+joff)
          sedgey(i,j) = sedgey(i,j) - dt*gamma/(2.d0*hx)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute \Gamma^{x-}
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          if (umac(i,j+joff) .gt. 0) then
             isign = 1.d0
             ioff = -1
          else
             isign = -1.d0
             ioff = 0
          endif

          v = 0.d0
          if (vmac(i,j)*vmac(i+ioff,j) .gt. 0) then
             v = vmac(i+ioff,j)
          endif

          p1(1) = isign*0.5d0*hx
          p1(2) = jsign*0.5d0*hy

          p2(1) = isign*0.5d0*hx
          p2(2) = jsign*0.5d0*hy - vmac(i,j)*dt

          p3(1) = isign*0.5d0*hx - umac(i,j+joff)*dt
          p3(2) = jsign*0.5d0*hy - v*dt

          do ll=1,2
             del(ll) = (p2(ll)+p3(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val1)

          do ll=1,2
             del(ll) = (p1(ll)+p3(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val2)

          do ll=1,2
             del(ll) = (p1(ll)+p2(ll))/2.d0
          end do
          call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val3)

          ! average these centroid values to get the average value
          gamma = (val1+val2+val3)/3.d0

          ! source term
          if (is_conservative) then
             gamma = gamma*(1.d0 - dt3*divu(i+ioff,j+joff))
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! correct sedgey with \Gamma^{x-}
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          gamma = gamma * umac(i,j+joff)
          sedgey(i,j) = sedgey(i,j) + dt*gamma/(2.d0*hx)

       end do
    end do

    deallocate(ux,vy,divu)

  end subroutine bdsconc_2d

  subroutine bdsconc_3d(lo,hi,s,ng_s,slope,ng_c,umac,vmac,wmac,ng_u,force,ng_f, &
                        sedgex,sedgey,sedgez,ng_se,dx,dt,is_conservative)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_c,ng_u,ng_f,ng_se
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::  slope(lo(1)-ng_c :,lo(2)-ng_c :,lo(3)-ng_c :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:)
    real(kind=dp_t), intent(inout) :: sedgez(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_conservative

    ! local variables
    integer i,j,k,ioff,joff,koff,ll

    real(kind=dp_t), allocatable ::   ux(:,:,:)
    real(kind=dp_t), allocatable ::   vy(:,:,:)
    real(kind=dp_t), allocatable ::   wz(:,:,:)
    real(kind=dp_t), allocatable :: divu(:,:,:)

    real(kind=dp_t) :: isign,jsign,ksign,hx,hy,hz
    real(kind=dp_t) :: del(3),p1(3),p2(3),p3(3),p4(3)
    real(kind=dp_t) :: val1,val2,val3,val4,val5
    real(kind=dp_t) :: u,v,w,uu,vv,ww,gamma,gamma2
    real(kind=dp_t) :: dt2,dt3,dt4,half,sixth

    allocate(  ux(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(  vy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(  wz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(divu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    dt2 = dt/2.d0
    dt3 = dt/3.d0
    dt4 = dt/4.d0

    half = 0.5d0
    sixth = 1.d0/6.d0

    ! compute cell-centered ux, vy, and wz
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             ux(i,j,k) = (umac(i+1,j,k) - umac(i,j,k)) / hx
             vy(i,j,k) = (vmac(i,j+1,k) - vmac(i,j,k)) / hy
             wz(i,j,k) = (wmac(i,j,k+1) - wmac(i,j,k)) / hz
             divu(i,j,k) = ux(i,j,k) + vy(i,j,k) + wz(i,j,k)
          end do
       end do
    end do

    ! compute sedgex on x-faces
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute sedgex without transverse corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i,j,k) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             ! centroid of rectangular volume
             del(1) = isign*0.5d0*hx - 0.5d0*umac(i,j,k)*dt
             del(2) = 0.d0
             del(3) = 0.d0
             call eval_3d(s(i+ioff,j,k),slope(i+ioff,j,k,:),del,sedgex(i,j,k))

             ! source term
             if (is_conservative) then
                sedgex(i,j,k) = sedgex(i,j,k)* &
                     (1.d0 - dt2*ux(i+ioff,j,k)) + dt2*force(i+ioff,j,k)
             else
                sedgex(i,j,k) = sedgex(i,j,k)* &
                     (1.d0 + dt2*(vy(i+ioff,j,k)+wz(i+ioff,j,k))) + dt2*force(i+ioff,j,k)
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i+ioff,j+1,k) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             u = 0.d0
             if (umac(i,j,k)*umac(i,j+joff,k) .gt. 0) then
                u = umac(i,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k)*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(ux(i+ioff,j+joff,k)+vy(i+ioff,j+joff,k)))
             else
                gamma = gamma*(1.d0 + dt3*wz(i+ioff,j+joff,k))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             uu = 0.d0
             if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vmac(i+ioff,j+1,k)*vmac(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             uu = 0.d0
             if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vmac(i+ioff,j+1,k)*vmac(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * wmac(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * vmac(i+ioff,j+1,k)
             sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i+ioff,j,k) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             u = 0.d0
             if (umac(i,j,k)*umac(i,j+joff,k) .gt. 0) then
                u = umac(i,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k)*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(ux(i+ioff,j+joff,k)+vy(i+ioff,j+joff,k)))
             else
                gamma = gamma*(1.d0 + dt3*wz(i+ioff,j+joff,k))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             uu = 0.d0
             if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vmac(i+ioff,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             uu = 0.d0
             if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vmac(i+ioff,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * wmac(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * vmac(i+ioff,j,k)
             sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i+ioff,j,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             u = 0.d0
             if (umac(i,j,k)*umac(i,j,k+koff) .gt. 0) then
                u = umac(i,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k+1)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(ux(i+ioff,j,k+koff)+wz(i+ioff,j,k+koff)))
             else
                gamma = gamma*(1.d0 + dt3*vy(i+ioff,j,k+koff))

             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             uu = 0.d0
             if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wmac(i+ioff,j,k+1)*wmac(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k+1)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             uu = 0.d0
             if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wmac(i+ioff,j,k+1)*wmac(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k+1)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * vmac(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wmac(i+ioff,j,k+1)
             sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i+ioff,j,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             u = 0.d0
             if (umac(i,j,k)*umac(i,j,k+koff) .gt. 0) then
                u = umac(i,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(ux(i+ioff,j,k+koff)+wz(i+ioff,j,k+koff)))
             else
                gamma = gamma*(1.d0 + dt3*vy(i+ioff,j,k+koff))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             uu = 0.d0
             if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wmac(i+ioff,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             uu = 0.d0
             if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wmac(i+ioff,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * vmac(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wmac(i+ioff,j,k)
             sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.d0*hz)

          enddo
       enddo
    enddo

    ! compute sedgey on y-faces    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute sedgey without transverse corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             ! centroid of rectangular volume
             if (vmac(i,j,k) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             del(1) = 0.d0
             del(2) = jsign*0.5d0*hy - 0.5d0*vmac(i,j,k)*dt
             del(3) = 0.d0
             call eval_3d(s(i,j+joff,k),slope(i,j+joff,k,:),del,sedgey(i,j,k))

             ! source term
             if (is_conservative) then
                sedgey(i,j,k) = sedgey(i,j,k)* &
                     (1.d0 - dt2*vy(i,j+joff,k)) + dt2*force(i,j+joff,k)
             else
                sedgey(i,j,k) = sedgey(i,j,k)* &
                     (1.d0 + dt2*(ux(i,j+joff,k)+wz(i,j+joff,k))) + dt2*force(i,j+joff,k)
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i+1,j+joff,k) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             v = 0.d0
             if (vmac(i,j,k)*vmac(i+ioff,j,k) .gt. 0) then
                v = vmac(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - umac(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(vy(i+ioff,j+joff,k)+ux(i+ioff,j+joff,k)))
             else
                gamma = gamma*(1.d0 + dt3*wz(i+ioff,j+joff,k))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             vv = 0.d0
             if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (umac(i+1,j+joff,k)*umac(i+1,j+joff,k+koff) .gt. 0) then
                uu = umac(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             vv = 0.d0
             if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (umac(i+1,j+joff,k)*umac(i+1,j+joff,k+koff) .gt. 0) then
                uu = umac(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * wmac(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
             gamma = gamma * umac(i+1,j+joff,k)
             sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i,j+joff,k) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             v = 0.d0
             if (vmac(i,j,k)*vmac(i+ioff,j,k) .gt. 0) then
                v = vmac(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - umac(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(vy(i+ioff,j+joff,k)+ux(i+ioff,j+joff,k)))
             else
                gamma = gamma*(1.d0 + dt3*wz(i+ioff,j+joff,k))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             vv = 0.d0
             if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (umac(i,j+joff,k)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             vv = 0.d0
             if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (umac(i,j+joff,k)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - umac(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * wmac(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * umac(i,j+joff,k)
             sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             v = 0.d0
             if (vmac(i,j,k)*vmac(i,j,k+koff) .gt. 0) then
                v = vmac(i,j,k+koff)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(vy(i,j+joff,k+koff)+wz(i,j+joff,k+koff)))
             else
                gamma = gamma*(1.d0 + dt3*ux(i,j+joff,k+koff))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             vv = 0.d0
             if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wmac(i,j+joff,k+1)*wmac(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k+1)*dt

             p4(1) = isign*0.5d0*hx - umac(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * umac(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             vv = 0.d0
             if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wmac(i,j+joff,k+1)*wmac(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k+1)*dt

             p4(1) = isign*0.5d0*hx - umac(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * umac(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wmac(i,j+joff,k+1)
             sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wmac(i,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             v = 0.d0
             if (vmac(i,j,k)*vmac(i,j,k+koff) .gt. 0) then
                v = vmac(i,j,k+koff)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(vy(i,j+joff,k+koff)+wz(i,j+joff,k+koff)))
             else
                gamma = gamma*(1.d0 + dt3*ux(i,j+joff,k+koff))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             vv = 0.d0
             if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wmac(i,j+joff,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k)*dt

             p4(1) = isign*0.5d0*hx - umac(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * umac(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             vv = 0.d0
             if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wmac(i,j+joff,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k)*dt

             p4(1) = isign*0.5d0*hx - umac(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * umac(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wmac(i,j+joff,k)
             sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.d0*hz)

          enddo
       enddo
    enddo

    ! compute sedgez on z-faces
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute sedgez without transverse corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             ! centroid of rectangular volume
             if (wmac(i,j,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             del(1) = 0.d0
             del(2) = 0.d0
             del(3) = ksign*0.5d0*hz - 0.5d0*wmac(i,j,k)*dt
             call eval_3d(s(i,j,k+koff),slope(i,j,k+koff,:),del,sedgez(i,j,k))

             ! source term
             if (is_conservative) then
                sedgez(i,j,k) = sedgez(i,j,k)* &
                     (1.d0 - dt2*wz(i,j,k+koff)) + dt2*force(i,j,k+koff)
             else
                sedgez(i,j,k) = sedgez(i,j,k)* &
                     (1.d0 + dt2*(ux(i,j,k+koff)+vy(i,j,k+koff))) + dt2*force(i,j,k+koff)
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i+1,j,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             w = 0.d0
             if (wmac(i,j,k)*wmac(i+ioff,j,k) .gt. 0) then
                w = wmac(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p3(1) = isign*0.5d0*hx - umac(i+1,j,k+koff)*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(wz(i+ioff,j,k+koff)+ux(i+ioff,j,k+koff)))
             else
                gamma = gamma*(1.d0 + dt3*vy(i+ioff,j,k+koff))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             ww = 0.d0
             if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (umac(i+1,j,k+koff)*umac(i+1,j+joff,k+koff) .gt. 0) then
                uu = umac(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - umac(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             ww = 0.d0
             if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (umac(i+1,j,k+koff)*umac(i+1,j+joff,k+koff) .gt. 0) then
                uu = umac(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - umac(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * vmac(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * umac(i+1,j,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) - dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i,j,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             w = 0.d0
             if (wmac(i,j,k)*wmac(i+ioff,j,k) .gt. 0) then
                w = wmac(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p3(1) = isign*0.5d0*hx - umac(i,j,k+koff)*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(wz(i+ioff,j,k+koff)+ux(i+ioff,j,k+koff)))
             else
                gamma = gamma*(1.d0 + dt3*vy(i+ioff,j,k+koff))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             ww = 0.d0
             if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (umac(i,j,k+koff)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - umac(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             ww = 0.d0
             if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (umac(i,j,k+koff)*umac(i,j+joff,k+koff) .gt. 0) then
                uu = umac(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - umac(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * vmac(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * umac(i,j,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) + dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             w = 0.d0
             if (wmac(i,j,k)*wmac(i,j+joff,k) .gt. 0) then
                w = wmac(i,j+joff,k)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - vmac(i,j+1,k+koff)*dt
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(wz(i,j+joff,k+koff)+vy(i,j+joff,k+koff)))
             else
                gamma = gamma*(1.d0 + dt3*ux(i,j+joff,k+koff))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             ww = 0.d0
             if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vmac(i,j+1,k+koff)*vmac(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - umac(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * umac(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             ww = 0.d0
             if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vmac(i,j+1,k+koff)*vmac(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - umac(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * umac(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
             gamma = gamma * vmac(i,j+1,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) - dt*gamma/(2.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vmac(i,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             w = 0.d0
             if (wmac(i,j,k)*wmac(i,j+joff,k) .gt. 0) then
                w = wmac(i,j+joff,k)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - vmac(i,j,k+koff)*dt
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             if (is_conservative) then
                gamma = gamma*(1.d0 - dt3*(wz(i,j+joff,k+koff)+vy(i,j+joff,k+koff)))
             else
                gamma = gamma*(1.d0 + dt3*ux(i,j+joff,k+koff))
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             ww = 0.d0
             if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vmac(i,j,k+koff)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - umac(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * umac(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (umac(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             ww = 0.d0
             if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
                ww = wmac(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vmac(i,j,k+koff)*vmac(i+ioff,j,k+koff) .gt. 0) then
                vv = vmac(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - umac(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             if (is_conservative) then
                gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
             end if

             gamma2 = gamma2 * umac(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * vmac(i,j,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) + dt*gamma/(2.d0*hy)

          enddo
       enddo
    enddo

    deallocate(ux,vy,wz,divu)

  end subroutine bdsconc_3d

  subroutine bds_velpred(u,umac,force,dx,dt,the_bc_level,mla)

    use ml_cc_restriction_module, only: ml_edge_restriction
    use create_umac_grown_module

    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: force(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(in   ) :: mla

    ! this will hold slx, sly, slxy, etc.
    type(multifab) :: slopeu(mla%nlevel)
    type(multifab) :: slopev(mla%nlevel)
    type(multifab) :: slopew(mla%nlevel)

    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: sup(:,:,:,:)
    real(kind=dp_t), pointer :: svp(:,:,:,:)
    real(kind=dp_t), pointer :: swp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)

    integer :: nlevs,dm,n,i,comp
    integer :: ng_u,ng_um,ng_f,ng_s
    integer :: lo(mla%dim),hi(mla%dim)

    nlevs = mla%nlevel
    dm = u(1)%dim

    if (dm .eq. 2) then
       ! 3 components and 2 ghost cells
       ! component 1 = slx
       ! component 2 = sly
       ! component 3 = slxy
       do n=1,nlevs
          call multifab_build(slopeu(n),mla%la(n),3,2)
          call multifab_build(slopev(n),mla%la(n),3,2)
       end do
    else if (dm .eq. 3) then
       ! 7 components and 2 ghost cells
       ! component 1 = slx
       ! component 2 = sly
       ! component 3 = slz
       ! component 4 = slxy
       ! component 5 = slxz
       ! component 6 = slyz
       ! component 7 = slxyz
       do n=1,nlevs
          call multifab_build(slopeu(n),mla%la(n),7,2)
          call multifab_build(slopev(n),mla%la(n),7,2)
          call multifab_build(slopew(n),mla%la(n),7,2)
       end do
    end if

    ng_u  = u(1)%ng
    ng_um = umac(1,1)%ng
    ng_f  = force(1)%ng
    ng_s  = slopeu(1)%ng

    do n=1,nlevs
       do i=1,nfabs(u(n))
          up  => dataptr(u(n), i)
          sup => dataptr(slopeu(n), i)
          svp => dataptr(slopev(n), i)
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)
          fp  => dataptr(force(n), i)
          lo =  lwb(get_box(u(n), i))
          hi =  upb(get_box(u(n), i))
          select case(dm)
          case (2)
             call bdsslope_2d(lo, hi, up(:,:,1,1), ng_u, sup(:,:,1,:), ng_s, dx(n,:))
             call bdsslope_2d(lo, hi, up(:,:,1,2), ng_u, svp(:,:,1,:), ng_s, dx(n,:))

             call bds_velpred_2d(lo, hi, dx(n,:), dt, up(:,:,1,:), ng_u, &
                                 sup(:,:,1,:), svp(:,:,1,:), ng_s, &
                                 ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                 fp(:,:,1,:), ng_f)

          case (3)
             swp => dataptr(slopew(n), i)
             wmp => dataptr(umac(n,3), i)
             call bdsslope_3d(lo, hi, up(:,:,:,1), ng_u, sup(:,:,:,:), ng_s, dx(n,:))
             call bdsslope_3d(lo, hi, up(:,:,:,2), ng_u, svp(:,:,:,:), ng_s, dx(n,:))
             call bdsslope_3d(lo, hi, up(:,:,:,3), ng_u, swp(:,:,:,:), ng_s, dx(n,:))

             call bl_error("bds_velpred_3d not written yet")

          end select
       end do
    end do

    if (nlevs .gt. 1) then
       do n=2,nlevs
          call create_umac_grown(umac(n,:),umac(n-1,:),the_bc_level(n-1),the_bc_level(n), &
               n.eq.nlevs)
       end do
    else
       do n=1,nlevs
          do comp=1,dm
             call multifab_fill_boundary(umac(n,comp))
          enddo
       end do
    end if

    do n = nlevs,2,-1
       do comp = 1,dm
          call ml_edge_restriction(umac(n-1,comp),umac(n,comp),mla%mba%rr(n-1,:),comp)
       end do
    end do

  end subroutine bds_velpred

  subroutine bds_velpred_2d(lo,hi,dx,dt,u,ng_u,uslope,vslope,ng_s,umac,vmac,ng_um,force,ng_f)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s,ng_um,ng_f
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,:)
    real(kind=dp_t), intent(in   ) :: uslope(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real(kind=dp_t), intent(in   ) :: vslope(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,:)

    ! local
    integer :: i,j,ioff,joff,ll

    real(kind=dp_t) :: del(2),p1(2),p2(2),p3(2)
    real(kind=dp_t) :: uavg,uu,vv,eps
    real(kind=dp_t) :: isign,jsign,hx,hy,dt2
    real(kind=dp_t) :: val1,val2,val3,gamma

    logical :: test

    real(kind=dp_t), allocatable ::  ulx(:,:,:)
    real(kind=dp_t), allocatable ::  urx(:,:,:)
    real(kind=dp_t), allocatable :: u1dx(:,:,:)
    real(kind=dp_t), allocatable ::  uly(:,:,:)
    real(kind=dp_t), allocatable ::  ury(:,:,:)
    real(kind=dp_t), allocatable :: u1dy(:,:,:)

    allocate( ulx(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1,2))
    allocate( urx(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1,2))
    allocate(u1dx(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1,2))
    allocate( uly(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2,2))
    allocate( ury(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2,2))
    allocate(u1dy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2,2))

    eps = 1.d-8

    hx = dx(1)
    hy = dx(2)

    dt2 = 0.5d0*dt

    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+2

          ! compute ul, vl at x-faces
          uu = max(u(i-1,j,1),0.d0)
          del(1) = 0.5d0*hx - 0.5d0*uu*dt
          del(2) = 0.d0

          call eval_2d(u(i-1,j,1),uslope(i-1,j,:),del,ulx(i,j,1))
          call eval_2d(u(i-1,j,2),vslope(i-1,j,:),del,ulx(i,j,2))

          ! compute ur, vr at x-faces
          uu = min(u(i,j,1),0.d0)
          del(1) = -0.5d0*hx - 0.5*uu*dt
          del(2) = 0.d0

          call eval_2d(u(i,j,1),uslope(i,j,:),del,urx(i,j,1))
          call eval_2d(u(i,j,2),vslope(i,j,:),del,urx(i,j,2))

          ! solve Riemann problems to get u1d at x-faces
          uavg = 0.5d0*(ulx(i,j,1)+urx(i,j,1))
          test = ((ulx(i,j,1) .le. 0.d0 .and. urx(i,j,1) .ge. 0.d0) .or. abs(uavg) .lt. eps)
          u1dx(i,j,1) = merge(ulx(i,j,1),urx(i,j,1),uavg .gt. 0.d0)
          u1dx(i,j,1) = merge(0.d0,u1dx(i,j,1),test)

          ! upwind to get v1d at x-faces
          u1dx(i,j,2) = merge(ulx(i,j,2),urx(i,j,2),u1dx(i,j,1) .gt. 0.d0)
          uavg = 0.5d0*(ulx(i,j,2)+urx(i,j,2))
          u1dx(i,j,2) = merge(uavg,u1dx(i,j,2),abs(u1dx(i,j,1)) .lt. eps)

       end do
    end do

    do j=lo(2)-1,hi(2)+2
       do i=lo(1)-1,hi(1)+1

          ! compute ul, vl at y-faces
          vv = max(u(i,j-1,2),0.d0)
          del(1) = 0.d0
          del(2) = 0.5d0*hy - 0.5d0*vv*dt

          call eval_2d(u(i,j-1,1),uslope(i,j-1,:),del,uly(i,j,1))
          call eval_2d(u(i,j-1,2),vslope(i,j-1,:),del,uly(i,j,2))

          ! compute ur, vr at y-faces
          vv = min(u(i,j,2),0.d0)
          del(1) = 0.d0
          del(2) = -0.5d0*hy - 0.5*vv*dt

          call eval_2d(u(i,j,1),uslope(i,j,:),del,ury(i,j,1))
          call eval_2d(u(i,j,2),vslope(i,j,:),del,ury(i,j,2))

          ! solve Riemann problems to get v1d at y-faces
          uavg = 0.5d0*(uly(i,j,2)+ury(i,j,2))
          test = ((uly(i,j,2) .le. 0.d0 .and. ury(i,j,2) .ge. 0.d0) .or. abs(uavg) .lt. eps)
          u1dy(i,j,2) = merge(uly(i,j,2),ury(i,j,2),uavg .gt. 0.d0)
          u1dy(i,j,2) = merge(0.d0,u1dy(i,j,2),test)

          ! upwind to get u1d at y-faces
          u1dy(i,j,1) = merge(uly(i,j,1),ury(i,j,1),u1dy(i,j,2) .gt. 0.d0)
          uavg = 0.5d0*(uly(i,j,1)+ury(i,j,1))
          u1dy(i,j,1) = merge(uavg,u1dy(i,j,1),abs(u1dy(i,j,2)) .lt. eps)

       end do
    end do

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1

          ! update ul at x-faces with source terms
          ulx(i,j,1) = ulx(i,j,1) - dt2*ulx(i,j,1)*(u1dy(i-1,j+1,2)-u1dy(i-1,j,2))/hy
          ulx(i,j,1) = ulx(i,j,1) + dt2*force(i-1,j,1)

          ! update ul at x-faces with transverse fluxes
          if (u(i-1,j,1) .gt. 0.d0) then

             ! flux on hi y-face
             if (u1dy(i-1,j+1,2) .ge. 0.d0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             end if

             uu = max(0.d0,u(i-1,j+joff,1))

             p1(1) = 0.5d0*hx
             p1(2) = jsign*0.5d0*hy

             p2(1) = 0.5d0*hx - u(i-1,j,1)*dt
             p2(2) = jsign*0.5d0*hy

             p3(1) = 0.5d0*hx - uu*dt
             p3(2) = jsign*0.5d0*hy - u1dy(i-1,j+1,2)*dt

             do ll=1,2
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i-1,j+joff,2),vslope(i-1,j+joff,:),del,val1)

             do ll=1,2
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i-1,j+joff,2),vslope(i-1,j+joff,:),del,val2)

             do ll=1,2
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_2d(u(i-1,j+joff,2),vslope(i-1,j+joff,:),del,val3)

             gamma = (val1+val2+val3)/3.d0

             ulx(i,j,1) = ulx(i,j,1) - dt2*gamma*u1dy(i-1,j+1,2)

             ! flux on lo y-face
             if (u1dy(i-1,j,2) .gt. 0.d0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             end if

             uu = max(0.d0,u(i-1,j+joff,1))

             p1(1) = 0.5d0*hx
             p1(2) = jsign*0.5d0*hy

             p2(1) = 0.5d0*hx - u(i-1,j,1)*dt
             p2(2) = jsign*0.5d0*hy

             p3(1) = 0.5d0*hx - uu*dt
             p3(2) = jsign*0.5d0*hy - u1dy(i-1,j,2)*dt

             do ll=1,2
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i-1,j+joff,2),vslope(i-1,j+joff,:),del,val1)

             do ll=1,2
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i-1,j+joff,2),vslope(i-1,j+joff,:),del,val2)

             do ll=1,2
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_2d(u(i-1,j+joff,2),vslope(i-1,j+joff,:),del,val3)

             gamma = (val1+val2+val3)/3.d0

             ulx(i,j,1) = ulx(i,j,1) + dt2*gamma*u1dy(i-1,j,2)             
            
          end if

          ! update ur at x-faces with source terms
          urx(i,j,1) = urx(i,j,1) - dt2*urx(i,j,1)*(u1dy(i,j+1,2)-u1dy(i,j,2))/hy
          urx(i,j,1) = urx(i,j,1) + dt2*force(i,j,1)

          ! update ur at x-faces with transverse fluxes
          if (u(i,j,1) .lt. 0.d0) then

             ! flux on hi y-face
             if (u1dy(i,j+1,2) .ge. 0.d0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             end if

             uu = min(0.d0,u(i,j+joff,1))

             p1(1) = -0.5d0*hx
             p1(2) = jsign*0.5d0*hy

             p2(1) = -0.5d0*hx - u(i,j,1)*dt
             p2(2) = jsign*0.5d0*hy

             p3(1) = -0.5d0*hx - uu*dt
             p3(2) = jsign*0.5d0*hy - u1dy(i,j+1,2)*dt

             do ll=1,2
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i,j+joff,2),vslope(i,j+joff,:),del,val1)

             do ll=1,2
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i,j+joff,2),vslope(i,j+joff,:),del,val2)

             do ll=1,2
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_2d(u(i,j+joff,2),vslope(i,j+joff,:),del,val3)

             gamma = (val1+val2+val3)/3.d0

             ulx(i,j,1) = ulx(i,j,1) - dt2*gamma*u1dy(i,j+1,2)

             ! flux on lo y-face
             if (u1dy(i,j,2) .ge. 0.d0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             end if

             uu = min(0.d0,u(i,j+joff,1))

             p1(1) = -0.5d0*hx
             p1(2) = jsign*0.5d0*hy

             p2(1) = -0.5d0*hx - u(i,j,1)*dt
             p2(2) = jsign*0.5d0*hy

             p3(1) = -0.5d0*hx - uu*dt
             p3(2) = jsign*0.5d0*hy - u1dy(i,j,2)*dt

             do ll=1,2
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i,j+joff,2),vslope(i,j+joff,:),del,val1)

             do ll=1,2
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i,j+joff,2),vslope(i,j+joff,:),del,val2)

             do ll=1,2
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_2d(u(i,j+joff,2),vslope(i,j+joff,:),del,val3)

             gamma = (val1+val2+val3)/3.d0

             ulx(i,j,1) = ulx(i,j,1) - dt2*gamma*u1dy(i,j,2)

          end if

          ! solve Riemann problem to get umac
          uavg = 0.5d0*(ulx(i,j,1)+urx(i,j,1))
          test = ((ulx(i,j,1) .le. 0.d0 .and. urx(i,j,1) .ge. 0.d0) .or. abs(uavg) .lt. eps)
          umac(i,j) = merge(ulx(i,j,1),urx(i,j,1),uavg .gt. 0.d0)
          umac(i,j) = merge(0.d0,umac(i,j),test)

       end do
    end do


    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)

          ! update vl at y-faces with source terms
          uly(i,j,2) = uly(i,j,2) - dt2*uly(i,j,2)*(u1dx(i+1,j-1,1)-u1dx(i,j-1,1))/hx
          uly(i,j,2) = uly(i,j,2) + dt2*force(i,j-1,2)

          ! update vl at y-faces with transverse fluxes
          if (u(i,j-1,2) .gt. 0.d0) then

             ! flux on hi x-face
             if (u1dx(i+1,j-1,1) .ge. 0.d0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             end if

             vv = max(0.d0,u(i+ioff,j-1,2))

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.5d0*hy

             p2(1) = isign*0.5d0*hx
             p2(2) = 0.5d0*hy - u(i,j-1,2)*dt

             p3(1) = isign*0.5d0*hx - u1dx(i+1,j-1,1)*dt
             p3(2) = 0.5d0*hy - vv*dt

             do ll=1,2
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j-1,1),uslope(i+ioff,j-1,:),del,val1)

             do ll=1,2
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j-1,1),uslope(i+ioff,j-1,:),del,val2)

             do ll=1,2
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j-1,1),uslope(i+ioff,j-1,:),del,val3)

             gamma = (val1+val2+val3)/3.d0

             uly(i,j,2) = uly(i,j,2) - dt2*gamma*u1dx(i+1,j-1,1)

             ! flux on lo x-face
             if (u1dx(i,j-1,1) .ge. 0.d0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             end if

             vv = max(0.d0,u(i+ioff,j-1,2))

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.5d0*hy

             p2(1) = isign*0.5d0*hx
             p2(2) = 0.5d0*hy - u(i,j-1,2)*dt

             p3(1) = isign*0.5d0*hx - u1dx(i,j-1,1)*dt
             p3(2) = 0.5d0*hy - vv*dt

             do ll=1,2
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j-1,1),uslope(i+ioff,j-1,:),del,val1)

             do ll=1,2
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j-1,1),uslope(i+ioff,j-1,:),del,val2)

             do ll=1,2
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j-1,1),uslope(i+ioff,j-1,:),del,val3)

             gamma = (val1+val2+val3)/3.d0

             uly(i,j,2) = uly(i,j,2) + dt2*gamma*u1dx(i,j-1,1)

          end if

          ! update vr at y-faces with source terms
          ury(i,j,2) = ury(i,j,2) - dt2*ury(i,j,2)*(u1dx(i+1,j,1)-u1dx(i,j,1))/hx
          ury(i,j,2) = ury(i,j,2) + dt2*force(i,j,2)

          ! update vr at y-faces with transverse fluxes
          if (u(i,j,2) .lt. 0.d0) then

             ! flux on hi x-face
             if (u1dx(i+1,j,1) .ge. 0.d0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             end if

             vv = max(0.d0,u(i+ioff,j,2))

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.5d0*hy

             p2(1) = isign*0.5d0*hx
             p2(2) = 0.5d0*hy - u(i,j,2)*dt

             p3(1) = isign*0.5d0*hx - u1dx(i+1,j,1)*dt
             p3(2) = 0.5d0*hy - vv*dt

             do ll=1,2
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j,1),uslope(i+ioff,j,:),del,val1)

             do ll=1,2
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j,1),uslope(i+ioff,j,:),del,val2)

             do ll=1,2
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j,1),uslope(i+ioff,j,:),del,val3)

             gamma = (val1+val2+val3)/3.d0

             uly(i,j,2) = uly(i,j,2) - dt2*gamma*u1dx(i+1,j,1)

             ! flux on lo x-face
             if (u1dx(i,j,1) .ge. 0.d0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             end if

             vv = max(0.d0,u(i+ioff,j,2))

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.5d0*hy

             p2(1) = isign*0.5d0*hx
             p2(2) = 0.5d0*hy - u(i,j,2)*dt

             p3(1) = isign*0.5d0*hx - u1dx(i,j,1)*dt
             p3(2) = 0.5d0*hy - vv*dt

             do ll=1,2
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j,1),uslope(i+ioff,j,:),del,val1)

             do ll=1,2
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j,1),uslope(i+ioff,j,:),del,val2)

             do ll=1,2
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval_2d(u(i+ioff,j,1),uslope(i+ioff,j,:),del,val3)

             gamma = (val1+val2+val3)/3.d0

             uly(i,j,2) = uly(i,j,2) + dt2*gamma*u1dx(i,j,1)
             
          end if

          ! solve Riemann problem to get vmac
          uavg = 0.5d0*(uly(i,j,2)+ury(i,j,2))
          test = ((uly(i,j,2) .le. 0.d0 .and. ury(i,j,2) .ge. 0.d0) .or. abs(uavg) .lt. eps)
          vmac(i,j) = merge(uly(i,j,2),ury(i,j,2),uavg .gt. 0.d0)
          vmac(i,j) = merge(0.d0,vmac(i,j),test)

       end do
    end do

    deallocate(ulx,urx,u1dx,uly,ury,u1dy)

  end subroutine bds_velpred_2d

  subroutine eval_2d(s,slope,del,val)

    real(kind=dp_t), intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: slope(:)
    real(kind=dp_t), intent(in   ) :: del(:)
    real(kind=dp_t), intent(  out) :: val

    val = s + del(1)*slope(1) + del(2)*slope(2) + del(1)*del(2)*slope(3)

  end subroutine eval_2d

  subroutine eval_3d(s,slope,del,val)

    real(kind=dp_t), intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: slope(:)
    real(kind=dp_t), intent(in   ) :: del(:)
    real(kind=dp_t), intent(  out) :: val

    val = s + del(1)*slope(1) + del(2)*slope(2) + del(3)*slope(3) &
         + del(1)*del(2)*slope(4) + del(1)*del(3)*slope(5) + del(2)*del(3)*slope(6) &
         + del(1)*del(2)*del(3)*slope(7)

  end subroutine eval_3d

end module bds_module
