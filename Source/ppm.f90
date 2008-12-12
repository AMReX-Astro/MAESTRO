module ppm_module

  use bl_types

  implicit none

  private

  public :: ppm_2d, ppm_fpu_2d, ppm_3d, ppm_fpu_3d

contains

  ! characteristics based on u
  subroutine ppm_2d(n,s,ng_s,u,ng_u,Ip,Im,w0,lo,hi,bc,dx,dt)

    use bc_module
    use bl_constants_module
    use geometry, only: nr

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_s,ng_u
    real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(in   ) ::    u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real(kind=dp_t), intent(inout) ::   Ip(lo(1)-1   :,lo(2)-1   :,:) 
    real(kind=dp_t), intent(inout) ::   Im(lo(1)-1   :,lo(2)-1   :,:) 
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: bc(:,:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i,j,edge_interp_type,spm_limiter_type

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, C
    real(kind=dp_t) :: sigma, s6, w0cen

    ! s_{\ib,+}, s_{\ib,-}
    real(kind=dp_t), allocatable :: sp(:,:,:)
    real(kind=dp_t), allocatable :: sm(:,:,:)

    ! \delta s_{\ib}^{vL}
    real(kind=dp_t), allocatable :: dsvl_x(:,:)
    real(kind=dp_t), allocatable :: dsvl_y(:,:)

    ! s_{i+\half}^{H.O.}
    real(kind=dp_t), allocatable :: sedgex(:,:)
    real(kind=dp_t), allocatable :: sedgey(:,:)

    ! cell-centered indexing
    allocate(sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

    ! cell-centered indexing
    allocate(dsvl_x(lo(1)-2:hi(1)+2,lo(2)-1:hi(2)+1))
    allocate(dsvl_y(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+2))

    ! edge-centered indexing
    allocate(sedgex(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1))
    allocate(sedgey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2))

    ! constant used in Colella 2008
    C = 1.25d0

    ! 1 = 4th order with van Leer limiting
    ! 2 = 4th order with Colella 2008 limiting
    edge_interp_type = 1

    ! 1 = "old" limiters described in Colella 2008
    ! 2 = "new" limiters described in Colella 2008
    spm_limiter_type = 1

    ! compute s at x-edges
    if (edge_interp_type .eq. 1) then

       ! compute van Leer slopes in x-direction
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+2
             dsc = HALF * (s(i+1,j) - s(i-1,j))
             dsl = TWO  * (s(i  ,j) - s(i-1,j))
             dsr = TWO  * (s(i+1,j) - s(i  ,j))
             dsvl_x(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          end do
       end do
       
       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+2
             sedgex(i,j) = HALF*(s(i,j)+s(i-1,j)) - SIXTH*(dsvl_x(i,j)-dsvl_x(i-1,j))
             ! make sure sedgex lies in between adjacent cell-centered values
             sedgex(i,j) = max(sedgex(i,j),min(s(i,j),s(i-1,j)))
             sedgex(i,j) = min(sedgex(i,j),max(s(i,j),s(i-1,j)))
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl_x
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+2
             dsvl_x(i,j) = HALF * (s(i+1,j) - s(i-1,j))
          end do
       end do
       
       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+2
             sedgex(i,j) = HALF*(s(i,j)+s(i-1,j)) - SIXTH*(dsvl_x(i,j)-dsvl_x(i-1,j))
             ! if sedgex is not in between the neighboring s values, we limit
             if (sedgex(i,j) .lt. min(s(i,j),s(i-1,j)) .or. &
                 sedgex(i,j) .gt. max(s(i,j),s(i-1,j))) then
                D2  = (THREE/dx(1)**2)*(s(i-1,j)-TWO*sedgex(i,j)+s(i,j))
                D2L = (ONE/dx(1)**2)*(s(i-2,j)-TWO*s(i-1,j)+s(i,j))
                D2R = (ONE/dx(1)**2)*(s(i-1,j)-TWO*s(i,j)+s(i+1,j))
                if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R)) then
                   sedgex(i,j) = HALF*(s(i-1,j)+s(i,j)) - (dx(1)**2/THREE) &
                        *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                else
                   sedgex(i,j) = HALF*(s(i-1,j)+s(i,j))
                end if
             end if
          end do
       end do

    end if

    ! fill x-component of sp and sm
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sp(i,j,1) = sedgex(i+1,j)
          sm(i,j,1) = sedgex(i  ,j)
       end do
    end do

    ! different stencil needed for EXT_DIR and HOEXTRAP bc's
    !
    !
    !

    ! limit sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify sp and sm using quadratic limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j,1)-s(i,j))*(s(i,j)-sm(i,j,1)) .le. ZERO) then
                sp(i,j,1) = s(i,j)
                sm(i,j,1) = s(i,j)
             else if (abs(sp(i,j,1)-s(i,j)) .ge. TWO*abs(sm(i,j,1)-s(i,j))) then
                sp(i,j,1) = THREE*s(i,j) - TWO*sm(i,j,1)
             else if (abs(sm(i,j,1)-s(i,j)) .ge. TWO*abs(sp(i,j,1)-s(i,j))) then
                sm(i,j,1) = THREE*s(i,j) - TWO*sp(i,j,1)
             end if
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

    end if

    ! compute Ip and Im
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sigma = abs(u(i,j,1))*dt/dx(1)
          s6 = SIX*s(i,j) - THREE*(sm(i,j,1)+sp(i,j,1))
          if (u(i,j,1) .gt. ZERO) then
             Ip(i,j,1) = sp(i,j,1) - (sigma/TWO)*(sp(i,j,1)-sm(i,j,1)-(ONE-TWO3RD*sigma)*s6)
             Im(i,j,1) = s(i,j)
          else if (u(i,j,1) .lt. ZERO) then
             Ip(i,j,1) = s(i,j)
             Im(i,j,1) = sm(i,j,1) + (sigma/TWO)*(sp(i,j,1)-sm(i,j,1)+(ONE-TWO3RD*sigma)*s6)
          else
             Ip(i,j,1) = s(i,j)
             Im(i,j,1) = s(i,j)
          end if
       end do
    end do

    ! compute s at y-edges
    if (edge_interp_type .eq. 1) then

       ! compute van Leer slopes in y-direction
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             dsc = HALF * (s(i,j+1) - s(i,j-1))
             dsl = TWO  * (s(i,j  ) - s(i,j-1))
             dsr = TWO  * (s(i,j+1) - s(i,j  ))
             dsvl_y(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          end do
       end do
       
       ! interpolate s to y-edges
       do j=lo(2)-1,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             sedgey(i,j) = HALF*(s(i,j)+s(i,j-1)) - SIXTH*(dsvl_y(i,j)-dsvl_y(i,j-1))
             ! make sure sedgey lies in between adjacent cell-centered values
             sedgey(i,j) = max(sedgey(i,j),min(s(i,j),s(i,j-1)))
             sedgey(i,j) = min(sedgey(i,j),max(s(i,j),s(i,j-1)))
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl_y
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             dsvl_y(i,j) = HALF * (s(i,j+1) - s(i,j-1))
          end do
       end do
       
       ! interpolate s to y-edges
       do j=lo(2)-1,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             sedgey(i,j) = HALF*(s(i,j)+s(i,j-1)) - SIXTH*(dsvl_y(i,j)-dsvl_y(i,j-1))
             ! if sedgey is not in between the neighboring s values, we limit
             if (sedgey(i,j) .lt. min(s(i,j),s(i,j-1)) .or. &
                 sedgey(i,j) .gt. max(s(i,j),s(i,j-1))) then
                D2  = (THREE/dx(2)**2)*(s(i,j-1)-TWO*sedgey(i,j)+s(i,j))
                D2L = (ONE/dx(2)**2)*(s(i,j-2)-TWO*s(i,j-1)+s(i,j))
                D2R = (ONE/dx(2)**2)*(s(i,j-1)-TWO*s(i,j)+s(i,j+1))
                if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R)) then
                   sedgey(i,j) = HALF*(s(i,j-1)+s(i,j)) - (dx(2)**2/THREE) &
                        *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                else
                   sedgey(i,j) = HALF*(s(i,j-1)+s(i,j))
                end if
             end if
          end do
       end do

    end if

    ! ! fill y-component of sp and sm
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sp(i,j,2) = sedgey(i,j+1)
          sm(i,j,2) = sedgey(i,j  )
       end do
    end do

    ! different stencil needed for EXT_DIR and HOEXTRAP bc's
    !
    !
    !

    if (spm_limiter_type .eq. 1) then

       ! modify sp and sm using quadratic limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j,2)-s(i,j))*(s(i,j)-sm(i,j,2)) .le. ZERO) then
                sp(i,j,2) = s(i,j)
                sm(i,j,2) = s(i,j)
             else if (abs(sp(i,j,2)-s(i,j)) .ge. TWO*abs(sm(i,j,2)-s(i,j))) then
                sp(i,j,2) = THREE*s(i,j) - TWO*sm(i,j,2)
             else if (abs(sm(i,j,2)-s(i,j)) .ge. TWO*abs(sp(i,j,2)-s(i,j))) then
                sm(i,j,2) = THREE*s(i,j) - TWO*sp(i,j,2)
             end if
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

    end if

    ! compute Ip and Im
    do j=lo(2)-1,hi(2)+1
       ! compute effect of w0
       if (j .le. 0) then
          w0cen = w0(0)
       else if (j .ge. nr(n)) then
          w0cen = w0(nr(n))
       else
          w0cen = HALF*(w0(j)+w0(j+1))
       end if
       do i=lo(1)-1,hi(1)+1
          sigma = abs(u(i,j,2)+w0cen)*dt/dx(2)
          s6 = SIX*s(i,j) - THREE*(sm(i,j,2)+sp(i,j,2))
          if (u(i,j,2)+w0cen .gt. ZERO) then
             Ip(i,j,2) = sp(i,j,2) - (sigma/TWO)*(sp(i,j,2)-sm(i,j,2)-(ONE-TWO3RD*sigma)*s6)
             Im(i,j,2) = s(i,j)
          else if (u(i,j,2)+w0cen .lt. ZERO) then
             Ip(i,j,2) = s(i,j)
             Im(i,j,2) = sm(i,j,2) + (sigma/TWO)*(sp(i,j,2)-sm(i,j,2)+(ONE-TWO3RD*sigma)*s6)
          else
             Ip(i,j,2) = s(i,j)
             Im(i,j,2) = s(i,j)
          end if
       end do
    end do

    deallocate(sp,sm,dsvl_x,dsvl_y,sedgex,sedgey)

  end subroutine ppm_2d

  ! characteristics based on umac
  subroutine ppm_fpu_2d(n,s,ng_s,umac,vmac,ng_u,Ip,Im,w0,lo,hi,bc,dx,dt)

    use bc_module
    use bl_constants_module
    use geometry, only: nr

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_s,ng_u
    real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) ::   Ip(lo(1)-1   :,lo(2)-1   :,:) 
    real(kind=dp_t), intent(inout) ::   Im(lo(1)-1   :,lo(2)-1   :,:) 
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: bc(:,:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i,j,edge_interp_type,spm_limiter_type

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, C
    real(kind=dp_t) :: sigmam, sigmap, s6, w0lo, w0hi

    ! s_{\ib,+}, s_{\ib,-}
    real(kind=dp_t), allocatable :: sp(:,:,:)
    real(kind=dp_t), allocatable :: sm(:,:,:)

    ! \delta s_{\ib}^{vL}
    real(kind=dp_t), allocatable :: dsvl_x(:,:)
    real(kind=dp_t), allocatable :: dsvl_y(:,:)

    ! s_{i+\half}^{H.O.}
    real(kind=dp_t), allocatable :: sedgex(:,:)
    real(kind=dp_t), allocatable :: sedgey(:,:)

    ! cell-centered indexing
    allocate(sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

    ! cell-centered indexing
    allocate(dsvl_x(lo(1)-2:hi(1)+2,lo(2)-1:hi(2)+1))
    allocate(dsvl_y(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+2))

    ! edge-centered indexing
    allocate(sedgex(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1))
    allocate(sedgey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2))

    ! constant used in Colella 2008
    C = 1.25d0

    ! 1 = 4th order with van Leer limiting
    ! 2 = 4th order with Colella 2008 limiting
    edge_interp_type = 1

    ! 1 = "old" limiters described in Colella 2008
    ! 2 = "new" limiters described in Colella 2008
    spm_limiter_type = 1

    ! compute s at x-edges
    if (edge_interp_type .eq. 1) then

       ! compute van Leer slopes in x-direction
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+2
             dsc = HALF * (s(i+1,j) - s(i-1,j))
             dsl = TWO  * (s(i  ,j) - s(i-1,j))
             dsr = TWO  * (s(i+1,j) - s(i  ,j))
             dsvl_x(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          end do
       end do
       
       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+2
             sedgex(i,j) = HALF*(s(i,j)+s(i-1,j)) - SIXTH*(dsvl_x(i,j)-dsvl_x(i-1,j))
             ! make sure sedgex lies in between adjacent cell-centered values
             sedgex(i,j) = max(sedgex(i,j),min(s(i,j),s(i-1,j)))
             sedgex(i,j) = min(sedgex(i,j),max(s(i,j),s(i-1,j)))
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl_x
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+2
             dsvl_x(i,j) = HALF * (s(i+1,j) - s(i-1,j))
          end do
       end do
       
       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+2
             sedgex(i,j) = HALF*(s(i,j)+s(i-1,j)) - SIXTH*(dsvl_x(i,j)-dsvl_x(i-1,j))
             ! if sedgex is not in between the neighboring s values, we limit
             if (sedgex(i,j) .lt. min(s(i,j),s(i-1,j)) .or. &
                 sedgex(i,j) .gt. max(s(i,j),s(i-1,j))) then
                D2  = (THREE/dx(1)**2)*(s(i-1,j)-TWO*sedgex(i,j)+s(i,j))
                D2L = (ONE/dx(1)**2)*(s(i-2,j)-TWO*s(i-1,j)+s(i,j))
                D2R = (ONE/dx(1)**2)*(s(i-1,j)-TWO*s(i,j)+s(i+1,j))
                if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R)) then
                   sedgex(i,j) = HALF*(s(i-1,j)+s(i,j)) - (dx(1)**2/THREE) &
                        *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                else
                   sedgex(i,j) = HALF*(s(i-1,j)+s(i,j))
                end if
             end if
          end do
       end do

    end if

    ! fill x-component of sp and sm
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sp(i,j,1) = sedgex(i+1,j)
          sm(i,j,1) = sedgex(i  ,j)
       end do
    end do

    ! different stencil needed for EXT_DIR and HOEXTRAP bc's
    !
    !
    !

    ! limit sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify sp and sm using quadratic limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j,1)-s(i,j))*(s(i,j)-sm(i,j,1)) .le. ZERO) then
                sp(i,j,1) = s(i,j)
                sm(i,j,1) = s(i,j)
             else if (abs(sp(i,j,1)-s(i,j)) .ge. TWO*abs(sm(i,j,1)-s(i,j))) then
                sp(i,j,1) = THREE*s(i,j) - TWO*sm(i,j,1)
             else if (abs(sm(i,j,1)-s(i,j)) .ge. TWO*abs(sp(i,j,1)-s(i,j))) then
                sm(i,j,1) = THREE*s(i,j) - TWO*sp(i,j,1)
             end if
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

    end if

    ! compute Ip and Im
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sigmap = abs(umac(i+1,j))*dt/dx(1)
          sigmam = abs(umac(i,j))*dt/dx(1)
          s6 = SIX*s(i,j) - THREE*(sm(i,j,1)+sp(i,j,1))
          if (umac(i+1,j) .gt. ZERO) then
             Ip(i,j,1) = sp(i,j,1) - (sigmap/TWO)*(sp(i,j,1)-sm(i,j,1)-(ONE-TWO3RD*sigmap)*s6)
          else
             Ip(i,j,1) = s(i,j)
          end if
          if (umac(i,j) .lt. ZERO) then
             Im(i,j,1) = sm(i,j,1) + (sigmam/TWO)*(sp(i,j,1)-sm(i,j,1)+(ONE-TWO3RD*sigmam)*s6)
          else
             Im(i,j,1) = s(i,j)
          end if
       end do
    end do

    ! compute s at y-edges
    if (edge_interp_type .eq. 1) then

       ! compute van Leer slopes in y-direction
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             dsc = HALF * (s(i,j+1) - s(i,j-1))
             dsl = TWO  * (s(i,j  ) - s(i,j-1))
             dsr = TWO  * (s(i,j+1) - s(i,j  ))
             dsvl_y(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          end do
       end do
       
       ! interpolate s to y-edges
       do j=lo(2)-1,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             sedgey(i,j) = HALF*(s(i,j)+s(i,j-1)) - SIXTH*(dsvl_y(i,j)-dsvl_y(i,j-1))
             ! make sure sedgey lies in between adjacent cell-centered values
             sedgey(i,j) = max(sedgey(i,j),min(s(i,j),s(i,j-1)))
             sedgey(i,j) = min(sedgey(i,j),max(s(i,j),s(i,j-1)))
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl_y
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             dsvl_y(i,j) = HALF * (s(i,j+1) - s(i,j-1))
          end do
       end do
       
       ! interpolate s to y-edges
       do j=lo(2)-1,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             sedgey(i,j) = HALF*(s(i,j)+s(i,j-1)) - SIXTH*(dsvl_y(i,j)-dsvl_y(i,j-1))
             ! if sedgey is not in between the neighboring s values, we limit
             if (sedgey(i,j) .lt. min(s(i,j),s(i,j-1)) .or. &
                 sedgey(i,j) .gt. max(s(i,j),s(i,j-1))) then
                D2  = (THREE/dx(2)**2)*(s(i,j-1)-TWO*sedgey(i,j)+s(i,j))
                D2L = (ONE/dx(2)**2)*(s(i,j-2)-TWO*s(i,j-1)+s(i,j))
                D2R = (ONE/dx(2)**2)*(s(i,j-1)-TWO*s(i,j)+s(i,j+1))
                if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R)) then
                   sedgey(i,j) = HALF*(s(i,j-1)+s(i,j)) - (dx(2)**2/THREE) &
                        *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                else
                   sedgey(i,j) = HALF*(s(i,j-1)+s(i,j))
                end if
             end if
          end do
       end do

    end if

    ! ! fill y-component of sp and sm
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sp(i,j,2) = sedgey(i,j+1)
          sm(i,j,2) = sedgey(i,j  )
       end do
    end do

    ! different stencil needed for EXT_DIR and HOEXTRAP bc's
    !
    !
    !

    if (spm_limiter_type .eq. 1) then

       ! modify sp and sm using quadratic limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j,2)-s(i,j))*(s(i,j)-sm(i,j,2)) .le. ZERO) then
                sp(i,j,2) = s(i,j)
                sm(i,j,2) = s(i,j)
             else if (abs(sp(i,j,2)-s(i,j)) .ge. TWO*abs(sm(i,j,2)-s(i,j))) then
                sp(i,j,2) = THREE*s(i,j) - TWO*sm(i,j,2)
             else if (abs(sm(i,j,2)-s(i,j)) .ge. TWO*abs(sp(i,j,2)-s(i,j))) then
                sm(i,j,2) = THREE*s(i,j) - TWO*sp(i,j,2)
             end if
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

    end if

    ! compute Ip and Im
    do j=lo(2)-1,hi(2)+1
       ! compute effect of w0
       if (j .lt. 0) then
          w0lo = w0(0)
          w0hi = w0(0)
       else if (j .gt. nr(n)-1) then
          w0lo = w0(nr(n))
          w0hi = w0(nr(n))
       else
          w0lo = w0(j)
          w0hi = w0(j+1)
       end if
       do i=lo(1)-1,hi(1)+1
          sigmap = abs(vmac(i,j+1)+w0hi)*dt/dx(2)
          sigmam = abs(vmac(i,j  )+w0lo)*dt/dx(2)
          s6 = SIX*s(i,j) - THREE*(sm(i,j,2)+sp(i,j,2))
          if (vmac(i,j+1)+w0hi .gt. ZERO) then
             Ip(i,j,2) = sp(i,j,2) - (sigmap/TWO)*(sp(i,j,2)-sm(i,j,2)-(ONE-TWO3RD*sigmap)*s6)
          else
             Ip(i,j,2) = s(i,j)
          end if
          if (vmac(i,j)+w0lo .lt. ZERO) then
             Im(i,j,2) = sm(i,j,2) + (sigmam/TWO)*(sp(i,j,2)-sm(i,j,2)+(ONE-TWO3RD*sigmam)*s6)
          else
             Im(i,j,2) = s(i,j)
          end if
       end do
    end do

    deallocate(sp,sm,dsvl_x,dsvl_y,sedgex,sedgey)

  end subroutine ppm_fpu_2d

  ! characteristics based on u
  subroutine ppm_3d(s,ng_s,u,ng_u,Ip,Im,lo,hi,bc,dx,dt)

    use bc_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_u
    real(kind=dp_t), intent(in   ) ::   s(lo(1)-ng_s:,lo(2)-ng_s:,hi(3)-ng_s:)
    real(kind=dp_t), intent(in   ) ::   u(lo(1)-ng_u:,lo(2)-ng_u:,hi(3)-ng_u:,:)
    real(kind=dp_t), intent(inout) ::  Ip(lo(1)-1   :,lo(2)-1   :,hi(3)-1   :,:) 
    real(kind=dp_t), intent(inout) ::  Im(lo(1)-1   :,lo(2)-1   :,hi(3)-1   :,:)  
    integer        , intent(in   ) :: bc(:,:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local



  end subroutine ppm_3d

  ! characteristics based on umac
  subroutine ppm_fpu_3d(s,ng_s,umac,vmac,wmac,ng_u,Ip,Im,lo,hi,bc,dx,dt)

    use bc_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_u
    real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng_s:,lo(2)-ng_s:,hi(3)-ng_s:)
    real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:,lo(2)-ng_u:,hi(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_u:,lo(2)-ng_u:,hi(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: wmac(lo(1)-ng_u:,lo(2)-ng_u:,hi(3)-ng_u:)
    real(kind=dp_t), intent(inout) ::   Ip(lo(1)-1   :,lo(2)-1   :,hi(3)-1   :,:) 
    real(kind=dp_t), intent(inout) ::   Im(lo(1)-1   :,lo(2)-1   :,hi(3)-1   :,:)  
    integer        , intent(in   ) :: bc(:,:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local



  end subroutine ppm_fpu_3d


end module ppm_module
