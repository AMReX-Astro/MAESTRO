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
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i,j,edge_interp_type,spm_limiter_type

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham, ds
    real(kind=dp_t) :: dI, sgn
    real(kind=dp_t) :: sigma, s6, w0cc

    ! s_{\ib,+}, s_{\ib,-}
    real(kind=dp_t), allocatable :: sp(:,:)
    real(kind=dp_t), allocatable :: sm(:,:)

    ! \delta s_{\ib}^{vL}
    real(kind=dp_t), allocatable :: dsvl(:,:)

    ! s_{i+\half}^{H.O.}
    real(kind=dp_t), allocatable :: sedge(:,:)

    ! cell-centered indexing
    allocate(sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))

    ! constant used in Colella 2008
    C = 1.25d0

    ! 1 = 4th order with van Leer limiting
    ! 2 = 4th order with Colella 2008 limiting
    edge_interp_type = 1

    ! 1 = "old" limiters described in Colella 2008
    ! 2 = "new" limiters described in Colella 2008
    spm_limiter_type = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2,lo(2)-1:hi(2)+1))

    ! edge-centered indexing for x-faces
    allocate(sedge(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1))

    ! compute s at x-edges
    if (edge_interp_type .eq. 1) then

       ! compute van Leer slopes in x-direction
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+2
             dsc = HALF * (s(i+1,j) - s(i-1,j))
             dsl = TWO  * (s(i  ,j) - s(i-1,j))
             dsr = TWO  * (s(i+1,j) - s(i  ,j))
             dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          end do
       end do
       
       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+2
             sedge(i,j) = HALF*(s(i,j)+s(i-1,j)) - SIXTH*(dsvl(i,j)-dsvl(i-1,j))
             ! make sure sedge lies in between adjacent cell-centered values
             sedge(i,j) = max(sedge(i,j),min(s(i,j),s(i-1,j)))
             sedge(i,j) = min(sedge(i,j),max(s(i,j),s(i-1,j)))
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+2
             dsvl(i,j) = HALF * (s(i+1,j) - s(i-1,j))
          end do
       end do
       
       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+2
             sedge(i,j) = HALF*(s(i,j)+s(i-1,j)) - SIXTH*(dsvl(i,j)-dsvl(i-1,j))
             ! if sedge is not in between the neighboring s values, we limit
             if (sedge(i,j) .lt. min(s(i,j),s(i-1,j)) .or. &
                 sedge(i,j) .gt. max(s(i,j),s(i-1,j))) then
                D2  = (THREE/dx(1)**2)*(s(i-1,j)-TWO*sedge(i,j)+s(i,j))
                D2L = (ONE/dx(1)**2)*(s(i-2,j)-TWO*s(i-1,j)+s(i,j))
                D2R = (ONE/dx(1)**2)*(s(i-1,j)-TWO*s(i,j)+s(i+1,j))
                if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R)) then
                   sedge(i,j) = HALF*(s(i-1,j)+s(i,j)) - (dx(1)**2/THREE) &
                        *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                else
                   sedge(i,j) = HALF*(s(i-1,j)+s(i,j))
                end if
             end if
          end do
       end do

    end if

    ! fill x-component of sp and sm
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sp(i,j) = sedge(i+1,j)
          sm(i,j) = sedge(i  ,j)
       end do
    end do

    ! limit x-component of sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify using quadratic limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

       ! modify using Colella 2008 limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO .or. &
                 (s(i+1,j)-s(i,j))*(s(i,j)-s(i-1,j)) .le. ZERO ) then
                s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
                D2  = -TWO*s6/dx(1)**2
                D2C = (ONE/dx(1)**2)*(s(i-1,j)-TWO*s(i,j)+s(i+1,j))
                D2L = (ONE/dx(1)**2)*(s(i-2,j)-TWO*s(i-1,j)+s(i,j))
                D2R = (ONE/dx(1)**2)*(s(i,j)-TWO*s(i+1,j)+s(i+2,j))
                if (sign(ONE,D2) .eq. sign(ONE,D2C) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R) .and. &
                    D2 .ne. ZERO) then
                   D2LIM = sign(ONE,D2)*min(C*abs(D2C),C*abs(D2L),C*abs(D2R),abs(D2))
                   sp(i,j) = s(i,j) + (sp(i,j)-s(i,j))*(D2LIM/D2)
                   sm(i,j) = s(i,j) + (sm(i,j)-s(i,j))*(D2LIM/D2)
                else
                   sp(i,j) = s(i,j)
                   sm(i,j) = s(i,j)
                end if
             else
                alphap = sp(i,j)-s(i,j)
                alpham = sm(i,j)-s(i,j)
                if (abs(alphap) .ge. TWO*abs(alpham)) then
                   dI = -alphap**2 / (FOUR*(alphap+alpham))
                   ds = s(i+1,j)-s(i,j)
                   sgn = sign(ONE,s(i+1,j)-s(i-1,j))
                   if (sgn*dI .ge. sgn*ds) then
                      sp(i,j) = s(i,j) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alpham))
                   end if
                else if (abs(alpham) .ge. TWO*abs(alphap)) then
                   dI = -alpham**2 / (FOUR*(alphap+alpham))
                   ds = s(i-1,j)-s(i,j)
                   sgn = sign(ONE,s(i+1,j)-s(i-1,j))
                   if (sgn*dI .ge. sgn*ds) then
                      sm(i,j) = s(i,j) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alphap))
                   end if
                end if
             end if
          end do
       end do

    end if

    ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
    if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sm(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)+1,lo(2)-1:hi(2)+1) = -FIFTH        *s(lo(1)-1,lo(2)-1:hi(2)+1) &
                                        + (THREE/FOUR)*s(lo(1)  ,lo(2)-1:hi(2)+1) &
                                        + HALF        *s(lo(1)+1,lo(2)-1:hi(2)+1) &
                                        - (ONE/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do j=lo(2)-1,hi(2)+1
          sedge(lo(1)+1,j) = max(sedge(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
          sedge(lo(1)+1,j) = min(sedge(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
       end do

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          sp(lo(1)  ,j) = sedge(lo(1)+1,j)
          sm(lo(1)+1,j) = sedge(lo(1)+1,j)
       end do
    end if

    if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sp(hi(1),lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(hi(1),lo(2)-1:hi(2)+1) = -FIFTH        *s(hi(1)+1,lo(2)-1:hi(2)+1) &
                                      + (THREE/FOUR)*s(hi(1)  ,lo(2)-1:hi(2)+1) &
                                      + HALF        *s(hi(1)-1,lo(2)-1:hi(2)+1) &
                                      - (ONE/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do j=lo(2)-1,hi(2)+1
          sedge(hi(1),j) = max(sedge(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
          sedge(hi(1),j) = min(sedge(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
       end do

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          sp(hi(1)-1,j) = sedge(hi(1),j)
          sm(hi(1)  ,j) = sedge(hi(1),j)
       end do
    end if

    ! compute x-component of Ip and Im
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sigma = abs(u(i,j,1))*dt/dx(1)
          s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
          if (u(i,j,1) .gt. ZERO) then
             Ip(i,j,1) = sp(i,j) - (sigma/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
             Im(i,j,1) = s(i,j)
          else if (u(i,j,1) .lt. ZERO) then
             Ip(i,j,1) = s(i,j)
             Im(i,j,1) = sm(i,j) + (sigma/TWO)*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          else
             Ip(i,j,1) = s(i,j)
             Im(i,j,1) = s(i,j)
          end if
       end do
    end do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+2))

    ! edge-centered indexing for y-faces
    allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2))

    ! compute s at y-edges
    if (edge_interp_type .eq. 1) then

       ! compute van Leer slopes in y-direction
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             dsc = HALF * (s(i,j+1) - s(i,j-1))
             dsl = TWO  * (s(i,j  ) - s(i,j-1))
             dsr = TWO  * (s(i,j+1) - s(i,j  ))
             dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          end do
       end do
       
       ! interpolate s to y-edges
       do j=lo(2)-1,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             sedge(i,j) = HALF*(s(i,j)+s(i,j-1)) - SIXTH*(dsvl(i,j)-dsvl(i,j-1))
             ! make sure sedge lies in between adjacent cell-centered values
             sedge(i,j) = max(sedge(i,j),min(s(i,j),s(i,j-1)))
             sedge(i,j) = min(sedge(i,j),max(s(i,j),s(i,j-1)))
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             dsvl(i,j) = HALF * (s(i,j+1) - s(i,j-1))
          end do
       end do
       
       ! interpolate s to y-edges
       do j=lo(2)-1,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             sedge(i,j) = HALF*(s(i,j)+s(i,j-1)) - SIXTH*(dsvl(i,j)-dsvl(i,j-1))
             ! if sedge is not in between the neighboring s values, we limit
             if (sedge(i,j) .lt. min(s(i,j),s(i,j-1)) .or. &
                 sedge(i,j) .gt. max(s(i,j),s(i,j-1))) then
                D2  = (THREE/dx(2)**2)*(s(i,j-1)-TWO*sedge(i,j)+s(i,j))
                D2L = (ONE/dx(2)**2)*(s(i,j-2)-TWO*s(i,j-1)+s(i,j))
                D2R = (ONE/dx(2)**2)*(s(i,j-1)-TWO*s(i,j)+s(i,j+1))
                if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R)) then
                   sedge(i,j) = HALF*(s(i,j-1)+s(i,j)) - (dx(2)**2/THREE) &
                        *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                else
                   sedge(i,j) = HALF*(s(i,j-1)+s(i,j))
                end if
             end if
          end do
       end do

    end if

    ! fill y-component of sp and sm
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sp(i,j) = sedge(i,j+1)
          sm(i,j) = sedge(i,j  )
       end do
    end do

    ! limit y-component of sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify using quadratic limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

       ! modify using Colella 2008 limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO .or. &
                 (s(i,j+1)-s(i,j))*(s(i,j)-s(i,j-1)) .le. ZERO ) then
                s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
                D2  = -TWO*s6/dx(2)**2
                D2C = (ONE/dx(2)**2)*(s(i,j-1)-TWO*s(i,j)+s(i,j+1))
                D2L = (ONE/dx(2)**2)*(s(i,j-2)-TWO*s(i,j-1)+s(i,j))
                D2R = (ONE/dx(2)**2)*(s(i,j)-TWO*s(i,j+1)+s(i,j+2))
                if (sign(ONE,D2) .eq. sign(ONE,D2C) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R) .and. &
                    D2 .ne. ZERO) then
                   D2LIM = sign(ONE,D2)*min(C*abs(D2C),C*abs(D2L),C*abs(D2R),abs(D2))
                   sp(i,j) = s(i,j) + (sp(i,j)-s(i,j))*(D2LIM/D2)
                   sm(i,j) = s(i,j) + (sm(i,j)-s(i,j))*(D2LIM/D2)
                else
                   sp(i,j) = s(i,j)
                   sm(i,j) = s(i,j)
                end if
             else
                alphap = sp(i,j)-s(i,j)
                alpham = sm(i,j)-s(i,j)
                if (abs(alphap) .ge. TWO*abs(alpham)) then
                   dI = -alphap**2 / (FOUR*(alphap+alpham))
                   ds = s(i,j+1)-s(i,j)
                   sgn = sign(ONE,s(i,j+1)-s(i,j-1))
                   if (sgn*dI .ge. sgn*ds) then
                      sp(i,j) = s(i,j) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alpham))
                   end if
                else if (abs(alpham) .ge. TWO*abs(alphap)) then
                   dI = -alpham**2 / (FOUR*(alphap+alpham))
                   ds = s(i,j-1)-s(i,j)
                   sgn = sign(ONE,s(i,j+1)-s(i,j-1))
                   if (sgn*dI .ge. sgn*ds) then
                      sm(i,j) = s(i,j) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alphap))
                   end if
                end if
             end if
          end do
       end do

    end if

    ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
    if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sm(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,lo(2)+1) = -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1) &
                                        + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)  ) &
                                        + HALF        *s(lo(1)-1:hi(1)+1,lo(2)+1) &
                                        - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

       ! make sure sedge lies in between adjacent cell-centered values
       do i=lo(1)-1,hi(1)+1
          sedge(i,lo(2)+1) = max(sedge(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
          sedge(i,lo(2)+1) = min(sedge(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
       end do

       ! copy sedge into sp and sm
       do i=lo(1)-1,hi(1)+1
          sp(i,lo(2)  ) = sedge(i,lo(2)+1)
          sm(i,lo(2)+1) = sedge(i,lo(2)+1)
       end do
    end if

    if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sp(lo(1)-1:hi(1)+1,hi(2)) = s(lo(1)-1:hi(1)+1,hi(2)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,hi(2)) = -FIFTH        *s(lo(1)-1:hi(1)+1,hi(2)+1) &
                                      + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,hi(2)  ) &
                                      + HALF        *s(lo(1)-1:hi(1)+1,hi(2)-1) &
                                      - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

       ! make sure sedge lies in between adjacent cell-centered values
       do i=lo(1)-1,hi(1)+1
          sedge(i,hi(2)) = max(sedge(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
          sedge(i,hi(2)) = min(sedge(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
       end do

       ! copy sedge into sp and sm
       do i=lo(1)-1,hi(1)+1
          sp(i,hi(2)-1) = sedge(i,hi(2))
          sm(i,hi(2)  ) = sedge(i,hi(2))
       end do
    end if

    ! compute y-component of Ip and Im
    do j=lo(2)-1,hi(2)+1
       ! compute effect of w0
       if (j .le. 0) then
          w0cc = w0(0)
       else if (j .ge. nr(n)) then
          w0cc = w0(nr(n))
       else
          w0cc = HALF*(w0(j)+w0(j+1))
       end if
       do i=lo(1)-1,hi(1)+1
          sigma = abs(u(i,j,2)+w0cc)*dt/dx(2)
          s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
          if (u(i,j,2)+w0cc .gt. ZERO) then
             Ip(i,j,2) = sp(i,j) - (sigma/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
             Im(i,j,2) = s(i,j)
          else if (u(i,j,2)+w0cc .lt. ZERO) then
             Ip(i,j,2) = s(i,j)
             Im(i,j,2) = sm(i,j) + (sigma/TWO)*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          else
             Ip(i,j,2) = s(i,j)
             Im(i,j,2) = s(i,j)
          end if
       end do
    end do

    deallocate(sp,sm,dsvl,sedge)

  end subroutine ppm_2d

  ! characteristics based on umac
  subroutine ppm_fpu_2d(n,s,ng_s,umac,vmac,ng_um,Ip,Im,w0,lo,hi,bc,dx,dt)

    use bc_module
    use bl_constants_module
    use geometry, only: nr

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_s,ng_um
    real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng_s :,lo(2)-ng_s :)
    real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(inout) ::   Ip(lo(1)-1    :,lo(2)-1    :,:)
    real(kind=dp_t), intent(inout) ::   Im(lo(1)-1    :,lo(2)-1    :,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i,j,edge_interp_type,spm_limiter_type

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham, ds
    real(kind=dp_t) :: dI, sgn
    real(kind=dp_t) :: sigmam, sigmap, s6, w0lo, w0hi

    ! s_{\ib,+}, s_{\ib,-}
    real(kind=dp_t), allocatable :: sp(:,:)
    real(kind=dp_t), allocatable :: sm(:,:)

    ! \delta s_{\ib}^{vL}
    real(kind=dp_t), allocatable :: dsvl(:,:)

    ! s_{i+\half}^{H.O.}
    real(kind=dp_t), allocatable :: sedge(:,:)

    ! cell-centered indexing
    allocate(sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))

    ! constant used in Colella 2008
    C = 1.25d0

    ! 1 = 4th order with van Leer limiting
    ! 2 = 4th order with Colella 2008 limiting
    edge_interp_type = 1

    ! 1 = "old" limiters described in Colella 2008
    ! 2 = "new" limiters described in Colella 2008
    spm_limiter_type = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2,lo(2)-1:hi(2)+1))

    ! edge-centered indexing for x-faces
    allocate(sedge(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1))

    ! compute s at x-edges
    if (edge_interp_type .eq. 1) then

       ! compute van Leer slopes in x-direction
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+2
             dsc = HALF * (s(i+1,j) - s(i-1,j))
             dsl = TWO  * (s(i  ,j) - s(i-1,j))
             dsr = TWO  * (s(i+1,j) - s(i  ,j))
             dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          end do
       end do
       
       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+2
             sedge(i,j) = HALF*(s(i,j)+s(i-1,j)) - SIXTH*(dsvl(i,j)-dsvl(i-1,j))
             ! make sure sedge lies in between adjacent cell-centered values
             sedge(i,j) = max(sedge(i,j),min(s(i,j),s(i-1,j)))
             sedge(i,j) = min(sedge(i,j),max(s(i,j),s(i-1,j)))
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+2
             dsvl(i,j) = HALF * (s(i+1,j) - s(i-1,j))
          end do
       end do
       
       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+2
             sedge(i,j) = HALF*(s(i,j)+s(i-1,j)) - SIXTH*(dsvl(i,j)-dsvl(i-1,j))
             ! if sedge is not in between the neighboring s values, we limit
             if (sedge(i,j) .lt. min(s(i,j),s(i-1,j)) .or. &
                 sedge(i,j) .gt. max(s(i,j),s(i-1,j))) then
                D2  = (THREE/dx(1)**2)*(s(i-1,j)-TWO*sedge(i,j)+s(i,j))
                D2L = (ONE/dx(1)**2)*(s(i-2,j)-TWO*s(i-1,j)+s(i,j))
                D2R = (ONE/dx(1)**2)*(s(i-1,j)-TWO*s(i,j)+s(i+1,j))
                if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R)) then
                   sedge(i,j) = HALF*(s(i-1,j)+s(i,j)) - (dx(1)**2/THREE) &
                        *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                else
                   sedge(i,j) = HALF*(s(i-1,j)+s(i,j))
                end if
             end if
          end do
       end do

    end if

    ! fill x-component of sp and sm
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sp(i,j) = sedge(i+1,j)
          sm(i,j) = sedge(i  ,j)
       end do
    end do

    ! limit x-component of sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify using quadratic limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

       ! modify using Colella 2008 limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO .or. &
                 (s(i+1,j)-s(i,j))*(s(i,j)-s(i-1,j)) .le. ZERO ) then
                s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
                D2  = -TWO*s6/dx(1)**2
                D2C = (ONE/dx(1)**2)*(s(i-1,j)-TWO*s(i,j)+s(i+1,j))
                D2L = (ONE/dx(1)**2)*(s(i-2,j)-TWO*s(i-1,j)+s(i,j))
                D2R = (ONE/dx(1)**2)*(s(i,j)-TWO*s(i+1,j)+s(i+2,j))
                if (sign(ONE,D2) .eq. sign(ONE,D2C) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R) .and. &
                    D2 .ne. ZERO) then
                   D2LIM = sign(ONE,D2)*min(C*abs(D2C),C*abs(D2L),C*abs(D2R),abs(D2))
                   sp(i,j) = s(i,j) + (sp(i,j)-s(i,j))*(D2LIM/D2)
                   sm(i,j) = s(i,j) + (sm(i,j)-s(i,j))*(D2LIM/D2)
                else
                   sp(i,j) = s(i,j)
                   sm(i,j) = s(i,j)
                end if
             else
                alphap = sp(i,j)-s(i,j)
                alpham = sm(i,j)-s(i,j)
                if (abs(alphap) .ge. TWO*abs(alpham)) then
                   dI = -alphap**2 / (FOUR*(alphap+alpham))
                   ds = s(i+1,j)-s(i,j)
                   sgn = sign(ONE,s(i+1,j)-s(i-1,j))
                   if (sgn*dI .ge. sgn*ds) then
                      sp(i,j) = s(i,j) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alpham))
                   end if
                else if (abs(alpham) .ge. TWO*abs(alphap)) then
                   dI = -alpham**2 / (FOUR*(alphap+alpham))
                   ds = s(i-1,j)-s(i,j)
                   sgn = sign(ONE,s(i+1,j)-s(i-1,j))
                   if (sgn*dI .ge. sgn*ds) then
                      sm(i,j) = s(i,j) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alphap))
                   end if
                end if
             end if
          end do
       end do

    end if

    ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
    if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sm(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)+1,lo(2)-1:hi(2)+1) = -FIFTH        *s(lo(1)-1,lo(2)-1:hi(2)+1) &
                                        + (THREE/FOUR)*s(lo(1)  ,lo(2)-1:hi(2)+1) &
                                        + HALF        *s(lo(1)+1,lo(2)-1:hi(2)+1) &
                                        - (ONE/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do j=lo(2)-1,hi(2)+1
          sedge(lo(1)+1,j) = max(sedge(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
          sedge(lo(1)+1,j) = min(sedge(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
       end do

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          sp(lo(1)  ,j) = sedge(lo(1)+1,j)
          sm(lo(1)+1,j) = sedge(lo(1)+1,j)
       end do
    end if

    if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sp(hi(1),lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(hi(1),lo(2)-1:hi(2)+1) = -FIFTH        *s(hi(1)+1,lo(2)-1:hi(2)+1) &
                                      + (THREE/FOUR)*s(hi(1)  ,lo(2)-1:hi(2)+1) &
                                      + HALF        *s(hi(1)-1,lo(2)-1:hi(2)+1) &
                                      - (ONE/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do j=lo(2)-1,hi(2)+1
          sedge(hi(1),j) = max(sedge(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
          sedge(hi(1),j) = min(sedge(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
       end do

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          sp(hi(1)-1,j) = sedge(hi(1),j)
          sm(hi(1)  ,j) = sedge(hi(1),j)
       end do
    end if

    ! compute x-component of Ip and Im
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sigmap = abs(umac(i+1,j))*dt/dx(1)
          sigmam = abs(umac(i,j))*dt/dx(1)
          s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
          if (umac(i+1,j) .gt. ZERO) then
             Ip(i,j,1) = sp(i,j) - (sigmap/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigmap)*s6)
          else
             Ip(i,j,1) = s(i,j)
          end if
          if (umac(i,j) .lt. ZERO) then
             Im(i,j,1) = sm(i,j) + (sigmam/TWO)*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigmam)*s6)
          else
             Im(i,j,1) = s(i,j)
          end if
       end do
    end do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+2))

    ! edge-centered indexing for y-faces
    allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2))

    ! compute s at y-edges
    if (edge_interp_type .eq. 1) then

       ! compute van Leer slopes in y-direction
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             dsc = HALF * (s(i,j+1) - s(i,j-1))
             dsl = TWO  * (s(i,j  ) - s(i,j-1))
             dsr = TWO  * (s(i,j+1) - s(i,j  ))
             dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          end do
       end do
       
       ! interpolate s to y-edges
       do j=lo(2)-1,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             sedge(i,j) = HALF*(s(i,j)+s(i,j-1)) - SIXTH*(dsvl(i,j)-dsvl(i,j-1))
             ! make sure sedge lies in between adjacent cell-centered values
             sedge(i,j) = max(sedge(i,j),min(s(i,j),s(i,j-1)))
             sedge(i,j) = min(sedge(i,j),max(s(i,j),s(i,j-1)))
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             dsvl(i,j) = HALF * (s(i,j+1) - s(i,j-1))
          end do
       end do
       
       ! interpolate s to y-edges
       do j=lo(2)-1,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             sedge(i,j) = HALF*(s(i,j)+s(i,j-1)) - SIXTH*(dsvl(i,j)-dsvl(i,j-1))
             ! if sedge is not in between the neighboring s values, we limit
             if (sedge(i,j) .lt. min(s(i,j),s(i,j-1)) .or. &
                 sedge(i,j) .gt. max(s(i,j),s(i,j-1))) then
                D2  = (THREE/dx(2)**2)*(s(i,j-1)-TWO*sedge(i,j)+s(i,j))
                D2L = (ONE/dx(2)**2)*(s(i,j-2)-TWO*s(i,j-1)+s(i,j))
                D2R = (ONE/dx(2)**2)*(s(i,j-1)-TWO*s(i,j)+s(i,j+1))
                if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R)) then
                   sedge(i,j) = HALF*(s(i,j-1)+s(i,j)) - (dx(2)**2/THREE) &
                        *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                else
                   sedge(i,j) = HALF*(s(i,j-1)+s(i,j))
                end if
             end if
          end do
       end do

    end if

    ! fill y-component of sp and sm
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sp(i,j) = sedge(i,j+1)
          sm(i,j) = sedge(i,j  )
       end do
    end do

    ! limit y-component of sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify using quadratic limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

       ! modify using Colella 2008 limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO .or. &
                 (s(i,j+1)-s(i,j))*(s(i,j)-s(i,j-1)) .le. ZERO ) then
                s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
                D2  = -TWO*s6/dx(2)**2
                D2C = (ONE/dx(2)**2)*(s(i,j-1)-TWO*s(i,j)+s(i,j+1))
                D2L = (ONE/dx(2)**2)*(s(i,j-2)-TWO*s(i,j-1)+s(i,j))
                D2R = (ONE/dx(2)**2)*(s(i,j)-TWO*s(i,j+1)+s(i,j+2))
                if (sign(ONE,D2) .eq. sign(ONE,D2C) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                    sign(ONE,D2) .eq. sign(ONE,D2R) .and. &
                    D2 .ne. ZERO) then
                   D2LIM = sign(ONE,D2)*min(C*abs(D2C),C*abs(D2L),C*abs(D2R),abs(D2))
                   sp(i,j) = s(i,j) + (sp(i,j)-s(i,j))*(D2LIM/D2)
                   sm(i,j) = s(i,j) + (sm(i,j)-s(i,j))*(D2LIM/D2)
                else
                   sp(i,j) = s(i,j)
                   sm(i,j) = s(i,j)
                end if
             else
                alphap = sp(i,j)-s(i,j)
                alpham = sm(i,j)-s(i,j)
                if (abs(alphap) .ge. TWO*abs(alpham)) then
                   dI = -alphap**2 / (FOUR*(alphap+alpham))
                   ds = s(i,j+1)-s(i,j)
                   sgn = sign(ONE,s(i,j+1)-s(i,j-1))
                   if (sgn*dI .ge. sgn*ds) then
                      sp(i,j) = s(i,j) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alpham))
                   end if
                else if (abs(alpham) .ge. TWO*abs(alphap)) then
                   dI = -alpham**2 / (FOUR*(alphap+alpham))
                   ds = s(i,j-1)-s(i,j)
                   sgn = sign(ONE,s(i,j+1)-s(i,j-1))
                   if (sgn*dI .ge. sgn*ds) then
                      sm(i,j) = s(i,j) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alphap))
                   end if
                end if
             end if
          end do
       end do

    end if

    ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
    if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sm(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,lo(2)+1) = -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1) &
                                        + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)  ) &
                                        + HALF        *s(lo(1)-1:hi(1)+1,lo(2)+1) &
                                        - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

       ! make sure sedge lies in between adjacent cell-centered values
       do i=lo(1)-1,hi(1)+1
          sedge(i,lo(2)+1) = max(sedge(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
          sedge(i,lo(2)+1) = min(sedge(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
       end do

       ! copy sedge into sp and sm
       do i=lo(1)-1,hi(1)+1
          sp(i,lo(2)  ) = sedge(i,lo(2)+1)
          sm(i,lo(2)+1) = sedge(i,lo(2)+1)
       end do
    end if

    if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sp(lo(1)-1:hi(1)+1,hi(2)) = s(lo(1)-1:hi(1)+1,hi(2)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,hi(2)) = -FIFTH        *s(lo(1)-1:hi(1)+1,hi(2)+1) &
                                      + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,hi(2)  ) &
                                      + HALF        *s(lo(1)-1:hi(1)+1,hi(2)-1) &
                                      - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

       ! make sure sedge lies in between adjacent cell-centered values
       do i=lo(1)-1,hi(1)+1
          sedge(i,hi(2)) = max(sedge(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
          sedge(i,hi(2)) = min(sedge(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
       end do

       ! copy sedge into sp and sm
       do i=lo(1)-1,hi(1)+1
          sp(i,hi(2)-1) = sedge(i,hi(2))
          sm(i,hi(2)  ) = sedge(i,hi(2))
       end do
    end if

    ! compute y-component of Ip and Im
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
          s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
          if (vmac(i,j+1)+w0hi .gt. ZERO) then
             Ip(i,j,2) = sp(i,j) - (sigmap/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigmap)*s6)
          else
             Ip(i,j,2) = s(i,j)
          end if
          if (vmac(i,j)+w0lo .lt. ZERO) then
             Im(i,j,2) = sm(i,j) + (sigmam/TWO)*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigmam)*s6)
          else
             Im(i,j,2) = s(i,j)
          end if
       end do
    end do

    deallocate(sp,sm,dsvl,sedge)

  end subroutine ppm_fpu_2d

  ! characteristics based on u
  subroutine ppm_3d(n,s,ng_s,u,ng_u,Ip,Im,w0,w0macx,w0macy,w0macz, &
                    ng_w0,lo,hi,bc,dx,dt)

    use bc_module
    use bl_constants_module
    use geometry, only: nr, spherical

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_s,ng_u,ng_w0
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,hi(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,hi(3)-ng_u :,:)
    real(kind=dp_t), intent(inout) ::     Ip(lo(1)-1    :,lo(2)-1    :,hi(3)-1    :,:)
    real(kind=dp_t), intent(inout) ::     Im(lo(1)-1    :,lo(2)-1    :,hi(3)-1    :,:)
    real(kind=dp_t), intent(in   ) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,hi(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,hi(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,hi(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i,j,k,edge_interp_type,spm_limiter_type

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham, ds
    real(kind=dp_t) :: dI, sgn
    real(kind=dp_t) :: sigma, s6, w0cc, velcc

    ! s_{\ib,+}, s_{\ib,-}
    real(kind=dp_t), allocatable :: sp(:,:,:)
    real(kind=dp_t), allocatable :: sm(:,:,:)

    ! \delta s_{\ib}^{vL}
    real(kind=dp_t), allocatable :: dsvl(:,:,:)

    ! s_{i+\half}^{H.O.}
    real(kind=dp_t), allocatable :: sedge(:,:,:)

    ! cell-centered indexing
    allocate(sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    ! constant used in Colella 2008
    C = 1.25d0

    ! 1 = 4th order with van Leer limiting
    ! 2 = 4th order with Colella 2008 limiting
    edge_interp_type = 1

    ! 1 = "old" limiters described in Colella 2008
    ! 2 = "new" limiters described in Colella 2008
    spm_limiter_type = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    ! edge-centered indexing for x-faces
    allocate(sedge(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    ! compute s at x-edges
    if (edge_interp_type .eq. 1) then
       
       ! compute van Leer slopes in x-direction
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-2,hi(1)+2
                dsc = HALF * (s(i+1,j,k) - s(i-1,j,k))
                dsl = TWO  * (s(i  ,j,k) - s(i-1,j,k))
                dsr = TWO  * (s(i+1,j,k) - s(i  ,j,k))
                dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             end do
          end do
       end do
       
       ! interpolate s to x-edges
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+2
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i-1,j,k)) - SIXTH*(dsvl(i,j,k)-dsvl(i-1,j,k))
                ! make sure sedge lies in between adjacent cell-centered values
                sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k),s(i-1,j,k)))
                sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k),s(i-1,j,k)))
             end do
          end do
       end do

    else if (edge_interp_type .eq. 2) then
       
       ! store centered differences in dsvl
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-2,hi(1)+2
                dsvl(i,j,k) = HALF * (s(i+1,j,k) - s(i-1,j,k))
             end do
          end do
       end do

       ! interpolate s to x-edges
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+2
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i-1,j,k)) - SIXTH*(dsvl(i,j,k)-dsvl(i-1,j,k))
                ! if sedge is not in between the neighboring s values, we limit
                if (sedge(i,j,k) .lt. min(s(i,j,k),s(i-1,j,k)) .or. &
                     sedge(i,j,k) .gt. max(s(i,j,k),s(i-1,j,k))) then
                   D2  = (THREE/dx(1)**2)*(s(i-1,j,k)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = (ONE/dx(1)**2)*(s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k))
                   D2R = (ONE/dx(1)**2)*(s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k))
                   if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R)) then
                      sedge(i,j,k) = HALF*(s(i-1,j,k)+s(i,j,k)) - (dx(1)**2/THREE) &
                           *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                   else
                      sedge(i,j,k) = HALF*(s(i-1,j,k)+s(i,j,k))
                   end if
                end if
             end do
          end do
       end do

    end if

    ! fill x-component of sp and sm
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j,k) = sedge(i+1,j,k)
             sm(i,j,k) = sedge(i  ,j,k)
          end do
       end do
    end do
    
    ! limit x-component of sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify using quadratic limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                   sp(i,j,k) = s(i,j,k)
                   sm(i,j,k) = s(i,j,k)
                else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                   sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                   sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                end if
             end do
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

       ! modify using Colella 2008 limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO .or. &
                     (s(i+1,j,k)-s(i,j,k))*(s(i,j,k)-s(i-1,j,k)) .le. ZERO ) then
                   s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                   D2  = -TWO*s6/dx(1)**2
                   D2C = (ONE/dx(1)**2)*(s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k))
                   D2L = (ONE/dx(1)**2)*(s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k))
                   D2R = (ONE/dx(1)**2)*(s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k))
                   if (sign(ONE,D2) .eq. sign(ONE,D2C) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R) .and. &
                        D2 .ne. ZERO) then
                      D2LIM = sign(ONE,D2)*min(C*abs(D2C),C*abs(D2L),C*abs(D2R),abs(D2))
                      sp(i,j,k) = s(i,j,k) + (sp(i,j,k)-s(i,j,k))*(D2LIM/D2)
                      sm(i,j,k) = s(i,j,k) + (sm(i,j,k)-s(i,j,k))*(D2LIM/D2)
                   else
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   end if
                else
                   alphap = sp(i,j,k)-s(i,j,k)
                   alpham = sm(i,j,k)-s(i,j,k)
                   if (abs(alphap) .ge. TWO*abs(alpham)) then
                      dI = -alphap**2 / (FOUR*(alphap+alpham))
                      ds = s(i+1,j,k)-s(i,j,k)
                      sgn = sign(ONE,s(i+1,j,k)-s(i-1,j,k))
                      if (sgn*dI .ge. sgn*ds) then
                         sp(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alpham))
                      end if
                   else if (abs(alpham) .ge. TWO*abs(alphap)) then
                      dI = -alpham**2 / (FOUR*(alphap+alpham))
                      ds = s(i-1,j,k)-s(i,j,k)
                      sgn = sign(ONE,s(i+1,j,k)-s(i-1,j,k))
                      if (sgn*dI .ge. sgn*ds) then
                         sm(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alphap))
                      end if
                   end if
                end if
             end do
          end do
       end do

    end if
    
    ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
    if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sm(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
            -FIFTH        *s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            + (THREE/FOUR)*s(lo(1)  ,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            + HALF        *s(lo(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            - (ONE/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             sedge(lo(1)+1,j,k) = max(sedge(lo(1)+1,j,k),min(s(lo(1)+1,j,k),s(lo(1),j,k)))
             sedge(lo(1)+1,j,k) = min(sedge(lo(1)+1,j,k),max(s(lo(1)+1,j,k),s(lo(1),j,k)))
          end do
       end do

       ! copy sedge into sp and sm
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             sp(lo(1)  ,j,k) = sedge(lo(1)+1,j,k)
             sm(lo(1)+1,j,k) = sedge(lo(1)+1,j,k)
          end do
       end do
    end if

    if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sp(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
            -FIFTH        *s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            + (THREE/FOUR)*s(hi(1)  ,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            + HALF        *s(hi(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            - (ONE/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             sedge(hi(1),j,k) = max(sedge(hi(1),j,k),min(s(hi(1)-1,j,k),s(hi(1),j,k)))
             sedge(hi(1),j,k) = min(sedge(hi(1),j,k),max(s(hi(1)-1,j,k),s(hi(1),j,k)))
          end do
       end do

       ! copy sedge into sp and sm
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             sp(hi(1)-1,j,k) = sedge(hi(1),j,k)
             sm(hi(1)  ,j,k) = sedge(hi(1),j,k)
          end do
       end do
    end if
    
    ! compute x-component of Ip and Im
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if (spherical .eq. 1) then
                velcc = u(i,j,k,1)+HALF*(w0macx(i+1,j,k)+w0macx(i,j,k))
             else
                velcc = u(i,j,k,1)
             end if
             sigma = abs(velcc)*dt/dx(1)
             s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
             if (velcc .gt. ZERO) then
                Ip(i,j,k,1) = sp(i,j,k) - (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,1) = s(i,j,k)
             else if (velcc .lt. ZERO) then
                Ip(i,j,k,1) = s(i,j,k)
                Im(i,j,k,1) = sm(i,j,k) + (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,1) = s(i,j,k)
                Im(i,j,k,1) = s(i,j,k)
             end if
          end do
       end do
    end do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+2,lo(3)-1:hi(3)+1))

    ! edge-centered indexing for y-faces
    allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+1))

    ! compute s at y-edges
    if (edge_interp_type .eq. 1) then
       
       ! compute van Leer slopes in y-direction
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-2,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                dsc = HALF * (s(i,j+1,k) - s(i,j-1,k))
                dsl = TWO  * (s(i,j  ,k) - s(i,j-1,k))
                dsr = TWO  * (s(i,j+1,k) - s(i,j  ,k))
                dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             end do
          end do
       end do

       ! interpolate s to y-edges
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i,j-1,k)) - SIXTH*(dsvl(i,j,k)-dsvl(i,j-1,k))
                ! make sure sedge lies in between adjacent cell-centered values
                sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k),s(i,j-1,k)))
                sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k),s(i,j-1,k)))
             end do
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-2,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                dsvl(i,j,k) = HALF * (s(i,j+1,k) - s(i,j-1,k))
             end do
          end do
       end do

       ! interpolate s to y-edges
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i,j-1,k)) - SIXTH*(dsvl(i,j,k)-dsvl(i,j-1,k))
                ! if sedge is not in between the neighboring s values, we limit
                if (sedge(i,j,k) .lt. min(s(i,j,k),s(i,j-1,k)) .or. &
                     sedge(i,j,k) .gt. max(s(i,j,k),s(i,j-1,k))) then
                   D2  = (THREE/dx(2)**2)*(s(i,j-1,k)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = (ONE/dx(2)**2)*(s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k))
                   D2R = (ONE/dx(2)**2)*(s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k))
                   if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R)) then
                      sedge(i,j,k) = HALF*(s(i,j-1,k)+s(i,j,k)) - (dx(2)**2/THREE) &
                           *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                   else
                      sedge(i,j,k) = HALF*(s(i,j-1,k)+s(i,j,k))
                   end if
                end if
             end do
          end do
       end do

    end if

    ! fill y-component of sp and sm
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j,k) = sedge(i,j+1,k)
             sm(i,j,k) = sedge(i,j  ,k)
          end do
       end do
    end do

    ! limit y-component of sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify using quadratic limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                   sp(i,j,k) = s(i,j,k)
                   sm(i,j,k) = s(i,j,k)
                else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                   sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                   sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                end if
             end do
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

       ! modify using Colella 2008 limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO .or. &
                     (s(i,j+1,k)-s(i,j,k))*(s(i,j,k)-s(i,j-1,k)) .le. ZERO ) then
                   s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                   D2  = -TWO*s6/dx(2)**2
                   D2C = (ONE/dx(2)**2)*(s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k))
                   D2L = (ONE/dx(2)**2)*(s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k))
                   D2R = (ONE/dx(2)**2)*(s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k))
                   if (sign(ONE,D2) .eq. sign(ONE,D2C) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R) .and. &
                        D2 .ne. ZERO) then
                      D2LIM = sign(ONE,D2)*min(C*abs(D2C),C*abs(D2L),C*abs(D2R),abs(D2))
                      sp(i,j,k) = s(i,j,k) + (sp(i,j,k)-s(i,j,k))*(D2LIM/D2)
                      sm(i,j,k) = s(i,j,k) + (sm(i,j,k)-s(i,j,k))*(D2LIM/D2)
                   else
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   end if
                else
                   alphap = sp(i,j,k)-s(i,j,k)
                   alpham = sm(i,j,k)-s(i,j,k)
                   if (abs(alphap) .ge. TWO*abs(alpham)) then
                      dI = -alphap**2 / (FOUR*(alphap+alpham))
                      ds = s(i,j+1,k)-s(i,j,k)
                      sgn = sign(ONE,s(i,j+1,k)-s(i,j-1,k))
                      if (sgn*dI .ge. sgn*ds) then
                         sp(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alpham))
                      end if
                   else if (abs(alpham) .ge. TWO*abs(alphap)) then
                      dI = -alpham**2 / (FOUR*(alphap+alpham))
                      ds = s(i,j-1,k)-s(i,j,k)
                      sgn = sign(ONE,s(i,j+1,k)-s(i,j-1,k))
                      if (sgn*dI .ge. sgn*ds) then
                         sm(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alphap))
                      end if
                   end if
                end if
             end do
          end do
       end do

    end if

    ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
    if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sm(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,lo(2)+1,lo(3)-1:hi(3)+1) = &
            -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1) &
            + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)  ,lo(3)-1:hi(3)+1) &
            + HALF        *s(lo(1)-1:hi(1)+1,lo(2)+1,lo(3)-1:hi(3)+1) &
            - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2,lo(3)-1:hi(3)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do k=lo(3)-1,hi(3)+1
          do i=lo(1)-1,hi(1)+1
             sedge(i,lo(2)+1,k) = max(sedge(i,lo(2)+1,k),min(s(i,lo(2)+1,k),s(i,lo(2),k)))
             sedge(i,lo(2)+1,k) = min(sedge(i,lo(2)+1,k),max(s(i,lo(2)+1,k),s(i,lo(2),k)))
          end do
       end do

       ! copy sedge into sp and sm
       do k=lo(3)-1,hi(3)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,lo(2)  ,k) = sedge(i,lo(2)+1,k)
             sm(i,lo(2)+1,k) = sedge(i,lo(2)+1,k)
          end do
       end do
    end if

    if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sp(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = &
            -FIFTH        *s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1) &
            + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,hi(2)  ,lo(3)-1:hi(3)+1) &
            + HALF        *s(lo(1)-1:hi(1)+1,hi(2)-1,lo(3)-1:hi(3)+1) &
            - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2,lo(3)-1:hi(3)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do k=lo(3)-1,hi(3)+1
          do i=lo(1)-1,hi(1)+1
             sedge(i,hi(2),k) = max(sedge(i,hi(2),k),min(s(i,hi(2)-1,k),s(i,hi(2),k)))
             sedge(i,hi(2),k) = min(sedge(i,hi(2),k),max(s(i,hi(2)-1,k),s(i,hi(2),k)))
          end do
       end do

       ! copy sedge into sp and sm
       do k=lo(3)-1,hi(3)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,hi(2)-1,k) = sedge(i,hi(2),k)
             sm(i,hi(2)  ,k) = sedge(i,hi(2),k)
          end do
       end do
    end if

    ! compute y-component of Ip and Im
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if (spherical .eq. 1) then
                velcc = u(i,j,k,2)+HALF*(w0macy(i,j+1,k)+w0macy(i,j,k))
             else
                velcc = u(i,j,k,2)
             end if
             sigma = abs(velcc)*dt/dx(2)
             s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
             if (velcc .gt. ZERO) then
                Ip(i,j,k,2) = sp(i,j,k) - (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,2) = s(i,j,k)
             else if (velcc .lt. ZERO) then
                Ip(i,j,k,2) = s(i,j,k)
                Im(i,j,k,2) = sm(i,j,k) + (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,2) = s(i,j,k)
                Im(i,j,k,2) = s(i,j,k)
             end if
          end do
       end do
    end do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra z-ghost cell
    allocate( dsvl(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-2:hi(3)+2))

    ! edge-centered indexing for z-faces
    allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+2))

    ! compute s at z-edges
    if (edge_interp_type .eq. 1) then
       
       ! compute van Leer slopes in z-direction
       do k=lo(3)-2,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                dsc = HALF * (s(i,j,k+1) - s(i,j,k-1))
                dsl = TWO  * (s(i,j,k  ) - s(i,j,k-1))
                dsr = TWO  * (s(i,j,k+1) - s(i,j,k  ))
                dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             end do
          end do
       end do

       ! interpolate s to z-edges
       do k=lo(3)-1,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i,j,k-1)) - SIXTH*(dsvl(i,j,k)-dsvl(i,j,k-1))
                ! make sure sedge lies in between adjacent cell-centered values
                sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k),s(i,j,k-1)))
                sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k),s(i,j,k-1)))
             end do
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl
       do k=lo(3)-2,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                dsvl(i,j,k) = HALF * (s(i,j,k+1) - s(i,j,k-1))
             end do
          end do
       end do

       ! interpolate s to z-edges
       do k=lo(3)-1,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i,j,k-1)) - SIXTH*(dsvl(i,j,k)-dsvl(i,j,k-1))
                ! if sedge is not in between the neighboring s values, we limit
                if (sedge(i,j,k) .lt. min(s(i,j,k),s(i,j,k-1)) .or. &
                     sedge(i,j,k) .gt. max(s(i,j,k),s(i,j,k-1))) then
                   D2  = (THREE/dx(3)**2)*(s(i,j,k-1)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = (ONE/dx(3)**2)*(s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k))
                   D2R = (ONE/dx(3)**2)*(s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1))
                   if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R)) then
                      sedge(i,j,k) = HALF*(s(i,j,k-1)+s(i,j,k)) - (dx(3)**2/THREE) &
                           *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                   else
                      sedge(i,j,k) = HALF*(s(i,j,k-1)+s(i,j,k))
                   end if
                end if
             end do
          end do
       end do

    end if

    ! fill z-component of sp and sm
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j,k) = sedge(i,j,k+1)
             sm(i,j,k) = sedge(i,j,k  )
          end do
       end do
    end do

    ! limit z-component of sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify using quadratic limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                   sp(i,j,k) = s(i,j,k)
                   sm(i,j,k) = s(i,j,k)
                else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                   sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                   sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                end if
             end do
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

       ! modify using Colella 2008 limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO .or. &
                     (s(i,j,k+1)-s(i,j,k))*(s(i,j,k)-s(i,j,k-1)) .le. ZERO ) then
                   s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                   D2  = -TWO*s6/dx(3)**2
                   D2C = (ONE/dx(3)**2)*(s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1))
                   D2L = (ONE/dx(3)**2)*(s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k))
                   D2R = (ONE/dx(3)**2)*(s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2))
                   if (sign(ONE,D2) .eq. sign(ONE,D2C) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R) .and. &
                        D2 .ne. ZERO) then
                      D2LIM = sign(ONE,D2)*min(C*abs(D2C),C*abs(D2L),C*abs(D2R),abs(D2))
                      sp(i,j,k) = s(i,j,k) + (sp(i,j,k)-s(i,j,k))*(D2LIM/D2)
                      sm(i,j,k) = s(i,j,k) + (sm(i,j,k)-s(i,j,k))*(D2LIM/D2)
                   else
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   end if
                else
                   alphap = sp(i,j,k)-s(i,j,k)
                   alpham = sm(i,j,k)-s(i,j,k)
                   if (abs(alphap) .ge. TWO*abs(alpham)) then
                      dI = -alphap**2 / (FOUR*(alphap+alpham))
                      ds = s(i,j,k+1)-s(i,j,k)
                      sgn = sign(ONE,s(i,j,k+1)-s(i,j,k-1))
                      if (sgn*dI .ge. sgn*ds) then
                         sp(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alpham))
                      end if
                   else if (abs(alpham) .ge. TWO*abs(alphap)) then
                      dI = -alpham**2 / (FOUR*(alphap+alpham))
                      ds = s(i,j,k-1)-s(i,j,k)
                      sgn = sign(ONE,s(i,j,k+1)-s(i,j,k-1))
                      if (sgn*dI .ge. sgn*ds) then
                         sm(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alphap))
                      end if
                   end if
                end if
             end do
          end do
       end do

    end if

    ! different stencil needed for z-component of EXT_DIR and HOEXTRAP bc's
    if (bc(3,1) .eq. EXT_DIR  .or. bc(3,1) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)+1) = &
            -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1) &
            + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)  ) &
            + HALF        *s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)+1) &
            - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)+2)

       ! make sure sedge lies in between adjacent cell-centered values
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sedge(i,j,lo(3)+1) = max(sedge(i,j,lo(3)+1),min(s(i,j,lo(3)+1),s(i,j,lo(3))))
             sedge(i,j,lo(3)+1) = min(sedge(i,j,lo(3)+1),max(s(i,j,lo(3)+1),s(i,j,lo(3))))
          end do
       end do

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j,lo(3)  ) = sedge(i,j,lo(3)+1)
             sm(i,j,lo(3)+1) = sedge(i,j,lo(3)+1)
          end do
       end do
    end if

    if (bc(3,2) .eq. EXT_DIR  .or. bc(3,2) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = &
            -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1) &
            + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)  ) &
            + HALF        *s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)-1) &
            - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)-2)

       ! make sure sedge lies in between adjacent cell-centered values
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sedge(i,j,hi(3)) = max(sedge(i,j,hi(3)),min(s(i,j,hi(3)-1),s(i,j,hi(3))))
             sedge(i,j,hi(3)) = min(sedge(i,j,hi(3)),max(s(i,j,hi(3)-1),s(i,j,hi(3))))
          end do
       end do

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j,hi(3)-1) = sedge(i,j,hi(3))
             sm(i,j,hi(3)  ) = sedge(i,j,hi(3))
          end do
       end do
    end if

    ! compute z-component of Ip and Im
    do k=lo(3)-1,hi(3)+1
       ! compute effect of w0 in plane-parallel
       if (spherical .eq. 0) then
          if (k .le. 0) then
             w0cc = w0(0)
          else if (k .ge. nr(n)) then
             w0cc = w0(nr(n))
          else
             w0cc = HALF*(w0(k)+w0(k+1))
          end if
       end if
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if (spherical .eq. 1) then
                velcc = u(i,j,k,3)+HALF*(w0macz(i,j,k+1)+w0macz(i,j,k))
             else
                velcc = u(i,j,k,3) + w0cc
             end if
             sigma = abs(velcc)*dt/dx(3)
             s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
             if (velcc .gt. ZERO) then
                Ip(i,j,k,3) = sp(i,j,k) - (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,3) = s(i,j,k)
             else if (velcc .lt. ZERO) then
                Ip(i,j,k,3) = s(i,j,k)
                Im(i,j,k,3) = sm(i,j,k) + (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,3) = s(i,j,k)
                Im(i,j,k,3) = s(i,j,k)
             end if
          end do
       end do
    end do

    deallocate(sp,sm,dsvl,sedge)

  end subroutine ppm_3d

  ! characteristics based on umac
  subroutine ppm_fpu_3d(n,s,ng_s,umac,vmac,wmac,ng_um,Ip,Im,w0,w0macx,w0macy,w0macz, &
                        ng_w0,lo,hi,bc,dx,dt)

    use bc_module
    use bl_constants_module
    use geometry, only: nr, spherical

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_s,ng_um,ng_w0
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,hi(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:,hi(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:,hi(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)-ng_um:,lo(2)-ng_um:,hi(3)-ng_um:)
    real(kind=dp_t), intent(inout) ::     Ip(lo(1)-1    :,lo(2)-1    :,hi(3)-1    :,:)
    real(kind=dp_t), intent(inout) ::     Im(lo(1)-1    :,lo(2)-1    :,hi(3)-1    :,:)
    real(kind=dp_t), intent(in   ) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,hi(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,hi(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,hi(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i,j,k,edge_interp_type,spm_limiter_type

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham, ds
    real(kind=dp_t) :: dI, sgn
    real(kind=dp_t) :: sigmam, sigmap, s6, w0lo, w0hi, vello, velhi

    ! s_{\ib,+}, s_{\ib,-}
    real(kind=dp_t), allocatable :: sp(:,:,:)
    real(kind=dp_t), allocatable :: sm(:,:,:)

    ! \delta s_{\ib}^{vL}
    real(kind=dp_t), allocatable :: dsvl(:,:,:)

    ! s_{i+\half}^{H.O.}
    real(kind=dp_t), allocatable :: sedge(:,:,:)

    ! cell-centered indexing
    allocate(sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    ! constant used in Colella 2008
    C = 1.25d0

    ! 1 = 4th order with van Leer limiting
    ! 2 = 4th order with Colella 2008 limiting
    edge_interp_type = 1

    ! 1 = "old" limiters described in Colella 2008
    ! 2 = "new" limiters described in Colella 2008
    spm_limiter_type = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    ! edge-centered indexing for x-faces
    allocate(sedge(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    ! compute s at x-edges
    if (edge_interp_type .eq. 1) then
       
       ! compute van Leer slopes in x-direction
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-2,hi(1)+2
                dsc = HALF * (s(i+1,j,k) - s(i-1,j,k))
                dsl = TWO  * (s(i  ,j,k) - s(i-1,j,k))
                dsr = TWO  * (s(i+1,j,k) - s(i  ,j,k))
                dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             end do
          end do
       end do
       
       ! interpolate s to x-edges
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+2
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i-1,j,k)) - SIXTH*(dsvl(i,j,k)-dsvl(i-1,j,k))
                ! make sure sedge lies in between adjacent cell-centered values
                sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k),s(i-1,j,k)))
                sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k),s(i-1,j,k)))
             end do
          end do
       end do

    else if (edge_interp_type .eq. 2) then
       
       ! store centered differences in dsvl
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-2,hi(1)+2
                dsvl(i,j,k) = HALF * (s(i+1,j,k) - s(i-1,j,k))
             end do
          end do
       end do

       ! interpolate s to x-edges
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+2
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i-1,j,k)) - SIXTH*(dsvl(i,j,k)-dsvl(i-1,j,k))
                ! if sedge is not in between the neighboring s values, we limit
                if (sedge(i,j,k) .lt. min(s(i,j,k),s(i-1,j,k)) .or. &
                     sedge(i,j,k) .gt. max(s(i,j,k),s(i-1,j,k))) then
                   D2  = (THREE/dx(1)**2)*(s(i-1,j,k)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = (ONE/dx(1)**2)*(s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k))
                   D2R = (ONE/dx(1)**2)*(s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k))
                   if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R)) then
                      sedge(i,j,k) = HALF*(s(i-1,j,k)+s(i,j,k)) - (dx(1)**2/THREE) &
                           *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                   else
                      sedge(i,j,k) = HALF*(s(i-1,j,k)+s(i,j,k))
                   end if
                end if
             end do
          end do
       end do

    end if

    ! fill x-component of sp and sm
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j,k) = sedge(i+1,j,k)
             sm(i,j,k) = sedge(i  ,j,k)
          end do
       end do
    end do
    
    ! limit x-component of sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify using quadratic limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                   sp(i,j,k) = s(i,j,k)
                   sm(i,j,k) = s(i,j,k)
                else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                   sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                   sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                end if
             end do
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

       ! modify using Colella 2008 limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO .or. &
                     (s(i+1,j,k)-s(i,j,k))*(s(i,j,k)-s(i-1,j,k)) .le. ZERO ) then
                   s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                   D2  = -TWO*s6/dx(1)**2
                   D2C = (ONE/dx(1)**2)*(s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k))
                   D2L = (ONE/dx(1)**2)*(s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k))
                   D2R = (ONE/dx(1)**2)*(s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k))
                   if (sign(ONE,D2) .eq. sign(ONE,D2C) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R) .and. &
                        D2 .ne. ZERO) then
                      D2LIM = sign(ONE,D2)*min(C*abs(D2C),C*abs(D2L),C*abs(D2R),abs(D2))
                      sp(i,j,k) = s(i,j,k) + (sp(i,j,k)-s(i,j,k))*(D2LIM/D2)
                      sm(i,j,k) = s(i,j,k) + (sm(i,j,k)-s(i,j,k))*(D2LIM/D2)
                   else
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   end if
                else
                   alphap = sp(i,j,k)-s(i,j,k)
                   alpham = sm(i,j,k)-s(i,j,k)
                   if (abs(alphap) .ge. TWO*abs(alpham)) then
                      dI = -alphap**2 / (FOUR*(alphap+alpham))
                      ds = s(i+1,j,k)-s(i,j,k)
                      sgn = sign(ONE,s(i+1,j,k)-s(i-1,j,k))
                      if (sgn*dI .ge. sgn*ds) then
                         sp(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alpham))
                      end if
                   else if (abs(alpham) .ge. TWO*abs(alphap)) then
                      dI = -alpham**2 / (FOUR*(alphap+alpham))
                      ds = s(i-1,j,k)-s(i,j,k)
                      sgn = sign(ONE,s(i+1,j,k)-s(i-1,j,k))
                      if (sgn*dI .ge. sgn*ds) then
                         sm(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alphap))
                      end if
                   end if
                end if
             end do
          end do
       end do

    end if
    
    ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
    if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sm(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
            -FIFTH        *s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            + (THREE/FOUR)*s(lo(1)  ,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            + HALF        *s(lo(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            - (ONE/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             sedge(lo(1)+1,j,k) = max(sedge(lo(1)+1,j,k),min(s(lo(1)+1,j,k),s(lo(1),j,k)))
             sedge(lo(1)+1,j,k) = min(sedge(lo(1)+1,j,k),max(s(lo(1)+1,j,k),s(lo(1),j,k)))
          end do
       end do

       ! copy sedge into sp and sm
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             sp(lo(1)  ,j,k) = sedge(lo(1)+1,j,k)
             sm(lo(1)+1,j,k) = sedge(lo(1)+1,j,k)
          end do
       end do
    end if

    if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sp(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
            -FIFTH        *s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            + (THREE/FOUR)*s(hi(1)  ,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            + HALF        *s(hi(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) &
            - (ONE/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             sedge(hi(1),j,k) = max(sedge(hi(1),j,k),min(s(hi(1)-1,j,k),s(hi(1),j,k)))
             sedge(hi(1),j,k) = min(sedge(hi(1),j,k),max(s(hi(1)-1,j,k),s(hi(1),j,k)))
          end do
       end do

       ! copy sedge into sp and sm
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             sp(hi(1)-1,j,k) = sedge(hi(1),j,k)
             sm(hi(1)  ,j,k) = sedge(hi(1),j,k)
          end do
       end do
    end if
    
    ! compute x-component of Ip and Im
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if (spherical .eq. 1) then
                velhi = umac(i+1,j,k) + w0macx(i+1,j,k)
                vello = umac(i  ,j,k) + w0macx(i  ,j,k)
             else
                velhi = umac(i+1,j,k)
                vello = umac(i  ,j,k)
             end if
             sigmap = abs(velhi)*dt/dx(1)
             sigmam = abs(vello)*dt/dx(1)
             s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
             if (velhi .gt. ZERO) then
                Ip(i,j,k,1) = sp(i,j,k) - (sigmap/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigmap)*s6)
             else
                Ip(i,j,k,1) = s(i,j,k)
             end if
             if (vello .lt. ZERO) then
                Im(i,j,k,1) = sm(i,j,k) + (sigmam/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigmam)*s6)
             else
                Im(i,j,k,1) = s(i,j,k)
             end if
          end do
       end do
    end do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+2,lo(3)-1:hi(3)+1))

    ! edge-centered indexing for y-faces
    allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+1))

    ! compute s at y-edges
    if (edge_interp_type .eq. 1) then
       
       ! compute van Leer slopes in y-direction
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-2,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                dsc = HALF * (s(i,j+1,k) - s(i,j-1,k))
                dsl = TWO  * (s(i,j  ,k) - s(i,j-1,k))
                dsr = TWO  * (s(i,j+1,k) - s(i,j  ,k))
                dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             end do
          end do
       end do

       ! interpolate s to y-edges
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i,j-1,k)) - SIXTH*(dsvl(i,j,k)-dsvl(i,j-1,k))
                ! make sure sedge lies in between adjacent cell-centered values
                sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k),s(i,j-1,k)))
                sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k),s(i,j-1,k)))
             end do
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-2,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                dsvl(i,j,k) = HALF * (s(i,j+1,k) - s(i,j-1,k))
             end do
          end do
       end do

       ! interpolate s to y-edges
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i,j-1,k)) - SIXTH*(dsvl(i,j,k)-dsvl(i,j-1,k))
                ! if sedge is not in between the neighboring s values, we limit
                if (sedge(i,j,k) .lt. min(s(i,j,k),s(i,j-1,k)) .or. &
                     sedge(i,j,k) .gt. max(s(i,j,k),s(i,j-1,k))) then
                   D2  = (THREE/dx(2)**2)*(s(i,j-1,k)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = (ONE/dx(2)**2)*(s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k))
                   D2R = (ONE/dx(2)**2)*(s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k))
                   if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R)) then
                      sedge(i,j,k) = HALF*(s(i,j-1,k)+s(i,j,k)) - (dx(2)**2/THREE) &
                           *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                   else
                      sedge(i,j,k) = HALF*(s(i,j-1,k)+s(i,j,k))
                   end if
                end if
             end do
          end do
       end do

    end if

    ! fill y-component of sp and sm
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j,k) = sedge(i,j+1,k)
             sm(i,j,k) = sedge(i,j  ,k)
          end do
       end do
    end do

    ! limit y-component of sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify using quadratic limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                   sp(i,j,k) = s(i,j,k)
                   sm(i,j,k) = s(i,j,k)
                else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                   sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                   sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                end if
             end do
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

       ! modify using Colella 2008 limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO .or. &
                     (s(i,j+1,k)-s(i,j,k))*(s(i,j,k)-s(i,j-1,k)) .le. ZERO ) then
                   s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                   D2  = -TWO*s6/dx(2)**2
                   D2C = (ONE/dx(2)**2)*(s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k))
                   D2L = (ONE/dx(2)**2)*(s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k))
                   D2R = (ONE/dx(2)**2)*(s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k))
                   if (sign(ONE,D2) .eq. sign(ONE,D2C) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R) .and. &
                        D2 .ne. ZERO) then
                      D2LIM = sign(ONE,D2)*min(C*abs(D2C),C*abs(D2L),C*abs(D2R),abs(D2))
                      sp(i,j,k) = s(i,j,k) + (sp(i,j,k)-s(i,j,k))*(D2LIM/D2)
                      sm(i,j,k) = s(i,j,k) + (sm(i,j,k)-s(i,j,k))*(D2LIM/D2)
                   else
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   end if
                else
                   alphap = sp(i,j,k)-s(i,j,k)
                   alpham = sm(i,j,k)-s(i,j,k)
                   if (abs(alphap) .ge. TWO*abs(alpham)) then
                      dI = -alphap**2 / (FOUR*(alphap+alpham))
                      ds = s(i,j+1,k)-s(i,j,k)
                      sgn = sign(ONE,s(i,j+1,k)-s(i,j-1,k))
                      if (sgn*dI .ge. sgn*ds) then
                         sp(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alpham))
                      end if
                   else if (abs(alpham) .ge. TWO*abs(alphap)) then
                      dI = -alpham**2 / (FOUR*(alphap+alpham))
                      ds = s(i,j-1,k)-s(i,j,k)
                      sgn = sign(ONE,s(i,j+1,k)-s(i,j-1,k))
                      if (sgn*dI .ge. sgn*ds) then
                         sm(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alphap))
                      end if
                   end if
                end if
             end do
          end do
       end do

    end if

    ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
    if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sm(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,lo(2)+1,lo(3)-1:hi(3)+1) = &
            -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1) &
            + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)  ,lo(3)-1:hi(3)+1) &
            + HALF        *s(lo(1)-1:hi(1)+1,lo(2)+1,lo(3)-1:hi(3)+1) &
            - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2,lo(3)-1:hi(3)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do k=lo(3)-1,hi(3)+1
          do i=lo(1)-1,hi(1)+1
             sedge(i,lo(2)+1,k) = max(sedge(i,lo(2)+1,k),min(s(i,lo(2)+1,k),s(i,lo(2),k)))
             sedge(i,lo(2)+1,k) = min(sedge(i,lo(2)+1,k),max(s(i,lo(2)+1,k),s(i,lo(2),k)))
          end do
       end do

       ! copy sedge into sp and sm
       do k=lo(3)-1,hi(3)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,lo(2)  ,k) = sedge(i,lo(2)+1,k)
             sm(i,lo(2)+1,k) = sedge(i,lo(2)+1,k)
          end do
       end do
    end if

    if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sp(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = &
            -FIFTH        *s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1) &
            + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,hi(2)  ,lo(3)-1:hi(3)+1) &
            + HALF        *s(lo(1)-1:hi(1)+1,hi(2)-1,lo(3)-1:hi(3)+1) &
            - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2,lo(3)-1:hi(3)+1)

       ! make sure sedge lies in between adjacent cell-centered values
       do k=lo(3)-1,hi(3)+1
          do i=lo(1)-1,hi(1)+1
             sedge(i,hi(2),k) = max(sedge(i,hi(2),k),min(s(i,hi(2)-1,k),s(i,hi(2),k)))
             sedge(i,hi(2),k) = min(sedge(i,hi(2),k),max(s(i,hi(2)-1,k),s(i,hi(2),k)))
          end do
       end do

       ! copy sedge into sp and sm
       do k=lo(3)-1,hi(3)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,hi(2)-1,k) = sedge(i,hi(2),k)
             sm(i,hi(2)  ,k) = sedge(i,hi(2),k)
          end do
       end do
    end if

    ! compute y-component of Ip and Im
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if (spherical .eq. 1) then
                velhi = vmac(i,j+1,k) + w0macy(i,j+1,k)
                vello = vmac(i,j  ,k) + w0macy(i,j  ,k)
             else
                velhi = vmac(i,j+1,k)
                vello = vmac(i,j  ,k)
             end if
             sigmap = abs(velhi)*dt/dx(2)
             sigmam = abs(vello)*dt/dx(2)
             s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
             if (velhi .gt. ZERO) then
                Ip(i,j,k,2) = sp(i,j,k) - (sigmap/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigmap)*s6)
             else
                Ip(i,j,k,2) = s(i,j,k)
             end if
             if (vello .lt. ZERO) then
                Im(i,j,k,2) = sm(i,j,k) + (sigmam/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigmam)*s6)
             else
                Im(i,j,k,2) = s(i,j,k)
             end if
          end do
       end do
    end do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra z-ghost cell
    allocate( dsvl(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-2:hi(3)+2))

    ! edge-centered indexing for z-faces
    allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+2))

    ! compute s at z-edges
    if (edge_interp_type .eq. 1) then
       
       ! compute van Leer slopes in z-direction
       do k=lo(3)-2,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                dsc = HALF * (s(i,j,k+1) - s(i,j,k-1))
                dsl = TWO  * (s(i,j,k  ) - s(i,j,k-1))
                dsr = TWO  * (s(i,j,k+1) - s(i,j,k  ))
                dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             end do
          end do
       end do

       ! interpolate s to z-edges
       do k=lo(3)-1,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i,j,k-1)) - SIXTH*(dsvl(i,j,k)-dsvl(i,j,k-1))
                ! make sure sedge lies in between adjacent cell-centered values
                sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k),s(i,j,k-1)))
                sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k),s(i,j,k-1)))
             end do
          end do
       end do

    else if (edge_interp_type .eq. 2) then

       ! store centered differences in dsvl
       do k=lo(3)-2,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                dsvl(i,j,k) = HALF * (s(i,j,k+1) - s(i,j,k-1))
             end do
          end do
       end do

       ! interpolate s to z-edges
       do k=lo(3)-1,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i,j,k-1)) - SIXTH*(dsvl(i,j,k)-dsvl(i,j,k-1))
                ! if sedge is not in between the neighboring s values, we limit
                if (sedge(i,j,k) .lt. min(s(i,j,k),s(i,j,k-1)) .or. &
                     sedge(i,j,k) .gt. max(s(i,j,k),s(i,j,k-1))) then
                   D2  = (THREE/dx(3)**2)*(s(i,j,k-1)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = (ONE/dx(3)**2)*(s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k))
                   D2R = (ONE/dx(3)**2)*(s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1))
                   if (sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R)) then
                      sedge(i,j,k) = HALF*(s(i,j,k-1)+s(i,j,k)) - (dx(3)**2/THREE) &
                           *sign(ONE,D2)*min(C*abs(D2L),C*abs(D2R),abs(D2))
                   else
                      sedge(i,j,k) = HALF*(s(i,j,k-1)+s(i,j,k))
                   end if
                end if
             end do
          end do
       end do

    end if

    ! fill z-component of sp and sm
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j,k) = sedge(i,j,k+1)
             sm(i,j,k) = sedge(i,j,k  )
          end do
       end do
    end do

    ! limit z-component of sp and sm
    if (spm_limiter_type .eq. 1) then

       ! modify using quadratic limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                   sp(i,j,k) = s(i,j,k)
                   sm(i,j,k) = s(i,j,k)
                else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                   sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                   sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                end if
             end do
          end do
       end do

    else if (spm_limiter_type .eq. 2) then

       ! modify using Colella 2008 limiters
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO .or. &
                     (s(i,j,k+1)-s(i,j,k))*(s(i,j,k)-s(i,j,k-1)) .le. ZERO ) then
                   s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                   D2  = -TWO*s6/dx(3)**2
                   D2C = (ONE/dx(3)**2)*(s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1))
                   D2L = (ONE/dx(3)**2)*(s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k))
                   D2R = (ONE/dx(3)**2)*(s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2))
                   if (sign(ONE,D2) .eq. sign(ONE,D2C) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2L) .and. &
                        sign(ONE,D2) .eq. sign(ONE,D2R) .and. &
                        D2 .ne. ZERO) then
                      D2LIM = sign(ONE,D2)*min(C*abs(D2C),C*abs(D2L),C*abs(D2R),abs(D2))
                      sp(i,j,k) = s(i,j,k) + (sp(i,j,k)-s(i,j,k))*(D2LIM/D2)
                      sm(i,j,k) = s(i,j,k) + (sm(i,j,k)-s(i,j,k))*(D2LIM/D2)
                   else
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   end if
                else
                   alphap = sp(i,j,k)-s(i,j,k)
                   alpham = sm(i,j,k)-s(i,j,k)
                   if (abs(alphap) .ge. TWO*abs(alpham)) then
                      dI = -alphap**2 / (FOUR*(alphap+alpham))
                      ds = s(i,j,k+1)-s(i,j,k)
                      sgn = sign(ONE,s(i,j,k+1)-s(i,j,k-1))
                      if (sgn*dI .ge. sgn*ds) then
                         sp(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alpham))
                      end if
                   else if (abs(alpham) .ge. TWO*abs(alphap)) then
                      dI = -alpham**2 / (FOUR*(alphap+alpham))
                      ds = s(i,j,k-1)-s(i,j,k)
                      sgn = sign(ONE,s(i,j,k+1)-s(i,j,k-1))
                      if (sgn*dI .ge. sgn*ds) then
                         sm(i,j,k) = s(i,j,k) - (TWO*ds + TWO*sgn*sqrt(ds**2 - ds*alphap))
                      end if
                   end if
                end if
             end do
          end do
       end do

    end if

    ! different stencil needed for z-component of EXT_DIR and HOEXTRAP bc's
    if (bc(3,1) .eq. EXT_DIR  .or. bc(3,1) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)+1) = &
            -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1) &
            + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)  ) &
            + HALF        *s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)+1) &
            - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)+2)

       ! make sure sedge lies in between adjacent cell-centered values
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sedge(i,j,lo(3)+1) = max(sedge(i,j,lo(3)+1),min(s(i,j,lo(3)+1),s(i,j,lo(3))))
             sedge(i,j,lo(3)+1) = min(sedge(i,j,lo(3)+1),max(s(i,j,lo(3)+1),s(i,j,lo(3))))
          end do
       end do

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j,lo(3)  ) = sedge(i,j,lo(3)+1)
             sm(i,j,lo(3)+1) = sedge(i,j,lo(3)+1)
          end do
       end do
    end if

    if (bc(3,2) .eq. EXT_DIR  .or. bc(3,2) .eq. HOEXTRAP) then
       ! the value in the first cc ghost cell represents the edge value
       sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)

       ! use a modified stencil to get sedge on the first interior edge
       sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = &
            -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1) &
            + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)  ) &
            + HALF        *s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)-1) &
            - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)-2)

       ! make sure sedge lies in between adjacent cell-centered values
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sedge(i,j,hi(3)) = max(sedge(i,j,hi(3)),min(s(i,j,hi(3)-1),s(i,j,hi(3))))
             sedge(i,j,hi(3)) = min(sedge(i,j,hi(3)),max(s(i,j,hi(3)-1),s(i,j,hi(3))))
          end do
       end do

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j,hi(3)-1) = sedge(i,j,hi(3))
             sm(i,j,hi(3)  ) = sedge(i,j,hi(3))
          end do
       end do
    end if

    ! compute z-component of Ip and Im
    do k=lo(3)-1,hi(3)+1
       ! compute effect of w0 in plane-parallel
       if (spherical .eq. 0) then
          if (k .lt. 0) then
             w0lo = w0(0)
             w0hi = w0(0)
          else if (k .gt. nr(n)-1) then
             w0lo = w0(nr(n))
             w0hi = w0(nr(n))
          else
             w0lo = w0(k)
             w0hi = w0(k+1)
          end if
       end if
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if (spherical .eq. 1) then
                velhi = vmac(i,j+1,k) + w0macy(i,j+1,k)
                vello = vmac(i,j  ,k) + w0macy(i,j  ,k)
             else
                velhi = vmac(i,j+1,k) + w0hi
                vello = vmac(i,j  ,k) + w0lo
             end if
             sigmap = abs(velhi)*dt/dx(2)
             sigmam = abs(vello)*dt/dx(2)
             s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
             if (velhi .gt. ZERO) then
                Ip(i,j,k,3) = sp(i,j,k) - (sigmap/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigmap)*s6)
             else
                Ip(i,j,k,3) = s(i,j,k)
             end if
             if (vello .lt. ZERO) then
                Im(i,j,k,3) = sm(i,j,k) + (sigmam/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigmam)*s6)
             else
                Im(i,j,k,3) = s(i,j,k)
             end if
          end do
       end do
    end do

    deallocate(sp,sm,dsvl,sedge)

  end subroutine ppm_fpu_3d

end module ppm_module
