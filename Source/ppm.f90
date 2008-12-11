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
    real(kind=dp_t), intent(in   ) ::  s(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(in   ) ::  u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real(kind=dp_t), intent(inout) :: Ip(lo(1)-1   :,lo(2)-1   :,:) 
    real(kind=dp_t), intent(inout) :: Im(lo(1)-1   :,lo(2)-1   :,:) 
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: bc(:,:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i,j

    real(kind=dp_t) :: dsl, dsr, dsc, sigma, s6, v_from_w0

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

    ! modify sp and sm using quadratic limiters
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          if ((sp(i,j,1)-s(i,j))*(s(i,j)-sm(i,j,1)) .le. ZERO) then
             sp(i,j,1) = s(i,j)
             sm(i,j,1) = s(i,j)
          end if
          if (abs(sp(i,j,1)-s(i,j)) .ge. TWO*abs(sm(i,j,1)-s(i,j))) then
             sp(i,j,1) = THREE*s(i,j) - TWO*sm(i,j,1)
          end if
          if (abs(sm(i,j,1)-s(i,j)) .ge. TWO*abs(sp(i,j,1)-s(i,j))) then
             sm(i,j,1) = THREE*s(i,j) - TWO*sp(i,j,1)
          end if
       end do
    end do

    ! compute Ip and Im
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sigma = abs(u(i,j,1))*dt/dx(1)
          s6 = SIX*s(i,j) - THREE*(sm(i,j,1)+sp(i,j,1))
          if (u(i,j,1) .gt. ZERO) then
             Ip(i,j,1) = sp(i,j,1) - (sigma/TWO)*(sp(i,j,1)-sm(i,j,1)-(ONE-TWO3RD*sigma)*s6)
             Im(i,j,1) = ZERO
          end if
          if (u(i,j,1) .lt. ZERO) then
             Ip(i,j,1) = ZERO
             Im(i,j,1) = sm(i,j,1) + (sigma/TWO)*(sp(i,j,1)-sm(i,j,1)+(ONE-TWO3RD*sigma)*s6)
          end if
       end do
    end do

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

    ! modify sp and sm using quadratic limiters
    ! modify sp and sm using quadratic limiters
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          if ((sp(i,j,2)-s(i,j))*(s(i,j)-sm(i,j,2)) .le. ZERO) then
             sp(i,j,2) = s(i,j)
             sm(i,j,2) = s(i,j)
          end if
          if (abs(sp(i,j,2)-s(i,j)) .ge. TWO*abs(sm(i,j,2)-s(i,j))) then
             sp(i,j,2) = THREE*s(i,j) - TWO*sm(i,j,2)
          end if
          if (abs(sm(i,j,2)-s(i,j)) .ge. TWO*abs(sp(i,j,2)-s(i,j))) then
             sm(i,j,2) = THREE*s(i,j) - TWO*sp(i,j,2)
          end if
       end do
    end do

    ! compute Ip and Im
    do j=lo(2)-1,hi(2)+1
       if (j .le. 0) then
          v_from_w0 = w0(0)
       else if (j .ge. nr(n)) then
          v_from_w0 = w0(nr(n))
       else
          v_from_w0 = HALF*(w0(j)+w0(j+1))
       end if
       do i=lo(1)-1,hi(1)+1
          sigma = abs(u(i,j,2)+v_from_w0)*dt/dx(2)
          s6 = SIX*s(i,j) - THREE*(sm(i,j,2)+sp(i,j,2))
          if (u(i,j,2)+v_from_w0 .gt. ZERO) then
             Ip(i,j,2) = sp(i,j,2) - (sigma/TWO)*(sp(i,j,2)-sm(i,j,2)-(ONE-TWO3RD*sigma)*s6)
             Im(i,j,2) = ZERO
          end if
          if (u(i,j,2)+v_from_w0 .lt. ZERO) then
             Ip(i,j,2) = ZERO
             Im(i,j,2) = sm(i,j,2) + (sigma/TWO)*(sp(i,j,2)-sm(i,j,2)+(ONE-TWO3RD*sigma)*s6)
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
    integer :: i,j

    real(kind=dp_t) :: dsl, dsr, dsc, sigmam, sigmap, s6, vlo, vhi

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

    ! modify sp and sm using quadratic limiters
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          if ((sp(i,j,1)-s(i,j))*(s(i,j)-sm(i,j,1)) .le. ZERO) then
             sp(i,j,1) = s(i,j)
             sm(i,j,1) = s(i,j)
          end if
          if (abs(sp(i,j,1)-s(i,j)) .ge. TWO*abs(sm(i,j,1)-s(i,j))) then
             sp(i,j,1) = THREE*s(i,j) - TWO*sm(i,j,1)
          end if
          if (abs(sm(i,j,1)-s(i,j)) .ge. TWO*abs(sp(i,j,1)-s(i,j))) then
             sm(i,j,1) = THREE*s(i,j) - TWO*sp(i,j,1)
          end if
       end do
    end do

    ! compute Ip and Im
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sigmam = abs(umac(i,j))*dt/dx(1)
          sigmap = abs(umac(i+1,j))*dt/dx(1)
          s6 = SIX*s(i,j) - THREE*(sm(i,j,1)+sp(i,j,1))
          if(umac(i+1,j) .gt. ZERO) then
             Ip(i,j,1) = sp(i,j,1) - (sigmap/TWO)*(sp(i,j,1)-sm(i,j,1)-(ONE-TWO3RD*sigmap)*s6)
          else
             Ip(i,j,1) = ZERO
          end if
          if(umac(i,j) .lt. ZERO) then
             Im(i,j,1) = sm(i,j,1) + (sigmam/TWO)*(sp(i,j,1)-sm(i,j,1)+(ONE-TWO3RD*sigmam)*s6)
          else
             Im(i,j,1) = ZERO
          end if
       end do
    end do

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

    ! modify sp and sm using quadratic limiters
    ! modify sp and sm using quadratic limiters
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          if ((sp(i,j,2)-s(i,j))*(s(i,j)-sm(i,j,2)) .le. ZERO) then
             sp(i,j,2) = s(i,j)
             sm(i,j,2) = s(i,j)
          end if
          if (abs(sp(i,j,2)-s(i,j)) .ge. TWO*abs(sm(i,j,2)-s(i,j))) then
             sp(i,j,2) = THREE*s(i,j) - TWO*sm(i,j,2)
          end if
          if (abs(sm(i,j,2)-s(i,j)) .ge. TWO*abs(sp(i,j,2)-s(i,j))) then
             sm(i,j,2) = THREE*s(i,j) - TWO*sp(i,j,2)
          end if
       end do
    end do

    ! compute Ip and Im
    do j=lo(2)-1,hi(2)+1
       ! compute effect of w0
       if (j .lt. 0) then
          vlo = w0(0)
          vhi = w0(0)
       else if (j .gt. nr(n)-1) then
          vlo = w0(nr(n))
          vhi = w0(nr(n))
       else
          vlo = w0(j)
          vhi = w0(j+1)
       end if
       do i=lo(1)-1,hi(1)+1
          sigmap = abs(vmac(i,j+1)+vhi)*dt/dx(2)
          sigmam = abs(vmac(i,j  )+vlo)*dt/dx(2)
          s6 = SIX*s(i,j) - THREE*(sm(i,j,2)+sp(i,j,2))
          if(vmac(i,j+1)+vhi .gt. ZERO) then
             Ip(i,j,2) = sp(i,j,2) - (sigmap/TWO)*(sp(i,j,2)-sm(i,j,2)-(ONE-TWO3RD*sigmap)*s6)
          else
             Ip(i,j,2) = ZERO
          end if
          if(vmac(i,j)+vlo .lt. ZERO) then
             Im(i,j,2) = sm(i,j,2) + (sigmam/TWO)*(sp(i,j,2)-sm(i,j,2)+(ONE-TWO3RD*sigmam)*s6)
          else
             Im(i,j,2) = ZERO
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
