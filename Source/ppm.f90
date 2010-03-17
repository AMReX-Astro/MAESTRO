module ppm_module

  use bl_types
  use bl_error_module
  use probin_module, only: ppm_type
  use variables, only: rel_eps

  implicit none

  private

  public :: ppm_1d, ppm_fpu_1d, ppm_2d, ppm_fpu_2d, ppm_3d, ppm_fpu_3d

contains

  ! characteristics based on u
  subroutine ppm_1d(n,s,ng_s,u,ng_u,Ip,Im,w0,lo,hi,bc,dx,dt)

    use bc_module
    use bl_constants_module
    use geometry, only: nr

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_s,ng_u
    real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng_s:)
    real(kind=dp_t), intent(in   ) ::    u(lo(1)-ng_u:)
    real(kind=dp_t), intent(inout) ::   Ip(lo(1)-1   :)
    real(kind=dp_t), intent(inout) ::   Im(lo(1)-1   :)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i
    logical :: extremum, bigp, bigm

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
    real(kind=dp_t) :: sgn, sigma, s6, w0cc, amax, delam, delap
    real(kind=dp_t) :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

    ! s_{\ib,+}, s_{\ib,-}
    real(kind=dp_t), allocatable :: sp(:)
    real(kind=dp_t), allocatable :: sm(:)

    ! \delta s_{\ib}^{vL}
    real(kind=dp_t), allocatable :: dsvl(:)

    ! s_{i+\half}^{H.O.}
    real(kind=dp_t), allocatable :: sedge(:)

    ! cell-centered indexing
    allocate(sp(lo(1)-1:hi(1)+1))
    allocate(sm(lo(1)-1:hi(1)+1))

    ! constant used in Colella 2008
    C = 1.25d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2))

    ! edge-centered indexing for x-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+2))
    else
       allocate(sedge(lo(1)-2:hi(1)+3))
    end if

    ! compute s at x-edges
    if (ppm_type .eq. 1) then

       ! compute van Leer slopes in x-direction
       do i=lo(1)-2,hi(1)+2
          dsc = HALF * (s(i+1) - s(i-1))
          dsl = TWO  * (s(i  ) - s(i-1))
          dsr = TWO  * (s(i+1) - s(i  ))
          dsvl(i) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do
       
       ! interpolate s to x-edges
       do i=lo(1)-1,hi(1)+2
          sedge(i) = HALF*(s(i)+s(i-1)) - SIXTH*(dsvl(i)-dsvl(i-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i) = max(sedge(i),min(s(i),s(i-1)))
          sedge(i) = min(sedge(i),max(s(i),s(i-1)))
       end do

       ! copy sedge into sp and sm
       do i=lo(1)-1,hi(1)+1
          sp(i) = sedge(i+1)
          sm(i) = sedge(i  )
       end do

       ! modify using quadratic limiters
       do i=lo(1)-1,hi(1)+1
          if ((sp(i)-s(i))*(s(i)-sm(i)) .le. ZERO) then
             sp(i) = s(i)
             sm(i) = s(i)
          else if (abs(sp(i)-s(i)) .ge. TWO*abs(sm(i)-s(i))) then
             sp(i) = THREE*s(i) - TWO*sm(i)
          else if (abs(sm(i)-s(i)) .ge. TWO*abs(sp(i)-s(i))) then
             sm(i) = THREE*s(i) - TWO*sp(i)
          end if
       end do
       
       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)) = s(lo(1)-1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)+1) = &
               -FIFTH        *s(lo(1)-1) &
               + (THREE/FOUR)*s(lo(1)  ) &
               + HALF        *s(lo(1)+1) &
               - (ONE/20.0d0)*s(lo(1)+2)

          ! make sure sedge lies in between adjacent cell-centered values
          sedge(lo(1)+1) = max(sedge(lo(1)+1),min(s(lo(1)+1),s(lo(1))))
          sedge(lo(1)+1) = min(sedge(lo(1)+1),max(s(lo(1)+1),s(lo(1))))

          ! copy sedge into sp and sm
          sp(lo(1)  ) = sedge(lo(1)+1)
          sm(lo(1)+1) = sedge(lo(1)+1)

          ! reset sp on second interior edge
          sp(lo(1)+1) = sedge(lo(1)+2)

          ! modify using quadratic limiters
          i = lo(1)+1
          if ((sp(i)-s(i))*(s(i)-sm(i)) .le. ZERO) then
             sp(i) = s(i)
             sm(i) = s(i)
          else if (abs(sp(i)-s(i)) .ge. TWO*abs(sm(i)-s(i))) then
             sp(i) = THREE*s(i) - TWO*sm(i)
          else if (abs(sm(i)-s(i)) .ge. TWO*abs(sp(i)-s(i))) then
             sm(i) = THREE*s(i) - TWO*sp(i)
          end if
       end if
       
       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1)) = s(hi(1)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(hi(1)) = &
               -FIFTH        *s(hi(1)+1) &
               + (THREE/FOUR)*s(hi(1)  ) &
               + HALF        *s(hi(1)-1) &
               - (ONE/20.0d0)*s(hi(1)-2)

          ! make sure sedge lies in between adjacent cell-centered values
          sedge(hi(1)) = max(sedge(hi(1)),min(s(hi(1)-1),s(hi(1))))
          sedge(hi(1)) = min(sedge(hi(1)),max(s(hi(1)-1),s(hi(1))))

          ! copy sedge into sp and sm
          sp(hi(1)-1) = sedge(hi(1))
          sm(hi(1)  ) = sedge(hi(1))

          ! reset sm on second interior edge
          sm(hi(1)-1) = sedge(hi(1)-1)

          ! modify using quadratic limiters
          i = hi(1)-1
          if ((sp(i)-s(i))*(s(i)-sm(i)) .le. ZERO) then
             sp(i) = s(i)
             sm(i) = s(i)
          else if (abs(sp(i)-s(i)) .ge. TWO*abs(sm(i)-s(i))) then
             sp(i) = THREE*s(i) - TWO*sm(i)
          else if (abs(sm(i)-s(i)) .ge. TWO*abs(sp(i)-s(i))) then
             sm(i) = THREE*s(i) - TWO*sp(i)
          end if
       end if

    else if (ppm_type .eq. 2) then
       
       if (ng_s .lt. 4) then
          call bl_error("Need 4 ghost cells for ppm_type=2")
       end if

       ! interpolate s to x-edges
       do i=lo(1)-2,hi(1)+3
          sedge(i) = (7.d0/12.d0)*(s(i-1)+s(i)) - (1.d0/12.d0)*(s(i-2)+s(i+1))
          ! limit sedge
          if ((sedge(i)-s(i-1))*(s(i)-sedge(i)) .lt. ZERO) then
             D2  = THREE*(s(i-1)-TWO*sedge(i)+s(i))
             D2L = s(i-2)-TWO*s(i-1)+s(i)
             D2R = s(i-1)-TWO*s(i)+s(i+1)
             sgn = sign(ONE,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
             sedge(i) = HALF*(s(i-1)+s(i)) - SIXTH*D2LIM
          end if
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
       do i=lo(1)-1,hi(1)+1

          alphap = sedge(i+1)-s(i)
          alpham = sedge(i  )-s(i)
          bigp = abs(alphap).gt.TWO*abs(alpham)
          bigm = abs(alpham).gt.TWO*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. ZERO) then
             extremum = .true.
          else if (bigp .or. bigm) then
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             dafacem = sedge(i) - sedge(i-1)
             dafacep = sedge(i+2) - sedge(i+1)
             dabarm = s(i) - s(i-1)
             dabarp = s(i+1) - s(i)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin= min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2  = SIX*(alpham + alphap)
             D2L = s(i-2)-TWO*s(i-1)+s(i)
             D2R = s(i)-TWO*s(i+1)+s(i+2)
             D2C = s(i-1)-TWO*s(i)+s(i+1)
             sgn = sign(ONE,D2)
             D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn = sign(ONE,alpham)
                amax = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1) - s(i)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -TWO*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn = sign(ONE,alphap)
                amax = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1) - s(i)
               if (sgn*amax .ge. sgn*delap) then
                  if (sgn*(delap - alphap).ge.1.d-10) then
                     alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                  else
                     alpham = -TWO*alphap
                  endif
               endif
             end if
          end if

          sm(i) = s(i) + alpham
          sp(i) = s(i) + alphap

       end do

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1))    = s(lo(1)-1)
          sedge(lo(1)) = s(lo(1)-1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)+1) = &
               -FIFTH        *s(lo(1)-1) &
               + (THREE/FOUR)*s(lo(1)  ) &
               + HALF        *s(lo(1)+1) &
               - (ONE/20.0d0)*s(lo(1)+2)

          ! make sure sedge lies in between adjacent cell-centered values
          sedge(lo(1)+1) = max(sedge(lo(1)+1),min(s(lo(1)+1),s(lo(1))))
          sedge(lo(1)+1) = min(sedge(lo(1)+1),max(s(lo(1)+1),s(lo(1))))

          ! copy sedge into sp
          sp(lo(1)  ) = sedge(lo(1)+1)

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do i=lo(1)+1,lo(1)+2

             alphap = sedge(i+1)-s(i)
             alpham = sedge(i  )-s(i)
             bigp = abs(alphap).gt.TWO*abs(alpham)
             bigm = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                dafacem = sedge(i) - sedge(i-1)
                dafacep = sedge(i+2) - sedge(i+1)
                dabarm = s(i) - s(i-1)
                dabarp = s(i+1) - s(i)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin= min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. 0.d0)
             end if

             if (extremum) then
                D2  = SIX*(alpham + alphap)
                D2L = s(i-2)-TWO*s(i-1)+s(i)
                D2R = s(i)-TWO*s(i+1)+s(i+2)
                D2C = s(i-1)-TWO*s(i)+s(i+1)
                sgn = sign(ONE,D2)
                D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn = sign(ONE,alpham)
                   amax = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i-1) - s(i)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn = sign(ONE,alphap)
                   amax = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i+1) - s(i)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge.1.d-10) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm(i) = s(i) + alpham
             sp(i) = s(i) + alphap

          end do
       end if

       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1)     ) = s(hi(1)+1)
          sedge(hi(1)+1) = s(hi(1)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(hi(1)) = &
               -FIFTH        *s(hi(1)+1) &
               + (THREE/FOUR)*s(hi(1)  ) &
               + HALF        *s(hi(1)-1) &
               - (ONE/20.0d0)*s(hi(1)-2)

          ! make sure sedge lies in between adjacent cell-centered values
          sedge(hi(1)) = max(sedge(hi(1)),min(s(hi(1)-1),s(hi(1))))
          sedge(hi(1)) = min(sedge(hi(1)),max(s(hi(1)-1),s(hi(1))))

          ! copy sedge into sm
          sm(hi(1)  ) = sedge(hi(1))

          ! reset sm on second interior edge
          sm(hi(1)-1) = sedge(hi(1)-1)

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do i=hi(1)-2,hi(1)-1

             alphap = sedge(i+1)-s(i)
             alpham = sedge(i  )-s(i)
             bigp = abs(alphap).gt.TWO*abs(alpham)
             bigm = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                dafacem = sedge(i) - sedge(i-1)
                dafacep = sedge(i+2) - sedge(i+1)
                dabarm = s(i) - s(i-1)
                dabarp = s(i+1) - s(i)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin= min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. 0.d0)
             end if

             if (extremum) then
                D2  = SIX*(alpham + alphap)
                D2L = s(i-2)-TWO*s(i-1)+s(i)
                D2R = s(i)-TWO*s(i+1)+s(i+2)
                D2C = s(i-1)-TWO*s(i)+s(i+1)
                sgn = sign(ONE,D2)
                D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn = sign(ONE,alpham)
                   amax = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i-1) - s(i)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn = sign(ONE,alphap)
                   amax = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i+1) - s(i)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge.1.d-10) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm(i) = s(i) + alpham
             sp(i) = s(i) + alphap

          end do
       end if

    end if

    ! compute x-component of Ip and Im
    do i=lo(1)-1,hi(1)+1
       if (i .le. 0) then
          w0cc = w0(0)
       else if (i .ge. nr(n)) then
          w0cc = w0(nr(n))
       else
          w0cc = HALF*(w0(i)+w0(I+1))
       end if
       sigma = abs(u(i)+w0cc)*dt/dx(1)
       s6 = SIX*s(i) - THREE*(sm(i)+sp(i))
       if (u(i)+w0cc .gt. rel_eps) then
          Ip(i) = sp(i) - (sigma/TWO)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma)*s6)
          Im(i) = s(i)
       else if (u(i)+w0cc .lt. -rel_eps) then
          Ip(i) = s(i)
          Im(i) = sm(i) + (sigma/TWO)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
       else
          Ip(i) = s(i)
          Im(i) = s(i)
       end if
    end do

    deallocate(sp,sm,dsvl,sedge)

  end subroutine ppm_1d

  ! characteristics based on u
  subroutine ppm_fpu_1d(n,s,ng_s,umac,ng_um,Ip,Im,w0,lo,hi,bc,dx,dt)

    use bc_module
    use bl_constants_module
    use geometry, only: nr

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_s,ng_um
    real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng_s:)
    real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_um:)
    real(kind=dp_t), intent(inout) ::   Ip(lo(1)-1   :)
    real(kind=dp_t), intent(inout) ::   Im(lo(1)-1   :)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i
    logical :: extremum, bigp, bigm

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
    real(kind=dp_t) :: sgn, sigmam, sigmap, s6, w0lo, w0hi, amax, delam, delap
    real(kind=dp_t) :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

    ! s_{\ib,+}, s_{\ib,-}
    real(kind=dp_t), allocatable :: sp(:)
    real(kind=dp_t), allocatable :: sm(:)

    ! \delta s_{\ib}^{vL}
    real(kind=dp_t), allocatable :: dsvl(:)

    ! s_{i+\half}^{H.O.}
    real(kind=dp_t), allocatable :: sedge(:)

    ! cell-centered indexing
    allocate(sp(lo(1)-1:hi(1)+1))
    allocate(sm(lo(1)-1:hi(1)+1))

    ! constant used in Colella 2008
    C = 1.25d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2))

    ! edge-centered indexing for x-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+2))
    else
       allocate(sedge(lo(1)-2:hi(1)+3))
    end if

    ! compute s at x-edges
    if (ppm_type .eq. 1) then

       ! compute van Leer slopes in x-direction
       do i=lo(1)-2,hi(1)+2
          dsc = HALF * (s(i+1) - s(i-1))
          dsl = TWO  * (s(i  ) - s(i-1))
          dsr = TWO  * (s(i+1) - s(i  ))
          dsvl(i) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do
       
       ! interpolate s to x-edges
       do i=lo(1)-1,hi(1)+2
          sedge(i) = HALF*(s(i)+s(i-1)) - SIXTH*(dsvl(i)-dsvl(i-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i) = max(sedge(i),min(s(i),s(i-1)))
          sedge(i) = min(sedge(i),max(s(i),s(i-1)))
       end do

       ! copy sedge into sp and sm
       do i=lo(1)-1,hi(1)+1
          sp(i) = sedge(i+1)
          sm(i) = sedge(i  )
       end do

       ! modify using quadratic limiters
       do i=lo(1)-1,hi(1)+1
          if ((sp(i)-s(i))*(s(i)-sm(i)) .le. ZERO) then
             sp(i) = s(i)
             sm(i) = s(i)
          else if (abs(sp(i)-s(i)) .ge. TWO*abs(sm(i)-s(i))) then
             sp(i) = THREE*s(i) - TWO*sm(i)
          else if (abs(sm(i)-s(i)) .ge. TWO*abs(sp(i)-s(i))) then
             sm(i) = THREE*s(i) - TWO*sp(i)
          end if
       end do
       
       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)) = s(lo(1)-1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)+1) = &
               -FIFTH        *s(lo(1)-1) &
               + (THREE/FOUR)*s(lo(1)  ) &
               + HALF        *s(lo(1)+1) &
               - (ONE/20.0d0)*s(lo(1)+2)

          ! make sure sedge lies in between adjacent cell-centered values
          sedge(lo(1)+1) = max(sedge(lo(1)+1),min(s(lo(1)+1),s(lo(1))))
          sedge(lo(1)+1) = min(sedge(lo(1)+1),max(s(lo(1)+1),s(lo(1))))

          ! copy sedge into sp and sm
          sp(lo(1)  ) = sedge(lo(1)+1)
          sm(lo(1)+1) = sedge(lo(1)+1)

          ! reset sp on second interior edge
          sp(lo(1)+1) = sedge(lo(1)+2)

          ! modify using quadratic limiters
          i = lo(1)+1
          if ((sp(i)-s(i))*(s(i)-sm(i)) .le. ZERO) then
             sp(i) = s(i)
             sm(i) = s(i)
          else if (abs(sp(i)-s(i)) .ge. TWO*abs(sm(i)-s(i))) then
             sp(i) = THREE*s(i) - TWO*sm(i)
          else if (abs(sm(i)-s(i)) .ge. TWO*abs(sp(i)-s(i))) then
             sm(i) = THREE*s(i) - TWO*sp(i)
          end if
       end if
       
       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1)) = s(hi(1)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(hi(1)) = &
               -FIFTH        *s(hi(1)+1) &
               + (THREE/FOUR)*s(hi(1)  ) &
               + HALF        *s(hi(1)-1) &
               - (ONE/20.0d0)*s(hi(1)-2)

          ! make sure sedge lies in between adjacent cell-centered values
          sedge(hi(1)) = max(sedge(hi(1)),min(s(hi(1)-1),s(hi(1))))
          sedge(hi(1)) = min(sedge(hi(1)),max(s(hi(1)-1),s(hi(1))))

          ! copy sedge into sp and sm
          sp(hi(1)-1) = sedge(hi(1))
          sm(hi(1)  ) = sedge(hi(1))

          ! reset sm on second interior edge
          sm(hi(1)-1) = sedge(hi(1)-1)

          ! modify using quadratic limiters
          i = hi(1)-1
          if ((sp(i)-s(i))*(s(i)-sm(i)) .le. ZERO) then
             sp(i) = s(i)
             sm(i) = s(i)
          else if (abs(sp(i)-s(i)) .ge. TWO*abs(sm(i)-s(i))) then
             sp(i) = THREE*s(i) - TWO*sm(i)
          else if (abs(sm(i)-s(i)) .ge. TWO*abs(sp(i)-s(i))) then
             sm(i) = THREE*s(i) - TWO*sp(i)
          end if
       end if

    else if (ppm_type .eq. 2) then
       
       if (ng_s .lt. 4) then
          call bl_error("Need 4 ghost cells for ppm_type=2")
       end if

       ! interpolate s to x-edges
       do i=lo(1)-2,hi(1)+3
          sedge(i) = (7.d0/12.d0)*(s(i-1)+s(i)) - (1.d0/12.d0)*(s(i-2)+s(i+1))
          ! limit sedge
          if ((sedge(i)-s(i-1))*(s(i)-sedge(i)) .lt. ZERO) then
             D2  = THREE*(s(i-1)-TWO*sedge(i)+s(i))
             D2L = s(i-2)-TWO*s(i-1)+s(i)
             D2R = s(i-1)-TWO*s(i)+s(i+1)
             sgn = sign(ONE,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
             sedge(i) = HALF*(s(i-1)+s(i)) - SIXTH*D2LIM
          end if
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
       do i=lo(1)-1,hi(1)+1

          alphap = sedge(i+1)-s(i)
          alpham = sedge(i  )-s(i)
          bigp = abs(alphap).gt.TWO*abs(alpham)
          bigm = abs(alpham).gt.TWO*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. ZERO) then
             extremum = .true.
          else if (bigp .or. bigm) then
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             dafacem = sedge(i) - sedge(i-1)
             dafacep = sedge(i+2) - sedge(i+1)
             dabarm = s(i) - s(i-1)
             dabarp = s(i+1) - s(i)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin= min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2  = SIX*(alpham + alphap)
             D2L = s(i-2)-TWO*s(i-1)+s(i)
             D2R = s(i)-TWO*s(i+1)+s(i+2)
             D2C = s(i-1)-TWO*s(i)+s(i+1)
             sgn = sign(ONE,D2)
             D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn = sign(ONE,alpham)
                amax = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1) - s(i)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -TWO*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn = sign(ONE,alphap)
                amax = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1) - s(i)
               if (sgn*amax .ge. sgn*delap) then
                  if (sgn*(delap - alphap).ge.1.d-10) then
                     alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                  else
                     alpham = -TWO*alphap
                  endif
               endif
             end if
          end if

          sm(i) = s(i) + alpham
          sp(i) = s(i) + alphap

       end do

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1))    = s(lo(1)-1)
          sedge(lo(1)) = s(lo(1)-1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)+1) = &
               -FIFTH        *s(lo(1)-1) &
               + (THREE/FOUR)*s(lo(1)  ) &
               + HALF        *s(lo(1)+1) &
               - (ONE/20.0d0)*s(lo(1)+2)

          ! make sure sedge lies in between adjacent cell-centered values
          sedge(lo(1)+1) = max(sedge(lo(1)+1),min(s(lo(1)+1),s(lo(1))))
          sedge(lo(1)+1) = min(sedge(lo(1)+1),max(s(lo(1)+1),s(lo(1))))

          ! copy sedge into sp
          sp(lo(1)  ) = sedge(lo(1)+1)

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do i=lo(1)+1,lo(1)+2

             alphap = sedge(i+1)-s(i)
             alpham = sedge(i  )-s(i)
             bigp = abs(alphap).gt.TWO*abs(alpham)
             bigm = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                dafacem = sedge(i) - sedge(i-1)
                dafacep = sedge(i+2) - sedge(i+1)
                dabarm = s(i) - s(i-1)
                dabarp = s(i+1) - s(i)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin= min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. 0.d0)
             end if

             if (extremum) then
                D2  = SIX*(alpham + alphap)
                D2L = s(i-2)-TWO*s(i-1)+s(i)
                D2R = s(i)-TWO*s(i+1)+s(i+2)
                D2C = s(i-1)-TWO*s(i)+s(i+1)
                sgn = sign(ONE,D2)
                D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn = sign(ONE,alpham)
                   amax = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i-1) - s(i)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn = sign(ONE,alphap)
                   amax = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i+1) - s(i)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge.1.d-10) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm(i) = s(i) + alpham
             sp(i) = s(i) + alphap

          end do
       end if

       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1)     ) = s(hi(1)+1)
          sedge(hi(1)+1) = s(hi(1)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(hi(1)) = &
               -FIFTH        *s(hi(1)+1) &
               + (THREE/FOUR)*s(hi(1)  ) &
               + HALF        *s(hi(1)-1) &
               - (ONE/20.0d0)*s(hi(1)-2)

          ! make sure sedge lies in between adjacent cell-centered values
          sedge(hi(1)) = max(sedge(hi(1)),min(s(hi(1)-1),s(hi(1))))
          sedge(hi(1)) = min(sedge(hi(1)),max(s(hi(1)-1),s(hi(1))))

          ! copy sedge into sm
          sm(hi(1)  ) = sedge(hi(1))

          ! reset sm on second interior edge
          sm(hi(1)-1) = sedge(hi(1)-1)

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do i=hi(1)-2,hi(1)-1

             alphap = sedge(i+1)-s(i)
             alpham = sedge(i  )-s(i)
             bigp = abs(alphap).gt.TWO*abs(alpham)
             bigm = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                dafacem = sedge(i) - sedge(i-1)
                dafacep = sedge(i+2) - sedge(i+1)
                dabarm = s(i) - s(i-1)
                dabarp = s(i+1) - s(i)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin= min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. 0.d0)
             end if

             if (extremum) then
                D2  = SIX*(alpham + alphap)
                D2L = s(i-2)-TWO*s(i-1)+s(i)
                D2R = s(i)-TWO*s(i+1)+s(i+2)
                D2C = s(i-1)-TWO*s(i)+s(i+1)
                sgn = sign(ONE,D2)
                D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn = sign(ONE,alpham)
                   amax = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i-1) - s(i)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn = sign(ONE,alphap)
                   amax = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i+1) - s(i)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge.1.d-10) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm(i) = s(i) + alpham
             sp(i) = s(i) + alphap

          end do
       end if

    end if

    ! compute x-component of Ip and Im
    do i=lo(1)-1,hi(1)+1
       if (i .lt. 0) then
          w0lo = w0(0)
          w0hi = w0(0)
       else if (i .gt. nr(n)-1) then
          w0lo = w0(nr(n))
          w0hi = w0(nr(n))
       else
          w0lo = w0(i)
          w0hi = w0(i+1)
       end if
       sigmap = abs(umac(i+1)+w0hi)*dt/dx(1)
       sigmam = abs(umac(i  )+w0lo)*dt/dx(1)
       s6 = SIX*s(i) - THREE*(sm(i)+sp(i))
       if (umac(i+1)+w0hi .gt. rel_eps) then
          Ip(i) = sp(i) - (sigmap/TWO)*(sp(i)-sm(i)-(ONE-TWO3RD*sigmap)*s6)
       else 
          Ip(i) = s(i)
       end if
       if (umac(i)+w0lo .lt. -rel_eps) then
          Im(i) = sm(i) + (sigmam/TWO)*(sp(i)-sm(i)+(ONE-TWO3RD*sigmam)*s6)
       else
          Im(i) = s(i)
       end if
    end do

    deallocate(sm,sp,dsvl,sedge)

  end subroutine ppm_fpu_1d

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
    integer :: i,j

    logical :: extremum, bigp, bigm

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
    real(kind=dp_t) :: sgn, sigma, s6, w0cc, amax, delam, delap
    real(kind=dp_t) :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2,lo(2)-1:hi(2)+1))

    ! edge-centered indexing for x-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1))
    else
       allocate(sedge(lo(1)-2:hi(1)+3,lo(2)-1:hi(2)+1))
    end if

    ! compute s at x-edges
    if (ppm_type .eq. 1) then

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

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j) = sedge(i+1,j)
             sm(i,j) = sedge(i  ,j)
          end do
       end do

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
       
       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)+1,lo(2)-1:hi(2)+1) = &
               -FIFTH        *s(lo(1)-1,lo(2)-1:hi(2)+1) &
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

          ! reset sp on second interior edge
          do j=lo(2)-1,hi(2)+1
             sp(lo(1)+1,j) = sedge(lo(1)+2,j)
          end do

          ! modify using quadratic limiters
          do j=lo(2)-1,hi(2)+1
             i = lo(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end if
       
       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1),lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(hi(1),lo(2)-1:hi(2)+1) = &
               -FIFTH        *s(hi(1)+1,lo(2)-1:hi(2)+1) &
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

          ! reset sm on second interior edge
          do j=lo(2)-1,hi(2)+1
             sm(hi(1)-1,j) = sedge(hi(1)-1,j)
          end do

          ! modify using quadratic limiters
          do j=lo(2)-1,hi(2)+1
             i = hi(1)-1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end if

    else if (ppm_type .eq. 2) then
       
       if (ng_s .lt. 4) then
          call bl_error("Need 4 ghost cells for ppm_type=2")
       end if

       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+3
             sedge(i,j) = (7.d0/12.d0)*(s(i-1,j)+s(i,j)) - (1.d0/12.d0)*(s(i-2,j)+s(i+1,j))
             ! limit sedge
             if ((sedge(i,j)-s(i-1,j))*(s(i,j)-sedge(i,j)) .lt. ZERO) then
                D2  = THREE*(s(i-1,j)-TWO*sedge(i,j)+s(i,j))
                D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                D2R = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedge(i,j) = HALF*(s(i-1,j)+s(i,j)) - SIXTH*D2LIM
             end if
          end do
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

             alphap = sedge(i+1,j)-s(i,j)
             alpham = sedge(i  ,j)-s(i,j)
             bigp = abs(alphap).gt.TWO*abs(alpham)
             bigm = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                dafacem = sedge(i,j) - sedge(i-1,j)
                dafacep = sedge(i+2,j) - sedge(i+1,j)
                dabarm = s(i,j) - s(i-1,j)
                dabarp = s(i+1,j) - s(i,j)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin= min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. 0.d0)
             end if

             if (extremum) then
                D2  = SIX*(alpham + alphap)
                D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                D2R = s(i,j)-TWO*s(i+1,j)+s(i+2,j)
                D2C = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                sgn = sign(ONE,D2)
                D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn = sign(ONE,alpham)
                   amax = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i-1,j) - s(i,j)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn = sign(ONE,alphap)
                   amax = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i+1,j) - s(i,j)
                  if (sgn*amax .ge. sgn*delap) then
                     if (sgn*(delap - alphap).ge.1.d-10) then
                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                     else
                        alpham = -TWO*alphap
                     endif
                  endif
                end if
             end if

             sm(i,j) = s(i,j) + alpham
             sp(i,j) = s(i,j) + alphap

          end do
       end do

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1),lo(2)-1:hi(2)+1)    = s(lo(1)-1,lo(2)-1:hi(2)+1)
          sedge(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)+1,lo(2)-1:hi(2)+1) = &
               -FIFTH        *s(lo(1)-1,lo(2)-1:hi(2)+1) &
               + (THREE/FOUR)*s(lo(1)  ,lo(2)-1:hi(2)+1) &
               + HALF        *s(lo(1)+1,lo(2)-1:hi(2)+1) &
               - (ONE/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

          ! make sure sedge lies in between adjacent cell-centered values
          do j=lo(2)-1,hi(2)+1
             sedge(lo(1)+1,j) = max(sedge(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
             sedge(lo(1)+1,j) = min(sedge(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
          end do

          ! copy sedge into sp
          do j=lo(2)-1,hi(2)+1
             sp(lo(1)  ,j) = sedge(lo(1)+1,j)
          end do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)+1,lo(1)+2

                alphap = sedge(i+1,j)-s(i,j)
                alpham = sedge(i  ,j)-s(i,j)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j) - sedge(i-1,j)
                   dafacep = sedge(i+2,j) - sedge(i+1,j)
                   dabarm = s(i,j) - s(i-1,j)
                   dabarp = s(i+1,j) - s(i,j)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                   D2R = s(i,j)-TWO*s(i+1,j)+s(i+2,j)
                   D2C = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i-1,j) - s(i,j)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i+1,j) - s(i,j)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j) = s(i,j) + alpham
                sp(i,j) = s(i,j) + alphap

             end do
          end do
       end if

       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1),lo(2)-1:hi(2)+1)      = s(hi(1)+1,lo(2)-1:hi(2)+1)
          sedge(hi(1)+1,lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(hi(1),lo(2)-1:hi(2)+1) = &
               -FIFTH        *s(hi(1)+1,lo(2)-1:hi(2)+1) &
               + (THREE/FOUR)*s(hi(1)  ,lo(2)-1:hi(2)+1) &
               + HALF        *s(hi(1)-1,lo(2)-1:hi(2)+1) &
               - (ONE/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

          ! make sure sedge lies in between adjacent cell-centered values
          do j=lo(2)-1,hi(2)+1
             sedge(hi(1),j) = max(sedge(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
             sedge(hi(1),j) = min(sedge(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
          end do

          ! copy sedge into sm
          do j=lo(2)-1,hi(2)+1
             sm(hi(1)  ,j) = sedge(hi(1),j)
          end do

          ! reset sm on second interior edge
          do j=lo(2)-1,hi(2)+1
             sm(hi(1)-1,j) = sedge(hi(1)-1,j)
          end do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do j=lo(2)-1,hi(2)+1
             do i=hi(1)-2,hi(1)-1

                alphap = sedge(i+1,j)-s(i,j)
                alpham = sedge(i  ,j)-s(i,j)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j) - sedge(i-1,j)
                   dafacep = sedge(i+2,j) - sedge(i+1,j)
                   dabarm = s(i,j) - s(i-1,j)
                   dabarp = s(i+1,j) - s(i,j)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                   D2R = s(i,j)-TWO*s(i+1,j)+s(i+2,j)
                   D2C = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i-1,j) - s(i,j)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i+1,j) - s(i,j)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j) = s(i,j) + alpham
                sp(i,j) = s(i,j) + alphap

             end do
          end do
       end if

    end if

    ! compute x-component of Ip and Im
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sigma = abs(u(i,j,1))*dt/dx(1)
          s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
          if (u(i,j,1) .gt. rel_eps) then
             Ip(i,j,1) = sp(i,j) - (sigma/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
             Im(i,j,1) = s(i,j)
          else if (u(i,j,1) .lt. -rel_eps) then
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
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2))
    else
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+3))
    end if

    ! compute s at y-edges
    if (ppm_type .eq. 1) then

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

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j) = sedge(i,j+1)
             sm(i,j) = sedge(i,j  )
          end do
       end do

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
       
       ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
       if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)-1:hi(1)+1,lo(2)+1) = &
               -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1) &
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

          ! reset sp on second interior edge
          do i=lo(1)-1,hi(1)+1
             sp(i,lo(2)+1) = sedge(i,lo(2)+2)
          end do

          ! modify using quadratic limiters
          do i=lo(1)-1,hi(1)+1
             j = lo(2)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end if

       if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,hi(2)) = s(lo(1)-1:hi(1)+1,hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)-1:hi(1)+1,hi(2)) = &
               -FIFTH        *s(lo(1)-1:hi(1)+1,hi(2)+1) &
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

          ! reset sm on second interior edge
          do i=lo(1)-1,hi(1)+1
             sm(i,hi(2)-1) = sedge(i,hi(2)-1)
          end do

          ! modify using quadratic limiters
          do i=lo(1)-1,hi(1)+1
             j = hi(2)-1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end if

    else if (ppm_type .eq. 2) then
       
       ! interpolate s to y-edges
       do j=lo(2)-2,hi(2)+3
          do i=lo(1)-1,hi(1)+1
             sedge(i,j) = (7.d0/12.d0)*(s(i,j-1)+s(i,j)) - (1.d0/12.d0)*(s(i,j-2)+s(i,j+1))
             ! limit sedge
             if ((sedge(i,j)-s(i,j-1))*(s(i,j)-sedge(i,j)) .lt. ZERO) then
                D2  = THREE*(s(i,j-1)-TWO*sedge(i,j)+s(i,j))
                D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                D2R = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedge(i,j) = HALF*(s(i,j-1)+s(i,j)) - SIXTH*D2LIM
             end if
          end do
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

             alphap = sedge(i,j+1)-s(i,j)
             alpham = sedge(i,j  )-s(i,j)
             bigp = abs(alphap).gt.TWO*abs(alpham)
             bigm = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                dafacem = sedge(i,j) - sedge(i,j-1)
                dafacep = sedge(i,j+2) - sedge(i,j+1)
                dabarm = s(i,j) - s(i,j-1)
                dabarp = s(i,j+1) - s(i,j)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin= min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. 0.d0)
             end if

             if (extremum) then
                D2  = SIX*(alpham + alphap)
                D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                D2R = s(i,j)-TWO*s(i,j+1)+s(i,j+2)
                D2C = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                sgn = sign(ONE,D2)
                D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn = sign(ONE,alpham)
                   amax = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i,j-1) - s(i,j)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn = sign(ONE,alphap)
                   amax = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i,j+1) - s(i,j)
                  if (sgn*amax .ge. sgn*delap) then
                     if (sgn*(delap - alphap).ge.1.d-10) then
                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                     else
                        alpham = -TWO*alphap
                     endif
                  endif
                end if
             end if

             sm(i,j) = s(i,j) + alpham
             sp(i,j) = s(i,j) + alphap

          end do
       end do

       ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
       if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2))    = s(lo(1)-1:hi(1)+1,lo(2)-1)
          sedge(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)-1:hi(1)+1,lo(2)+1) = &
               -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1) &
               + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)  ) &
               + HALF        *s(lo(1)-1:hi(1)+1,lo(2)+1) &
               - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

          ! make sure sedge lies in between adjacent cell-centered values
          do i=lo(1)-1,hi(1)+1
             sedge(i,lo(2)+1) = max(sedge(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
             sedge(i,lo(2)+1) = min(sedge(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
          end do

          ! copy sedge into sp
          do i=lo(1)-1,hi(1)+1
             sp(i,lo(2)  ) = sedge(i,lo(2)+1)
          end do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do j=lo(2)+1,lo(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i,j+1)-s(i,j)
                alpham = sedge(i,j  )-s(i,j)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j) - sedge(i,j-1)
                   dafacep = sedge(i,j+2) - sedge(i,j+1)
                   dabarm = s(i,j) - s(i,j-1)
                   dabarp = s(i,j+1) - s(i,j)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                   D2R = s(i,j)-TWO*s(i,j+1)+s(i,j+2)
                   D2C = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i,j-1) - s(i,j)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i,j+1) - s(i,j)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j) = s(i,j) + alpham
                sp(i,j) = s(i,j) + alphap

             end do
          end do
       end if

       if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,hi(2))      = s(lo(1)-1:hi(1)+1,hi(2)+1)
          sedge(lo(1)-1:hi(1)+1,hi(2)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)-1:hi(1)+1,hi(2)) = &
               -FIFTH        *s(lo(1)-1:hi(1)+1,hi(2)+1) &
               + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,hi(2)  ) &
               + HALF        *s(lo(1)-1:hi(1)+1,hi(2)-1) &
               - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

          ! make sure sedge lies in between adjacent cell-centered values
          do i=lo(1)-1,hi(1)+1
             sedge(i,hi(2)) = max(sedge(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
             sedge(i,hi(2)) = min(sedge(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
          end do

          ! copy sedge into sm
          do i=lo(1)-1,hi(1)+1
             sm(i,hi(2)  ) = sedge(i,hi(2))
          end do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do j=hi(2)-2,hi(2)-1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i,j+1)-s(i,j)
                alpham = sedge(i,j  )-s(i,j)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j) - sedge(i,j-1)
                   dafacep = sedge(i,j+2) - sedge(i,j+1)
                   dabarm = s(i,j) - s(i,j-1)
                   dabarp = s(i,j+1) - s(i,j)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                   D2R = s(i,j)-TWO*s(i,j+1)+s(i,j+2)
                   D2C = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i,j-1) - s(i,j)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i,j+1) - s(i,j)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j) = s(i,j) + alpham
                sp(i,j) = s(i,j) + alphap

             end do
          end do
       end if

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
          if (u(i,j,2)+w0cc .gt. rel_eps) then
             Ip(i,j,2) = sp(i,j) - (sigma/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
             Im(i,j,2) = s(i,j)
          else if (u(i,j,2)+w0cc .lt. -rel_eps) then
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
    integer :: i,j

    logical :: extremum, bigp, bigm

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
    real(kind=dp_t) :: sgn, sigmam, sigmap, s6, w0lo, w0hi, amax, delam, delap
    real(kind=dp_t) :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2,lo(2)-1:hi(2)+1))

    ! edge-centered indexing for x-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1))
    else
       allocate(sedge(lo(1)-2:hi(1)+3,lo(2)-1:hi(2)+1))
    end if

    ! compute s at x-edges
    if (ppm_type .eq. 1) then

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

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j) = sedge(i+1,j)
             sm(i,j) = sedge(i  ,j)
          end do
       end do

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
       
       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)+1,lo(2)-1:hi(2)+1) = &
               -FIFTH        *s(lo(1)-1,lo(2)-1:hi(2)+1) &
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

          ! reset sp on second interior edge
          do j=lo(2)-1,hi(2)+1
             sp(lo(1)+1,j) = sedge(lo(1)+2,j)
          end do

          ! modify using quadratic limiters
          do j=lo(2)-1,hi(2)+1
             i = lo(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end if
       
       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1),lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(hi(1),lo(2)-1:hi(2)+1) = &
               -FIFTH        *s(hi(1)+1,lo(2)-1:hi(2)+1) &
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

          ! reset sm on second interior edge
          do j=lo(2)-1,hi(2)+1
             sm(hi(1)-1,j) = sedge(hi(1)-1,j)
          end do

          ! modify using quadratic limiters
          do j=lo(2)-1,hi(2)+1
             i = hi(1)-1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end if

    else if (ppm_type .eq. 2) then
       
       if (ng_s .lt. 4) then
          call bl_error("Need 4 ghost cells for ppm_type=2")
       end if

       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+3
             sedge(i,j) = (7.d0/12.d0)*(s(i-1,j)+s(i,j)) - (1.d0/12.d0)*(s(i-2,j)+s(i+1,j))
             ! limit sedge
             if ((sedge(i,j)-s(i-1,j))*(s(i,j)-sedge(i,j)) .lt. ZERO) then
                D2  = THREE*(s(i-1,j)-TWO*sedge(i,j)+s(i,j))
                D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                D2R = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedge(i,j) = HALF*(s(i-1,j)+s(i,j)) - SIXTH*D2LIM
             end if
          end do
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

             alphap = sedge(i+1,j)-s(i,j)
             alpham = sedge(i  ,j)-s(i,j)
             bigp = abs(alphap).gt.TWO*abs(alpham)
             bigm = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                dafacem = sedge(i,j) - sedge(i-1,j)
                dafacep = sedge(i+2,j) - sedge(i+1,j)
                dabarm = s(i,j) - s(i-1,j)
                dabarp = s(i+1,j) - s(i,j)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin= min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. 0.d0)
             end if

             if (extremum) then
                D2  = SIX*(alpham + alphap)
                D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                D2R = s(i,j)-TWO*s(i+1,j)+s(i+2,j)
                D2C = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                sgn = sign(ONE,D2)
                D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn = sign(ONE,alpham)
                   amax = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i-1,j) - s(i,j)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn = sign(ONE,alphap)
                   amax = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i+1,j) - s(i,j)
                  if (sgn*amax .ge. sgn*delap) then
                     if (sgn*(delap - alphap).ge.1.d-10) then
                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                     else
                        alpham = -TWO*alphap
                     endif
                  endif
                end if
             end if

             sm(i,j) = s(i,j) + alpham
             sp(i,j) = s(i,j) + alphap

          end do
       end do

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1),lo(2)-1:hi(2)+1)    = s(lo(1)-1,lo(2)-1:hi(2)+1)
          sedge(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)+1,lo(2)-1:hi(2)+1) = &
               -FIFTH        *s(lo(1)-1,lo(2)-1:hi(2)+1) &
               + (THREE/FOUR)*s(lo(1)  ,lo(2)-1:hi(2)+1) &
               + HALF        *s(lo(1)+1,lo(2)-1:hi(2)+1) &
               - (ONE/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

          ! make sure sedge lies in between adjacent cell-centered values
          do j=lo(2)-1,hi(2)+1
             sedge(lo(1)+1,j) = max(sedge(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
             sedge(lo(1)+1,j) = min(sedge(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
          end do

          ! copy sedge into sp
          do j=lo(2)-1,hi(2)+1
             sp(lo(1)  ,j) = sedge(lo(1)+1,j)
          end do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)+1,lo(1)+2

                alphap = sedge(i+1,j)-s(i,j)
                alpham = sedge(i  ,j)-s(i,j)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j) - sedge(i-1,j)
                   dafacep = sedge(i+2,j) - sedge(i+1,j)
                   dabarm = s(i,j) - s(i-1,j)
                   dabarp = s(i+1,j) - s(i,j)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                   D2R = s(i,j)-TWO*s(i+1,j)+s(i+2,j)
                   D2C = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i-1,j) - s(i,j)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i+1,j) - s(i,j)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j) = s(i,j) + alpham
                sp(i,j) = s(i,j) + alphap

             end do
          end do
       end if

       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1),lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)
          sedge(hi(1)+1,lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(hi(1),lo(2)-1:hi(2)+1) = &
               -FIFTH        *s(hi(1)+1,lo(2)-1:hi(2)+1) &
               + (THREE/FOUR)*s(hi(1)  ,lo(2)-1:hi(2)+1) &
               + HALF        *s(hi(1)-1,lo(2)-1:hi(2)+1) &
               - (ONE/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

          ! make sure sedge lies in between adjacent cell-centered values
          do j=lo(2)-1,hi(2)+1
             sedge(hi(1),j) = max(sedge(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
             sedge(hi(1),j) = min(sedge(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
          end do

          ! copy sedge into sm
          do j=lo(2)-1,hi(2)+1
             sm(hi(1)  ,j) = sedge(hi(1),j)
          end do

          ! reset sm on second interior edge
          do j=lo(2)-1,hi(2)+1
             sm(hi(1)-1,j) = sedge(hi(1)-1,j)
          end do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do j=lo(2)-1,hi(2)+1
             do i=hi(1)-2,hi(1)-1

                alphap = sedge(i+1,j)-s(i,j)
                alpham = sedge(i  ,j)-s(i,j)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j) - sedge(i-1,j)
                   dafacep = sedge(i+2,j) - sedge(i+1,j)
                   dabarm = s(i,j) - s(i-1,j)
                   dabarp = s(i+1,j) - s(i,j)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                   D2R = s(i,j)-TWO*s(i+1,j)+s(i+2,j)
                   D2C = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i-1,j) - s(i,j)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i+1,j) - s(i,j)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j) = s(i,j) + alpham
                sp(i,j) = s(i,j) + alphap

             end do
          end do
       end if

    end if

    ! compute x-component of Ip and Im
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          sigmap = abs(umac(i+1,j))*dt/dx(1)
          sigmam = abs(umac(i,j))*dt/dx(1)
          s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
          if (umac(i+1,j) .gt. rel_eps) then
             Ip(i,j,1) = sp(i,j) - (sigmap/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigmap)*s6)
          else
             Ip(i,j,1) = s(i,j)
          end if
          if (umac(i,j) .lt. -rel_eps) then
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
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2))
    else
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+3))
    end if

    ! compute s at y-edges
    if (ppm_type .eq. 1) then

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

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j) = sedge(i,j+1)
             sm(i,j) = sedge(i,j  )
          end do
       end do

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
       
       ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
       if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)-1:hi(1)+1,lo(2)+1) = &
               -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1) &
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

          ! reset sp on second interior edge
          do i=lo(1)-1,hi(1)+1
             sp(i,lo(2)+1) = sedge(i,lo(2)+2)
          end do

          ! modify using quadratic limiters
          do i=lo(1)-1,hi(1)+1
             j = lo(2)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end if

       if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,hi(2)) = s(lo(1)-1:hi(1)+1,hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)-1:hi(1)+1,hi(2)) = &
               -FIFTH        *s(lo(1)-1:hi(1)+1,hi(2)+1) &
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

          ! reset sm on second interior edge
          do i=lo(1)-1,hi(1)+1
             sm(i,hi(2)-1) = sedge(i,hi(2)-1)
          end do

          ! modify using quadratic limiters
          do i=lo(1)-1,hi(1)+1
             j = hi(2)-1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end if

    else if (ppm_type .eq. 2) then
       
       ! interpolate s to y-edges
       do j=lo(2)-2,hi(2)+3
          do i=lo(1)-1,hi(1)+1
             sedge(i,j) = (7.d0/12.d0)*(s(i,j-1)+s(i,j)) - (1.d0/12.d0)*(s(i,j-2)+s(i,j+1))
             ! limit sedge
             if ((sedge(i,j)-s(i,j-1))*(s(i,j)-sedge(i,j)) .lt. ZERO) then
                D2  = THREE*(s(i,j-1)-TWO*sedge(i,j)+s(i,j))
                D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                D2R = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedge(i,j) = HALF*(s(i,j-1)+s(i,j)) - SIXTH*D2LIM
             end if
          end do
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

             alphap = sedge(i,j+1)-s(i,j)
             alpham = sedge(i,j  )-s(i,j)
             bigp = abs(alphap).gt.TWO*abs(alpham)
             bigm = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                dafacem = sedge(i,j) - sedge(i,j-1)
                dafacep = sedge(i,j+2) - sedge(i,j+1)
                dabarm = s(i,j) - s(i,j-1)
                dabarp = s(i,j+1) - s(i,j)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin= min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. 0.d0)
             end if

             if (extremum) then
                D2  = SIX*(alpham + alphap)
                D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                D2R = s(i,j)-TWO*s(i,j+1)+s(i,j+2)
                D2C = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                sgn = sign(ONE,D2)
                D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn = sign(ONE,alpham)
                   amax = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i,j-1) - s(i,j)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn = sign(ONE,alphap)
                   amax = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i,j+1) - s(i,j)
                  if (sgn*amax .ge. sgn*delap) then
                     if (sgn*(delap - alphap).ge.1.d-10) then
                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                     else
                        alpham = -TWO*alphap
                     endif
                  endif
                end if
             end if

             sm(i,j) = s(i,j) + alpham
             sp(i,j) = s(i,j) + alphap

          end do
       end do

       ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
       if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)
          sedge(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)-1:hi(1)+1,lo(2)+1) = &
               -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1) &
               + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)  ) &
               + HALF        *s(lo(1)-1:hi(1)+1,lo(2)+1) &
               - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

          ! make sure sedge lies in between adjacent cell-centered values
          do i=lo(1)-1,hi(1)+1
             sedge(i,lo(2)+1) = max(sedge(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
             sedge(i,lo(2)+1) = min(sedge(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
          end do

          ! copy sedge into sp
          do i=lo(1)-1,hi(1)+1
             sp(i,lo(2)  ) = sedge(i,lo(2)+1)
          end do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do j=lo(2)+1,lo(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i,j+1)-s(i,j)
                alpham = sedge(i,j  )-s(i,j)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j) - sedge(i,j-1)
                   dafacep = sedge(i,j+2) - sedge(i,j+1)
                   dabarm = s(i,j) - s(i,j-1)
                   dabarp = s(i,j+1) - s(i,j)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                   D2R = s(i,j)-TWO*s(i,j+1)+s(i,j+2)
                   D2C = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i,j-1) - s(i,j)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i,j+1) - s(i,j)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j) = s(i,j) + alpham
                sp(i,j) = s(i,j) + alphap

             end do
          end do
       end if

       if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,hi(2)) = s(lo(1)-1:hi(1)+1,hi(2)+1)
          sedge(lo(1)-1:hi(1)+1,hi(2)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1)

          ! use a modified stencil to get sedge on the first interior edge
          sedge(lo(1)-1:hi(1)+1,hi(2)) = &
               -FIFTH        *s(lo(1)-1:hi(1)+1,hi(2)+1) &
               + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,hi(2)  ) &
               + HALF        *s(lo(1)-1:hi(1)+1,hi(2)-1) &
               - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

          ! make sure sedge lies in between adjacent cell-centered values
          do i=lo(1)-1,hi(1)+1
             sedge(i,hi(2)) = max(sedge(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
             sedge(i,hi(2)) = min(sedge(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
          end do

          ! copy sedge into sm
          do i=lo(1)-1,hi(1)+1
             sm(i,hi(2)  ) = sedge(i,hi(2))
          end do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
          do j=hi(2)-2,hi(2)-1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i,j+1)-s(i,j)
                alpham = sedge(i,j  )-s(i,j)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j) - sedge(i,j-1)
                   dafacep = sedge(i,j+2) - sedge(i,j+1)
                   dabarm = s(i,j) - s(i,j-1)
                   dabarp = s(i,j+1) - s(i,j)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                   D2R = s(i,j)-TWO*s(i,j+1)+s(i,j+2)
                   D2C = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i,j-1) - s(i,j)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i,j+1) - s(i,j)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j) = s(i,j) + alpham
                sp(i,j) = s(i,j) + alphap

             end do
          end do
       end if

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
          if (vmac(i,j+1)+w0hi .gt. rel_eps) then
             Ip(i,j,2) = sp(i,j) - (sigmap/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigmap)*s6)
          else
             Ip(i,j,2) = s(i,j)
          end if
          if (vmac(i,j)+w0lo .lt. -rel_eps) then
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
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)
    real(kind=dp_t), intent(inout) ::     Ip(lo(1)-1    :,lo(2)-1    :,lo(3)-1    :,:)
    real(kind=dp_t), intent(inout) ::     Im(lo(1)-1    :,lo(2)-1    :,lo(3)-1    :,:)
    real(kind=dp_t), intent(in   ) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i,j,k

    logical :: extremum, bigp, bigm

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
    real(kind=dp_t) :: sgn, sigma, s6, w0cc, velcc
    real(kind=dp_t) :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp
    real(kind=dp_t) :: amax, delam, delap

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    ! edge-centered indexing for x-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    else
       allocate(sedge(lo(1)-2:hi(1)+3,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    end if

    ! compute s at x-edges
    if (ppm_type .eq. 1) then
       
       ! compute van Leer slopes in x-direction
!$omp parallel do private(i,j,k,dsc,dsl,dsr)
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
!$omp end parallel do
       
       ! interpolate s to x-edges
!$omp parallel do private(i,j,k)
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
!$omp end parallel do

       ! copy sedge into sp and sm
!$omp parallel do private(i,j,k)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sp(i,j,k) = sedge(i+1,j,k)
                sm(i,j,k) = sedge(i  ,j,k)
             end do
          end do
       end do
!$omp end parallel do

       ! modify using quadratic limiters
!$omp parallel do private(i,j,k)
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
!$omp end parallel do
       
       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
               s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(lo(1)+1,j,k) = -FIFTH        *s(lo(1)-1,j,k) &
                                     + (THREE/FOUR)*s(lo(1)  ,j,k) &
                                     + HALF        *s(lo(1)+1,j,k) &
                                     - (ONE/20.0d0)*s(lo(1)+2,j,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(lo(1)+1,j,k) = max(sedge(lo(1)+1,j,k),min(s(lo(1)+1,j,k),s(lo(1),j,k)))
                sedge(lo(1)+1,j,k) = min(sedge(lo(1)+1,j,k),max(s(lo(1)+1,j,k),s(lo(1),j,k)))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                ! copy sedge into sp and sm
                sp(lo(1)  ,j,k) = sedge(lo(1)+1,j,k)
                sm(lo(1)+1,j,k) = sedge(lo(1)+1,j,k)
                ! reset sp on second interior edge
                sp(lo(1)+1,j,k) = sedge(lo(1)+2,j,k)
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          i = lo(1)+1
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
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
!$omp end parallel do
       end if

       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
               s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(hi(1),j,k) = -FIFTH        *s(hi(1)+1,j,k) &
                                   + (THREE/FOUR)*s(hi(1)  ,j,k) &
                                   + HALF        *s(hi(1)-1,j,k) &
                                   - (ONE/20.0d0)*s(hi(1)-2,j,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(hi(1),j,k) = max(sedge(hi(1),j,k),min(s(hi(1)-1,j,k),s(hi(1),j,k)))
                sedge(hi(1),j,k) = min(sedge(hi(1),j,k),max(s(hi(1)-1,j,k),s(hi(1),j,k)))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                ! copy sedge into sp and sm
                sp(hi(1)-1,j,k) = sedge(hi(1),j,k)
                sm(hi(1)  ,j,k) = sedge(hi(1),j,k)
                ! reset sm on second interior edge
                sm(hi(1)-1,j,k) = sedge(hi(1)-1,j,k)
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          i = hi(1)-1
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
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
!$omp end parallel do
       end if

    else if (ppm_type .eq. 2) then

       if (ng_s .lt. 4) then
          call bl_error("Need 4 ghost cells for ppm_type=2")
       end if

       ! interpolate s to x-edges
!$omp parallel do private(i,j,k,D2,D2L,D2R,sgn,D2LIM)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-2,hi(1)+3
                sedge(i,j,k) = (7.d0/12.d0)*(s(i-1,j,k)+s(i,j,k)) &
                     - (1.d0/12.d0)*(s(i-2,j,k)+s(i+1,j,k))
                ! limit sedge
                if ((sedge(i,j,k)-s(i-1,j,k))*(s(i,j,k)-sedge(i,j,k)) .lt. ZERO) then
                   D2  = THREE*(s(i-1,j,k)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                   D2R = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                   sgn = sign(ONE,D2)
                   D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                   sedge(i,j,k) = HALF*(s(i-1,j,k)+s(i,j,k)) - SIXTH*D2LIM
                end if
             end do
          end do
       end do
!$omp end parallel do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i+1,j,k)-s(i,j,k)
                alpham = sedge(i  ,j,k)-s(i,j,k)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j,k) - sedge(i-1,j,k)
                   dafacep = sedge(i+2,j,k) - sedge(i+1,j,k)
                   dabarm = s(i,j,k) - s(i-1,j,k)
                   dabarp = s(i+1,j,k) - s(i,j,k)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                   D2R = s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k)
                   D2C = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i-1,j,k) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i+1,j,k) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if
                
                sm(i,j,k) = s(i,j,k) + alpham
                sp(i,j,k) = s(i,j,k) + alphap

             end do
          end do
       end do
!$omp end parallel do

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)    = &
               s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
          sedge(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
               s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(lo(1)+1,j,k) = -FIFTH        *s(lo(1)-1,j,k) &
                                     + (THREE/FOUR)*s(lo(1)  ,j,k) &
                                     + HALF        *s(lo(1)+1,j,k) &
                                     - (ONE/20.0d0)*s(lo(1)+2,j,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(lo(1)+1,j,k) = max(sedge(lo(1)+1,j,k),min(s(lo(1)+1,j,k),s(lo(1),j,k)))
                sedge(lo(1)+1,j,k) = min(sedge(lo(1)+1,j,k),max(s(lo(1)+1,j,k),s(lo(1),j,k)))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sp
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sp(lo(1)  ,j,k) = sedge(lo(1)+1,j,k)
             end do
          end do
!$omp end parallel do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                do i=lo(1)+1,lo(1)+2

                   alphap = sedge(i+1,j,k)-s(i,j,k)
                   alpham = sedge(i  ,j,k)-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is 
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i-1,j,k)
                      dafacep = sedge(i+2,j,k) - sedge(i+1,j,k)
                      dabarm = s(i,j,k) - s(i-1,j,k)
                      dabarp = s(i+1,j,k) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k)
                      D2C = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i-1,j,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i+1,j,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
               s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(hi(1),j,k) = -FIFTH        *s(hi(1)+1,j,k) &
                                   + (THREE/FOUR)*s(hi(1)  ,j,k) &
                                   + HALF        *s(hi(1)-1,j,k) &
                                   - (ONE/20.0d0)*s(hi(1)-2,j,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(hi(1),j,k) = max(sedge(hi(1),j,k),min(s(hi(1)-1,j,k),s(hi(1),j,k)))
                sedge(hi(1),j,k) = min(sedge(hi(1),j,k),max(s(hi(1)-1,j,k),s(hi(1),j,k)))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sm
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sm(hi(1)  ,j,k) = sedge(hi(1),j,k)
             end do
          end do
!$omp end parallel do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                do i=hi(1)-2,hi(1)-1

                   alphap = sedge(i+1,j,k)-s(i,j,k)
                   alpham = sedge(i  ,j,k)-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is 
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i-1,j,k)
                      dafacep = sedge(i+2,j,k) - sedge(i+1,j,k)
                      dabarm = s(i,j,k) - s(i-1,j,k)
                      dabarp = s(i+1,j,k) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k)
                      D2C = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i-1,j,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i+1,j,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

    end if
    
    ! compute x-component of Ip and Im
!$omp parallel do private(i,j,k,velcc,sigma,s6)
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
             if (velcc .gt. rel_eps) then
                Ip(i,j,k,1) = sp(i,j,k) - &
                     (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,1) = s(i,j,k)
             else if (velcc .lt. -rel_eps) then
                Ip(i,j,k,1) = s(i,j,k)
                Im(i,j,k,1) = sm(i,j,k) + &
                     (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,1) = s(i,j,k)
                Im(i,j,k,1) = s(i,j,k)
             end if
          end do
       end do
    end do
!$omp end parallel do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+2,lo(3)-1:hi(3)+1))

    ! edge-centered indexing for y-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+1))
    else
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+3,lo(3)-1:hi(3)+1))
    end if

    ! compute s at y-edges
    if (ppm_type .eq. 1) then
       
       ! compute van Leer slopes in y-direction
!$omp parallel do private(i,j,k,dsc,dsl,dsr)
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
!$omp end parallel do

       ! interpolate s to y-edges
!$omp parallel do private(i,j,k)
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
!$omp end parallel do

       ! copy sedge into sp and sm
!$omp parallel do private(i,j,k)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sp(i,j,k) = sedge(i,j+1,k)
                sm(i,j,k) = sedge(i,j  ,k)
             end do
          end do
       end do
!$omp end parallel do

       ! modify using quadratic limiters
!$omp parallel do private(i,j,k)
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
!$omp end parallel do
       
       ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
       if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,lo(2)+1,k) = -FIFTH        *s(i,lo(2)-1,k) &
                                     + (THREE/FOUR)*s(i,lo(2)  ,k) &
                                     + HALF        *s(i,lo(2)+1,k) &
                                     - (ONE/20.0d0)*s(i,lo(2)+2,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,lo(2)+1,k) = max(sedge(i,lo(2)+1,k),min(s(i,lo(2)+1,k),s(i,lo(2),k)))
                sedge(i,lo(2)+1,k) = min(sedge(i,lo(2)+1,k),max(s(i,lo(2)+1,k),s(i,lo(2),k)))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                ! copy sedge into sp and sm
                sp(i,lo(2)  ,k) = sedge(i,lo(2)+1,k)
                sm(i,lo(2)+1,k) = sedge(i,lo(2)+1,k)
                ! reset sp on second interior edge
                sp(i,lo(2)+1,k) = sedge(i,lo(2)+2,k)
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          j = lo(2)+1
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
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
!$omp end parallel do
       end if

       if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,hi(2),k) = -FIFTH        *s(i,hi(2)+1,k) &
                                   + (THREE/FOUR)*s(i,hi(2)  ,k) &
                                   + HALF        *s(i,hi(2)-1,k) &
                                   - (ONE/20.0d0)*s(i,hi(2)-2,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,hi(2),k) = max(sedge(i,hi(2),k),min(s(i,hi(2)-1,k),s(i,hi(2),k)))
                sedge(i,hi(2),k) = min(sedge(i,hi(2),k),max(s(i,hi(2)-1,k),s(i,hi(2),k)))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                ! copy sedge into sp and sm
                sp(i,hi(2)-1,k) = sedge(i,hi(2),k)
                sm(i,hi(2)  ,k) = sedge(i,hi(2),k)
                ! reset sm on second interior edge
                sm(i,hi(2)-1,k) = sedge(i,hi(2)-1,k)
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          j = hi(2)-1
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
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
!$omp end parallel do
       end if

    else if (ppm_type .eq. 2) then

       ! interpolate s to y-edges
!$omp parallel do private(i,j,k,D2,D2L,D2R,sgn,D2LIM)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-2,hi(2)+3
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = (7.d0/12.d0)*(s(i,j-1,k)+s(i,j,k)) &
                     - (1.d0/12.d0)*(s(i,j-2,k)+s(i,j+1,k))
                ! limit sedge
                if ((sedge(i,j,k)-s(i,j-1,k))*(s(i,j,k)-sedge(i,j,k)) .lt. ZERO) then
                   D2  = THREE*(s(i,j-1,k)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                   D2R = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                   sgn = sign(ONE,D2)
                   D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                   sedge(i,j,k) = HALF*(s(i,j-1,k)+s(i,j,k)) - SIXTH*D2LIM
                end if
             end do
          end do
       end do
!$omp end parallel do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i,j+1,k)-s(i,j,k)
                alpham = sedge(i,j  ,k)-s(i,j,k)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j,k) - sedge(i,j-1,k)
                   dafacep = sedge(i,j+2,k) - sedge(i,j+1,k)
                   dabarm = s(i,j,k) - s(i,j-1,k)
                   dabarp = s(i,j+1,k) - s(i,j,k)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                   D2R = s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k)
                   D2C = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i,j-1,k) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i,j+1,k) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if
                
                sm(i,j,k) = s(i,j,k) + alpham
                sp(i,j,k) = s(i,j,k) + alphap

             end do
          end do
       end do
!$omp end parallel do

       ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
       if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)
          sedge(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,lo(2)+1,k) = -FIFTH        *s(i,lo(2)-1,k) &
                                     + (THREE/FOUR)*s(i,lo(2)  ,k) &
                                     + HALF        *s(i,lo(2)+1,k) &
                                     - (ONE/20.0d0)*s(i,lo(2)+2,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,lo(2)+1,k) = max(sedge(i,lo(2)+1,k),min(s(i,lo(2)+1,k),s(i,lo(2),k)))
                sedge(i,lo(2)+1,k) = min(sedge(i,lo(2)+1,k),max(s(i,lo(2)+1,k),s(i,lo(2),k)))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sp
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sp(i,lo(2)  ,k) = sedge(i,lo(2)+1,k)
             end do
          end do
!$omp end parallel do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)+1,lo(2)+2
                do i=lo(1)-1,hi(1)+1

                   alphap = sedge(i,j+1,k)-s(i,j,k)
                   alpham = sedge(i,j  ,k)-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i,j-1,k)
                      dafacep = sedge(i,j+2,k) - sedge(i,j+1,k)
                      dabarm = s(i,j,k) - s(i,j-1,k)
                      dabarp = s(i,j+1,k) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k)
                      D2C = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j-1,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j+1,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

       if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)
          sedge(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,hi(2),k) = -FIFTH        *s(i,hi(2)+1,k) &
                                   + (THREE/FOUR)*s(i,hi(2)  ,k) &
                                   + HALF        *s(i,hi(2)-1,k) &
                                   - (ONE/20.0d0)*s(i,hi(2)-2,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,hi(2),k) = max(sedge(i,hi(2),k),min(s(i,hi(2)-1,k),s(i,hi(2),k)))
                sedge(i,hi(2),k) = min(sedge(i,hi(2),k),max(s(i,hi(2)-1,k),s(i,hi(2),k)))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sm
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sm(i,hi(2)  ,k) = sedge(i,hi(2),k)
             end do
          end do
!$omp end parallel do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do k=lo(3)-1,hi(3)+1
             do j=hi(2)-2,hi(2)-1
                do i=lo(1)-1,hi(1)+1

                   alphap = sedge(i,j+1,k)-s(i,j,k)
                   alpham = sedge(i,j  ,k)-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i,j-1,k)
                      dafacep = sedge(i,j+2,k) - sedge(i,j+1,k)
                      dabarm = s(i,j,k) - s(i,j-1,k)
                      dabarp = s(i,j+1,k) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k)
                      D2C = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j-1,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j+1,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

    end if

    ! compute y-component of Ip and Im
!$omp parallel do private(i,j,k,velcc,sigma,s6)
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
             if (velcc .gt. rel_eps) then
                Ip(i,j,k,2) = sp(i,j,k) - &
                     (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,2) = s(i,j,k)
             else if (velcc .lt. -rel_eps) then
                Ip(i,j,k,2) = s(i,j,k)
                Im(i,j,k,2) = sm(i,j,k) + &
                     (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,2) = s(i,j,k)
                Im(i,j,k,2) = s(i,j,k)
             end if
          end do
       end do
    end do
!$omp end parallel do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra z-ghost cell
    allocate( dsvl(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-2:hi(3)+2))

    ! edge-centered indexing for z-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+2))
    else
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-2:hi(3)+3))
    end if

    ! compute s at z-edges
    if (ppm_type .eq. 1) then
       
       ! compute van Leer slopes in z-direction
!$omp parallel do private(i,j,k,dsc,dsl,dsr)
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
!$omp end parallel do

       ! interpolate s to z-edges
!$omp parallel do private(i,j,k)
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
!$omp end parallel do

       ! copy sedge into sp and sm
!$omp parallel do private(i,j,k)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sp(i,j,k) = sedge(i,j,k+1)
                sm(i,j,k) = sedge(i,j,k  )
             end do
          end do
       end do
!$omp end parallel do

       ! modify using quadratic limiters
!$omp parallel do private(i,j,k)
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
!$omp end parallel do
       
       ! different stencil needed for z-component of EXT_DIR and HOEXTRAP bc's
       if (bc(3,1) .eq. EXT_DIR  .or. bc(3,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,lo(3)+1) = -FIFTH        *s(i,j,lo(3)-1) &
                                     + (THREE/FOUR)*s(i,j,lo(3)  ) &
                                     + HALF        *s(i,j,lo(3)+1) &
                                     - (ONE/20.0d0)*s(i,j,lo(3)+2)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,lo(3)+1) = max(sedge(i,j,lo(3)+1),min(s(i,j,lo(3)+1),s(i,j,lo(3))))
                sedge(i,j,lo(3)+1) = min(sedge(i,j,lo(3)+1),max(s(i,j,lo(3)+1),s(i,j,lo(3))))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                ! copy sedge into sp and sm
                sp(i,j,lo(3)  ) = sedge(i,j,lo(3)+1)
                sm(i,j,lo(3)+1) = sedge(i,j,lo(3)+1)
                ! reset sp on second interior edge
                sp(i,j,lo(3)+1) = sedge(i,j,lo(3)+2)
             end do
          end do
!$omp end parallel do

!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          k = lo(3)+1
!$omp parallel do private(i,j)
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
!$omp end parallel do
       end if

       if (bc(3,2) .eq. EXT_DIR  .or. bc(3,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,hi(3)) = -FIFTH        *s(i,j,hi(3)+1) &
                                   + (THREE/FOUR)*s(i,j,hi(3)  ) &
                                   + HALF        *s(i,j,hi(3)-1) &
                                   - (ONE/20.0d0)*s(i,j,hi(3)-2)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,hi(3)) = max(sedge(i,j,hi(3)),min(s(i,j,hi(3)-1),s(i,j,hi(3))))
                sedge(i,j,hi(3)) = min(sedge(i,j,hi(3)),max(s(i,j,hi(3)-1),s(i,j,hi(3))))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                ! copy sedge into sp and sm
                sp(i,j,hi(3)-1) = sedge(i,j,hi(3))
                sm(i,j,hi(3)  ) = sedge(i,j,hi(3))
                ! reset sm on second interior edge
                sm(i,j,hi(3)-1) = sedge(i,j,hi(3)-1)
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          k = hi(3)-1
!$omp parallel do private(i,j)
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
!$omp end parallel do
       end if

    else if (ppm_type .eq. 2) then

       ! interpolate s to z-edges
!$omp parallel do private(i,j,k,D2,D2L,D2R,sgn,D2LIM)
       do k=lo(3)-2,hi(3)+3
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = (7.d0/12.d0)*(s(i,j,k-1)+s(i,j,k)) &
                     - (1.d0/12.d0)*(s(i,j,k-2)+s(i,j,k+1))
                ! limit sedge
                if ((sedge(i,j,k)-s(i,j,k-1))*(s(i,j,k)-sedge(i,j,k)) .lt. ZERO) then
                   D2  = THREE*(s(i,j,k-1)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                   D2R = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                   sgn = sign(ONE,D2)
                   D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                   sedge(i,j,k) = HALF*(s(i,j,k-1)+s(i,j,k)) - SIXTH*D2LIM
                end if
             end do
          end do
       end do
!$omp end parallel do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i,j,k+1)-s(i,j,k)
                alpham = sedge(i,j,k  )-s(i,j,k)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j,k) - sedge(i,j,k-1)
                   dafacep = sedge(i,j,k+2) - sedge(i,j,k+1)
                   dabarm = s(i,j,k) - s(i,j,k-1)
                   dabarp = s(i,j,k+1) - s(i,j,k)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                   D2R = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
                   D2C = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i,j,k-1) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i,j,k+1) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if
                
                sm(i,j,k) = s(i,j,k) + alpham
                sp(i,j,k) = s(i,j,k) + alphap

             end do
          end do
       end do
!$omp end parallel do

       ! different stencil needed for z-component of EXT_DIR and HOEXTRAP bc's
       if (bc(3,1) .eq. EXT_DIR  .or. bc(3,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)
          sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,lo(3)+1) = -FIFTH        *s(i,j,lo(3)-1) &
                                     + (THREE/FOUR)*s(i,j,lo(3)  ) &
                                     + HALF        *s(i,j,lo(3)+1) &
                                     - (ONE/20.0d0)*s(i,j,lo(3)+2)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,lo(3)+1) = max(sedge(i,j,lo(3)+1),min(s(i,j,lo(3)+1),s(i,j,lo(3))))
                sedge(i,j,lo(3)+1) = min(sedge(i,j,lo(3)+1),max(s(i,j,lo(3)+1),s(i,j,lo(3))))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sp
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sp(i,j,lo(3)  ) = sedge(i,j,lo(3)+1)
             end do
          end do
!$omp end parallel do
          
          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do j=lo(2)-1,hi(2)+1
             do k=lo(3)+1,lo(3)+2
                do i=lo(1)-1,hi(1)+1

                   alphap = sedge(i,j,k+1)-s(i,j,k)
                   alpham = sedge(i,j,k  )-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i,j,k-1)
                      dafacep = sedge(i,j,k+2) - sedge(i,j,k+1)
                      dabarm = s(i,j,k) - s(i,j,k-1)
                      dabarp = s(i,j,k+1) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
                      D2C = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j,k-1) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j,k+1) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

       if (bc(3,2) .eq. EXT_DIR  .or. bc(3,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)
          sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,hi(3)) = -FIFTH        *s(i,j,hi(3)+1) &
                                   + (THREE/FOUR)*s(i,j,hi(3)  ) &
                                   + HALF        *s(i,j,hi(3)-1) &
                                   - (ONE/20.0d0)*s(i,j,hi(3)-2)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,hi(3)) = max(sedge(i,j,hi(3)),min(s(i,j,hi(3)-1),s(i,j,hi(3))))
                sedge(i,j,hi(3)) = min(sedge(i,j,hi(3)),max(s(i,j,hi(3)-1),s(i,j,hi(3))))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sm
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sm(i,j,hi(3)  ) = sedge(i,j,hi(3))
             end do
          end do
!$omp end parallel do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do j=lo(2)-1,hi(2)+1
             do k=hi(3)-2,hi(3)-1
                do i=lo(1)-1,hi(1)+1

                   alphap = sedge(i,j,k+1)-s(i,j,k)
                   alpham = sedge(i,j,k  )-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i,j,k-1)
                      dafacep = sedge(i,j,k+2) - sedge(i,j,k+1)
                      dabarm = s(i,j,k) - s(i,j,k-1)
                      dabarp = s(i,j,k+1) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
                      D2C = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j,k-1) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j,k+1) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

    end if

    ! compute z-component of Ip and Im
!$omp parallel do private(i,j,k,w0cc,velcc,sigma,s6)
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
             if (velcc .gt. rel_eps) then
                Ip(i,j,k,3) = sp(i,j,k) - &
                     (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,3) = s(i,j,k)
             else if (velcc .lt. -rel_eps) then
                Ip(i,j,k,3) = s(i,j,k)
                Im(i,j,k,3) = sm(i,j,k) + &
                     (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,3) = s(i,j,k)
                Im(i,j,k,3) = s(i,j,k)
             end if
          end do
       end do
    end do
!$omp end parallel do

    deallocate(sp,sm,dsvl,sedge)

  end subroutine ppm_3d

  ! characteristics based on umac
  subroutine ppm_fpu_3d(n,s,ng_s,umac,vmac,wmac,ng_um,Ip,Im,w0,w0macx,w0macy,w0macz, &
                        ng_w0,lo,hi,bc,dx,dt)

    use bc_module
    use bl_constants_module
    use geometry, only: nr, spherical

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_s,ng_um,ng_w0
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(inout) ::     Ip(lo(1)-1    :,lo(2)-1    :,lo(3)-1    :,:)
    real(kind=dp_t), intent(inout) ::     Im(lo(1)-1    :,lo(2)-1    :,lo(3)-1    :,:)
    real(kind=dp_t), intent(in   ) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt

    ! local
    integer :: i,j,k

    logical :: extremum, bigp, bigm

    real(kind=dp_t) :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
    real(kind=dp_t) :: sgn, sigmam, sigmap, s6, w0lo, w0hi, vello, velhi
    real(kind=dp_t) :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp
    real(kind=dp_t) :: amax, delam, delap

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    ! edge-centered indexing for x-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    else
       allocate(sedge(lo(1)-2:hi(1)+3,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    end if

    ! compute s at x-edges
    if (ppm_type .eq. 1) then
       
       ! compute van Leer slopes in x-direction
!$omp parallel do private(i,j,k,dsc,dsl,dsr)
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
!$omp end parallel do
       
       ! interpolate s to x-edges
!$omp parallel do private(i,j,k)
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
!$omp end parallel do

       ! copy sedge into sp and sm
!$omp parallel do private(i,j,k)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sp(i,j,k) = sedge(i+1,j,k)
                sm(i,j,k) = sedge(i  ,j,k)
             end do
          end do
       end do
!$omp end parallel do

       ! modify using quadratic limiters
!$omp parallel do private(i,j,k)
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
!$omp end parallel do
       
       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
               s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(lo(1)+1,j,k) = -FIFTH        *s(lo(1)-1,j,k) &
                                     + (THREE/FOUR)*s(lo(1)  ,j,k) &
                                     + HALF        *s(lo(1)+1,j,k) &
                                     - (ONE/20.0d0)*s(lo(1)+2,j,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(lo(1)+1,j,k) = max(sedge(lo(1)+1,j,k),min(s(lo(1)+1,j,k),s(lo(1),j,k)))
                sedge(lo(1)+1,j,k) = min(sedge(lo(1)+1,j,k),max(s(lo(1)+1,j,k),s(lo(1),j,k)))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                ! copy sedge into sp and sm
                sp(lo(1)  ,j,k) = sedge(lo(1)+1,j,k)
                sm(lo(1)+1,j,k) = sedge(lo(1)+1,j,k)
                ! reset sp on second interior edge
                sp(lo(1)+1,j,k) = sedge(lo(1)+2,j,k)
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          i = lo(1)+1
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
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
!$omp end parallel do
       end if

       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
               s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(hi(1),j,k) = -FIFTH        *s(hi(1)+1,j,k) &
                                   + (THREE/FOUR)*s(hi(1)  ,j,k) &
                                   + HALF        *s(hi(1)-1,j,k) &
                                   - (ONE/20.0d0)*s(hi(1)-2,j,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(hi(1),j,k) = max(sedge(hi(1),j,k),min(s(hi(1)-1,j,k),s(hi(1),j,k)))
                sedge(hi(1),j,k) = min(sedge(hi(1),j,k),max(s(hi(1)-1,j,k),s(hi(1),j,k)))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                ! copy sedge into sp and sm
                sp(hi(1)-1,j,k) = sedge(hi(1),j,k)
                sm(hi(1)  ,j,k) = sedge(hi(1),j,k)
                ! reset sm on second interior edge
                sm(hi(1)-1,j,k) = sedge(hi(1)-1,j,k)
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          i = hi(1)-1
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
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
!$omp end parallel do
       end if

    else if (ppm_type .eq. 2) then

       if (ng_s .lt. 4) then
          call bl_error("Need 4 ghost cells for ppm_type=2")
       end if

       ! interpolate s to x-edges
!$omp parallel do private(i,j,k,D2,D2L,D2R,sgn,D2LIM)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-2,hi(1)+3
                sedge(i,j,k) = (7.d0/12.d0)*(s(i-1,j,k)+s(i,j,k)) &
                     - (1.d0/12.d0)*(s(i-2,j,k)+s(i+1,j,k))
                ! limit sedge
                if ((sedge(i,j,k)-s(i-1,j,k))*(s(i,j,k)-sedge(i,j,k)) .lt. ZERO) then
                   D2  = THREE*(s(i-1,j,k)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                   D2R = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                   sgn = sign(ONE,D2)
                   D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                   sedge(i,j,k) = HALF*(s(i-1,j,k)+s(i,j,k)) - SIXTH*D2LIM
                end if
             end do
          end do
       end do
!$omp end parallel do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i+1,j,k)-s(i,j,k)
                alpham = sedge(i  ,j,k)-s(i,j,k)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j,k) - sedge(i-1,j,k)
                   dafacep = sedge(i+2,j,k) - sedge(i+1,j,k)
                   dabarm = s(i,j,k) - s(i-1,j,k)
                   dabarp = s(i+1,j,k) - s(i,j,k)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                   D2R = s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k)
                   D2C = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i-1,j,k) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i+1,j,k) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if
                
                sm(i,j,k) = s(i,j,k) + alpham
                sp(i,j,k) = s(i,j,k) + alphap

             end do
          end do
       end do
!$omp end parallel do

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
       if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)    = &
               s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
          sedge(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
               s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(lo(1)+1,j,k) = -FIFTH        *s(lo(1)-1,j,k) &
                                     + (THREE/FOUR)*s(lo(1)  ,j,k) &
                                     + HALF        *s(lo(1)+1,j,k) &
                                     - (ONE/20.0d0)*s(lo(1)+2,j,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(lo(1)+1,j,k) = max(sedge(lo(1)+1,j,k),min(s(lo(1)+1,j,k),s(lo(1),j,k)))
                sedge(lo(1)+1,j,k) = min(sedge(lo(1)+1,j,k),max(s(lo(1)+1,j,k),s(lo(1),j,k)))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sp
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sp(lo(1)  ,j,k) = sedge(lo(1)+1,j,k)
             end do
          end do
!$omp end parallel do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                do i=lo(1)+1,lo(1)+2

                   alphap = sedge(i+1,j,k)-s(i,j,k)
                   alpham = sedge(i  ,j,k)-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is 
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i-1,j,k)
                      dafacep = sedge(i+2,j,k) - sedge(i+1,j,k)
                      dabarm = s(i,j,k) - s(i-1,j,k)
                      dabarp = s(i+1,j,k) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k)
                      D2C = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i-1,j,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i+1,j,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

       if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
               s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(hi(1),j,k) = -FIFTH        *s(hi(1)+1,j,k) &
                                   + (THREE/FOUR)*s(hi(1)  ,j,k) &
                                   + HALF        *s(hi(1)-1,j,k) &
                                   - (ONE/20.0d0)*s(hi(1)-2,j,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sedge(hi(1),j,k) = max(sedge(hi(1),j,k),min(s(hi(1)-1,j,k),s(hi(1),j,k)))
                sedge(hi(1),j,k) = min(sedge(hi(1),j,k),max(s(hi(1)-1,j,k),s(hi(1),j,k)))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sm
!$omp parallel do private(j,k)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                sm(hi(1)  ,j,k) = sedge(hi(1),j,k)
             end do
          end do
!$omp end parallel do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)-1,hi(2)+1
                do i=hi(1)-2,hi(1)-1

                   alphap = sedge(i+1,j,k)-s(i,j,k)
                   alpham = sedge(i  ,j,k)-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is 
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i-1,j,k)
                      dafacep = sedge(i+2,j,k) - sedge(i+1,j,k)
                      dabarm = s(i,j,k) - s(i-1,j,k)
                      dabarp = s(i+1,j,k) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k)
                      D2C = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i-1,j,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i+1,j,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

    end if
    
    ! compute x-component of Ip and Im
!$omp parallel do private(i,j,k,velhi,vello,sigmap,sigmam,s6)
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
             if (velhi .gt. rel_eps) then
                Ip(i,j,k,1) = sp(i,j,k) - (sigmap/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigmap)*s6)
             else
                Ip(i,j,k,1) = s(i,j,k)
             end if
             if (vello .lt. -rel_eps) then
                Im(i,j,k,1) = sm(i,j,k) + (sigmam/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigmam)*s6)
             else
                Im(i,j,k,1) = s(i,j,k)
             end if
          end do
       end do
    end do
!$omp end parallel do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+2,lo(3)-1:hi(3)+1))

    ! edge-centered indexing for y-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+1))
    else
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+3,lo(3)-1:hi(3)+1))
    end if

    ! compute s at y-edges
    if (ppm_type .eq. 1) then
       
       ! compute van Leer slopes in y-direction
!$omp parallel do private(i,j,k,dsc,dsl,dsr)
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
!$omp end parallel do

       ! interpolate s to y-edges
!$omp parallel do private(i,j,k)
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
!$omp end parallel do

       ! copy sedge into sp and sm
!$omp parallel do private(i,j,k)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sp(i,j,k) = sedge(i,j+1,k)
                sm(i,j,k) = sedge(i,j  ,k)
             end do
          end do
       end do
!$omp end parallel do

       ! modify using quadratic limiters
!$omp parallel do private(i,j,k)
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
!$omp end parallel do
       
       ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
       if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,lo(2)+1,k) = -FIFTH        *s(i,lo(2)-1,k) &
                                     + (THREE/FOUR)*s(i,lo(2)  ,k) &
                                     + HALF        *s(i,lo(2)+1,k) &
                                     - (ONE/20.0d0)*s(i,lo(2)+2,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,lo(2)+1,k) = max(sedge(i,lo(2)+1,k),min(s(i,lo(2)+1,k),s(i,lo(2),k)))
                sedge(i,lo(2)+1,k) = min(sedge(i,lo(2)+1,k),max(s(i,lo(2)+1,k),s(i,lo(2),k)))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                ! copy sedge into sp and sm
                sp(i,lo(2)  ,k) = sedge(i,lo(2)+1,k)
                sm(i,lo(2)+1,k) = sedge(i,lo(2)+1,k)
                ! reset sp on second interior edge
                sp(i,lo(2)+1,k) = sedge(i,lo(2)+2,k)
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          j = lo(2)+1
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
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
!$omp end parallel do
       end if

       if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,hi(2),k) = -FIFTH        *s(i,hi(2)+1,k) &
                                   + (THREE/FOUR)*s(i,hi(2)  ,k) &
                                   + HALF        *s(i,hi(2)-1,k) &
                                   - (ONE/20.0d0)*s(i,hi(2)-2,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,hi(2),k) = max(sedge(i,hi(2),k),min(s(i,hi(2)-1,k),s(i,hi(2),k)))
                sedge(i,hi(2),k) = min(sedge(i,hi(2),k),max(s(i,hi(2)-1,k),s(i,hi(2),k)))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                ! copy sedge into sp and sm
                sp(i,hi(2)-1,k) = sedge(i,hi(2),k)
                sm(i,hi(2)  ,k) = sedge(i,hi(2),k)
                ! reset sm on second interior edge
                sm(i,hi(2)-1,k) = sedge(i,hi(2)-1,k)
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          j = hi(2)-1
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
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
!$omp end parallel do
       end if

    else if (ppm_type .eq. 2) then

       ! interpolate s to y-edges
!$omp parallel do private(i,j,k,D2,D2L,D2R,sgn,D2LIM)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-2,hi(2)+3
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = (7.d0/12.d0)*(s(i,j-1,k)+s(i,j,k)) &
                     - (1.d0/12.d0)*(s(i,j-2,k)+s(i,j+1,k))
                ! limit sedge
                if ((sedge(i,j,k)-s(i,j-1,k))*(s(i,j,k)-sedge(i,j,k)) .lt. ZERO) then
                   D2  = THREE*(s(i,j-1,k)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                   D2R = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                   sgn = sign(ONE,D2)
                   D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                   sedge(i,j,k) = HALF*(s(i,j-1,k)+s(i,j,k)) - SIXTH*D2LIM
                end if
             end do
          end do
       end do
!$omp end parallel do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i,j+1,k)-s(i,j,k)
                alpham = sedge(i,j  ,k)-s(i,j,k)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j,k) - sedge(i,j-1,k)
                   dafacep = sedge(i,j+2,k) - sedge(i,j+1,k)
                   dabarm = s(i,j,k) - s(i,j-1,k)
                   dabarp = s(i,j+1,k) - s(i,j,k)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                   D2R = s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k)
                   D2C = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i,j-1,k) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i,j+1,k) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if
                
                sm(i,j,k) = s(i,j,k) + alpham
                sp(i,j,k) = s(i,j,k) + alphap

             end do
          end do
       end do
!$omp end parallel do

       ! different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
       if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)
          sedge(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,lo(2)+1,k) = -FIFTH        *s(i,lo(2)-1,k) &
                                     + (THREE/FOUR)*s(i,lo(2)  ,k) &
                                     + HALF        *s(i,lo(2)+1,k) &
                                     - (ONE/20.0d0)*s(i,lo(2)+2,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,lo(2)+1,k) = max(sedge(i,lo(2)+1,k),min(s(i,lo(2)+1,k),s(i,lo(2),k)))
                sedge(i,lo(2)+1,k) = min(sedge(i,lo(2)+1,k),max(s(i,lo(2)+1,k),s(i,lo(2),k)))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sp
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sp(i,lo(2)  ,k) = sedge(i,lo(2)+1,k)
             end do
          end do
!$omp end parallel do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do k=lo(3)-1,hi(3)+1
             do j=lo(2)+1,lo(2)+2
                do i=lo(1)-1,hi(1)+1

                   alphap = sedge(i,j+1,k)-s(i,j,k)
                   alpham = sedge(i,j  ,k)-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i,j-1,k)
                      dafacep = sedge(i,j+2,k) - sedge(i,j+1,k)
                      dabarm = s(i,j,k) - s(i,j-1,k)
                      dabarp = s(i,j+1,k) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k)
                      D2C = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j-1,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j+1,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

       if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)
          sedge(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,hi(2),k) = -FIFTH        *s(i,hi(2)+1,k) &
                                   + (THREE/FOUR)*s(i,hi(2)  ,k) &
                                   + HALF        *s(i,hi(2)-1,k) &
                                   - (ONE/20.0d0)*s(i,hi(2)-2,k)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,hi(2),k) = max(sedge(i,hi(2),k),min(s(i,hi(2)-1,k),s(i,hi(2),k)))
                sedge(i,hi(2),k) = min(sedge(i,hi(2),k),max(s(i,hi(2)-1,k),s(i,hi(2),k)))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sm
!$omp parallel do private(i,k)
          do k=lo(3)-1,hi(3)+1
             do i=lo(1)-1,hi(1)+1
                sm(i,hi(2)  ,k) = sedge(i,hi(2),k)
             end do
          end do
!$omp end parallel do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do k=lo(3)-1,hi(3)+1
             do j=hi(2)-2,hi(2)-1
                do i=lo(1)-1,hi(1)+1

                   alphap = sedge(i,j+1,k)-s(i,j,k)
                   alpham = sedge(i,j  ,k)-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i,j-1,k)
                      dafacep = sedge(i,j+2,k) - sedge(i,j+1,k)
                      dabarm = s(i,j,k) - s(i,j-1,k)
                      dabarp = s(i,j+1,k) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k)
                      D2C = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j-1,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j+1,k) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

    end if

    ! compute y-component of Ip and Im
!$omp parallel do private(i,j,k,velhi,vello,sigmap,sigmam,s6)
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
             if (velhi .gt. rel_eps) then
                Ip(i,j,k,2) = sp(i,j,k) - (sigmap/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigmap)*s6)
             else
                Ip(i,j,k,2) = s(i,j,k)
             end if
             if (vello .lt. -rel_eps) then
                Im(i,j,k,2) = sm(i,j,k) + (sigmam/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigmam)*s6)
             else
                Im(i,j,k,2) = s(i,j,k)
             end if
          end do
       end do
    end do
!$omp end parallel do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra z-ghost cell
    allocate( dsvl(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-2:hi(3)+2))

    ! edge-centered indexing for z-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+2))
    else
       allocate(sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-2:hi(3)+3))
    end if

    ! compute s at z-edges
    if (ppm_type .eq. 1) then
       
       ! compute van Leer slopes in z-direction
!$omp parallel do private(i,j,k,dsc,dsl,dsr)
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
!$omp end parallel do

       ! interpolate s to z-edges
!$omp parallel do private(i,j,k)
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
!$omp end parallel do

       ! copy sedge into sp and sm
!$omp parallel do private(i,j,k)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sp(i,j,k) = sedge(i,j,k+1)
                sm(i,j,k) = sedge(i,j,k  )
             end do
          end do
       end do
!$omp end parallel do

       ! modify using quadratic limiters
!$omp parallel do private(i,j,k)
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
!$omp end parallel do
       
       ! different stencil needed for z-component of EXT_DIR and HOEXTRAP bc's
       if (bc(3,1) .eq. EXT_DIR  .or. bc(3,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,lo(3)+1) = -FIFTH        *s(i,j,lo(3)-1) &
                                     + (THREE/FOUR)*s(i,j,lo(3)  ) &
                                     + HALF        *s(i,j,lo(3)+1) &
                                     - (ONE/20.0d0)*s(i,j,lo(3)+2)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,lo(3)+1) = max(sedge(i,j,lo(3)+1),min(s(i,j,lo(3)+1),s(i,j,lo(3))))
                sedge(i,j,lo(3)+1) = min(sedge(i,j,lo(3)+1),max(s(i,j,lo(3)+1),s(i,j,lo(3))))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                ! copy sedge into sp and sm
                sp(i,j,lo(3)  ) = sedge(i,j,lo(3)+1)
                sm(i,j,lo(3)+1) = sedge(i,j,lo(3)+1)
                ! reset sp on second interior edge
                sp(i,j,lo(3)+1) = sedge(i,j,lo(3)+2)
             end do
          end do
!$omp end parallel do

!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          k = lo(3)+1
!$omp parallel do private(i,j)
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
!$omp end parallel do
       end if

       if (bc(3,2) .eq. EXT_DIR  .or. bc(3,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,hi(3)) = -FIFTH        *s(i,j,hi(3)+1) &
                                   + (THREE/FOUR)*s(i,j,hi(3)  ) &
                                   + HALF        *s(i,j,hi(3)-1) &
                                   - (ONE/20.0d0)*s(i,j,hi(3)-2)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,hi(3)) = max(sedge(i,j,hi(3)),min(s(i,j,hi(3)-1),s(i,j,hi(3))))
                sedge(i,j,hi(3)) = min(sedge(i,j,hi(3)),max(s(i,j,hi(3)-1),s(i,j,hi(3))))
             end do
          end do
!$omp end parallel do

!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                ! copy sedge into sp and sm
                sp(i,j,hi(3)-1) = sedge(i,j,hi(3))
                sm(i,j,hi(3)  ) = sedge(i,j,hi(3))
                ! reset sm on second interior edge
                sm(i,j,hi(3)-1) = sedge(i,j,hi(3)-1)
             end do
          end do
!$omp end parallel do

          ! modify using quadratic limiters
          k = hi(3)-1
!$omp parallel do private(i,j)
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
!$omp end parallel do
       end if

    else if (ppm_type .eq. 2) then

       ! interpolate s to z-edges
!$omp parallel do private(i,j,k,D2,D2L,D2R,sgn,D2LIM)
       do k=lo(3)-2,hi(3)+3
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,k) = (7.d0/12.d0)*(s(i,j,k-1)+s(i,j,k)) &
                     - (1.d0/12.d0)*(s(i,j,k-2)+s(i,j,k+1))
                ! limit sedge
                if ((sedge(i,j,k)-s(i,j,k-1))*(s(i,j,k)-sedge(i,j,k)) .lt. ZERO) then
                   D2  = THREE*(s(i,j,k-1)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                   D2R = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                   sgn = sign(ONE,D2)
                   D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                   sedge(i,j,k) = HALF*(s(i,j,k-1)+s(i,j,k)) - SIXTH*D2LIM
                end if
             end do
          end do
       end do
!$omp end parallel do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i,j,k+1)-s(i,j,k)
                alpham = sedge(i,j,k  )-s(i,j,k)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i,j,k) - sedge(i,j,k-1)
                   dafacep = sedge(i,j,k+2) - sedge(i,j,k+1)
                   dabarm = s(i,j,k) - s(i,j,k-1)
                   dabarp = s(i,j,k+1) - s(i,j,k)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                   D2R = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
                   D2C = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                   alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i,j,k-1) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else 
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i,j,k+1) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if
                
                sm(i,j,k) = s(i,j,k) + alpham
                sp(i,j,k) = s(i,j,k) + alphap

             end do
          end do
       end do
!$omp end parallel do

       ! different stencil needed for z-component of EXT_DIR and HOEXTRAP bc's
       if (bc(3,1) .eq. EXT_DIR  .or. bc(3,1) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)
          sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,lo(3)+1) = -FIFTH        *s(i,j,lo(3)-1) &
                                     + (THREE/FOUR)*s(i,j,lo(3)  ) &
                                     + HALF        *s(i,j,lo(3)+1) &
                                     - (ONE/20.0d0)*s(i,j,lo(3)+2)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,lo(3)+1) = max(sedge(i,j,lo(3)+1),min(s(i,j,lo(3)+1),s(i,j,lo(3))))
                sedge(i,j,lo(3)+1) = min(sedge(i,j,lo(3)+1),max(s(i,j,lo(3)+1),s(i,j,lo(3))))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sp
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sp(i,j,lo(3)  ) = sedge(i,j,lo(3)+1)
             end do
          end do
!$omp end parallel do
          
          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do j=lo(2)-1,hi(2)+1
             do k=lo(3)+1,lo(3)+2
                do i=lo(1)-1,hi(1)+1

                   alphap = sedge(i,j,k+1)-s(i,j,k)
                   alpham = sedge(i,j,k  )-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i,j,k-1)
                      dafacep = sedge(i,j,k+2) - sedge(i,j,k+1)
                      dabarm = s(i,j,k) - s(i,j,k-1)
                      dabarp = s(i,j,k+1) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
                      D2C = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j,k-1) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j,k+1) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

       if (bc(3,2) .eq. EXT_DIR  .or. bc(3,2) .eq. HOEXTRAP) then
          ! the value in the first cc ghost cell represents the edge value
          sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)
          sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1) = &
               s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)

          ! use a modified stencil to get sedge on the first interior edge
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,hi(3)) = -FIFTH        *s(i,j,hi(3)+1) &
                                   + (THREE/FOUR)*s(i,j,hi(3)  ) &
                                   + HALF        *s(i,j,hi(3)-1) &
                                   - (ONE/20.0d0)*s(i,j,hi(3)-2)
             end do
          end do
!$omp end parallel do

          ! make sure sedge lies in between adjacent cell-centered values
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sedge(i,j,hi(3)) = max(sedge(i,j,hi(3)),min(s(i,j,hi(3)-1),s(i,j,hi(3))))
                sedge(i,j,hi(3)) = min(sedge(i,j,hi(3)),max(s(i,j,hi(3)-1),s(i,j,hi(3))))
             end do
          end do
!$omp end parallel do

          ! copy sedge into sm
!$omp parallel do private(i,j)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sm(i,j,hi(3)  ) = sedge(i,j,hi(3))
             end do
          end do
!$omp end parallel do

          ! apply Colella 2008 limiters to compute sm and sp in the second
          ! and third inner cells
!$omp parallel do private(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep, &
!$omp dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax, &
!$omp delam,delap)
          do j=lo(2)-1,hi(2)+1
             do k=hi(3)-2,hi(3)-1
                do i=lo(1)-1,hi(1)+1

                   alphap = sedge(i,j,k+1)-s(i,j,k)
                   alpham = sedge(i,j,k  )-s(i,j,k)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is
                      ! the largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge(i,j,k) - sedge(i,j,k-1)
                      dafacep = sedge(i,j,k+2) - sedge(i,j,k+1)
                      dabarm = s(i,j,k) - s(i,j,k-1)
                      dabarp = s(i,j,k+1) - s(i,j,k)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                      D2R = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
                      D2C = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                      alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j,k-1) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else 
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j,k+1) - s(i,j,k)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j,k) = s(i,j,k) + alpham
                   sp(i,j,k) = s(i,j,k) + alphap

                end do
             end do
          end do
!$omp end parallel do
       end if

    end if

    ! compute z-component of Ip and Im
!$omp parallel do private(i,j,k,w0lo,w0hi,velhi,vello,sigmap,sigmam,s6)
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
                velhi = wmac(i,j,k+1) + w0macz(i,j,k+1)
                vello = wmac(i,j,k  ) + w0macz(i,j,k  )
             else
                velhi = wmac(i,j,k+1) + w0hi
                vello = wmac(i,j,k  ) + w0lo
             end if
             sigmap = abs(velhi)*dt/dx(2)
             sigmam = abs(vello)*dt/dx(2)
             s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
             if (velhi .gt. rel_eps) then
                Ip(i,j,k,3) = sp(i,j,k) - (sigmap/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigmap)*s6)
             else
                Ip(i,j,k,3) = s(i,j,k)
             end if
             if (vello .lt. -rel_eps) then
                Im(i,j,k,3) = sm(i,j,k) + (sigmam/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigmam)*s6)
             else
                Im(i,j,k,3) = s(i,j,k)
             end if
          end do
       end do
    end do
!$omp end parallel do

    deallocate(sp,sm,dsvl,sedge)

  end subroutine ppm_fpu_3d

end module ppm_module
