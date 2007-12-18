module vort_module

  use bl_types
  use bc_module
  use bl_constants_module
  use define_bc_module
  use multifab_module
  use variables

  implicit none

  private
  public :: make_vorticity, make_magvel

contains

   subroutine make_vorticity (vort,comp,u,dx,bc)

      integer        , intent(in   ) :: comp
      type(multifab) , intent(in   ) :: vort
      type(multifab) , intent(inout) :: u
      real(kind=dp_t), intent(in   ) :: dx(:)
      type(bc_level) , intent(in   ) :: bc

      real(kind=dp_t), pointer:: up(:,:,:,:)
      real(kind=dp_t), pointer:: vp(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      integer :: i

      ng = u%ng
      dm = u%dim
      call multifab_fill_boundary(u)

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         up => dataptr(u, i)
         vp => dataptr(vort, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call makevort_2d(vp(:,:,1,comp),up(:,:,1,:), lo, hi, ng, dx, bc%phys_bc_level_array(i,:,:))
            case (3)
              call makevort_3d(vp(:,:,:,comp),up(:,:,:,:), lo, hi, ng, dx, bc%phys_bc_level_array(i,:,:))
         end select
      end do

   end subroutine make_vorticity

   subroutine makevort_2d (vort,u,lo,hi,ng,dx,bc)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: vort(lo(1):,lo(2):)  
      real (kind = dp_t), intent(in   ) ::    u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      integer           , intent(in   ) :: bc(:,:)

!     Local variables
      integer :: i, j
      real (kind = dp_t) :: vx,uy
   
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           vx = (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
           uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
           vort(i,j) = vx - uy
        enddo
      enddo

      if (bc(1,1) .eq. INLET .or. bc(1,1) .eq. SLIP_WALL .or. &
          bc(1,1) .eq. NO_SLIP_WALL) then
        i = lo(1)
        do j = lo(2), hi(2)
           vx = (u(i+1,j,2) + 3.d0*u(i,j,2) - 4.d0*u(i-1,j,2)) / dx(1)
           uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
           vort(i,j) = vx - uy
        end do
      end if

      if (bc(1,2) .eq. INLET .or. bc(1,2) .eq. SLIP_WALL .or. &
          bc(1,2) .eq. NO_SLIP_WALL) then
        i = hi(1)
        do j = lo(2), hi(2)
           vx = -(u(i-1,j,2) + 3.d0*u(i,j,2) - 4.d0*u(i+1,j,2)) / dx(1)
           uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
           vort(i,j) = vx - uy
        end do
      end if

      if (bc(2,1) .eq. INLET .or. bc(2,1) .eq. SLIP_WALL .or. &
          bc(2,1) .eq. NO_SLIP_WALL) then
        j = lo(2)
        do i = lo(1), hi(1)
           vx = (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
           uy = (u(i,j+1,1) + 3.d0*u(i,j,1) - 4.d0*u(i,j-1,1)) / dx(2)
           vort(i,j) = vx - uy
        end do
      end if

      if (bc(2,2) .eq. INLET .or. bc(2,2) .eq. SLIP_WALL .or. &
          bc(2,2) .eq. NO_SLIP_WALL) then
        j = hi(2)
        do i = lo(1), hi(1)
           vx =  (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
           uy = -(u(i,j-1,1) + 3.d0*u(i,j,1) - 4.d0*u(i,j+1,1)) / dx(2)
           vort(i,j) = vx - uy
        end do
      end if

   end subroutine makevort_2d

   subroutine makevort_3d (vort,u,lo,hi,ng,dx,bc)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: vort(lo(1):,lo(2):,lo(3):)
      real (kind = dp_t), intent(in   ) ::    u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      integer           , intent(in   ) :: bc(:,:)

!     Local variables
      integer :: i, j, k
      logical :: fix_lo_x,fix_hi_x,fix_lo_y,fix_hi_y,fix_lo_z,fix_hi_z
      real (kind = dp_t) :: wy,vz,uz,wx,vx,uy

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         uy = uycen(i,j,k)
         uz = uzcen(i,j,k)
         vx = vxcen(i,j,k)
         vz = vzcen(i,j,k)
         wx = wxcen(i,j,k)
         wy = wycen(i,j,k)
         vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
      enddo
      enddo
      enddo

      fix_lo_x = ( bc(1,1) .eq. INLET .or. bc(1,1) .eq. NO_SLIP_WALL )
      fix_hi_x = ( bc(1,2) .eq. INLET .or. bc(1,2) .eq. NO_SLIP_WALL )
                   
      fix_lo_y = ( bc(2,1) .eq. INLET .or. bc(2,1) .eq. NO_SLIP_WALL )
      fix_hi_y = ( bc(2,2) .eq. INLET .or. bc(2,2) .eq. NO_SLIP_WALL )

      fix_lo_z = ( bc(3,1) .eq. INLET .or. bc(3,1) .eq. NO_SLIP_WALL )
      fix_hi_z = ( bc(3,2) .eq. INLET .or. bc(3,2) .eq. NO_SLIP_WALL )

!
!     First do all the faces
!
      if (fix_lo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = vxlo(i,j,k)
               wx = wxlo(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fix_hi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = vxhi(i,j,k)
               wx = wxhi(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fix_lo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uylo(i,j,k)
               wy = wylo(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fix_hi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uyhi(i,j,k)
               wy = wyhi(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fix_lo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzlo(i,j,k)
               vz = vzlo(i,j,k)
               vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fix_hi_z) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzhi(i,j,k)
               vz = vzhi(i,j,k)
               vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if
!
!     Next do all the edges
!
      if (fix_lo_x .and. fix_lo_y) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = vxlo(i,j,k)
            wx = wxlo(i,j,k)
            uy = uylo(i,j,k)
            wy = wylo(i,j,k)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_hi_x .and. fix_lo_y) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = vxhi(i,j,k)
            wx = wxhi(i,j,k)
            uy = uylo(i,j,k)
            wy = wylo(i,j,k)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_lo_x .and. fix_hi_y) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = vxlo(i,j,k)
            wx = wxlo(i,j,k)
            uy = uyhi(i,j,k)
            wy = wyhi(i,j,k)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_hi_x .and. fix_hi_y) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = vxhi(i,j,k)
            wx = wxhi(i,j,k)
            uy = uyhi(i,j,k)
            wy = wyhi(i,j,k)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_lo_x .and. fix_lo_z) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = vxlo(i,j,k)
            wx = wxlo(i,j,k)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = uzlo(i,j,k)
            vz = vzlo(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_hi_x .and. fix_lo_z) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = vxhi(i,j,k)
            wx = wxhi(i,j,k)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = uzlo(i,j,k)
            vz = vzlo(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_lo_x .and. fix_hi_z) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = vxlo(i,j,k)
            wx = wxlo(i,j,k)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = uzhi(i,j,k)
            vz = vzhi(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_hi_x .and. fix_hi_z) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = vxhi(i,j,k)
            wx = wxhi(i,j,k)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = uzhi(i,j,k)
            vz = vzhi(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_lo_y .and. fix_lo_z) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = uylo(i,j,k)
            wy = wylo(i,j,k)
            uz = uzlo(i,j,k)
            vz = vzlo(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_hi_y .and. fix_lo_z) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = uyhi(i,j,k)
            wy = wyhi(i,j,k)
            uz = uzlo(i,j,k)
            vz = vzlo(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_lo_y .and. fix_hi_z) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = uylo(i,j,k)
            wy = wylo(i,j,k)
            uz = uzhi(i,j,k)
            vz = vzhi(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if (fix_hi_y .and. fix_hi_z) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = uyhi(i,j,k)
            wy = wyhi(i,j,k)
            uz = uzhi(i,j,k)
            vz = vzhi(i,j,k)
            vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if
!
!     Finally do all the corners
!
      if (fix_lo_x .and. fix_lo_y .and. fix_lo_z) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         vx = vxlo(i,j,k)
         wx = wxlo(i,j,k)
         uy = uylo(i,j,k)
         wy = wylo(i,j,k)
         uz = uzlo(i,j,k)
         vz = vzlo(i,j,k)
         vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if (fix_hi_x .and. fix_lo_y .and. fix_lo_z) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         vx = vxhi(i,j,k)
         wx = wxhi(i,j,k)
         uy = uylo(i,j,k)
         wy = wylo(i,j,k)
         uz = uzlo(i,j,k)
         vz = vzlo(i,j,k)
         vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if (fix_lo_x .and. fix_hi_y .and. fix_lo_z) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         vx = vxlo(i,j,k)
         wx = wxlo(i,j,k)
         uy = uyhi(i,j,k)
         wy = wyhi(i,j,k)
         uz = uzlo(i,j,k)
         vz = vzlo(i,j,k)
         vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if (fix_hi_x .and. fix_hi_y .and. fix_lo_z) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         vx = vxhi(i,j,k)
         wx = wxhi(i,j,k)
         uy = uyhi(i,j,k)
         wy = wyhi(i,j,k)
         uz = uzlo(i,j,k)
         vz = vzlo(i,j,k)
         vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if (fix_lo_x .and. fix_lo_y .and. fix_hi_z) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         vx = vxlo(i,j,k)
         wx = wxlo(i,j,k)
         uy = uylo(i,j,k)
         wy = wylo(i,j,k)
         uz = uzhi(i,j,k)
         vz = vzhi(i,j,k)
         vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if (fix_hi_x .and. fix_lo_y .and. fix_hi_z) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         vx = vxhi(i,j,k)
         wx = wxhi(i,j,k)
         uy = uylo(i,j,k)
         wy = wylo(i,j,k)
         uz = uzhi(i,j,k)
         vz = vzhi(i,j,k)
         vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if (fix_lo_x .and. fix_hi_y .and. fix_hi_z) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         vx = vxlo(i,j,k)
         wx = wxlo(i,j,k)
         uy = uyhi(i,j,k)
         wy = wyhi(i,j,k)
         uz = uzhi(i,j,k)
         vz = vzhi(i,j,k)
         vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if (fix_hi_x .and. fix_hi_y .and. fix_hi_z) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         vx = vxhi(i,j,k)
         wx = wxhi(i,j,k)
         uy = uyhi(i,j,k)
         wy = wyhi(i,j,k)
         uz = uzhi(i,j,k)
         vz = vzhi(i,j,k)
         vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

    contains

      function uycen(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r = HALF*(u(i,j+1,k,1)-u(i,j-1,k,1))/dx(2)
      end function uycen

      function uylo(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r = (u(i,j+1,k,1)+THREE*u(i,j,k,1)-FOUR*u(i,j-1,k,1))/(THREE*dx(2))
      end function uylo

      function uyhi(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r = -(u(i,j-1,k,1)+THREE*u(i,j,k,1)-FOUR*u(i,j+1,k,1))/(THREE*dx(2))
      end function uyhi

      function uzcen(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r = HALF*(u(i,j,k+1,1)-u(i,j,k-1,1))/dx(3)
      end function uzcen

      function uzlo(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r = (u(i,j,k+1,1)+THREE*u(i,j,k,1)-FOUR*u(i,j,k-1,1))/(THREE*dx(3))
      end function uzlo

      function uzhi(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r =-(u(i,j,k-1,1)+THREE*u(i,j,k,1)-FOUR*u(i,j,k+1,1))/(THREE*dx(3))
      end function uzhi

      function vxcen(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r = HALF*(u(i+1,j,k,2)-u(i-1,j,k,2))/dx(1)
      end function vxcen

      function vxlo(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r = (u(i+1,j,k,2)+THREE*u(i,j,k,2)-FOUR*u(i-1,j,k,2))/(THREE*dx(1))
      end function vxlo

      function vxhi(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r =-(u(i-1,j,k,2)+THREE*u(i,j,k,2)-FOUR*u(i+1,j,k,2))/(THREE*dx(1))
      end function vxhi

      function vzcen(i,j,k) result(r) 
        integer :: i,j,k
        real(dp_t) :: r
        r = HALF*(u(i,j,k+1,2)-u(i,j,k-1,2))/dx(3)
      end function vzcen

      function vzlo(i,j,k) result(r) 
        integer :: i,j,k
        real(dp_t) :: r
        r = (u(i,j,k+1,2)+THREE*u(i,j,k,2)-FOUR*u(i,j,k-1,2))/(THREE*dx(3))
      end function vzlo

      function vzhi(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r =-(u(i,j,k-1,2)+THREE*u(i,j,k,2)-FOUR*u(i,j,k+1,2))/(THREE*dx(3))
      end function vzhi

      function wxcen(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r = HALF*(u(i+1,j,k,3)-u(i-1,j,k,3))/dx(1)
      end function wxcen

      function wxlo(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r = (u(i+1,j,k,3)+THREE*u(i,j,k,3)-FOUR*u(i-1,j,k,3))/(THREE*dx(1))
      end function wxlo

      function wxhi(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r =-(u(i-1,j,k,3)+THREE*u(i,j,k,3)-FOUR*u(i+1,j,k,3))/(THREE*dx(1))
      end function wxhi

      function wycen(i,j,k) result(r) 
        integer :: i,j,k
        real(dp_t) :: r
        r = HALF*(u(i,j+1,k,3)-u(i,j-1,k,3))/dx(2)
      end function wycen

      function wylo(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r = (u(i,j+1,k,3)+THREE*u(i,j,k,3)-FOUR*u(i,j-1,k,3))/(THREE*dx(2))
      end function wylo

      function wyhi(i,j,k) result(r)
        integer :: i,j,k
        real(dp_t) :: r
        r =-(u(i,j-1,k,3)+THREE*u(i,j,k,3)-FOUR*u(i,j+1,k,3))/(THREE*dx(2))
      end function wyhi

      function vorfun(uy,uz,vx,vz,wx,wy) result(r)
        real(dp_t) :: uy,uz,vx,vz,wx,wy
        real(dp_t) :: r
        r = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
      end function vorfun

   end subroutine makevort_3d

   subroutine make_magvel (plotdata,comp_magvel,comp_mom,u,s)

      integer        , intent(in   ) :: comp_magvel,comp_mom
      type(multifab) , intent(in   ) :: plotdata
      type(multifab) , intent(inout) :: u,s

      real(kind=dp_t), pointer:: pp(:,:,:,:)
      real(kind=dp_t), pointer:: up(:,:,:,:)
      real(kind=dp_t), pointer:: sp(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      integer :: i

      ng = u%ng
      dm = u%dim
      call multifab_fill_boundary(u)
      call multifab_fill_boundary(s)

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         pp => dataptr(plotdata, i)
         up => dataptr(u, i)
         sp => dataptr(s, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call makemagvel_2d(pp(:,:,1,comp_magvel),pp(:,:,1,comp_mom),up(:,:,1,:), sp(:,:,1,rho_comp), lo, hi, ng)
            case (3)
              call makemagvel_3d(pp(:,:,:,comp_magvel),pp(:,:,:,comp_mom),up(:,:,:,:), sp(:,:,:,rho_comp), lo, hi, ng)
         end select
      end do

   end subroutine make_magvel

   subroutine makemagvel_2d (magvel,mom,u,rho,lo,hi,ng)

      integer           , intent(in   ) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: magvel(lo(1):,lo(2):)  
      real (kind = dp_t), intent(  out) ::    mom(lo(1):,lo(2):)  
      real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng:,lo(2)-ng:)  

!     Local variables
      integer :: i, j
   
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           magvel(i,j) = sqrt(u(i,j,1)**2 + u(i,j,2)**2)
              mom(i,j) = rho(i,j) * magvel(i,j)
        enddo
      enddo

   end subroutine makemagvel_2d

   subroutine makemagvel_3d (magvel,mom,u,rho,lo,hi,ng)

      integer           , intent(in   ) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: magvel(lo(1):,lo(2):,lo(3):)
      real (kind = dp_t), intent(  out) ::    mom(lo(1):,lo(2):,lo(3):)
      real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  

!     Local variables
      integer :: i, j, k

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           magvel(i,j,k) = sqrt(u(i,j,k,1)**2 + u(i,j,k,2)**2 + u(i,j,k,3)**2)
              mom(i,j,k) = rho(i,j,k) * magvel(i,j,k)
        enddo
      enddo
      enddo

   end subroutine makemagvel_3d

end module vort_module
