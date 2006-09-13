module machno_module

  use bl_types
  use bc_module
  use multifab_module
  use eos_module
  use network

  implicit none

contains

   subroutine make_machno (machno,deltap,u,s,p0,temp0)

      type(multifab) , intent(inout) :: machno
      type(multifab) , intent(inout) :: deltap
      type(multifab) , intent(in   ) :: u,s
      real(kind=dp_t), intent(in   ) :: p0(:),temp0(:)

      real(kind=dp_t), pointer:: up(:,:,:,:)
      real(kind=dp_t), pointer:: sp(:,:,:,:)
      real(kind=dp_t), pointer:: mp(:,:,:,:)
      real(kind=dp_t), pointer:: dp(:,:,:,:)
      integer :: lo(s%dim),hi(s%dim),ng,dm
      integer :: i

      ng = s%ng
      dm = s%dim

      do i = 1, s%nboxes
         if ( multifab_remote(s, i) ) cycle
         up => dataptr(u, i)
         sp => dataptr(s, i)
         mp => dataptr(machno, i)
         dp => dataptr(deltap, i)
         lo =  lwb(get_box(s, i))
         hi =  upb(get_box(s, i))
         select case (dm)
            case (2)
              call makemachno_2d(mp(:,:,1,1),dp(:,:,1,1),up(:,:,1,:),sp(:,:,1,:), lo, hi, ng, p0, temp0)
            case (3)
              call makemachno_3d(mp(:,:,:,1),dp(:,:,:,1),up(:,:,:,:),sp(:,:,:,:), lo, hi, ng, p0, temp0)
         end select
      end do

   end subroutine make_machno

   subroutine makemachno_2d (machno,deltap,u,s,lo,hi,ng,p0,temp0)

      implicit none
      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: machno(lo(1):,lo(2):)  
      real (kind = dp_t), intent(  out) :: deltap(lo(1):,lo(2):)  
      real (kind = dp_t), intent(in   ) ::       u(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(in   ) ::      p0(lo(2):)
      real (kind = dp_t), intent(in   ) ::   temp0(lo(2):)

!     Local variables
      integer :: i, j
      real (kind = dp_t) :: vel

      integer :: input_flag
      integer :: nx
      
      do_diag = .false.

      do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             den_row(1) = s(i,j,1)
               p_row(1) = p0(j)

             ! (rho,h) --> T,p, etc
!            h_row(1) = s(i,j,2) / s(i,j,1)
!            input_flag = 2

             ! (rho,P) --> T,h, etc
             input_flag = 4

             call eos(input_flag, den_row, temp_row, npts, nspec, &
                  xmass, aion, zion, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, &
                  pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
                  s_row, do_diag)
             vel = sqrt(u(i,j,1)*u(i,j,1) + u(i,j,2)*u(i,j,2))
             machno(i,j) = vel / cs_row(1)
             deltap(i,j) = (p_row(1)-p0(j))/ p0(j)
          enddo
!       end if
      enddo

   end subroutine makemachno_2d

   subroutine makemachno_3d (machno,deltap,u,s,lo,hi,ng,p0,temp0)

      implicit none
      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: machno(lo(1):,lo(2):,lo(3):)  
      real (kind = dp_t), intent(  out) :: deltap(lo(1):,lo(2):,lo(3):)  
      real (kind = dp_t), intent(in   ) ::       u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::      p0(lo(3):)
      real (kind = dp_t), intent(in   ) ::   temp0(lo(3):)

!     Local variables
      integer :: i, j, k
      real (kind = dp_t) :: vel

      integer :: input_flag
      integer :: nx
      
      do_diag = .false.

      do k = lo(3), hi(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             den_row(1) = s(i,j,k,1)
               p_row(1) = p0(k)

             ! (rho,h) --> T,p, etc
!            h_row(1) = s(i,j,k,2) / s(i,j,k,1)
!            input_flag = 2

             ! (rho,P) --> T,h, etc
             input_flag = 4

             call eos(input_flag, den_row, temp_row, npts, nspec, &
                  xmass, aion, zion, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, &
                  pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
                  s_row, do_diag)
             vel = sqrt(u(i,j,k,1)*u(i,j,k,1) + u(i,j,k,2)*u(i,j,k,2) + u(i,j,k,3)*u(i,j,k,3))
             machno(i,j,k) = vel / cs_row(1)
             deltap(i,j,k) = (p_row(1)-p0(k))/ p0(k)
          enddo
          enddo
!       end if
      enddo

   end subroutine makemachno_3d

end module machno_module
