module macrhs_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_module
  use heating_module
  use eos_module
  use network
  use variables

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_macrhs (macrhs,s,u,div_coeff,p0,t0,gam1,dx,time)

      type(multifab) , intent(inout) :: macrhs
      type(multifab) , intent(in   ) :: s,u
      real(kind=dp_t), intent(in   ) :: div_coeff(:),p0(:),t0(:),gam1(:)
      real(kind=dp_t), intent(in   ) :: dx(:), time

      real(kind=dp_t), pointer:: mp(:,:,:,:),sp(:,:,:,:),up(:,:,:,:)
      integer :: lo(s%dim),hi(s%dim),ng,dm
      integer :: i

      ng = s%ng
      dm = s%dim

      do i = 1, s%nboxes
         if ( multifab_remote(s, i) ) cycle
         mp => dataptr(macrhs, i)
         sp => dataptr(s, i)
         up => dataptr(u, i)
         lo =  lwb(get_box(s, i))
         hi =  upb(get_box(s, i))
         select case (dm)
            case (2)
              call make_macrhs_2d(lo,hi,mp(:,:,1,1),sp(:,:,1,:), up(:,:,1,:), &
                                  ng, div_coeff, p0, t0, gam1, dx, time)
            case (3)
              call make_macrhs_3d(lo,hi,mp(:,:,:,1),sp(:,:,:,:), up(:,:,:,:), &
                                  ng, div_coeff, p0, t0, gam1, dx, time)
         end select
      end do

   end subroutine make_macrhs

   subroutine make_macrhs_2d (lo,hi,rhs,s,u,ng,div_coeff,p0,t0,gam1,dx,time)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(in   ) :: u(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(in   ) :: div_coeff(lo(2):)
      real (kind=dp_t), intent(in   ) ::        p0(lo(2):)
      real (kind=dp_t), intent(in   ) ::        t0(lo(2):)
      real (kind=dp_t), intent(in   ) ::      gam1(lo(2):)
      real (kind=dp_t), intent(in   ) :: dx(:), time

!     Local variables
      integer :: i, j
      integer :: imax, jmax

      real(kind=dp_t) :: x,y,rhs_max
      real(kind=dp_t) :: sigma_H,denom
      real(kind=dp_t), allocatable :: H(:,:)

      real(kind=dp_t) :: gradp0,dgam,dgam_max,gam1_save

      allocate(H(lo(1):hi(1),lo(2):hi(2)))

      rhs = zero
      rhs_max = ZERO

      denom = 1.0_dp_t / dble(hi(1)-lo(1)+1)

      call get_H_2d(H,lo,hi,dx,time)

      do_diag = .false.

      dgam_max = ZERO

      do j = lo(2), hi(2)

        sigma_H = 0.0_dp_t

        if (j.eq.lo(2)) then
          gradp0 =        (p0(j+1) - p0(j)) / dx(2)
        else if (j.eq.hi(2)) then
          gradp0 =        (p0(j) - p0(j-1)) / dx(2)
        else
          gradp0 = HALF * (p0(j+1) - p0(j-1)) / dx(2)
        end if

        do i = lo(1), hi(1)

           den_row(1) = s(i,j,rho_comp)
           temp_row(1) = t0(j)
           p_row(1) = p0(j)
           xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

           
           ! (rho, P) --> T
           input_flag = 4

           call eos(input_flag, den_row, temp_row, &
                    npts, nspec, &
                    xn_zone, aion, zion, &
                    p_row, h_row, e_row, & 
                    cv_row, cp_row, xne_row, eta_row, pele_row, &
                    dpdt_row, dpdr_row, dedt_row, dedr_row, &
                    dpdX_row, dhdX_row, &
                    gam1_row, cs_row, s_row, &
                    do_diag)

!          dgam = abs(gam1_row(1) - gam1(j))
!          if (dgam .gt. dgam_max) then
!             dgam_max = dgam
!             gam1_save = gam1_row(1)
!             imax = i
!             jmax = j
!          end if
!          rhs(i,j) = (gam1_row(1) - gam1(j)) / (gam1_row(1) * gam1(j)) *  u(i,j,2) * gradp0 / p0(j)
!          rhs(i,j) = div_coeff(j) * rhs(i,j)

           H(i,j) = H(i,j) * dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))
           sigma_H = sigma_H + H(i,j)
        enddo
        sigma_H = sigma_H * denom

        do i = lo(1), hi(1)
           rhs(i,j) = div_coeff(j) * (H(i,j) - sigma_H)
           rhs_max = max(rhs_max, abs(rhs(i,j)))
        enddo
      enddo

      print *,'MACRHS: DIVU AT TIME ',time, rhs_max

      deallocate(H)
 
   end subroutine make_macrhs_2d

   subroutine make_macrhs_3d (lo,hi,rhs,s,u,ng,div_coeff,p0,t0,gam1,dx,time)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(in   ) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(in   ) :: div_coeff(lo(3):)
      real (kind=dp_t), intent(in   ) ::        p0(lo(3):)
      real (kind=dp_t), intent(in   ) ::        t0(lo(3):)
      real (kind=dp_t), intent(in   ) ::      gam1(lo(3):)
      real (kind=dp_t), intent(in   ) :: dx(:), time

!     Local variables
      integer :: i, j, k 
      integer :: imax, jmax, kmax

      real(kind=dp_t) :: x,y,z,rhs_max
      real(kind=dp_t) :: sigma_H,denom
      real(kind=dp_t), allocatable :: H(:,:,:)

      real(kind=dp_t) :: gradp0,dgam,dgam_max,gam1_save

      allocate(H(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      rhs = zero
      rhs_max = ZERO

      denom = 1.0_dp_t / dble(hi(1)-lo(1)+1)

      call get_H_3d(H,lo,hi,dx,time)

      do_diag = .false.

      dgam_max = ZERO

      do k = lo(3), hi(3)

        sigma_H = 0.0_dp_t

        if (k.eq.lo(3)) then
          gradp0 =        (p0(k+1) - p0(k)) / dx(2)
        else if (k.eq.hi(3)) then
          gradp0 =        (p0(j) - p0(k-1)) / dx(2)
        else
          gradp0 = HALF * (p0(k+1) - p0(k-1)) / dx(2)
        end if

        do j = lo(3), hi(2)
           do i = lo(1), hi(1)

              den_row(1) = s(i,j,k,rho_comp)
              temp_row(1) = t0(k)
              p_row(1) = p0(k)
              xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

              ! (rho, P) --> T
              input_flag = 4

              call eos(input_flag, den_row, temp_row, &
                       npts, nspec, &
                       xn_zone, aion, zion, &
                       p_row, h_row, e_row, & 
                       cv_row, cp_row, xne_row, eta_row, pele_row, &
                       dpdt_row, dpdr_row, dedt_row, dedr_row, &
                       dpdX_row, dhdX_row, &
                       gam1_row, cs_row, s_row, &
                       do_diag)

!          dgam = abs(gam1_row(1) - gam1(k))
!          if (dgam .gt. dgam_max) then
!             dgam_max = dgam
!             gam1_save = gam1_row(1)
!             imax = i
!             jmax = j
!             kmax = k
!          end if
!          rhs(i,j,k) = (gam1_row(1) - gam1(k)) / (gam1_row(1) * gam1(k)) *  u(i,j,k,3) * gradp0 / p0(k)
!          rhs(i,j,k) = div_coeff(k) * rhs(i,j,k)

              H(i,j,k) = H(i,j,k) * dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))
              sigma_H = sigma_H + H(i,j,k)
           enddo
        enddo
        sigma_H = sigma_H * denom

        do j = lo(1), hi(1)
           do i = lo(1), hi(1)
              rhs(i,j,k) = div_coeff(k) * (H(i,j,k) - sigma_H)
              rhs_max = max(rhs_max, abs(rhs(i,j,k)))
           enddo
        enddo
     enddo

     print *,'MACRHS: DIVU AT TIME ',time, rhs_max

     deallocate(H)
 
   end subroutine make_macrhs_3d

end module macrhs_module
