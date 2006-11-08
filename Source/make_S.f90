module make_S_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use heating_module
  use eos_module
  use fill_3d_module
  use network
  use geometry
  use variables

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_S (Source,state,p0,t0,gam1,dx,time)

      type(multifab) , intent(inout) :: Source
      type(multifab) , intent(in   ) :: state
      real(kind=dp_t), intent(in   ) :: p0(:),t0(:),gam1(:)
      real(kind=dp_t), intent(in   ) :: dx(:), time

      real(kind=dp_t), pointer:: srcp(:,:,:,:),sp(:,:,:,:),up(:,:,:,:)
      integer :: lo(state%dim),hi(state%dim),ng,dm
      integer :: i

      ng = state%ng
      dm = state%dim

      do i = 1, state%nboxes
         if ( multifab_remote(state, i) ) cycle
         srcp => dataptr(Source, i)
         sp => dataptr(state, i)
         lo =  lwb(get_box(state, i))
         hi =  upb(get_box(state, i))
         select case (dm)
            case (2)
              call make_S_2d(lo,hi,srcp(:,:,1,1),sp(:,:,1,:), &
                             ng, p0, t0, gam1, dx, time)
            case (3)
              call make_S_3d(lo,hi,srcp(:,:,:,1),sp(:,:,:,:), &
                             ng, p0, t0, gam1, dx, time)
         end select
      end do

   end subroutine make_S

   subroutine make_S_2d (lo,hi,Source,s,ng,p0,t0,gam1,dx,time)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) :: Source(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(in   ) ::        p0(lo(2):)
      real (kind=dp_t), intent(in   ) ::        t0(lo(2):)
      real (kind=dp_t), intent(in   ) ::      gam1(lo(2):)
      real (kind=dp_t), intent(in   ) :: dx(:), time

!     Local variables
      integer :: i, j
      integer :: imax, jmax

      real(kind=dp_t) :: x,y,Smax
      real(kind=dp_t), allocatable :: H(:,:)

      allocate(H(lo(1):hi(1),lo(2):hi(2)))

      Source = zero
      Smax = ZERO

      call get_H_2d(H,lo,hi,dx,time)

      do_diag = .false.

      do j = lo(2), hi(2)

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

           Source(i,j) = H(i,j) * dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))
           Smax = max(Smax, abs(Source(i,j)))
        enddo
      enddo

      print *,'new S at time ',time, Smax

      deallocate(H)
 
   end subroutine make_S_2d

   subroutine make_S_3d (lo,hi,Source,s,ng,p0,t0,gam1,dx,time)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) :: Source(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(in   ) ::        p0(lo(3):)
      real (kind=dp_t), intent(in   ) ::        t0(lo(3):)
      real (kind=dp_t), intent(in   ) ::      gam1(lo(3):)
      real (kind=dp_t), intent(in   ) :: dx(:), time

!     Local variables
      integer :: i, j, k 
      integer :: imax, jmax, kmax

      real(kind=dp_t) :: x,y,z,Smax
      real(kind=dp_t), allocatable :: p0_cart(:,:,:)
      real(kind=dp_t), allocatable :: t0_cart(:,:,:)
      real(kind=dp_t), allocatable :: H(:,:,:)

      allocate(H(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      if (spherical .eq. 1) then
        allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        allocate(t0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call fill_3d_data(p0_cart,p0,dx,0)
        call fill_3d_data(t0_cart,t0,dx,0)
      end if

      Source = zero
      Smax = zero

      call get_H_3d(H,lo,hi,dx,time)

      do_diag = .false.

      do k = lo(3), hi(3)

        do j = lo(3), hi(2)
           do i = lo(1), hi(1)

              den_row(1) = s(i,j,k,rho_comp)

              if (spherical .eq. 1) then
                temp_row(1) = t0(k)
                p_row(1) = p0(k)
              else
                temp_row(1) = t0_cart(i,j,k)
                p_row(1) = p0_cart(i,j,k)
              end if

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

              Source(i,j,k) = H(i,j,k) * dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))
              Smax = max(Smax, abs(Source(i,j,k)))
           enddo
        enddo
      enddo

      print *,'new S at time ',time, Smax

      deallocate(H)
      if (spherical .eq. 1) then
        deallocate(p0_cart,t0_cart)
      end if
 
   end subroutine make_S_3d

end module make_S_module
