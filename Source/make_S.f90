module make_S_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use eos_module
  use fill_3d_module
  use network
  use geometry
  use variables

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_S (Source,gamma1_term,state,u,rho_omegadot,rho_Hext,p0,t0,gam1,dx,time)

      type(multifab) , intent(inout) :: Source, gamma1_term
      type(multifab) , intent(in   ) :: state, u
      type(multifab) , intent(in   ) :: rho_omegadot
      type(multifab) , intent(in   ) :: rho_Hext
      real(kind=dp_t), intent(in   ) :: p0(:),t0(:),gam1(:)
      real(kind=dp_t), intent(in   ) :: dx(:), time

      real(kind=dp_t), pointer:: srcp(:,:,:,:),gp(:,:,:,:),sp(:,:,:,:),up(:,:,:,:)
      real(kind=dp_t), pointer:: omegap(:,:,:,:), hp(:,:,:,:)
      integer :: lo(state%dim),hi(state%dim),ng,dm
      integer :: i

      ng = state%ng
      dm = state%dim

      do i = 1, state%nboxes
         if ( multifab_remote(state, i) ) cycle
         srcp => dataptr(Source, i)
         gp => dataptr(gamma1_term, i)
         sp => dataptr(state, i)
         up => dataptr(u, i)
         omegap => dataptr(rho_omegadot, i)
         hp     => dataptr(rho_Hext, i)
         lo =  lwb(get_box(state, i))
         hi =  upb(get_box(state, i))
         select case (dm)
            case (2)
              call make_S_2d(lo,hi,srcp(:,:,1,1),gp(:,:,1,1),sp(:,:,1,:),up(:,:,1,:), &
                             omegap(:,:,1,:), hp(:,:,1,1), &
                             ng, p0, t0, gam1, dx, time)
            case (3)
              call make_S_3d(lo,hi,srcp(:,:,:,1),gp(:,:,:,1),sp(:,:,:,:),up(:,:,:,:), &
                             omegap(:,:,:,:), hp(:,:,:,1), &
                             ng, p0, t0, gam1, dx, time)
         end select
      end do

   end subroutine make_S

   subroutine make_S_2d (lo,hi,Source,gamma1_term,s,u, &
                         rho_omegadot,rho_Hext,ng,p0,t0,gam1,dx,time)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) :: Source(lo(1):,lo(2):)
      real (kind=dp_t), intent(  out) :: gamma1_term(lo(1):,lo(2):)
      real (kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(in   ) :: u(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1):,lo(2):,:)
      real (kind=dp_t), intent(in   ) ::     rho_Hext(lo(1):,lo(2):)
      real (kind=dp_t), intent(in   ) ::        p0(lo(2):)
      real (kind=dp_t), intent(in   ) ::        t0(lo(2):)
      real (kind=dp_t), intent(in   ) ::      gam1(lo(2):)
      real (kind=dp_t), intent(in   ) :: dx(:), time

!     Local variables
      integer :: i, j, n
      integer :: imax, jmax

      real(kind=dp_t) :: x,y,Smax
      real(kind=dp_t) :: sigma, react_term, pres_term, gradp0

      Source = zero
      Smax = ZERO

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
                    dsdt_row, dsdr_row, &
                    do_diag)

           sigma = dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))

           react_term = ZERO
           pres_term = ZERO
           do n = 1, nspec
              react_term = react_term - &
                   (dhdX_row(1,n) + ebin(n))*rho_omegadot(i,j,n)/den_row(1)

              pres_term = pres_term + &
                   dpdX_row(1,n)*rho_omegadot(i,j,n)/den_row(1)
           enddo

           Source(i,j) = sigma*(rho_Hext(i,j)/den_row(1) + react_term) + &
                pres_term/(den_row(1)*dpdr_row(1)) 

           if (j .eq. lo(2)) then
              gradp0 = (p0(j+1) - p0(j))/dx(2)
           else if (j .eq. hi(2)) then
              gradp0 = (p0(j) - p0(j-1))/dx(2)
           else
              gradp0 = HALF*(p0(j+1) - p0(j-1))/dx(2)
           endif

!           gamma1_term(i,j) = (gam1_row(1) - gam1(j))*u(i,j,2)*gradp0/(gam1(j)*gam1(j)*p0(j))
           gamma1_term(i,j) = 0.0_dp_t

           Smax = max(Smax, abs(Source(i,j)))
        enddo
      enddo

      print *,'new S at time ',time, Smax

 
   end subroutine make_S_2d

   subroutine make_S_3d (lo,hi,Source,gamma1_term,s,u, &
                         rho_omegadot,rho_Hext,ng,p0,t0,gam1,dx,time)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) :: Source(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(  out) :: gamma1_term(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(in   ) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1):,lo(2):,lo(3):,:)
      real (kind=dp_t), intent(in   ) ::     rho_Hext(lo(1):,lo(2):,lo(3):)
      real (kind=dp_t), intent(in   ) ::        p0(lo(3):)
      real (kind=dp_t), intent(in   ) ::        t0(lo(3):)
      real (kind=dp_t), intent(in   ) ::      gam1(lo(3):)
      real (kind=dp_t), intent(in   ) :: dx(:), time

!     Local variables
      integer :: i, j, k , n
      integer :: imax, jmax, kmax

      real(kind=dp_t) :: x,y,z,Smax
      real(kind=dp_t), allocatable :: p0_cart(:,:,:)
      real(kind=dp_t), allocatable :: t0_cart(:,:,:)
      real(kind=dp_t) :: sigma, react_term, pres_term

      if (spherical .eq. 1) then
        allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        allocate(t0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call fill_3d_data(p0_cart,p0,lo,hi,dx,0)
        call fill_3d_data(t0_cart,t0,lo,hi,dx,0)
      end if

      Source = zero
      Smax = zero

      do_diag = .false.

      do k = lo(3), hi(3)

        do j = lo(3), hi(2)
           do i = lo(1), hi(1)

              den_row(1) = s(i,j,k,rho_comp)

              if (spherical .eq. 0) then
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
                       dsdt_row, dsdr_row, &
                       do_diag)

              sigma = dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))

              react_term = ZERO
              pres_term = ZERO
              do n = 1, nspec
                 react_term = react_term - &
                      (dhdX_row(1,n) + ebin(n))*rho_omegadot(i,j,k,n)/den_row(1)

                 pres_term = pres_term + &
                      dpdX_row(1,n)*rho_omegadot(i,j,k,n)/den_row(1)
              enddo
  
              Source(i,j,k) = sigma*(rho_Hext(i,j,k)/den_row(1) + react_term) + &
                   pres_term/(den_row(1)*dpdr_row(1))

              Smax = max(Smax, abs(Source(i,j,k)))

              gamma1_term(i,j,k) = 0.0_dp_t
           enddo
        enddo
      enddo

      print *,'new S at time ',time, Smax

      if (spherical .eq. 1) then
        deallocate(p0_cart,t0_cart)
      end if
 
   end subroutine make_S_3d

end module make_S_module
