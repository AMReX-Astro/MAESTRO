! compute gamma1 for the full state.

module make_gamma_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: make_gamma

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_gamma(nlevs,gamma,s,p0,tempbar)

    use bl_prof_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: gamma(:)
    type(multifab) , intent(in   ) :: s(:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: tempbar(:,0:,:)
    
    real(kind=dp_t), pointer:: gamp(:,:,:,:),sp(:,:,:,:)
    integer :: lo(s(1)%dim),hi(s(1)%dim),dm
    integer :: i,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_gamma")
    
    dm = s(1)%dim

    do n = 1, nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          gamp => dataptr(gamma(n), i)
          sp   => dataptr(s(n), i)
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
          select case (dm)
          case (2)
             call make_gamma_2d(lo,hi,gamp(:,:,1,1),sp(:,:,1,:),p0(n,:),tempbar(n,:,1))
          case (3)
             call make_gamma_3d(lo,hi,gamp(:,:,:,1),sp(:,:,:,:),p0(n,:),tempbar(n,:,1))
          end select
       end do
    end do

    call destroy(bpt)

   end subroutine make_gamma

   subroutine make_gamma_2d(lo,hi,gamma,s,p0,tempbar)

      use eos_module
      use variables, only: rho_comp, rhoh_comp, spec_comp

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(  out) :: gamma(lo(1)  :,lo(2)  :)
      real (kind=dp_t), intent(in   ) ::     s(lo(1)-3:,lo(2)-3:,:)
      real (kind=dp_t), intent(in   ) :: p0(0:)
      real (kind=dp_t), intent(in   ) :: tempbar(0:)

      ! local variables
      integer :: i, j

      do_diag = .false.

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            
            den_eos(1) = s(i,j,rho_comp)
            p_eos(1) = p0(j)
            xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
            temp_eos(1) = tempbar(j)
            
            ! dens, pres, and xmass are inputs
            call eos(eos_input_rp, den_eos, temp_eos, &
                     npts, nspec, &
                     xn_eos, &
                     p_eos, h_eos, e_eos, & 
                     cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                     dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                     dpdX_eos, dhdX_eos, &
                     gam1_eos, cs_eos, s_eos, &
                     dsdt_eos, dsdr_eos, &
                     do_diag)
            
            gamma(i,j) = gam1_eos(1)

        end do
      end do
 
   end subroutine make_gamma_2d

    subroutine make_gamma_3d(lo,hi,gamma,s,p0,tempbar)

      use eos_module
      use variables, only: rho_comp, rhoh_comp, spec_comp

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(  out) :: gamma(lo(1)  :,lo(2)  :,lo(3)  :)
      real (kind=dp_t), intent(in   ) ::     s(lo(1)-3:,lo(2)-3:,lo(3)-3:,:)
      real (kind=dp_t), intent(in   ) :: p0(0:)
      real (kind=dp_t), intent(in   ) :: tempbar(0:)

      ! local variables
      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               
               den_eos(1) = s(i,j,k,rho_comp)
               p_eos(1) = p0(k)
               xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
               temp_eos(1) = tempbar(k)
               
               ! dens, pres, and xmass are inputs
               call eos(eos_input_rp, den_eos, temp_eos, &
                        npts, nspec, &
                        xn_eos, &
                        p_eos, h_eos, e_eos, & 
                        cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                        dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                        dpdX_eos, dhdX_eos, &
                        gam1_eos, cs_eos, s_eos, &
                        dsdt_eos, dsdr_eos, &
                        do_diag)

!               den_eos(1) = s(i,j,k,rho_comp)
!               h_eos(1) =  s(i,j,k,rhoh_comp)/s(i,j,k,rho_comp)
!               xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
!               
!               ! dens, enthalpy, and xmass are inputs
!               call eos(eos_input_rh, den_eos, temp_eos, &
!                        npts, nspec, &
!                        xn_eos, &
!                        p_eos, h_eos, e_eos, & 
!                        cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
!                        dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
!                        dpdX_eos, dhdX_eos, &
!                        gam1_eos, cs_eos, s_eos, &
!                        dsdt_eos, dsdr_eos, &
!                        do_diag)
               
               gamma(i,j,k) = gam1_eos(1)
               
            end do
         end do
      end do
 
   end subroutine make_gamma_3d

end module make_gamma_module
