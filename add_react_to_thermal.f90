module add_react_to_thermal_module

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

   subroutine add_react_to_thermal(thermal,rho_omegadot,s)

      type(multifab) , intent(inout) :: thermal
      type(multifab) , intent(in   ) :: rho_omegadot,s

      real(kind=dp_t), pointer:: sp(:,:,:,:)
      real(kind=dp_t), pointer:: thermalp(:,:,:,:)
      real(kind=dp_t), pointer:: rhowp(:,:,:,:)
      integer :: lo(thermal%dim),hi(thermal%dim),dm
      integer :: i

      dm = thermal%dim

      do i = 1,thermal%nboxes
         if ( multifab_remote(thermal, i) ) cycle
         thermalp => dataptr(thermal,i)
         rhowp => dataptr(rho_omegadot,i)
         sp => dataptr(s,i)
         lo =  lwb(get_box(thermal, i))
         hi =  upb(get_box(thermal, i))
         select case (dm)
            case (2)
              call add_react_to_thermal_2d(lo,hi,thermalp(:,:,1,1), &
                                           rhowp(:,:,1,:),sp(:,:,1,:))
            case (3)

         end select
      end do
      call multifab_fill_boundary(thermal)

    end subroutine add_react_to_thermal

   subroutine add_react_to_thermal_2d(lo,hi,thermal,rho_omegadot,s)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(inout) :: thermal(lo(1):,lo(2):)
      real (kind=dp_t), intent(in  ) :: rho_omegadot(lo(1):,lo(2):,:)
      real (kind=dp_t), intent(in   ) :: s(lo(1)-3:,lo(2)-3:,:)

!     Local variables
      integer :: i, j, n, nr
      real(kind=dp_t) :: react_term

      do_diag = .false.

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           den_eos(1) = s(i,j,rho_comp)
           temp_eos(1) = s(i,j,temp_comp)
           xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
           
           ! dens, temp, and xmass are inputs
           call eos(eos_input_rt, den_eos, temp_eos, &
                    npts, nspec, &
                    xn_eos, &
                    p_eos, h_eos, e_eos, & 
                    cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                    dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                    dpdX_eos, dhdX_eos, &
                    gam1_eos, cs_eos, s_eos, &
                    dsdt_eos, dsdr_eos, &
                    do_diag)

           react_term = ZERO
           do n = 1, nspec
              react_term = react_term - &
                   (dhdX_eos(1,n) + ebin(n))*rho_omegadot(i,j,n)
           enddo

           thermal(i,j) = thermal(i,j) + react_term
        enddo
      enddo
 
   end subroutine add_react_to_thermal_2d

end module add_react_to_thermal_module
