module add_react_to_thermal_module 

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: add_react_to_thermal

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_react_to_thermal(nlevs,thermal,rho_omegadot,s,the_bc_level,mla,dx)

    use bl_constants_module
    use eos_module
    use variables
    use fill_3d_module
    use network
    use geometry
    use ml_restriction_module, only : ml_cc_restriction
    use multifab_physbc_module
    use multifab_fill_ghost_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: thermal(:)
    type(multifab) , intent(in   ) :: rho_omegadot(:)
    type(multifab) , intent(in   ) :: s(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: tp(:,:,:,:)
    real(kind=dp_t), pointer :: rwp(:,:,:,:)
    integer                  :: lo(thermal(1)%dim),hi(thermal(1)%dim)
    integer                  :: dm,i,n
    
    dm = thermal(1)%dim
    
    do n=1,nlevs
       do i=1,thermal(n)%nboxes
          if ( multifab_remote(thermal(n),i) ) cycle
          tp  => dataptr(thermal(n),i)
          rwp => dataptr(rho_omegadot(n),i)
          sp  => dataptr(s(n),i)
          lo = lwb(get_box(thermal(n), i))
          hi = upb(get_box(thermal(n), i))
          select case (dm)
          case (2)
             call add_react_to_thermal_2d(lo,hi,tp(:,:,1,1),rwp(:,:,1,:),sp(:,:,1,:))
          case (3)
             call add_react_to_thermal_3d(lo,hi,tp(:,:,:,1),rwp(:,:,:,:),sp(:,:,:,:))
          end select
       end do

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary conditions
       call multifab_fill_boundary(thermal(n))
       
       call multifab_physbc(thermal(n),1,foextrap_comp,1,dx(n,:),the_bc_level(n))
    enddo

    do n=nlevs,2,-1
       ! make sure that coarse cells are the average of the fine cells covering it.
       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       call ml_cc_restriction(thermal(n-1),thermal(n),mla%mba%rr(n-1,:))

       ! fill fine ghost cells using interpolation from the underlying coarse data
       call multifab_fill_ghost_cells(thermal(n),thermal(n-1), &
                                      1,mla%mba%rr(n-1,:), &
                                      the_bc_level(n-1), the_bc_level(n), &
                                      1,foextrap_comp,1)
    enddo
       
  end subroutine add_react_to_thermal
  
  subroutine add_react_to_thermal_2d(lo,hi,thermal,rho_omegadot,s)

    use variables
    use eos_module
    use bl_constants_module
    
    implicit none
    
    integer         , intent(in   ) :: lo(:),hi(:)
    real (kind=dp_t), intent(inout) :: thermal(lo(1)-1:,lo(2)-1:)
    real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1):,lo(2):,:)
    real (kind=dp_t), intent(in   ) :: s(lo(1)-3:,lo(2)-3:,:)
    
    ! Local variables
    integer         :: i,j,comp
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
          do comp = 1, nspec
             react_term = react_term - &
                  (dhdX_eos(1,comp) + ebin(comp))*rho_omegadot(i,j,comp)
          enddo
          
          thermal(i,j) = thermal(i,j) + react_term
       enddo
    enddo
    
  end subroutine add_react_to_thermal_2d

  subroutine add_react_to_thermal_3d(lo,hi,thermal,rho_omegadot,s)

    use variables
    use eos_module
    use bl_constants_module
    
    implicit none
    
    integer         , intent(in   ) :: lo(:),hi(:)
    real (kind=dp_t), intent(inout) :: thermal(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1):,lo(2):,lo(3):,:)
    real (kind=dp_t), intent(in   ) :: s(lo(1)-3:,lo(2)-3:,lo(3):,:)
    
    ! Local variables
    integer :: i,j,k,comp
    real(kind=dp_t) :: react_term
    
    do_diag = .false.
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
          
             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = s(i,j,k,temp_comp)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
             
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
             do comp = 1, nspec
                react_term = react_term - &
                     (dhdX_eos(1,comp) + ebin(comp))*rho_omegadot(i,j,k,comp)
             enddo
             
             thermal(i,j,k) = thermal(i,j,k) + react_term
          enddo
       enddo
    enddo

  end subroutine add_react_to_thermal_3d
  
end module add_react_to_thermal_module
