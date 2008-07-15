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

  subroutine add_react_to_thermal(nlevs,thermal,rho_omegadot,s,rho_Hext,the_bc_level,mla, &
                                  dx,time)

    use bl_prof_module
    use bl_constants_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use variables, only: foextrap_comp
    use ml_restriction_module, only : ml_cc_restriction
    use heating_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: thermal(:)
    type(multifab) , intent(in   ) :: rho_omegadot(:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    real(kind=dp_t), intent(in   ) :: dx(:,:),time

    ! local
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: tp(:,:,:,:)
    real(kind=dp_t), pointer :: rwp(:,:,:,:)
    real(kind=dp_t), pointer :: hep(:,:,:,:)
    integer                  :: lo(thermal(1)%dim),hi(thermal(1)%dim)
    integer                  :: dm,i,n
    integer :: ng_t, ng_rw, ng_s, ng_he

    type(bl_prof_timer), save :: bpt

    call build(bpt, "add_react_to_thermal")
    
    dm = thermal(1)%dim

    ng_t = thermal(1)%ng
    ng_rw = rho_omegadot(1)%ng
    ng_s = s(1)%ng
    ng_he = rho_Hext(1)%ng

    call get_rho_Hext(nlevs,mla,s,rho_Hext,dx,time)
    
    do n=1,nlevs
       do i=1,thermal(n)%nboxes
          if ( multifab_remote(thermal(n),i) ) cycle
          tp  => dataptr(thermal(n),i)
          rwp => dataptr(rho_omegadot(n),i)
          sp  => dataptr(s(n),i)
          hep => dataptr(rho_Hext(n), i)
          lo = lwb(get_box(thermal(n), i))
          hi = upb(get_box(thermal(n), i))
          select case (dm)
          case (2)
             call add_react_to_thermal_2d(lo,hi,tp(:,:,1,1),ng_t,rwp(:,:,1,:),ng_rw, &
                                          sp(:,:,1,:),ng_s,hep(:,:,1,1),ng_he)
          case (3)
             call add_react_to_thermal_3d(lo,hi,tp(:,:,:,1),ng_t,rwp(:,:,:,:),ng_rw, &
                                          sp(:,:,:,:),ng_s,hep(:,:,:,1),ng_he)
          end select
       end do
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(thermal(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(thermal(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(thermal(n-1),thermal(n),mla%mba%rr(n-1,:))
          
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(thermal(n),thermal(n-1),1,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n),1, &
                                         foextrap_comp,1)
       enddo

    end if

    call destroy(bpt)
       
  end subroutine add_react_to_thermal
  
  subroutine add_react_to_thermal_2d(lo,hi,thermal,ng_t,rho_omegadot,ng_rw, &
                                     s,ng_s,rho_Hext,ng_he)

    use eos_module
    use bl_constants_module
    use variables, only: temp_comp, rho_comp, spec_comp
    
    integer         , intent(in   ) :: lo(:),hi(:),ng_t,ng_rw,ng_s,ng_he
    real (kind=dp_t), intent(inout) ::      thermal(lo(1)-ng_t :,lo(2)-ng_t :)
    real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,:)
    real (kind=dp_t), intent(in   ) ::            s(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real (kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:)
    
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
          
          thermal(i,j) = thermal(i,j) + react_term + rho_Hext(i,j)
       enddo
    enddo
    
  end subroutine add_react_to_thermal_2d

  subroutine add_react_to_thermal_3d(lo,hi,thermal,ng_t,rho_omegadot,ng_rw, &
                                     s,ng_s,rho_Hext,ng_he)

    use eos_module
    use bl_constants_module
    use variables, only: spec_comp, rho_comp, temp_comp
    
    integer         , intent(in   ) :: lo(:),hi(:),ng_t,ng_rw,ng_s,ng_he
    real (kind=dp_t), intent(inout) ::      thermal(lo(1)-ng_t :,lo(2)-ng_t :,lo(3)-ng_t :)
    real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,lo(3)-ng_rw:,:)
    real (kind=dp_t), intent(in   ) ::            s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real (kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:,lo(3)-ng_he:)
    
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
             
             thermal(i,j,k) = thermal(i,j,k) + react_term + rho_Hext(i,j,k)
          enddo
       enddo
    enddo

  end subroutine add_react_to_thermal_3d
  
end module add_react_to_thermal_module
