! compute gamma1 for the full state.

module make_gamma_module

  use bl_types
  use multifab_module
  use ml_layout_module
  implicit none

  private

  public :: make_gamma

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_gamma(nlevs,mla,gamma,s,p0,tempbar,dx,the_bc_level)

    use variables, only: foextrap_comp
    use bl_prof_module
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: gamma(:)
    type(multifab) , intent(in   ) :: s(:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: tempbar(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    real(kind=dp_t), pointer:: gamp(:,:,:,:),sp(:,:,:,:)
    integer :: lo(s(1)%dim),hi(s(1)%dim),dm
    integer :: i,n
    integer :: ng_g, ng_s

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_gamma")
    
    dm = s(1)%dim

    ng_g = gamma(1)%ng
    ng_s = s(1)%ng

    do n = 1, nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          gamp => dataptr(gamma(n), i)
          sp   => dataptr(s(n), i)
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
          select case (dm)
          case (2)
             call make_gamma_2d(lo,hi,gamp(:,:,1,1),ng_g,sp(:,:,1,:),ng_s,p0(n,:),tempbar(n,:))
          case (3)
             call make_gamma_3d(n,lo,hi,gamp(:,:,:,1),ng_g,sp(:,:,:,:),ng_s,p0(n,:),tempbar(n,:),dx(n,:))
          end select
       end do
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(gamma(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(gamma(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(gamma(n-1)    ,gamma(n)    ,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(gamma(n),gamma(n-1), &
                                         ng_g,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         1,foextrap_comp,1)
       enddo

    end if

    call destroy(bpt)

   end subroutine make_gamma

   subroutine make_gamma_2d(lo,hi,gamma,ng_g,s,ng_s,p0,tempbar)

      use eos_module
      use variables, only: rho_comp, rhoh_comp, spec_comp

      integer         , intent(in   ) :: lo(:), hi(:), ng_g, ng_s
      real (kind=dp_t), intent(  out) :: gamma(lo(1)-ng_g:,lo(2)-ng_g:)
      real (kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s:,lo(2)-ng_s:,:)
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

    subroutine make_gamma_3d(n,lo,hi,gamma,ng_g,s,ng_s,p0,tempbar,dx)

      use eos_module
      use variables, only: rho_comp, rhoh_comp, spec_comp
      use fill_3d_module
      use geometry, only: spherical

      integer         , intent(in   ) :: n
      integer         , intent(in   ) :: lo(:), hi(:), ng_g, ng_s
      real (kind=dp_t), intent(  out) :: gamma(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real (kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
      real (kind=dp_t), intent(in   ) :: p0(0:)
      real (kind=dp_t), intent(in   ) :: tempbar(0:)
      real (kind=dp_t), intent(in   ) :: dx(:)

      ! local variables
      integer :: i, j, k

      real (kind=dp_t), allocatable :: tempbar_cart(:,:,:,:)
      real (kind=dp_t), allocatable :: p0_cart(:,:,:,:)

      if (spherical .eq. 1) then
         allocate(tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
         call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,tempbar,tempbar_cart, &
                                           lo,hi,dx,0,0)
         
         allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
         call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,p0,p0_cart, &
                                           lo,hi,dx,0,0)
      endif

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               
               den_eos(1) = s(i,j,k,rho_comp)

               if (spherical .eq. 1) then
                  p_eos(1) = p0_cart(i,j,k,1)
                  temp_eos(1) = tempbar_cart(i,j,k,1)
               else
                  p_eos(1) = p0(k)
                  temp_eos(1) = tempbar(k)
               endif

               xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

               
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

      if (spherical .eq. 1) then
         deallocate(tempbar_cart, p0_cart)
      endif

   end subroutine make_gamma_3d

end module make_gamma_module
