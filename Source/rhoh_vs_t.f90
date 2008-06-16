module rhoh_vs_t_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: makeRhoHfromT, makeTfromRhoH, makeRhoHfromP, makePfromRhoH
  
contains
  
  subroutine makeRhoHfromT(nlevs,u,sedge,rho0_old,rhoh0_old,rho0_edge_old,rhoh0_edge_old, &
                           rho0_new,rhoh0_new,rho0_edge_new,rhoh0_edge_new,the_bc_level,dx)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical, nr_fine, r_end_coord
    use variables
    use network
    use fill_3d_module
    use multifab_physbc_module
    
    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:),      rhoh0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(:,0:), rhoh0_edge_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:),      rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(:,0:), rhoh0_edge_new(:,0:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    
    ! local
    integer :: i,r,dm,n
    integer :: lo(u(1)%dim),hi(u(1)%dim)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer ::   rp(:,:,:,:)
    real(kind=dp_t), pointer ::  rhp(:,:,:,:)

    real(kind=dp_t), allocatable ::  rho0_halftime(:)
    real(kind=dp_t), allocatable :: rhoh0_halftime(:)
    type(multifab)               ::  rho0_cart
    type(multifab)               :: rhoh0_cart

    type(bl_prof_timer), save :: bpt

    call build(bpt, "makeRhoHfromT")

    dm = u(1)%dim

    if (spherical .eq. 1) then
      allocate( rho0_halftime(0:nr_fine-1))
      allocate(rhoh0_halftime(0:nr_fine-1))
      do r=0,r_end_coord(nlevs)
         rho0_halftime(r)  = HALF * (rho0_old(nlevs,r) + &
                                     rho0_new(nlevs,r) )

         rhoh0_halftime(r) = HALF * (rhoh0_old(nlevs,r) + &
                                     rhoh0_new(nlevs,r) )
      end do
   endif

   do n=1,nlevs

      if (spherical .eq. 1) then

         call multifab_build( rho0_cart,u(n)%la,1,2)
         call multifab_build(rhoh0_cart,u(n)%la,1,2)

         do i=1,rho0_cart%nboxes
            if ( multifab_remote(u(n),i) ) cycle
            rp  => dataptr( rho0_cart, i)
            rhp => dataptr(rhoh0_cart, i)
            lo = lwb(get_box(rho0_cart,i))
            hi = upb(get_box(rho0_cart,i))

            call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,rho0_halftime,rp, &
                                              lo,hi,dx(n,:),2)
            call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,rhoh0_halftime,rhp, &
                                              lo,hi,dx(n,:),2)
         enddo

         ! fill ghost cells for two adjacent grids at the same level
         ! this includes periodic domain boundary ghost cells
         call multifab_fill_boundary( rho0_cart)
         call multifab_fill_boundary(rhoh0_cart)

         ! fill non-periodic domain boundary ghost cells
         call multifab_physbc( rho0_cart,1,dm+rho_comp ,1,the_bc_level(n))
         call multifab_physbc(rhoh0_cart,1,dm+rhoh_comp,1,the_bc_level(n))

      endif

       do i=1,u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          sepx => dataptr(sedge(n,1), i)
          sepy => dataptr(sedge(n,2), i)
          lo = lwb(get_box(u(n),i))
          hi = upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call makeRhoHfromT_2d(sepx(:,:,1,:), sepy(:,:,1,:), &
                                   rho0_old(n,:), rhoh0_old(n,:), &
                                   rho0_edge_old(n,:), rhoh0_edge_old(n,:), &
                                   rho0_new(n,:), rhoh0_new(n,:), &
                                   rho0_edge_new(n,:), rhoh0_edge_new(n,:), &
                                   lo, hi)
          case (3)
             sepz => dataptr(sedge(n,3),i)
             if (spherical .eq. 1) then
               rp   => dataptr( rho0_cart, i)
               rhp  => dataptr(rhoh0_cart, i)
               call makeRhoHfromT_3d_sphr(sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                          rp(:,:,:,1), rhp(:,:,:,1), lo, hi, rho0_cart%ng)
             else
               call makeRhoHfromT_3d_cart(sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                          rho0_old(n,:), rhoh0_old(n,:), &
                                          rho0_edge_old(n,:), rhoh0_edge_old(n,:), &
                                          rho0_new(n,:), rhoh0_new(n,:), &
                                          rho0_edge_new(n,:), rhoh0_edge_new(n,:), &
                                          lo, hi)
             end if
          end select
       end do

    end do

    if (spherical .eq. 1) then
      deallocate( rho0_halftime)
      deallocate(rhoh0_halftime)

      call destroy( rho0_cart)
      call destroy(rhoh0_cart)
    end if

    call destroy(bpt)
    
  end subroutine makeRhoHfromT

  subroutine makeRhoHfromT_2d(sx,sy,rho0_old,rhoh0_old,rho0_edge_old,rhoh0_edge_old, &
                              rho0_new,rhoh0_new,rho0_edge_new,rhoh0_edge_new,lo,hi)

    use bl_constants_module
    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use eos_module
    use probin_module, only: enthalpy_pred_type, small_temp, predict_rho
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:),      rhoh0_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(0:), rhoh0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:),      rhoh0_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(0:), rhoh0_edge_new(0:)
 
    integer :: i, j
    
    do_diag = .false.
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          
          temp_eos(1) = max(sx(i,j,temp_comp),small_temp)

          if (predict_rho) then
             den_eos(1)  = sx(i,j,rho_comp)
          else
             den_eos(1)  = sx(i,j,rho_comp) + HALF * (rho0_old(j) + rho0_new(j))
          end if

          ! sx(i,j,spec_comp:spec_comp+nspec-1) holds X
          xn_eos(1,:) = sx(i,j,spec_comp:spec_comp+nspec-1)
          
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
          
          if (enthalpy_pred_type .eq. predict_T_then_h) then
             sx(i,j,rhoh_comp) = h_eos(1)
          else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
             sx(i,j,rhoh_comp) = den_eos(1)*h_eos(1)
          end if
          
          if (enthalpy_pred_type .eq. predict_T_then_rhohprime) &
             sx(i,j,rhoh_comp) = sx(i,j,rhoh_comp) - &
                  HALF * (rhoh0_old(j) + rhoh0_new(j))
          
       enddo
    enddo
    
    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          
          temp_eos(1) = max(sy(i,j,temp_comp),small_temp)

          if (predict_rho) then
             den_eos(1)  = sy(i,j,rho_comp)
          else
             den_eos(1)  = sy(i,j,rho_comp) + HALF * (rho0_edge_old(j) + rho0_edge_new(j))
          end if
          
          ! sy(i,j,spec_comp:spec_comp+nspec-1) holds X
          xn_eos(1,:) = sy(i,j,spec_comp:spec_comp+nspec-1)
          
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
          
          if (enthalpy_pred_type .eq. predict_T_then_h) then
             sy(i,j,rhoh_comp) = h_eos(1) 
          else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
             sy(i,j,rhoh_comp) = den_eos(1)*h_eos(1) 
          end if
          
          if (enthalpy_pred_type .eq. predict_T_then_rhohprime) &
             sy(i,j,rhoh_comp) = sy(i,j,rhoh_comp) - &
                  HALF * (rhoh0_edge_old(j) + rhoh0_edge_new(j))
          
       enddo
    enddo
    
  end subroutine makeRhoHfromT_2d
  
  subroutine makeRhoHfromT_3d_cart(sx,sy,sz,rho0_old,rhoh0_old,rho0_edge_old, &
                                   rhoh0_edge_old,rho0_new,rhoh0_new,rho0_edge_new, &
                                   rhoh0_edge_new,lo,hi)

    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use eos_module
    use probin_module, only: enthalpy_pred_type, small_temp, predict_rho
    use pred_parameters
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sz(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:),      rhoh0_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(0:), rhoh0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:),      rhoh0_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(0:), rhoh0_edge_new(0:)
    
    integer :: i, j, k
    
    do_diag = .false.
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             
             temp_eos(1) = max(sx(i,j,k,temp_comp),small_temp)

             if (predict_rho) then
                den_eos(1) = sx(i,j,k,rho_comp)
             else
                den_eos(1) = sx(i,j,k,rho_comp) + &
                     HALF * (rho0_old(k) + rho0_new(k))
             end if

             ! sx(i,j,k,spec_comp:spec_comp+nspec-1) holds X
             xn_eos(1,:) = sx(i,j,k,spec_comp:spec_comp+nspec-1)
             
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
             
             if (enthalpy_pred_type .eq. predict_T_then_h) then
                sx(i,j,k,rhoh_comp) = h_eos(1)
             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                sx(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             end if
             
             if (enthalpy_pred_type .eq. predict_T_then_rhohprime) &
                sx(i,j,k,rhoh_comp) = sx(i,j,k,rhoh_comp) - &
                     HALF * (rhoh0_old(k) + rhoh0_new(k))
             
          enddo
       enddo
    enddo
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             
             temp_eos(1) = max(sy(i,j,k,temp_comp),small_temp)

             if (predict_rho) then
                den_eos(1)  = sy(i,j,k,rho_comp)
             else
                den_eos(1)  = sy(i,j,k,rho_comp) + &
                     HALF * (rho0_old(k) + rho0_new(k))
             end if

             ! sy(i,j,k,spec_comp:spec_comp+nspec-1) holds X
             xn_eos(1,:) = sy(i,j,k,spec_comp:spec_comp+nspec-1)
             
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
             
             if (enthalpy_pred_type .eq. predict_T_then_h) then
                sy(i,j,k,rhoh_comp) = h_eos(1)
             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                sy(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             end if
             
             if (enthalpy_pred_type .eq. predict_T_then_rhohprime) &
                sy(i,j,k,rhoh_comp) = sy(i,j,k,rhoh_comp) - &
                     HALF * (rhoh0_old(k) + rhoh0_new(k))
             
          enddo
       enddo
    enddo

    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             temp_eos(1) = max(sz(i,j,k,temp_comp),small_temp)

             if (predict_rho) then
                den_eos(1) = sz(i,j,k,rho_comp)
             else
                den_eos(1) = sz(i,j,k,rho_comp) + &
                     HALF * (rho0_edge_old(k) + rho0_edge_new(k))
             end if

             ! sz(i,j,k,spec_comp:spec_comp+nspec-1) X
             xn_eos(1,:) = sz(i,j,k,spec_comp:spec_comp+nspec-1)

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
             
             if (enthalpy_pred_type .eq. predict_T_then_h) then
                sz(i,j,k,rhoh_comp) = h_eos(1)
             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                sz(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             end if
             
             if (enthalpy_pred_type .eq. predict_T_then_rhohprime) &
                sz(i,j,k,rhoh_comp) = sz(i,j,k,rhoh_comp) - &
                     HALF * (rhoh0_edge_old(k) + rhoh0_edge_new(k))
             
          enddo
       enddo
    enddo
    
  end subroutine makeRhoHfromT_3d_cart

  subroutine makeRhoHfromT_3d_sphr(sx,sy,sz,rho0_cart,rhoh0_cart,lo,hi,ngc)

    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use geometry,      only: spherical
    use eos_module
    use probin_module, only: enthalpy_pred_type, small_temp, predict_rho
    use pred_parameters
    use bl_constants_module

    integer        , intent(in   ) :: ngc
    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sz(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(in   ) ::  rho0_cart(lo(1)-ngc:,lo(2)-ngc:,lo(3)-ngc:)
    real(kind=dp_t), intent(in   ) :: rhoh0_cart(lo(1)-ngc:,lo(2)-ngc:,lo(3)-ngc:)
    
    ! Local variables
    integer :: i, j, k
    real(kind=dp_t) rho0_edge, rho0min, rho0max
    real(kind=dp_t) rhoh0_edge, rhoh0min, rhoh0max
    
    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             
             temp_eos(1) = max(sx(i,j,k,temp_comp),small_temp)

             if (predict_rho) then
                den_eos(1) = sx(i,j,k,rho_comp)
             else
                rho0_edge = 7.d0/12.d0 * (rho0_cart(i  ,j,k) + rho0_cart(i-1,j,k)) &
                     -1.d0/12.d0 * (rho0_cart(i+1,j,k) + rho0_cart(i-2,j,k))
                
                rho0min = min(rho0_cart(i,j,k),rho0_cart(i-1,j,k))
                rho0max = max(rho0_cart(i,j,k),rho0_cart(i-1,j,k))
                
                rho0_edge = max(rho0_edge, rho0min)
                rho0_edge = min(rho0_edge, rho0max)

                den_eos(1) = sx(i,j,k,rho_comp) + rho0_edge
             end if

             ! sx(i,j,k,spec_comp:spec_comp+nspec-1) holds X
             xn_eos(1,:) = sx(i,j,k,spec_comp:spec_comp+nspec-1)

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
             
             if (enthalpy_pred_type .eq. predict_T_then_h) then
                sx(i,j,k,rhoh_comp) = h_eos(1)
             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                sx(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             end if
             
             if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                rhoh0_edge = 7.d0/12.d0 * (rhoh0_cart(i  ,j,k) + rhoh0_cart(i-1,j,k)) &
                            -1.d0/12.d0 * (rhoh0_cart(i+1,j,k) + rhoh0_cart(i-2,j,k))
                rhoh0min = min(rhoh0_cart(i,j,k),rhoh0_cart(i-1,j,k))
                rhoh0max = max(rhoh0_cart(i,j,k),rhoh0_cart(i-1,j,k))
                rhoh0_edge = max(rhoh0_edge,rhoh0min)
                rhoh0_edge = min(rhoh0_edge,rhoh0max)
             
                sx(i,j,k,rhoh_comp) = sx(i,j,k,rhoh_comp) - rhoh0_edge
             end if
             
          enddo
       enddo
    enddo

    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             
             temp_eos(1) = max(sy(i,j,k,temp_comp),small_temp)

             if (predict_rho) then
                den_eos(1) = sy(i,j,k,rho_comp)
             else
                rho0_edge = 7.d0/12.d0 * (rho0_cart(i,j  ,k) + rho0_cart(i,j-1,k)) &
                     -1.d0/12.d0 * (rho0_cart(i,j+1,k) + rho0_cart(i,j-2,k))

                rho0min = min(rho0_cart(i,j,k),rho0_cart(i,j-1,k))
                rho0max = max(rho0_cart(i,j,k),rho0_cart(i,j-1,k))
                
                rho0_edge = max(rho0_edge, rho0min)
                rho0_edge = min(rho0_edge, rho0max)
                
                den_eos(1) = sy(i,j,k,rho_comp) + rho0_edge
             end if

             ! sy(i,j,k,spec_comp:spec_comp+nspec-1) holds X
             xn_eos(1,:) = sy(i,j,k,spec_comp:spec_comp+nspec-1) 

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
             
             if (enthalpy_pred_type .eq. predict_T_then_h) then
                sy(i,j,k,rhoh_comp) = h_eos(1)
             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                sy(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             end if
             
             if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                rhoh0_edge = 7.d0/12.d0 * (rhoh0_cart(i,j  ,k) + rhoh0_cart(i,j-1,k)) &
                            -1.d0/12.d0 * (rhoh0_cart(i,j+1,k) + rhoh0_cart(i,j-2,k))
                rhoh0min = min(rhoh0_cart(i,j,k),rhoh0_cart(i,j-1,k))
                rhoh0max = max(rhoh0_cart(i,j,k),rhoh0_cart(i,j-1,k))
                rhoh0_edge = max(rhoh0_edge,rhoh0min)
                rhoh0_edge = min(rhoh0_edge,rhoh0max)

                sy(i,j,k,rhoh_comp) = sy(i,j,k,rhoh_comp) - rhoh0_edge
             end if
             
          enddo
       enddo
    enddo

    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             temp_eos(1) = max(sz(i,j,k,temp_comp),small_temp)

             if (predict_rho) then
                den_eos(1) = sz(i,j,k,rho_comp)
             else
                rho0_edge = 7.d0/12.d0 * (rho0_cart(i,j,k  ) + rho0_cart(i,j,k-1)) &
                     -1.d0/12.d0 * (rho0_cart(i,j,k+1) + rho0_cart(i,j,k-2))
                
                rho0min = min(rho0_cart(i,j,k),rho0_cart(i,j,k-1))
                rho0max = max(rho0_cart(i,j,k),rho0_cart(i,j,k-1))
                
                rho0_edge = max(rho0_edge, rho0min)
                rho0_edge = min(rho0_edge, rho0max)
                
                den_eos(1) = sz(i,j,k,rho_comp) + rho0_edge
             end if
             
             ! sz(i,j,k,spec_comp:spec_comp+nspec-1) holds X
             xn_eos(1,:) = sz(i,j,k,spec_comp:spec_comp+nspec-1) 

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
             
             if (enthalpy_pred_type .eq. predict_T_then_h) then
                sz(i,j,k,rhoh_comp) = h_eos(1)
             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                sz(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             end if
             
             if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                rhoh0_edge = 7.d0/12.d0 * (rhoh0_cart(i,j,k  ) + rhoh0_cart(i,j,k-1)) &
                            -1.d0/12.d0 * (rhoh0_cart(i,j,k+1) + rhoh0_cart(i,j,k-2))
                rhoh0min = min(rhoh0_cart(i,j,k),rhoh0_cart(i,j,k-1))
                rhoh0max = max(rhoh0_cart(i,j,k),rhoh0_cart(i,j,k-1))
                rhoh0_edge = max(rhoh0_edge,rhoh0min)
                rhoh0_edge = min(rhoh0_edge,rhoh0max)

                sz(i,j,k,rhoh_comp) = sz(i,j,k,rhoh_comp) - rhoh0_edge
             end if
             
          enddo
       enddo
    enddo
    
  end subroutine makeRhoHfromT_3d_sphr
  
  subroutine makeTfromRhoH(nlevs,s,p0,tempbar,mla,the_bc_level,dx)

    use variables,             only: temp_comp
    use bl_prof_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_physbc_module
    use multifab_fill_ghost_module

    integer           , intent(in   ) :: nlevs
    type(multifab)    , intent(inout) :: s(:)
    real (kind = dp_t), intent(in   ) :: p0(:,0:)
    real (kind = dp_t), intent(in   ) :: tempbar(:,0:)
    type(ml_layout)   , intent(inout) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t)   , intent(in   ) :: dx(:,:)

    ! local
    integer                  :: i,ng,dm,n
    integer                  :: lo(s(1)%dim),hi(s(1)%dim)
    real(kind=dp_t), pointer :: snp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "makeTfromRhoH")

    dm = s(1)%dim
    ng = s(1)%ng

    do n=1,nlevs

       do i=1,s(n)%nboxes
          if (multifab_remote(s(n),i)) cycle
          snp => dataptr(s(n),i)
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call makeTfromRhoH_2d(snp(:,:,1,:), lo, hi, ng, p0(n,:), tempbar(n,:))
          case (3)
             call makeTfromRhoH_3d(snp(:,:,:,:), lo, hi, ng, p0(n,:), tempbar(n,:), dx(n,:), n)
          end select
       end do

    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s(nlevs),temp_comp,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),temp_comp,dm+temp_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(s(n-1),temp_comp,s(n),temp_comp,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n  ), &
                                         temp_comp,dm+temp_comp,1)
       enddo

    end if

    call destroy(bpt)

  end subroutine makeTfromRhoH

  subroutine makeTfromRhoH_2d (state,lo,hi,ng,p0,tempbar)

    use variables,     only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use eos_module
    use probin_module, only: use_tfromp

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:)
    real (kind = dp_t), intent(in   ) ::  tempbar(0:)
    
    ! Local variables
    integer :: i, j
    
    do_diag = .false.


    if (.not. use_tfromp) then

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, H) --> T, p
          
             den_eos(1)  = state(i,j,rho_comp)
             temp_eos(1) = tempbar(j)
             xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
             
             h_eos(1) = state(i,j,rhoh_comp) / state(i,j,rho_comp)

             call eos(eos_input_rh, den_eos, temp_eos, &
                      npts, nspec, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)

             state(i,j,temp_comp) = temp_eos(1)

          enddo
       enddo

    else

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, p) --> T
          
             den_eos(1)  = state(i,j,rho_comp)
             temp_eos(1) = tempbar(j)
             xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
             
             p_eos(1) = p0(j)

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

             state(i,j,temp_comp) = temp_eos(1)

          enddo
       enddo

    endif

  end subroutine makeTfromRhoH_2d

  subroutine makeTfromRhoH_3d(state,lo,hi,ng,p0,tempbar,dx,n)

    use variables,      only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use eos_module
    use geometry,       only: spherical
    use fill_3d_module
    use probin_module, only: use_tfromp

    integer, intent(in) :: lo(:), hi(:), ng, n
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:)
    real (kind = dp_t), intent(in   ) ::  tempbar(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j, k
    real(kind=dp_t), allocatable :: p0_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: tempbar_cart(:,:,:,:)

    if (spherical .eq. 1) then
       allocate(tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,tempbar,tempbar_cart,lo,hi,dx,0)

       allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,p0,p0_cart,lo,hi,dx,0)
    endif

    do_diag = .false.


    if (.not. use_tfromp) then
    
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                
                ! (rho, H) --> T, p
             
                den_eos(1)  = state(i,j,k,rho_comp)

                if (spherical .eq. 1) then
                   temp_eos(1) = tempbar_cart(i,j,k,1)
                else
                   temp_eos(1) = tempbar(k)
                endif

                xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
                h_eos(1) = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
                
                call eos(eos_input_rh, den_eos, temp_eos, &
                         npts, nspec, &
                         xn_eos, &
                         p_eos, h_eos, e_eos, &
                         cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                         dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                         dpdX_eos, dhdX_eos, &
                         gam1_eos, cs_eos, s_eos, &
                         dsdt_eos, dsdr_eos, &
                         do_diag)
                
                state(i,j,k,temp_comp) = temp_eos(1)
             
             enddo
          enddo
       enddo

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                
                ! (rho, p) --> T
             
                den_eos(1)  = state(i,j,k,rho_comp)

                if (spherical .eq. 1) then
                   temp_eos(1) = tempbar_cart(i,j,k,1)
                   p_eos(1) = p0_cart(i,j,k,1)
                else
                   temp_eos(1) = tempbar(k)
                   p_eos(1) = p0(k)
                endif

                xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
                
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
                
                state(i,j,k,temp_comp) = temp_eos(1)
             
             enddo
          enddo
       enddo

    endif

       
    if (spherical .eq. 1) then
       deallocate(tempbar_cart)
       deallocate(p0_cart)
    endif

  end subroutine makeTfromRhoH_3d

  subroutine makePfromRhoH(nlevs,s,p,tempbar,mla,the_bc_level,dx)

    use variables,             only: foextrap_comp
    use bl_prof_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_physbc_module
    use multifab_fill_ghost_module

    integer           , intent(in   ) :: nlevs
    type(multifab)    , intent(in   ) :: s(:)
    type(multifab)    , intent(inout) :: p(:)
    real (kind = dp_t), intent(in   ) :: tempbar(:,0:)
    type(ml_layout)   , intent(inout) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t)   , intent(in   ) :: dx(:,:)

    ! local
    integer                  :: i,ng_s,ng_p,dm,n
    integer                  :: lo(s(1)%dim),hi(s(1)%dim)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: pnp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "makePfromRhoH")

    dm = s(1)%dim
    ng_s = s(1)%ng
    ng_p = p(1)%ng

    do n=1,nlevs

       do i=1,s(n)%nboxes
          if (multifab_remote(s(n),i)) cycle
          snp => dataptr(s(n),i)
          pnp => dataptr(p(n),i)
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call makePfromRhoH_2d(snp(:,:,1,:), pnp(:,:,1,1), lo, hi, ng_s, ng_p, &
                                   tempbar(n,:))
          case (3)
             call makePfromRhoH_3d(snp(:,:,:,:), pnp(:,:,:,1), lo, hi, ng_s, ng_p, &
                                   tempbar(n,:), dx(n,:), n)
          end select
       end do

    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(p(nlevs),1,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(p(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(p(n-1),1,p(n),1,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(p(n),p(n-1),ng_p,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n),1,foextrap_comp,1)
       enddo

    end if

    call destroy(bpt)

  end subroutine makePfromRhoH

  subroutine makePfromRhoH_2d(state,p,lo,hi,ng_s,ng_p,tempbar)

    use variables,     only: rho_comp, spec_comp, rhoh_comp
    use eos_module

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_p
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind = dp_t), intent(inout) ::     p(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(in   ) ::  tempbar(0:)
    
    ! Local variables
    integer :: i, j
    
    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! (rho, H) --> T, p
          
          den_eos(1)  = state(i,j,rho_comp)
          temp_eos(1) = tempbar(j)
          xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

          h_eos(1) = state(i,j,rhoh_comp) / state(i,j,rho_comp)

          call eos(eos_input_rh, den_eos, temp_eos, &
                   npts, nspec, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   do_diag)

          p(i,j) = p_eos(1)

       enddo
    enddo

  end subroutine makePfromRhoH_2d

  subroutine makePfromRhoH_3d(state,p,lo,hi,ng_s,ng_p,tempbar,dx,n)

    use variables,      only: rho_comp, spec_comp, rhoh_comp
    use eos_module
    use geometry,       only: spherical
    use fill_3d_module

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_p, n
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(inout) ::     p(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tempbar(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j, k
    real(kind=dp_t), allocatable :: tempbar_cart(:,:,:,:)

    if (spherical .eq. 1) then
       allocate(tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,tempbar,tempbar_cart,lo,hi,dx,0)
    endif

    do_diag = .false.
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             ! (rho, H) --> T, p
             
             den_eos(1)  = state(i,j,k,rho_comp)

             if (spherical .eq. 1) then
                temp_eos(1) = tempbar_cart(i,j,k,1)
             else
                temp_eos(1) = tempbar(k)
             endif

             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
             h_eos(1) = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
             
             call eos(eos_input_rh, den_eos, temp_eos, &
                      npts, nspec, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)
             
             p(i,j,k) = p_eos(1)
             
          enddo
       enddo
    enddo

    if (spherical .eq. 1) then
       deallocate(tempbar_cart)
    endif

  end subroutine makePfromRhoH_3d

  subroutine makeRhoHfromP(nlevs,u,sedge, &
                           rho0_old,rho0_edge_old,&
                           rho0_new,rho0_edge_new, &
                             p0_old, p0_new, the_bc_level,dx)

    use bl_prof_module
    use bl_constants_module
    use geometry
    use variables
    use network
    use fill_3d_module
    use multifab_physbc_module
    
    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(:,0:)
    real(kind=dp_t), intent(in   ) ::   p0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::   p0_new(:,0:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    
    ! local
    integer :: i,r,dm,n
    integer :: lo(u(1)%dim),hi(u(1)%dim)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer ::   rp(:,:,:,:)
    real(kind=dp_t), pointer ::  rhp(:,:,:,:)

    real(kind=dp_t), allocatable ::  rho0_halftime(:)
    real(kind=dp_t), allocatable :: rhoh0_halftime(:)
    type(multifab)               ::  rho0_cart
    type(multifab)               :: rhoh0_cart

    type(bl_prof_timer), save :: bpt

    call build(bpt, "makeRhoHfromP")

    dm = u(1)%dim

   do n=1,nlevs

       do i=1,u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          sepx => dataptr(sedge(n,1), i)
          sepy => dataptr(sedge(n,2), i)
          lo = lwb(get_box(u(n),i))
          hi = upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call makeRhoHfromP_2d(n, sepx(:,:,1,:), sepy(:,:,1,:), &
                                   rho0_old(n,:), rho0_edge_old(n,:), &
                                   rho0_new(n,:), rho0_edge_new(n,:), &
                                     p0_old(n,:), p0_new(n,:), lo, hi, dx(n,:))
          case (3)
             if (spherical .eq. 1) then
                print *,'NO MAKERHOHFROMP FOR SPHERICAL '
                stop
             else
                call makeRhoHfromP_3d(n, sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                      rho0_old(n,:), rho0_edge_old(n,:), &
                                      rho0_new(n,:), rho0_edge_new(n,:), &
                                        p0_old(n,:), p0_new(n,:), lo, hi, dx(n,:))
             end if
          end select
       end do

    end do

    if (spherical .eq. 1) then
      deallocate( rho0_halftime)
      deallocate(rhoh0_halftime)

      call destroy( rho0_cart)
      call destroy(rhoh0_cart)
    end if

    call destroy(bpt)
    
  end subroutine makeRhoHfromP

  subroutine makeRhoHfromP_2d(n,sx,sy,rho0_old,rho0_edge_old,&
                                      rho0_new,rho0_edge_new,p0_old,p0_new,lo,hi,dx)

    use bl_constants_module
    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use eos_module
    use probin_module, only: enthalpy_pred_type, small_temp, predict_rho, grav_const, &
         base_cutoff_density
    use geometry, only: nr

    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:),n
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) ::   p0_old(0:)
    real(kind=dp_t), intent(in   ) ::   p0_new(0:)
    real(kind=dp_t), intent(in   ) ::   dx(:)
 
    integer         :: i, j
    real(kind=dp_t) :: p0_old_edge, p0_new_edge

    do_diag = .false.
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          
          temp_eos(1) = max(sx(i,j,temp_comp),small_temp)
           den_eos(1) = sx(i,j,rho_comp) + HALF * (rho0_old(j) + rho0_new(j))
             p_eos(1) = HALF * (p0_old(j) + p0_new(j))

          ! sx(i,j,spec_comp:spec_comp+nspec-1) holds X
          xn_eos(1,:) = sx(i,j,spec_comp:spec_comp+nspec-1)
          
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
          
          sx(i,j,rhoh_comp) = h_eos(1)
          
       enddo
    enddo

    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          
          temp_eos(1) = max(sy(i,j,temp_comp),small_temp)
           den_eos(1) = sy(i,j,rho_comp) + HALF * (rho0_edge_old(j) + rho0_edge_new(j))

          if (j .eq. nr(n)+1) then
             p0_old_edge = p0_old(j-1)
             p0_new_edge = p0_new(j-1)
          else if (rho0_edge_old(j) .lt. base_cutoff_density) then
             p0_old_edge = p0_old(j)
             p0_new_edge = p0_new(j)
          else 
             p0_old_edge = p0_old(j) + HALF * rho0_old(j) * abs(grav_const) * dx(2)
             p0_new_edge = p0_new(j) + HALF * rho0_new(j) * abs(grav_const) * dx(2)
          end if

          p_eos(1) = HALF * (p0_old_edge + p0_new_edge)
          
          ! sy(i,j,spec_comp:spec_comp+nspec-1) holds X
          xn_eos(1,:) = sy(i,j,spec_comp:spec_comp+nspec-1)
          
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
          
          sy(i,j,rhoh_comp) = h_eos(1) 
          
       enddo
    enddo
    
  end subroutine makeRhoHfromP_2d

  subroutine makeRhoHfromP_3d(n,sx,sy,sz,rho0_old,rho0_edge_old,&
                                         rho0_new,rho0_edge_new,p0_old,p0_new,lo,hi,dx)

    use bl_constants_module
    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use eos_module
    use probin_module, only: enthalpy_pred_type, small_temp, predict_rho, grav_const, &
         base_cutoff_density
    use geometry, only: nr

    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:),n
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sz(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(0:)
    real(kind=dp_t), intent(in   ) ::   p0_old(0:)
    real(kind=dp_t), intent(in   ) ::   p0_new(0:)
    real(kind=dp_t), intent(in   ) ::   dx(:)
 
    integer         :: i, j, k
    real(kind=dp_t) :: p0_old_edge, p0_new_edge
    
    do_diag = .false.
    
    do k = lo(2), hi(2)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
          
          temp_eos(1) = max(sx(i,j,k,temp_comp),small_temp)
           den_eos(1) = sx(i,j,k,rho_comp) + HALF * (rho0_old(k) + rho0_new(k))
             p_eos(1) = HALF * (p0_old(k) + p0_new(k))

          ! sx(i,j,spec_comp:spec_comp+nspec-1) holds X
          xn_eos(1,:) = sx(i,j,k,spec_comp:spec_comp+nspec-1)
          
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
          
          sx(i,j,k,rhoh_comp) = h_eos(1)
          
        enddo
      enddo
    enddo
    
    do k = lo(2), hi(2)
      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
          
          temp_eos(1) = max(sy(i,j,k,temp_comp),small_temp)
           den_eos(1) = sy(i,j,k,rho_comp) + HALF * (rho0_old(k) + rho0_new(k))
             p_eos(1) = HALF * (p0_old(k) + p0_new(k))

          ! sy(i,j,spec_comp:spec_comp+nspec-1) holds X
          xn_eos(1,:) = sy(i,j,k,spec_comp:spec_comp+nspec-1)
          
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
          
          sy(i,j,k,rhoh_comp) = h_eos(1)
          
        enddo
      enddo
    enddo

    do k = lo(3), hi(3)+1
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          
          temp_eos(1) = max(sz(i,j,k,temp_comp),small_temp)
           den_eos(1) = sz(i,j,k,rho_comp) + HALF * (rho0_edge_old(k) + rho0_edge_new(k))

          if (k .eq. nr(n)+1) then
             p0_old_edge = p0_old(k-1)
             p0_new_edge = p0_new(k-1)
          else if (rho0_edge_old(k) .lt. base_cutoff_density) then
             p0_old_edge = p0_old(k)
             p0_new_edge = p0_new(k)
          else 
             p0_old_edge = p0_old(k) + HALF * rho0_old(k) * abs(grav_const) * dx(3)
             p0_new_edge = p0_new(k) + HALF * rho0_new(k) * abs(grav_const) * dx(3)
          end if

          p_eos(1) = HALF * (p0_old_edge + p0_new_edge)
          
          ! sz(i,j,spec_comp:spec_comp+nspec-1) holds X
          xn_eos(1,:) = sz(i,j,k,spec_comp:spec_comp+nspec-1)
          
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
          
          sz(i,j,k,rhoh_comp) = h_eos(1) 
          
        enddo
      enddo
    enddo
    
  end subroutine makeRhoHfromP_3d
  
end module rhoh_vs_t_module
