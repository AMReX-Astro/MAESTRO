module rhoh_vs_t_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: makeRhoHfromT, makeTfromRhoH
  
contains
  
  subroutine makeRhoHfromT(nlevs,u,sedge,s0_old,s0_edge_old,s0_new,s0_edge_new, &
                           the_bc_level,dx)

    use bl_prof_module
    use bl_constants_module
    use geometry
    use variables
    use network
    use fill_3d_module, only: fill_3d_data
    use multifab_physbc_module
    
    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,0:,:), s0_edge_old(:,0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(:,0:,:), s0_edge_new(:,0:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    
    ! local
    integer :: i,k,r,dm,n,comp
    integer :: lo(u(1)%dim),hi(u(1)%dim)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer ::  xnp(:,:,:,:)
    real(kind=dp_t), pointer ::  rhp(:,:,:,:)

    real(kind=dp_t), allocatable ::   xn0_halftime(:,:)
    real(kind=dp_t), allocatable :: rhoh0_halftime(:)
    type(multifab)               :: xn0_cart
    type(multifab)               :: rhoh0_cart

    type(bl_prof_timer), save :: bpt

    call build(bpt, "makeRhoHfromT")

    dm = u(1)%dim

    if (spherical .eq. 1) then
      allocate(  xn0_halftime(0:nr(nlevs),nspec))
      allocate(rhoh0_halftime(0:nr(nlevs)))
      do r = 0,nr(nlevs)-1
         xn0_halftime(r,1:nspec) = HALF * (s0_old(nlevs,r,spec_comp:spec_comp+nspec-1) &
                                         + s0_new(nlevs,r,spec_comp:spec_comp+nspec-1) )
         rhoh0_halftime(r) = HALF * (s0_old(nlevs,r,rhoh_comp) + s0_new(nlevs,r,rhoh_comp) )
      end do
   endif

   do n=1,nlevs

      if (spherical .eq. 1) then

         call multifab_build(  xn0_cart,u(n)%la,nspec,2)
         call multifab_build(rhoh0_cart,u(n)%la,1    ,2)

         do i=1,xn0_cart%nboxes
            if ( multifab_remote(u(n),i) ) cycle
            xnp => dataptr(  xn0_cart, i)
            rhp => dataptr(rhoh0_cart, i)
            lo = lwb(get_box(xn0_cart,i))
            hi = upb(get_box(xn0_cart,i))
            do comp = 1,nspec
               call fill_3d_data(n,xnp(:,:,:,comp),xn0_halftime(0:,comp),lo,hi,dx(n,:),xn0_cart%ng)
            end do
            call fill_3d_data(n,rhp(:,:,:,1),rhoh0_halftime(0:),lo,hi,dx(n,:),rhoh0_cart%ng)
         enddo

         call multifab_fill_boundary(  xn0_cart)
         call multifab_fill_boundary(rhoh0_cart)

         call multifab_physbc(xn0_cart,1,dm+spec_comp,nspec,the_bc_level(n))

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
                                   s0_old(n,:,:), s0_edge_old(n,:,:), &
                                   s0_new(n,:,:), s0_edge_new(n,:,:), lo, hi)
          case (3)
             sepz => dataptr(sedge(n,3),i)
             if (spherical .eq. 1) then
               xnp  => dataptr(  xn0_cart, i)
               rhp  => dataptr(rhoh0_cart, i)
               call makeRhoHfromT_3d_sphr(n,sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                          xnp(:,:,:,:), rhp(:,:,:,1), lo, hi, xn0_cart%ng)
             else
               call makeRhoHfromT_3d_cart(sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                          s0_old(n,:,:), s0_edge_old(n,:,:), &
                                          s0_new(n,:,:), s0_edge_new(n,:,:), lo, hi)
             end if
          end select
       end do

    end do

    if (spherical .eq. 1) then
      deallocate(xn0_halftime)
      deallocate(rhoh0_halftime)

      call destroy(xn0_cart)
      call destroy(rhoh0_cart)
    end if

    call destroy(bpt)
    
  end subroutine makeRhoHfromT

  subroutine makeRhoHfromT_2d (sx,sy,s0_old,s0_edge_old,s0_new,s0_edge_new,lo,hi)

    use bl_constants_module
    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use eos_module
    use probin_module, only: use_big_h

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    
    integer :: i, j, comp
    real(kind=dp_t) qreact
    
    do_diag = .false.

    sx(lo(1):hi(1)+1, lo(2):hi(2), rho_comp) = ZERO

    do comp = 1,nspec    
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             sx(i,j,rho_comp) = sx(i,j,rho_comp) + sx(i,j,spec_comp+comp-1)
          end do
       end do
    end do
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          
          temp_eos(1) = sx(i,j,temp_comp)
          den_eos(1)  = sx(i,j,rho_comp) + HALF * (s0_old(j,rho_comp) + s0_new(j,rho_comp))
          xn_eos(1,:) = (sx(i,j,spec_comp:spec_comp+nspec-1)  + &
               HALF * ( s0_old(j,spec_comp:spec_comp+nspec-1) + &
               s0_new(j,spec_comp:spec_comp+nspec-1) ) ) /den_eos(1) 
          
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
          
          sx(i,j,rhoh_comp) = den_eos(1)*h_eos(1)
          
          qreact = 0.0d0
          if (use_big_h) then
             do comp=1,nspec
                qreact = qreact + ebin(comp)*xn_eos(1,comp)
             enddo
             sx(i,j,rhoh_comp) = sx(i,j,rhoh_comp) + den_eos(1) * qreact
          endif
          
          sx(i,j,rhoh_comp) = sx(i,j,rhoh_comp) - &
               HALF * (s0_old(j,rhoh_comp) + s0_new(j,rhoh_comp))
          
       enddo
    enddo

    sy(lo(1):hi(1), lo(2):hi(2)+1, rho_comp) = ZERO

    do comp = 1,nspec    
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             sy(i,j,rho_comp) = sy(i,j,rho_comp) + sy(i,j,spec_comp+comp-1)
          end do
       end do
    end do
    
    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          
          temp_eos(1) = sy(i,j,temp_comp)
          den_eos(1)  = sy(i,j,rho_comp) + &
               HALF * (s0_edge_old(j,rho_comp) + s0_edge_new(j,rho_comp))
          xn_eos(1,:) = (sy(i,j,spec_comp:spec_comp+nspec-1)  + &
               HALF * ( s0_edge_old(j,spec_comp:spec_comp+nspec-1) + &
               s0_edge_new(j,spec_comp:spec_comp+nspec-1) ) ) /den_eos(1) 
          
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
          
          sy(i,j,rhoh_comp) = den_eos(1)*h_eos(1) 
          
          qreact = 0.0d0
          if (use_big_h) then
             do comp=1,nspec
                qreact = qreact + ebin(comp)*xn_eos(1,comp)
             enddo
             sy(i,j,rhoh_comp) = sy(i,j,rhoh_comp) + den_eos(1) * qreact
          endif
          
          sy(i,j,rhoh_comp) = sy(i,j,rhoh_comp) - &
               HALF * (s0_edge_old(j,rhoh_comp) + s0_edge_new(j,rhoh_comp))
          
       enddo
    enddo
    
  end subroutine makeRhoHfromT_2d
  
  subroutine makeRhoHfromT_3d_cart (sx,sy,sz,s0_old,s0_edge_old,s0_new,s0_edge_new,lo,hi)

    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use geometry,      only: spherical
    use eos_module
    use probin_module, only: use_big_h
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sz(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    
    integer :: i, j, k, comp
    real(kind=dp_t) qreact
    
    do_diag = .false.

    sx(lo(1):hi(1)+1, lo(2):hi(2), lo(3):hi(3), rho_comp) = ZERO

    do comp = 1,nspec    
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1
                sx(i,j,k,rho_comp) = sx(i,j,k,rho_comp) + sx(i,j,k,spec_comp+comp-1)
             end do
          end do
       end do
    end do
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             
             temp_eos(1) = sx(i,j,k,temp_comp)
             den_eos(1) = sx(i,j,k,rho_comp) + &
                  HALF * (s0_old(k,rho_comp) + s0_new(k,rho_comp))
             xn_eos(1,:) = (sx(i,j,k,spec_comp:spec_comp+nspec-1)  + &
                  HALF * ( s0_old(k,spec_comp:spec_comp+nspec-1) + &
                  s0_new(k,spec_comp:spec_comp+nspec-1) ) ) /den_eos(1) 
             
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
             
             sx(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             
             qreact = 0.0d0
             if(use_big_h) then
                do comp=1,nspec
                   qreact = qreact + ebin(comp)*xn_eos(1,comp)
                enddo
                sx(i,j,k,rhoh_comp) = sx(i,j,k,rhoh_comp) + den_eos(1) * qreact
             endif
             
             sx(i,j,k,rhoh_comp) = sx(i,j,k,rhoh_comp) - &
                  HALF * (s0_old(k,rhoh_comp) + s0_new(k,rhoh_comp))
             
          enddo
       enddo
    enddo

    sy(lo(1):hi(1), lo(2):hi(2)+1, lo(3):hi(3), rho_comp) = ZERO

    do comp = 1,nspec    
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)
                sy(i,j,k,rho_comp) = sy(i,j,k,rho_comp) + sy(i,j,k,spec_comp+comp-1)
             end do
          end do
       end do
    end do
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             
             temp_eos(1) = sy(i,j,k,temp_comp)
             den_eos(1)  = sy(i,j,k,rho_comp) + &
                  HALF * (s0_old(k,rho_comp) + s0_new(k,rho_comp))
             xn_eos(1,:) = (sy(i,j,k,spec_comp:spec_comp+nspec-1)  + &
                  HALF * ( s0_old(k,spec_comp:spec_comp+nspec-1) + &
                  s0_new(k,spec_comp:spec_comp+nspec-1) ) ) /den_eos(1) 
             
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
             
             sy(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             
             qreact = 0.0d0
             if (use_big_h) then
                do comp=1,nspec
                   qreact = qreact + ebin(comp)*xn_eos(1,comp)
                enddo
                sy(i,j,k,rhoh_comp) = sy(i,j,k,rhoh_comp) + den_eos(1) * qreact
             endif
             
             sy(i,j,k,rhoh_comp) = sy(i,j,k,rhoh_comp) - &
                  HALF * (s0_old(k,rhoh_comp) + s0_new(k,rhoh_comp))
             
          enddo
       enddo
    enddo

    sz(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)+1, rho_comp) = ZERO

    do comp = 1,nspec    
       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                sz(i,j,k,rho_comp) = sz(i,j,k,rho_comp) + sz(i,j,k,spec_comp+comp-1)
             end do
          end do
       end do
    end do
    
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             temp_eos(1) = sz(i,j,k,temp_comp)
             den_eos(1) = sz(i,j,k,rho_comp) + &
                  HALF * (s0_edge_old(k,rho_comp) + s0_edge_new(k,rho_comp))
             xn_eos(1,:) = (sz(i,j,k,spec_comp:spec_comp+nspec-1)  + &
                  HALF * ( s0_edge_old(k,spec_comp:spec_comp+nspec-1) + &
                  s0_edge_new(k,spec_comp:spec_comp+nspec-1) ) ) /den_eos(1)

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
             
             sz(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             
             qreact = 0.0d0
             if(use_big_h) then
                do comp=1,nspec
                   qreact = qreact + ebin(comp)*xn_eos(1,comp)
                enddo
                sz(i,j,k,rhoh_comp) = sz(i,j,k,rhoh_comp) + den_eos(1) * qreact
             endif
             
             sz(i,j,k,rhoh_comp) = sz(i,j,k,rhoh_comp) - &
                  HALF * (s0_edge_old(k,rhoh_comp) + s0_edge_new(k,rhoh_comp))
             
          enddo
       enddo
    enddo
    
  end subroutine makeRhoHfromT_3d_cart

  subroutine makeRhoHfromT_3d_sphr (n,sx,sy,sz,xn0_cart,rhoh0_cart,lo,hi,ngc)

    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use geometry,      only: spherical
    use eos_module
    use probin_module, only: use_big_h
    use bl_constants_module

    integer        , intent(in   ) :: n,ngc
    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sz(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(in   ) ::   xn0_cart(lo(1)-ngc:,lo(2)-ngc:,lo(3)-ngc:,:)
    real(kind=dp_t), intent(in   ) :: rhoh0_cart(lo(1)-ngc:,lo(2)-ngc:,lo(3)-ngc:)
    
    ! Local variables
    integer :: i, j, k, comp
    real(kind=dp_t) qreact
    real(kind=dp_t) xn0_edge(nspec), xn0min(nspec), xn0max(nspec)
    real(kind=dp_t) rhoh0_edge, rhoh0min, rhoh0max
    
    do_diag = .false.
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             
             temp_eos(1) = sx(i,j,k,temp_comp)

             xn0_edge(:) = 7.d0/12.d0 * (xn0_cart(i  ,j,k,:) + xn0_cart(i-1,j,k,:)) &
                          -1.d0/12.d0 * (xn0_cart(i+1,j,k,:) + xn0_cart(i-2,j,k,:))

             xn0min(:) = min(xn0_cart(i,j,k,:),xn0_cart(i-1,j,k,:))
             xn0max(:) = max(xn0_cart(i,j,k,:),xn0_cart(i-1,j,k,:))
             xn0_edge(:) = max(xn0_edge(:),xn0min(:))

             xn0_edge(:) = min(xn0_edge(:),xn0max(:))

             xn_eos(1,:) = sx(i,j,k,spec_comp:spec_comp+nspec-1) + xn0_edge(:)

             den_eos(1) = ZERO
             do comp = 1,nspec
               den_eos(1) = den_eos(1) + xn_eos(1,comp)
             end do

             xn_eos(1,:) = xn_eos(1,:) / den_eos(1)

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
             
             sx(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             
             qreact = 0.0d0
             if(use_big_h) then
                do comp=1,nspec
                   qreact = qreact + ebin(comp)*xn_eos(1,comp)
                enddo
                sx(i,j,k,rhoh_comp) = sx(i,j,k,rhoh_comp) + den_eos(1) * qreact
             endif

             rhoh0_edge = 7.d0/12.d0 * (rhoh0_cart(i  ,j,k) + rhoh0_cart(i-1,j,k)) &
                         -1.d0/12.d0 * (rhoh0_cart(i+1,j,k) + rhoh0_cart(i-2,j,k))
             rhoh0min = min(rhoh0_cart(i,j,k),rhoh0_cart(i-1,j,k))
             rhoh0max = max(rhoh0_cart(i,j,k),rhoh0_cart(i-1,j,k))
             rhoh0_edge = max(rhoh0_edge,rhoh0min)
             rhoh0_edge = min(rhoh0_edge,rhoh0max)
             
             sx(i,j,k,rhoh_comp) = sx(i,j,k,rhoh_comp) - rhoh0_edge
             
          enddo
       enddo
    enddo
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             
             temp_eos(1) = sy(i,j,k,temp_comp)

             xn0_edge(:) = 7.d0/12.d0 * (xn0_cart(i,j  ,k,:) + xn0_cart(i,j-1,k,:)) &
                          -1.d0/12.d0 * (xn0_cart(i,j+1,k,:) + xn0_cart(i,j-2,k,:))
             xn0min(:) = min(xn0_cart(i,j,k,:),xn0_cart(i,j-1,k,:))
             xn0max(:) = max(xn0_cart(i,j,k,:),xn0_cart(i,j-1,k,:))
             xn0_edge(:) = max(xn0_edge(:),xn0min(:))
             xn0_edge(:) = min(xn0_edge(:),xn0max(:))

             xn_eos(1,:) = sy(i,j,k,spec_comp:spec_comp+nspec-1)  + xn0_edge(:)

             den_eos(1) = ZERO
             do comp = 1,nspec
               den_eos(1) = den_eos(1) + xn_eos(1,comp)
             end do

             xn_eos(1,:) = xn_eos(1,:) / den_eos(1)
             
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
             
             sy(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             
             qreact = 0.0d0
             if(use_big_h) then
                do comp=1,nspec
                   qreact = qreact + ebin(comp)*xn_eos(1,comp)
                enddo
                sy(i,j,k,rhoh_comp) = sy(i,j,k,rhoh_comp) + den_eos(1) * qreact
             endif

             rhoh0_edge = 7.d0/12.d0 * (rhoh0_cart(i,j  ,k) + rhoh0_cart(i,j-1,k)) &
                         -1.d0/12.d0 * (rhoh0_cart(i,j+1,k) + rhoh0_cart(i,j-2,k))
             rhoh0min = min(rhoh0_cart(i,j,k),rhoh0_cart(i,j-1,k))
             rhoh0max = max(rhoh0_cart(i,j,k),rhoh0_cart(i,j-1,k))
             rhoh0_edge = max(rhoh0_edge,rhoh0min)
             rhoh0_edge = min(rhoh0_edge,rhoh0max)
             
             sy(i,j,k,rhoh_comp) = sy(i,j,k,rhoh_comp) - rhoh0_edge
             
          enddo
       enddo
    enddo
    
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             temp_eos(1) = sz(i,j,k,temp_comp)

             xn0_edge(:) = 7.d0/12.d0 * (xn0_cart(i,j,k  ,:) + xn0_cart(i,j,k-1,:)) &
                          -1.d0/12.d0 * (xn0_cart(i,j,k+1,:) + xn0_cart(i,j,k-2,:))
             xn0min(:) = min(xn0_cart(i,j,k,:),xn0_cart(i,j,k-1,:))
             xn0max(:) = max(xn0_cart(i,j,k,:),xn0_cart(i,j,k-1,:))
             xn0_edge(:) = max(xn0_edge(:),xn0min(:))
             xn0_edge(:) = min(xn0_edge(:),xn0max(:))

             xn_eos(1,:) = sz(i,j,k,spec_comp:spec_comp+nspec-1) + xn0_edge(:)

             den_eos(1) = ZERO
             do comp = 1,nspec
               den_eos(1) = den_eos(1) + xn_eos(1,comp)
             end do

             xn_eos(1,:) = xn_eos(1,:) / den_eos(1)

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
             
             sz(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             
             qreact = 0.0d0
             if(use_big_h) then
                do comp=1,nspec
                   qreact = qreact + ebin(comp)*xn_eos(1,comp)
                enddo
                sz(i,j,k,rhoh_comp) = sz(i,j,k,rhoh_comp) + den_eos(1) * qreact
             endif

             rhoh0_edge = 7.d0/12.d0 * (rhoh0_cart(i,j,k  ) + rhoh0_cart(i,j,k-1)) &
                         -1.d0/12.d0 * (rhoh0_cart(i,j,k+1) + rhoh0_cart(i,j,k-2))
             rhoh0min = min(rhoh0_cart(i,j,k),rhoh0_cart(i,j,k-1))
             rhoh0max = max(rhoh0_cart(i,j,k),rhoh0_cart(i,j,k-1))
             rhoh0_edge = max(rhoh0_edge,rhoh0min)
             rhoh0_edge = min(rhoh0_edge,rhoh0max)
             
             sz(i,j,k,rhoh_comp) = sz(i,j,k,rhoh_comp) - rhoh0_edge
             
          enddo
       enddo
    enddo
    
  end subroutine makeRhoHfromT_3d_sphr
  
  subroutine makeTfromRhoH(nlevs,s,t0,mla,the_bc_level,dx)

    use variables,             only: temp_comp
    use bl_prof_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_physbc_module
    use multifab_fill_ghost_module

    integer           , intent(in   ) :: nlevs
    type(multifab)    , intent(inout) :: s(:)
    real (kind = dp_t), intent(in   ) :: t0(:,0:)
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
             call makeTfromRhoH_2d(snp(:,:,1,:), lo, hi, ng, t0(n,:))
          case (3)
             call makeTfromRhoH_3d(snp(:,:,:,:), lo, hi, dx(n,:), ng, n, t0(n,:))
          end select
       end do

       call multifab_fill_boundary_c(s(n),temp_comp,1)
       call multifab_physbc(s(n),temp_comp,dm+temp_comp,1,the_bc_level(n))
    end do

    do n=nlevs,2,-1
       call ml_cc_restriction_c(s(n-1),temp_comp,s(n),temp_comp,mla%mba%rr(n-1,:),1)
       call multifab_fill_ghost_cells(s(n),s(n-1), &
            ng,mla%mba%rr(n-1,:), &
            the_bc_level(n-1), &
            the_bc_level(n  ), &
            temp_comp,dm+temp_comp,1)
    enddo

    call destroy(bpt)

  end subroutine makeTfromRhoH

  subroutine makeTfromRhoH_2d (state,lo,hi,ng,t0)

    use variables,     only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use eos_module
    use probin_module, only: use_big_h

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  t0(0:)
    
    ! Local variables
    integer :: i, j, comp
    real(kind=dp_t) qreact
    
    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! (rho, H) --> T, p
          
          den_eos(1)  = state(i,j,rho_comp)
          temp_eos(1) = t0(j)
          xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

          qreact = 0.0d0
          if(use_big_h) then
             do comp=1,nspec
                qreact = qreact + ebin(comp)*xn_eos(1,comp)
             enddo
             h_eos(1) = state(i,j,rhoh_comp) / state(i,j,rho_comp) - qreact
          else
             h_eos(1) = state(i,j,rhoh_comp) / state(i,j,rho_comp)
          endif

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

  end subroutine makeTfromRhoH_2d

  subroutine makeTfromRhoH_3d (state,lo,hi,dx,ng,n,t0)

    use variables,      only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use eos_module
    use probin_module,  only: use_big_h
    use geometry,       only: spherical
    use fill_3d_module, only: fill_3d_data

    integer, intent(in) :: lo(:), hi(:), ng, n
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  t0(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j, k, comp
    real(kind=dp_t) qreact
    real(kind=dp_t), allocatable :: t0_cart(:,:,:)

    if (spherical .eq. 1) then
       allocate(t0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
       call fill_3d_data(n,t0_cart,t0,lo,hi,dx,0)
    endif

    do_diag = .false.
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             ! (rho, H) --> T, p
             
             den_eos(1)  = state(i,j,k,rho_comp)

             if (spherical .eq. 1) then
                temp_eos(1) = t0_cart(i,j,k)
             else
                temp_eos(1) = t0(k)
             endif

             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
             
             qreact = 0.0d0
             if(use_big_h) then
                do comp=1,nspec
                   qreact = qreact + ebin(comp)*xn_eos(1,comp)
                enddo
                h_eos(1) = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - qreact
             else
                h_eos(1) = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
             endif
             
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

    if (spherical .eq. 1) then
       deallocate(t0_cart)
    endif

  end subroutine makeTfromRhoH_3d
  
end module rhoh_vs_t_module
