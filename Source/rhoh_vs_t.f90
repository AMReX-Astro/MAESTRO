module rhoh_vs_t_module

  use bl_types
  use bl_constants_module
  use variables
  use geometry
  use network
  use probin_module, ONLY: use_big_h
  use eos_module
  use multifab_module

  implicit none

  private
  public :: makeRhoHfromT
  public :: makeTfromRhoH
  
contains
  
  subroutine makeRhoHfromT(u,sedge,s0_old,s0_edge_old,s0_new,s0_edge_new)
    
    implicit none
    type(multifab) , intent(in   ) :: u
    type(multifab) , intent(inout) :: sedge(:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    
    ! local
    integer :: i,dm
    integer :: lo(u%dim),hi(u%dim)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    
    dm = u%dim
    
    do i=1,u%nboxes
       if ( multifab_remote(u,i) ) cycle
       sepx => dataptr(sedge(1), i)
       sepy => dataptr(sedge(2), i)
       lo = lwb(get_box(u,i))
       hi = upb(get_box(u,i))
       select case (dm)
       case (2)
          call makeRhoHfromT_2d(sepx(:,:,1,:), sepy(:,:,1,:), &
                                s0_old, s0_edge_old, s0_new, s0_edge_new, lo, hi)
       case (3)
          sepz => dataptr(sedge(3),i)
          call makeRhoHfromT_3d(sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                s0_old, s0_edge_old, s0_new, s0_edge_new, lo, hi)
       end select
    end do
    
  end subroutine makeRhoHfromT

  subroutine makeRhoHfromT_2d (sx,sy,s0_old,s0_edge_old,s0_new,s0_edge_new,lo,hi)
    
    implicit none
    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    
    ! Local variables
    integer :: i, j, n
    real(kind=dp_t) qreact
    
    do_diag = .false.
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          sx(i,j,rho_comp) = 0.d0
          do n = 1,nspec
             sx(i,j,rho_comp) = sx(i,j,rho_comp) + sx(i,j,spec_comp+n-1)
          end do
       end do
    end do
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          
          temp_eos(1) = sx(i,j,temp_comp)
          den_eos(1) = sx(i,j,rho_comp) + HALF * (s0_old(j,rho_comp) + s0_new(j,rho_comp))
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
          if(use_big_h) then
             do n=1,nspec
                qreact = qreact + ebin(n)*xn_eos(1,n)
             enddo
             sx(i,j,rhoh_comp) = sx(i,j,rhoh_comp) + den_eos(1) * qreact
          endif
          
          sx(i,j,rhoh_comp) = sx(i,j,rhoh_comp) - &
               HALF * (s0_old(j,rhoh_comp) + s0_new(j,rhoh_comp))
          
       enddo
    enddo
    
    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          sy(i,j,rho_comp) = 0.d0
          do n = 1,nspec
             sy(i,j,rho_comp) = sy(i,j,rho_comp) + sy(i,j,spec_comp+n-1)
          end do
       end do
    end do
    
    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          
          temp_eos(1) = sy(i,j,temp_comp)
          den_eos(1) = sy(i,j,rho_comp) + &
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
          if(use_big_h) then
             do n=1,nspec
                qreact = qreact + ebin(n)*xn_eos(1,n)
             enddo
             sy(i,j,rhoh_comp) = sy(i,j,rhoh_comp) + den_eos(1) * qreact
          endif
          
          sy(i,j,rhoh_comp) = sy(i,j,rhoh_comp) - &
               HALF * (s0_edge_old(j,rhoh_comp) + s0_edge_new(j,rhoh_comp))
          
       enddo
    enddo
    
  end subroutine makeRhoHfromT_2d
  
  subroutine makeRhoHfromT_3d (sx,sy,sz,s0_old,s0_edge_old,s0_new,s0_edge_new,lo,hi)
    
    implicit none
    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sz(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    
    ! Local variables
    integer :: i, j, k, n
    real(kind=dp_t) qreact
    
    do_diag = .false.
    
    if (spherical .eq. 1) then
       print *,'MAKERHOHFROMT_3D NOT YET SET UP FOR SPHERICAL '
       stop
    end if
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             sx(i,j,k,rho_comp) = 0.d0
             do n = 1,nspec
                sx(i,j,k,rho_comp) = sx(i,j,k,rho_comp) + sx(i,j,k,spec_comp+n-1)
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
                do n=1,nspec
                   qreact = qreact + ebin(n)*xn_eos(1,n)
                enddo
                sx(i,j,k,rhoh_comp) = sx(i,j,k,rhoh_comp) + sx(i,j,k,rho_comp) * qreact
             endif
             
             sx(i,j,k,rhoh_comp) = sx(i,j,k,rhoh_comp) - &
                  HALF * (s0_old(k,rhoh_comp) + s0_new(k,rhoh_comp))
             
          enddo
       enddo
    enddo
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             sy(i,j,k,rho_comp) = 0.d0
             do n = 1,nspec
                sy(i,j,k,rho_comp) = sy(i,j,k,rho_comp) + sy(i,j,k,spec_comp+n-1)
             end do
          end do
       end do
    end do
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             
             temp_eos(1) = sy(i,j,k,temp_comp)
             den_eos(1) = sy(i,j,k,rho_comp) + &
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
             if(use_big_h) then
                do n=1,nspec
                   qreact = qreact + ebin(n)*xn_eos(1,n)
                enddo
                sy(i,j,k,rhoh_comp) = sy(i,j,k,rhoh_comp) + sy(i,j,k,rho_comp) * qreact
             endif
             
             sy(i,j,k,rhoh_comp) = sy(i,j,k,rhoh_comp) - &
                  HALF * (s0_old(k,rhoh_comp) + s0_new(k,rhoh_comp))
             
          enddo
       enddo
    enddo
    
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             sz(i,j,k,rho_comp) = 0.d0
             do n = 1,nspec
                sz(i,j,k,rho_comp) = sz(i,j,k,rho_comp) + sz(i,j,k,spec_comp+n-1)
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
                do n=1,nspec
                   qreact = qreact + ebin(n)*xn_eos(1,n)
                enddo
                sz(i,j,k,rhoh_comp) = sz(i,j,k,rhoh_comp) + sz(i,j,k,rho_comp) * qreact
             endif
             
             sz(i,j,k,rhoh_comp) = sz(i,j,k,rhoh_comp) - &
                  HALF * (s0_edge_old(k,rhoh_comp) + s0_edge_new(k,rhoh_comp))
             
          enddo
       enddo
    enddo
    
  end subroutine makeRhoHfromT_3d
  
  subroutine makeTfromRhoH(s,t0)

    implicit none
    type(multifab)    , intent(inout) :: s
    real (kind = dp_t), intent(in   ) :: t0(0:)

    ! local
    integer                  :: i,ng,dm
    integer                  :: lo(s%dim),hi(s%dim)
    real(kind=dp_t), pointer :: snp(:,:,:,:)

    dm = s%dim
    ng = s%ng

    do i=1,s%nboxes
       if (multifab_remote(s,i)) cycle
       snp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          call makeTfromRhoH_2d(snp(:,:,1,:), lo, hi, 3, t0)
       case (3)
          call makeTfromRhoH_3d(snp(:,:,:,:), lo, hi, 3, t0)
       end select
    end do
    
  end subroutine makeTfromRhoH

  subroutine makeTfromRhoH_2d (state,lo,hi,ng,t0)
    
    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  t0(0:)
    
    ! Local variables
    integer :: i, j, n
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
             do n=1,nspec
                qreact = qreact + ebin(n)*xn_eos(1,n)
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

  subroutine makeTfromRhoH_3d (state,lo,hi,ng,t0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  t0(0:)

    ! Local variables
    integer :: i, j, k, n
    real(kind=dp_t) qreact
    
    do_diag = .false.
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             ! (rho, H) --> T, p
             
             den_eos(1)  = state(i,j,k,rho_comp)
             temp_eos(1) = t0(k)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
             
             qreact = 0.0d0
             if(use_big_h) then
                do n=1,nspec
                   qreact = qreact + ebin(n)*xn_eos(1,n)
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
    
  end subroutine makeTfromRhoH_3d
  
end module rhoh_vs_t_module
