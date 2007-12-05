module react_state_module

  use bl_types
  use multifab_module
  use define_bc_module
  use multifab_physbc_module
  use eos_module
  use network
  use burner_module
  use variables
  use geometry
  use fill_3d_module
  use heating_module
  use probin_module, ONLY: use_big_h
  use ml_restriction_module
  use ml_layout_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: react_state

contains

  subroutine react_state (nlevs,mla,s_in,s_out,rho_omegadot,rho_Hext,dt,dx,the_bc_level,time)

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(in   ) :: s_in(:)
    type(multifab) , intent(inout) :: s_out(:)
    type(multifab) , intent(inout) :: rho_omegadot(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:,:),time
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local
    real(kind=dp_t), pointer:: sinp(:,:,:,:)
    real(kind=dp_t), pointer:: sotp(:,:,:,:)
    real(kind=dp_t), pointer::   rp(:,:,:,:)
    real(kind=dp_t), pointer::   hp(:,:,:,:)

    integer :: lo(s_in(1)%dim),hi(s_in(1)%dim),ng,dm
    integer :: i,n,bc_comp,comp

    ng = s_in(1)%ng
    dm = s_in(1)%dim

    do n = 1, nlevs
       
       do i = 1, s_in(n)%nboxes
          if ( multifab_remote(s_in(n), i) ) cycle
          sinp => dataptr(s_in(n) , i)
          sotp => dataptr(s_out(n), i)
          rp => dataptr(rho_omegadot(n), i)
          hp => dataptr(rho_Hext(n), i)
          lo =  lwb(get_box(s_in(n), i))
          hi =  upb(get_box(s_in(n), i))
          select case (dm)
          case (2)
             call react_state_2d(sinp(:,:,1,:),sotp(:,:,1,:),rp(:,:,1,:), &
                                 hp(:,:,1,1),dt,dx(n,:),lo,hi,ng,time)
          case (3)
             call react_state_3d(sinp(:,:,:,:),sotp(:,:,:,:),rp(:,:,:,:), &
                                 hp(:,:,:,1),dt,dx(n,:),lo,hi,ng,time)
          end select
       end do

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary conditions
       call multifab_fill_boundary(s_out(n))
       
       ! fill physical boundary conditions at domain boundaries
       call multifab_physbc(s_out(n),rho_comp,dm+rho_comp,nscal,dx(n,:),the_bc_level(n))
       
    enddo

    do n=nlevs,2,-1
       ! make sure that coarse cells are the average of the fine cells covering it.
       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       call ml_cc_restriction(s_out(n-1),s_out(n),mla%mba%rr(n-1,:))

       ! fill fine ghost cells using interpolation from the underlying coarse data
       call multifab_fill_ghost_cells(s_out(n),s_out(n-1), &
                                      ng,mla%mba%rr(n-1,:), &
                                      the_bc_level(n-1), the_bc_level(n), &
                                      rho_comp,dm+rho_comp,nscal)
    enddo

  end subroutine react_state

  subroutine react_state_2d(s_in,s_out,rho_omegadot,rho_Hext, &
                            dt,dx,lo,hi,ng,time)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) :: s_in (lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(  out) :: s_out(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(  out) :: rho_omegadot(lo(1):,lo(2):,:)
    real (kind = dp_t), intent(  out) :: rho_Hext(lo(1):,lo(2):)
    real (kind = dp_t), intent(in   ) :: dt,dx(:),time

    !     Local variables
    integer :: i, j, comp
    real (kind = dp_t), allocatable :: x_in(:),x_out(:),rhowdot(:)
    real (kind = dp_t), allocatable :: H(:,:)
    real (kind = dp_t) :: rho,T_in,h_in,h_out,qreact

    allocate(x_in(nspec),x_out(nspec),rhowdot(nspec),H(lo(1):hi(1),lo(2):hi(2)))
    call get_H_2d(H,lo,hi,dx,time)

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          rho = s_in(i,j,rho_comp)
          x_in(:) = s_in(i,j,spec_comp:spec_comp+nspec-1) / rho

          qreact = 0.0d0
          if(use_big_h) then
             do comp = 1, nspec
                qreact = qreact + x_in(comp)*ebin(comp)
             enddo
             h_in = s_in(i,j,rhoh_comp) / rho - qreact
          else
             h_in = s_in(i,j,rhoh_comp) / rho
          endif
          
          T_in = s_in(i,j,temp_comp)
          
          call burner(rho, T_in, x_in, h_in, dt, x_out, h_out, rhowdot)
          
          s_out(i,j,rho_comp) = s_in(i,j,rho_comp)
          s_out(i,j,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho
          
          rho_Hext(i,j) = s_in(i,j,rho_comp) * H(i,j)
          rho_omegadot(i,j,1:nspec) = rhowdot(1:nspec)
          
          if(use_big_h) then
             s_out(i,j,rhoh_comp) = s_in(i,j,rhoh_comp) + dt * rho_Hext(i,j)
          else
             s_out(i,j,rhoh_comp) = rho * h_out + dt * rho_Hext(i,j)
          endif

          ! now compute temperature and put it into s_out
          ! dens, enthalpy, and xmass are inputs
!
!          den_eos(1) = rho
!          h_eos(1) = h_out
!          xn_eos(1,:) = x_out(1:nspec)
!          temp_eos(1) = T_in
!
!          call eos(eos_input_rh, den_eos, temp_eos, &
!               npts, nspec, &
!               xn_eos, &
!               p_eos, h_eos, e_eos, &
!               cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
!               dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
!               dpdX_eos, dhdX_eos, &
!               gam1_eos, cs_eos, s_eos, &
!               dsdt_eos, dsdr_eos, &
!               do_diag)
!
!          s_out(i,j,temp_comp) = temp_eos(1)
          s_out(i,j,temp_comp) = s_in(i,j,temp_comp)

       enddo
    enddo

    deallocate(x_in,x_out,rhowdot,H)

  end subroutine react_state_2d

  subroutine react_state_3d(s_in,s_out,rho_omegadot,rho_Hext,dt,dx,lo,hi,ng,time)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) :: s_in (lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) :: s_out(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) :: rho_omegadot(lo(1):,lo(2):,lo(3):,:)
    real (kind = dp_t), intent(  out) :: rho_Hext(lo(1):,lo(2):,lo(3):)
    real (kind = dp_t), intent(in   ) :: dt,dx(:),time

    !     Local variables
    integer :: i, j, k, comp
    real (kind = dp_t), allocatable :: x_in(:),x_out(:),rhowdot(:)
    real (kind = dp_t), allocatable :: H(:,:,:)
    real (kind = dp_t) :: rho,T_in,h_in,h_out,qreact

    allocate(x_in(nspec),x_out(nspec),rhowdot(nspec),H(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call get_H_3d(H,lo,hi,dx,time)

    do k = lo(3), hi(3)
     do j = lo(2), hi(2)

       do i = lo(1), hi(1)
          rho = s_in(i,j,k,rho_comp)
          x_in(:) = s_in(i,j,k,spec_comp:spec_comp+nspec-1) / rho

          qreact = 0.0d0
          if(use_big_h) then
             do comp = 1, nspec
                qreact = qreact + x_in(comp)*ebin(comp)
             enddo
             h_in = s_in(i,j,k,rhoh_comp) / rho - qreact
          else
             h_in = s_in(i,j,k,rhoh_comp) / rho
          endif
          
          T_in = s_in(i,j,k,temp_comp)
          
          call burner(rho, T_in, x_in, h_in, dt, x_out, h_out, rhowdot)
          
          s_out(i,j,k,rho_comp) = s_in(i,j,k,rho_comp)
          s_out(i,j,k,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho
          
          rho_Hext(i,j,k) = s_in(i,j,k,rho_comp) * H(i,j,k)
          rho_omegadot(i,j,k,1:nspec) = rhowdot(1:nspec)
          
          if(use_big_h) then
             s_out(i,j,k,rhoh_comp) = s_in(i,j,k,rhoh_comp) + dt * rho_Hext(i,j,k)
          else
             s_out(i,j,k,rhoh_comp) = rho * h_out + dt * rho_Hext(i,j,k)
          endif
  
          ! now compute temperature and put it into s_out
          ! dens, enthalpy, and xmass are inputs
!
!          den_eos(1) = rho
!          h_eos(1) = h_out
!          xn_eos(1,:) = x_out(1:nspec)
!          temp_eos(1) = T_in
!
!          call eos(eos_input_rh, den_eos, temp_eos, &
!               npts, nspec, &
!               xn_eos, &
!               p_eos, h_eos, e_eos, &
!               cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
!               dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
!               dpdX_eos, dhdX_eos, &
!               gam1_eos, cs_eos, s_eos, &
!               dsdt_eos, dsdr_eos, &
!               do_diag)
!
!          s_out(i,j,k,temp_comp) = temp_eos(1)
          s_out(i,j,k,temp_comp) = s_in(i,j,k,temp_comp)

       enddo
     enddo
    enddo

    deallocate(x_in,x_out,rhowdot,H)

  end subroutine react_state_3d

end module react_state_module
