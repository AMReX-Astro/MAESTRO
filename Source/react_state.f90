module react_state_module

  use bl_types
  use multifab_module
  use eos_module
  use network
  use variables

  implicit none
  
contains

  subroutine react_state (s_in,s_out,rho_omegadot,dt)

    type(multifab) , intent(in   ) :: s_in
    type(multifab) , intent(inout) :: s_out
    type(multifab) , intent(inout) :: rho_omegadot
    real(kind=dp_t), intent(in   ) :: dt

    real(kind=dp_t), pointer:: sinp(:,:,:,:)
    real(kind=dp_t), pointer:: sotp(:,:,:,:)
    real(kind=dp_t), pointer::   rp(:,:,:,:)

    integer :: lo(s_in%dim),hi(s_in%dim),ng,dm
    integer :: i

    ng = s_in%ng
    dm = s_in%dim

    print *,"<<< react state >>> " 

    do i = 1, s_in%nboxes
       if ( multifab_remote(s_in, i) ) cycle
       sinp => dataptr(s_in , i)
       sotp => dataptr(s_out, i)
         rp => dataptr(rho_omegadot, i)
       lo =  lwb(get_box(s_in, i))
       hi =  upb(get_box(s_in, i))
       select case (dm)
       case (2)
          call react_state_2d(sinp(:,:,1,:),sotp(:,:,1,:),rp(:,:,1,:),dt,lo,hi,ng)
       case (3)
          call react_state_3d(sinp(:,:,:,:),sotp(:,:,:,:),rp(:,:,:,:),dt,lo,hi,ng)
       end select
    end do

  end subroutine react_state

  subroutine react_state_2d (s_in,s_out,rho_omegadot,dt,lo,hi,ng)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) :: s_in (lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(  out) :: s_out(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(  out) :: rho_omegadot(lo(1):,lo(2):,:)
    real (kind = dp_t), intent(in   ) :: dt

    !     Local variables
    integer :: i, j
    real (kind = dp_t), allocatable :: x_in(:),x_out(:),rhowdot(:)
    real (kind = dp_t) :: rho,T_in,h_in,h_out

    allocate(x_in(nspec),x_out(nspec),rhowdot(nspec))

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          rho = s_in(i,j,rho_comp)
          x_in(:) = s_in(i,j,spec_comp:spec_comp+nspec-1) / rho
          h_in = s_in(i,j,rhoh_comp) / rho

          ! (rho, H) --> T, p
          input_flag = 2
          xn_zone(:) = x_in(:)

          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)

          T_in = temp_row(1)

          call burner(rho, T_in, x_in, h_in, dt, x_out, h_out, rhowdot)

          s_out(i,j,rho_comp) = s_in(i,j,rho_comp)
          s_out(i,j,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho
          s_out(i,j,rhoh_comp) = rho * h_out
          rho_omegadot(i,j,1:nspec) = rhowdot(1:nspec)

       enddo
    enddo

    deallocate(x_in,x_out,rhowdot)

  end subroutine react_state_2d

  subroutine react_state_3d (s_in,s_out,rho_omegadot,dt,lo,hi,ng)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) :: s_in (lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) :: s_out(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) :: rho_omegadot(lo(1):,lo(2):,lo(3):,:)
    real (kind = dp_t), intent(in   ) :: dt

    !     Local variables
    integer :: i, j, k
    real (kind = dp_t), allocatable :: x_in(:),x_out(:),rhowdot(:)
    real (kind = dp_t) :: rho,T_in,h_in,h_out

    allocate(x_in(nspec),x_out(nspec),rhowdot(nspec))

    do k = lo(3), hi(3)
     do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          rho = s_in(i,j,k,rho_comp)
          x_in(:) = s_in(i,j,k,spec_comp:spec_comp+nspec-1) / rho
          h_in = s_in(i,j,k,rhoh_comp) / rho

          ! (rho, H) --> T, p
          input_flag = 2
          xn_zone(:) = x_in(:)

          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)

          T_in = temp_row(1)

          call burner(rho, T_in, x_in, h_in, dt, x_out, h_out, rhowdot)

          s_out(i,j,k,rho_comp) = s_in(i,j,k,rho_comp)
          s_out(i,j,k,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho
          s_out(i,j,k,rhoh_comp) = rho * h_out
          rho_omegadot(i,j,k,1:nspec) = rhowdot(1:nspec)

       enddo
     enddo
    enddo

    deallocate(x_in,x_out,rhowdot)

  end subroutine react_state_3d

end module react_state_module
