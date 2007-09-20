module react_state_module

  use bl_types
  use multifab_module
  use define_bc_module
  use setbc_module
  use eos_module
  use network
  use burner_module
  use variables
  use geometry
  use fill_3d_module
  use heating_module
  use probin_module

  implicit none
  
contains

  subroutine react_state (s_in,s_out,rho_omegadot,rho_Hext,t0,dt,dx,the_bc_level,time)

    type(multifab) , intent(in   ) :: s_in
    type(multifab) , intent(inout) :: s_out
    type(multifab) , intent(inout) :: rho_omegadot
    type(multifab) , intent(inout) :: rho_Hext
    real(kind=dp_t), intent(in   ) :: t0(0:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:),time
    type(bc_level) , intent(in   ) :: the_bc_level

    real(kind=dp_t), pointer:: sinp(:,:,:,:)
    real(kind=dp_t), pointer:: sotp(:,:,:,:)
    real(kind=dp_t), pointer::   rp(:,:,:,:)
    real(kind=dp_t), pointer::   hp(:,:,:,:)

    integer :: lo(s_in%dim),hi(s_in%dim),ng,dm
    integer :: i,n,bc_comp

    ng = s_in%ng
    dm = s_in%dim

    do i = 1, s_in%nboxes
       if ( multifab_remote(s_in, i) ) cycle
       sinp => dataptr(s_in , i)
       sotp => dataptr(s_out, i)
         rp => dataptr(rho_omegadot, i)
         hp => dataptr(rho_Hext, i)
       lo =  lwb(get_box(s_in, i))
       hi =  upb(get_box(s_in, i))
       select case (dm)
       case (2)
          call react_state_2d(sinp(:,:,1,:),sotp(:,:,1,:),rp(:,:,1,:), &
               hp(:,:,1,1),t0,dt,dx,lo,hi,ng,time)
       case (3)
          call react_state_3d(sinp(:,:,:,:),sotp(:,:,:,:),rp(:,:,:,:), &
               hp(:,:,:,1),t0,dt,dx,lo,hi,ng,time)
       end select
    end do

    ! Fill ghost cells on periodic boundaries and in between patches
    call multifab_fill_boundary(s_out)

    do i = 1, s_in%nboxes
       if ( multifab_remote(s_in, i) ) cycle
       sotp => dataptr(s_out, i)
       lo =  lwb(get_box(s_in, i))
       hi =  upb(get_box(s_in, i))
       select case (dm)
       case (2)
          ! Impose bc's
          do n = rho_comp,rho_comp+nscal-1
             bc_comp = dm+n 
             call setbc_2d(sotp(:,:,1,n), lo, ng, &
                           the_bc_level%adv_bc_level_array(i,:,:,bc_comp), &
                           dx,bc_comp)
          enddo
       case (3)
          ! Impose bc's
          do n = rho_comp,rho_comp+nscal-1
             bc_comp = dm+n 
             call setbc_3d(sotp(:,:,:,n), lo, ng, &
                           the_bc_level%adv_bc_level_array(i,:,:,bc_comp), &
                           dx,bc_comp)
          enddo
       end select
    end do

  end subroutine react_state

  subroutine react_state_2d (s_in,s_out,rho_omegadot,rho_Hext,t0, &
       dt,dx,lo,hi,ng,time)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) :: s_in (lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(  out) :: s_out(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(  out) :: rho_omegadot(lo(1):,lo(2):,:)
    real (kind = dp_t), intent(  out) :: rho_Hext(lo(1):,lo(2):)
    real (kind = dp_t), intent(in   ) :: t0(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:),time

    !     Local variables
    integer :: i, j, n
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
             do n=1,nspec
                qreact = qreact + x_in(n)*ebin(n)
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
!          input_flag = 2
!
!          den_row(1) = rho
!          h_row(1) = h_out
!          xn_zone(:) = x_out(1:nspec)
!          temp_row(1) = T_in
!
!          call eos(input_flag, den_row, temp_row, &
!               npts, nspec, &
!               xn_zone, aion, zion, &
!               p_row, h_row, e_row, &
!               cv_row, cp_row, xne_row, eta_row, pele_row, &
!               dpdt_row, dpdr_row, dedt_row, dedr_row, &
!               dpdX_row, dhdX_row, &
!               gam1_row, cs_row, s_row, &
!               dsdt_row, dsdr_row, &
!               do_diag)
!
!          s_out(i,j,temp_comp) = temp_row(1)
          s_out(i,j,temp_comp) = s_in(i,j,temp_comp)

       enddo
    enddo

    deallocate(x_in,x_out,rhowdot,H)

  end subroutine react_state_2d

  subroutine react_state_3d (s_in,s_out,rho_omegadot,rho_Hext,t0,dt,dx,lo,hi,ng,time)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) :: s_in (lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) :: s_out(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) :: rho_omegadot(lo(1):,lo(2):,lo(3):,:)
    real (kind = dp_t), intent(  out) :: rho_Hext(lo(1):,lo(2):,lo(3):)
    real (kind = dp_t), intent(in   ) :: t0(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:),time

    real (kind = dp_t), allocatable :: t0_cart(:,:,:)

    !     Local variables
    integer :: i, j, k, n
    real (kind = dp_t), allocatable :: x_in(:),x_out(:),rhowdot(:)
    real (kind = dp_t), allocatable :: H(:,:,:)
    real (kind = dp_t) :: rho,T_in,h_in,h_out,qreact

    allocate(x_in(nspec),x_out(nspec),rhowdot(nspec),H(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call get_H_3d(H,lo,hi,dx,time)

    if (spherical == 1) then
       allocate(t0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
       call fill_3d_data(t0_cart,t0,lo,hi,dx,0)
    endif

    do k = lo(3), hi(3)
     do j = lo(2), hi(2)

       do i = lo(1), hi(1)
          rho = s_in(i,j,k,rho_comp)
          x_in(:) = s_in(i,j,k,spec_comp:spec_comp+nspec-1) / rho

          qreact = 0.0d0
          if(use_big_h) then
             do n=1,nspec
                qreact = qreact + x_in(n)*ebin(n)
             enddo
             h_in = s_in(i,j,k,rhoh_comp) / rho - qreact
          else
             h_in = s_in(i,j,k,rhoh_comp) / rho
          endif
          
          if(spherical == 0) then
             T_in = s_in(i,j,k,temp_comp)
          else
             temp_row(1) = t0_cart(i,j,k)
          endif
          
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
!          input_flag = 2
!
!          den_row(1) = rho
!           h_row(1) = h_out
!          xn_zone(:) = x_out(1:nspec)
!         temp_row(1) = T_in
!
!          call eos(input_flag, den_row, temp_row, &
!               npts, nspec, &
!               xn_zone, aion, zion, &
!               p_row, h_row, e_row, &
!               cv_row, cp_row, xne_row, eta_row, pele_row, &
!               dpdt_row, dpdr_row, dedt_row, dedr_row, &
!               dpdX_row, dhdX_row, &
!               gam1_row, cs_row, s_row, &
!               dsdt_row, dsdr_row, &
!               do_diag)
!
!          s_out(i,j,k,temp_comp) = temp_row(1)
          s_out(i,j,k,temp_comp) = s_in(i,j,k,temp_comp)

       enddo
     enddo
    enddo

    deallocate(x_in,x_out,rhowdot,H)

    if (spherical == 1) then
       deallocate(t0_cart)
    endif

  end subroutine react_state_3d

end module react_state_module
