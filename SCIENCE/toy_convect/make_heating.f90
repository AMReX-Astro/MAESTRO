! This routine returns the externally imposed (i.e. not reactions)
! heating source term to the enthalpy equation (actually rho * H
! is returned, where H has units of erg/g/s).

module heating_module

  use bl_types
  use define_bc_module
  use bl_error_module

  implicit none
  
  private
  public :: get_rho_Hext
  
contains

  subroutine get_rho_Hext(mla,s,rho_Hext,dx,time,dt)

    use multifab_module
    use ml_layout_module
    use ml_restriction_module
    use geometry, only: dm, nlevs
    use variables, only: foextrap_comp

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time,dt

    ! local
    integer                  :: n,i,ng_s,ng_h
    integer                  :: lo(dm),hi(dm)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: hp(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    call build(bpt, "get_rho_Hext")

    ng_s = s(1)%ng
    ng_h = rho_Hext(1)%ng

    do n=1,nlevs

       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          hp => dataptr(rho_Hext(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (1)
             call bl_error("1-d heating not implemented")
          case (2)
             call get_rho_Hext_2d(hp(:,:,1,1),ng_h,sp(:,:,1,:),ng_s,lo,hi,dx(n,:),time)
          case (3)
             call bl_error("3-d heating not implemented")
          end select
       end do

    end do

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       call ml_cc_restriction(rho_Hext(n-1), rho_Hext(n), mla%mba%rr(n-1,:))
    end do

    call destroy(bpt)

  end subroutine get_rho_Hext
  
  
  subroutine get_rho_Hext_2d(rho_Hext,ng_h,s,ng_s,lo,hi,dx,time)

    use bl_constants_module
    use variables
    use network

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    integer :: i, j
    integer, save :: ih1, ic12, in14, io16

    real(kind=dp_t) :: rho, T6, T613, X_CNO, X_1, g14, eps_CNO

    logical, save :: firstCall = .true.

    if (firstCall) then
       ih1 =  network_species_index("hydrogen-1")
       ic12 = network_species_index("carbon-12")
       in14 = network_species_index("nitrogen-14")
       io16 = network_species_index("oxygen-16")

       firstCall = .false.
    endif

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rho = s(i,j,rho_comp)
          T6 = s(i,j,temp_comp) / 1.0e6_dp_t
          T613 = T6**THIRD

          ! total CNO abundance
          X_CNO = (s(i,j,spec_comp-1+ic12) + &
                   s(i,j,spec_comp-1+in14) + &
                   s(i,j,spec_comp-1+io16)) / rho

          ! H abundance
          X_1 = s(i,j,spec_comp-1+ih1) / rho
          
          ! CNO heating from Kippenhahn & Weigert, Eq. 18.65
          g14 = 1.0_dp_t + 2.7d-3*T613 - 7.78d-3*T613**2 - 1.49d-4*T6
          eps_CNO = 8.67e27_dp_t * g14 * X_CNO * X_1 * rho * exp(-152.28_dp_t/T613) / T613**2
          
          rho_Hext(i,j) = rho * eps_CNO
       enddo
    enddo
    
  end subroutine get_rho_Hext_2d
  
  
end module heating_module
