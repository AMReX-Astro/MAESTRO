! This routine returns the externally imposed (i.e. not reactions)
! heating source term to the enthalpy equation (actually rho * H
! is returned, where H has units of erg/g/s).

module heating_module

  use bl_types
  use bl_error_module

  implicit none
  
  private
  public :: get_rho_Hext
  
contains

  subroutine get_rho_Hext(mla,tempbar_init,s,rho_Hext,the_bc_level,dx,dt)

    use multifab_module
    use define_bc_module
    use ml_layout_module
    use ml_cc_restriction_module, only : ml_cc_restriction
    use variables, only: foextrap_comp

    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: tempbar_init(:,0:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt

    ! local
    integer                  :: n,i,ng_s,ng_h,dm,nlevs
    integer                  :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: hp(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    call build(bpt, "get_rho_Hext")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s = nghost(s(1))
    ng_h = nghost(rho_Hext(1))

    do n=1,nlevs

       do i = 1, nfabs(s(n))
          sp => dataptr(s(n) , i)
          hp => dataptr(rho_Hext(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (1)
             call bl_error("1-d heating not implemented")
          case (2)
             call get_rho_Hext_2d(tempbar_init(n,:), &
                                  hp(:,:,1,1),ng_h,sp(:,:,1,:),ng_s, &
                                  lo,hi,dx(n,:))
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
  
  
  subroutine get_rho_Hext_2d(tempbar_init,rho_Hext,ng_h,s,ng_s,lo,hi,dx)

    use bl_constants_module
    use variables
    use network, only: network_species_index
    use probin_module, only: drive_initial_convection

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer :: i, j
    integer, save :: ilh1, ilc12, iln14, ilo16, ihe4

    real(kind=dp_t) :: rho, temp, T6, T613, X_CNO, X_1, g14, eps_CNO

    logical, save :: firstCall = .true.

    if (firstCall) then
       ilh1 =  network_species_index("hydrogen-1")
       ihe4 =  network_species_index("helium-4")
       ilc12 = network_species_index("carbon-12")
       iln14 = network_species_index("nitrogen-14")
       ilo16 = network_species_index("oxygen-16")

       firstCall = .false.
    endif
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rho = s(i,j,rho_comp)

          if (drive_initial_convection) then
             temp = tempbar_init(j)
          else
             temp = s(i,j,temp_comp)
          endif

          T6 = temp / 1.0e6_dp_t
          T613 = T6**THIRD

          ! total CNO abundance
          X_CNO = (s(i,j,spec_comp-1+ilc12) + &
               s(i,j,spec_comp-1+iln14) + &
               s(i,j,spec_comp-1+ilo16)) / rho

          ! H abundance
          X_1 = s(i,j,spec_comp-1+ilh1) / rho

          ! CNO heating from Kippenhahn & Weigert, Eq. 18.65
          g14 = 1.0_dp_t + 2.7d-3*T613 - 7.78d-3*T613**2 - 1.49d-4*T6
!          eps_CNO = 8.67e27_dp_t * g14 * X_CNO * X_1 * rho * exp(-152.28_dp_t/T613) / T613**2
          eps_CNO = HCNO_energy_generation(X_CNO)

          eps_CNO = eps_CNO + triple_alpha_energy_generation(s(i,j,spec_comp+ihe4-1)/rho,rho,temp)

          rho_Hext(i,j) = rho * eps_CNO
       enddo
    enddo
    
  end subroutine get_rho_Hext_2d


  ! the temperature-insensitive Hot CNO cycle energy generation
  function HCNO_energy_generation(XCNO) result(r)
    real(kind=dp_t), intent(in   ) :: XCNO
    real(kind=dp_t) :: r

    ! this factor is from Wallace & Woosley (1981) ApJS 45
    real(kind=dp_t), parameter :: factor = 5.86e15
    
    r = factor * XCNO
  end function HCNO_energy_generation

  ! the triple alpha reaction energy generation
  ! taken from Arnett's book pg 225
  ! needs more accurate screening factor; just setting it to unity now
  function triple_alpha_energy_generation(Y,rho,T) result(r)
    real(kind=dp_t), intent(in   ) :: Y, rho, T
    real(kind=dp_t) :: r

    real(kind=dp_t), parameter :: eps_0 = 3.9e11, &
                                  f3a = 1.0d0
    real(kind=dp_t) :: t8i, t8i3

    t8i = 1.e8/T
    t8i3 = t8i*t8i*t8i
      
    r = eps_0*f3a*rho*rho*Y*Y*Y*exp(-42.94d0*t8i)*t8i3

  end function triple_alpha_energy_generation
    
    
end module heating_module
