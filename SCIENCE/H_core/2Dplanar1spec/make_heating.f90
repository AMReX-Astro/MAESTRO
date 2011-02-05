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
    use ml_restriction_module
    use variables, only: foextrap_comp, rho_comp, temp_comp, rhoh_comp
    use geometry, only: nlevs_radial, nr_fine, nr, dr
    use average_module
    use fill_3d_module

    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: tempbar_init(:,0:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt

    ! local
    integer                  :: n,i,ng_s,ng_h, nlevs
    integer                  :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: hp(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    logical, parameter :: symmetric_heat = .true.
    real(kind=dp_t), allocatable :: rho0(:,:)
    real(kind=dp_t), allocatable :: T0(:,:)
    real(kind=dp_t), allocatable :: rhoh0(:,:)

    call build(bpt, "get_rho_Hext")

    nlevs = mla%nlevel

    ng_s = s(1)%ng
    ng_h = rho_Hext(1)%ng

    if (symmetric_heat) then
       
       allocate( rho0(nlevs_radial,0:nr_fine-1))
       allocate(rhoh0(nlevs_radial,0:nr_fine-1))
       allocate(   T0(nlevs_radial,0:nr_fine-1))

       ! average rho and T to get rho0 and T0
       ! rho0 is a 2D array with the first index being the level
       call average(mla,s,rho0,dx,rho_comp)
       call average(mla,s,T0,dx,temp_comp)
          
       ! compute heating using base state
       do n=1,nlevs_radial

          call get_symm_rho_Hext(rhoh0(n,:),rho0(n,:),T0(n,:),nr(n),dr(n))
       
       enddo

       ! average back onto cartesian grid
       call put_1d_array_on_cart(rhoh0,rho_Hext,(mla%dim)+rhoh_comp,.false., &
                                  .false.,dx,the_bc_level,mla)

       deallocate( rho0)
       deallocate(rhoh0)
       deallocate(   T0)


    else

       do n=1,nlevs

          do i = 1, s(n)%nboxes
             if ( multifab_remote(s(n), i) ) cycle
             sp => dataptr(s(n) , i)
             hp => dataptr(rho_Hext(n) , i)
             lo =  lwb(get_box(s(n), i))
             hi =  upb(get_box(s(n), i))
             select case (mla%dim)
             case (1)
                call get_rho_Hext_1d(hp(:,1,1,1),ng_h,sp(:,1,1,:),ng_s,lo,hi,dx(n,:))
             case (2)
                call get_rho_Hext_2d(hp(:,:,1,1),ng_h,sp(:,:,1,:),ng_s,lo,hi,dx(n,:))
             case (3)
                call get_rho_Hext_3d(hp(:,:,:,1),ng_h,sp(:,:,:,:),ng_s,lo,hi,dx(n,:))
             end select
          end do

       end do
 
    end if

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       call ml_cc_restriction(rho_Hext(n-1), rho_Hext(n), mla%mba%rr(n-1,:))
    end do

    call destroy(bpt)

  end subroutine get_rho_Hext
  
  subroutine get_symm_rho_Hext(rho_Hext,rho0,T0,nr,dr)
    use bl_constants_module
    use variables
    use network
    
    integer        , intent(in   ) :: nr
    real(kind=dp_t), intent(inout) :: rho_Hext(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: T0(0:)
    real(kind=dp_t), intent(in   ) :: dr

    integer                        :: i, ibegin
    integer                        :: h1_comp
    integer                        :: c12_comp
    integer                        :: n14_comp
    integer                        :: o16_comp
    real(kind=dp_t)                :: rho, T_6_third, X_CNO, X_1, g14
    real(kind=dp_t)                :: tmp1, tmp2, tmp3
    real(kind=dp_t)                :: intgrl


    rho_Hext = 0.0_dp_t

!..............................................................................
!     h1_comp = spec_comp - 1 + network_species_index("hydrogen-1")
!     c12_comp = spec_comp - 1 + network_species_index("carbon-12")
!     n14_comp = spec_comp - 1 + network_species_index("nitrogen-14")
!     o16_comp = spec_comp - 1 + network_species_index("oxygen-16")

! FIXME! HACK for anelastic comparison
    X_CNO = 0.00952433417 !0.9524325640064698d-2
    X_1   = 0.40137771765 !0.4013777831131446d0

    ibegin = -1
    intgrl = 0.d0

    do i = 0, nr-1
       rho = rho0(i)

       T_6_third = (T0(i) / 1.0d6) ** THIRD
! FIXME! HACK
!          tmp1 = s(i,j,c12_comp)
!          tmp2 = s(i,j,n14_comp)
!          tmp3 = s(i,j,o16_comp)
!          X_CNO = (tmp1 + tmp2 + tmp3) / rho
!          X_1 = s(i,j,h1_comp) / rho
       tmp1 =   2.7d-3 * T_6_third
       tmp2 = -7.78d-3 * T_6_third**2
       tmp3 = -1.49d-4 * T_6_third**3
       g14 = 1.0_dp_t + tmp1 + tmp2 + tmp3
       tmp1 = 8.67d27 * g14 * X_CNO * X_1 * rho / T_6_third**2
       tmp2 = dexp(-1.5228d2 / T_6_third)
       rho_Hext(i) = rho * tmp1 * tmp2
          
       if (rho .gt. 5.944) then

          intgrl = intgrl + rho_Hext(i)*(dble(i)+HALF)*dr*(dble(i)+HALF)*dr*dr

          rho_Hext(i) = 0.d0

       else

          if (ibegin .lt. 0) then
             ! add in a "surface flux"
             rho_Hext(i) = rho_Hext(i) + 3.d0*intgrl / &
                           (dr*dr*dr*(i*i*i + 3.d0*i*i + 3.d0*i +1))

             ibegin = i
          end if

       end if
    end do
!.............................................................................. 

  end subroutine get_symm_rho_Hext

  subroutine get_rho_Hext_1d(rho_Hext,ng_h,s,ng_s,lo,hi,dx)
    
    use bl_constants_module
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    call bl_error('1d external heating not written yet')
    
  end subroutine get_rho_Hext_1d
  
  subroutine get_rho_Hext_2d(rho_Hext,ng_h,s,ng_s,lo,hi,dx)

    use bl_constants_module
    use variables
    use network
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer                        :: i, j
    integer                        :: h1_comp
    integer                        :: c12_comp
    integer                        :: n14_comp
    integer                        :: o16_comp
    real(kind=dp_t)                :: rho, T_6_third, X_CNO, X_1, g14
    real(kind=dp_t)                :: tmp1, tmp2, tmp3

!    call bl_error('2d external heating not written yet')
    rho_Hext = 0.0_dp_t

!..............................................................................
!     h1_comp = spec_comp - 1 + network_species_index("hydrogen-1")
!     c12_comp = spec_comp - 1 + network_species_index("carbon-12")
!     n14_comp = spec_comp - 1 + network_species_index("nitrogen-14")
!     o16_comp = spec_comp - 1 + network_species_index("oxygen-16")

! FIXME! HACK for anelastic comparison
    X_CNO = 0.00952433417 !0.9524325640064698d-2
    X_1   = 0.40137771765 !0.4013777831131446d0


    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rho = s(i,j,rho_comp)
          T_6_third = (s(i,j,temp_comp) / 1.0d6) ** THIRD
! FIXME! HACK
!          tmp1 = s(i,j,c12_comp)
!          tmp2 = s(i,j,n14_comp)
!          tmp3 = s(i,j,o16_comp)
!          X_CNO = (tmp1 + tmp2 + tmp3) / rho
!          X_1 = s(i,j,h1_comp) / rho
          tmp1 =   2.7d-3 * T_6_third
          tmp2 = -7.78d-3 * T_6_third**2
          tmp3 = -1.49d-4 * T_6_third**3
          g14 = 1.0_dp_t + tmp1 + tmp2 + tmp3
          tmp1 = 8.67d27 * g14 * X_CNO * X_1 * rho / T_6_third**2
          tmp2 = dexp(-1.5228d2 / T_6_third)
          rho_Hext(i,j) = rho * tmp1 * tmp2
       enddo
    enddo
!.............................................................................. 

  end subroutine get_rho_Hext_2d
  
  subroutine get_rho_Hext_3d(rho_Hext,ng_h,s,ng_s,lo,hi,dx)
    
    use bl_constants_module
    use variables
    use network

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

!..............................................................................
    integer                        :: i, j, k
    integer                        :: h1_comp
    integer                        :: c12_comp
    integer                        :: n14_comp
    integer                        :: o16_comp
    real(kind=dp_t)                :: rho, T_6_third, X_CNO, X_1, g14
    real(kind=dp_t)                :: tmp1, tmp2, tmp3
!..............................................................................
    
    rho_Hext = 0.0_dp_t

!..............................................................................
!     h1_comp = spec_comp - 1 + network_species_index("hydrogen-1")
!     c12_comp = spec_comp - 1 + network_species_index("carbon-12")
!     n14_comp = spec_comp - 1 + network_species_index("nitrogen-14")
!     o16_comp = spec_comp - 1 + network_species_index("oxygen-16")
!     do k = lo(3), hi(3)
!        do j = lo(2), hi(2)
!           do i = lo(1), hi(1)
!              rho = s(i,j,k,rho_comp)
!              T_6_third = (s(i,j,k,temp_comp) / 1.0d6) ** THIRD
!              tmp1 = s(i,j,k,c12_comp)
!              tmp2 = s(i,j,k,n14_comp)
!              tmp3 = s(i,j,k,o16_comp)
!              X_CNO = (tmp1 + tmp2 + tmp3) / rho
!              X_1 = s(i,j,k,h1_comp) / rho
!              tmp1 =   2.7d-3 * T_6_third
!              tmp2 = -7.78d-3 * T_6_third**2
!              tmp3 = -1.49d-4 * T_6_third**3
!              g14 = 1.0_dp_t + tmp1 + tmp2 + tmp3
!              tmp1 = 8.67d27 * g14 * X_CNO * X_1 * rho / T_6_third**2
!              tmp2 = dexp(-1.5228d2 / T_6_third)
!              rho_Hext(i,j,k) = rho * tmp1 * tmp2
!           enddo
!       enddo
!     enddo
!..............................................................................
    
  end subroutine get_rho_Hext_3d
  
end module heating_module
