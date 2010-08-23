! This routine returns the externally imposed (i.e. not reactions)
! heating source term to the enthalpy equation (actually rho * H
! is returned, where H has units of erg/g/s).

module heating_module

  use bl_types
  use define_bc_module

  implicit none
  
  private
  public :: get_rho_Hext
  
contains

  subroutine get_rho_Hext(mla,s,rho_Hext,dx,time,dt)

    use geometry, only: nlevs
    use multifab_module
    use ml_layout_module
    use ml_restriction_module
    use geometry, only: dm
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
          case (2)
             call get_rho_Hext_2d(hp(:,:,1,1),ng_h,sp(:,:,1,:),ng_s,lo,hi,dx(n,:),time)
          case (3)
             call get_rho_Hext_3d(hp(:,:,:,1),ng_h,sp(:,:,:,:),ng_s,lo,hi,dx(n,:),time,dt)
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
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    rho_Hext = 0.0_dp_t
    
  end subroutine get_rho_Hext_2d
  
  subroutine get_rho_Hext_3d(rho_Hext,ng_h,s,ng_s,lo,hi,dx,time,dt)
    
    use bl_constants_module
    use geometry, only: center
    use variables, only: rho_comp, temp_comp, spec_comp
    use probin_module, only: heating_time, heating_rad, heating_peak, &
         heating_sigma, prob_type
    use network
    use bl_error_module

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time,dt

    integer :: i, j, k
    real(kind=dp_t) :: x, y, z, r, fac

    integer                        :: h1_comp
    integer                        :: c12_comp
    integer                        :: n14_comp
    integer                        :: o16_comp
    real(kind=dp_t)                :: rho, T_6_third, X_CNO, X_1, g14
    real(kind=dp_t)                :: tmp1, tmp2, tmp3


    rho_Hext = 0.d0

    if (prob_type .eq. 1) then
       
       if (time .lt. heating_time) then

          if ( (time+dt) .gt. heating_time ) then
             fac = (heating_time - time) / dt
          else
             fac = 1.d0
          end if

          do k = lo(3), hi(3)
             z = (dble(k) + HALF)*dx(3) - center(3)

             do j = lo(2), hi(2)
                y = (dble(j) + HALF)*dx(2) - center(2)

                do i = lo(1), hi(1)
                   x = (dble(i) + HALF)*dx(1) - center(1)

                   r = sqrt(x**2 + y**2 + z**2)

                   rho_Hext(i,j,k) = fac * s(i,j,k,rho_comp) * heating_peak &
                        * exp(-(r-heating_rad)**2/heating_sigma)

                enddo
             enddo
          enddo

       end if

    elseif (prob_type .eq. 2) then

       h1_comp = spec_comp - 1 + network_species_index("hydrogen-1")
       c12_comp = spec_comp - 1 + network_species_index("carbon-12")
       n14_comp = spec_comp - 1 + network_species_index("nitrogen-14")
       o16_comp = spec_comp - 1 + network_species_index("oxygen-16")

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rho = s(i,j,k,rho_comp)
                T_6_third = (s(i,j,k,temp_comp) / 1.0d6) ** THIRD
                tmp1 = s(i,j,k,c12_comp)
                tmp2 = s(i,j,k,n14_comp)
                tmp3 = s(i,j,k,o16_comp)
                X_CNO = (tmp1 + tmp2 + tmp3) / rho
                X_1 = s(i,j,k,h1_comp) / rho
                tmp1 =   2.7d-3 * T_6_third
                tmp2 = -7.78d-3 * T_6_third**2
                tmp3 = -1.49d-4 * T_6_third**3
                g14 = 1.0_dp_t + tmp1 + tmp2 + tmp3
                tmp1 = 8.67d27 * g14 * X_CNO * X_1 * rho / T_6_third**2
                tmp2 = dexp(-1.5228d2 / T_6_third)
                rho_Hext(i,j,k) = rho * tmp1 * tmp2
             enddo
          enddo
       enddo

    else
       call bl_error("prob_type not yet supported.")       
    endif

  end subroutine get_rho_Hext_3d
  
end module heating_module
