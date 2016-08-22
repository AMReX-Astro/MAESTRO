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

  subroutine get_rho_Hext(mla,tempbar_init,s,rho_Hext,the_bc_level,dx,dt)

    use multifab_module
    use ml_layout_module
    use ml_cc_restriction_module
    use variables, only: foextrap_comp

    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: tempbar_init(:,0:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt

    ! local
    integer                  :: n,i,ng_s,ng_h
    integer                  :: lo(mla%dim),hi(mla%dim),dm,nlevs
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: hp(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    call build(bpt, "get_rho_Hext")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s = s(1)%ng
    ng_h = rho_Hext(1)%ng

    do n=1,nlevs

       do i = 1, nfabs(s(n))
          sp => dataptr(s(n) , i)
          hp => dataptr(rho_Hext(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          call get_rho_Hext_3d(hp(:,:,:,1), ng_h, sp(:,:,:,:), ng_s, lo, hi, dx(n,:), &
		tempbar_init(n,:))
       end do

    end do

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       call ml_cc_restriction(rho_Hext(n-1), rho_Hext(n), mla%mba%rr(n-1,:))
    end do

    call destroy(bpt)

  end subroutine get_rho_Hext

  subroutine get_rho_Hext_3d(rho_Hext,ng_h,s,ng_s,lo,hi,dx,tempbar_init)

    use bl_constants_module
    use variables, only: rho_comp, spec_comp, temp_comp, press_comp
    use probin_module, only: prob_lo, drive_initial_convection
    use time_module, only: time
    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(0:)


    integer :: i,j,k
    real(kind=dp_t) :: T_in, T9, rho 
    real(kind=dp_t) :: z,z_layer
    real(kind=dp_t) :: ez,Hmax
    real(kind=dp_t) :: pi,L_x,L_y
    real(kind=dp_t) :: xspec(nspec), xcno
    real(kind=dp_t) :: epp, ecno 
    real(kind=dp_t) :: gpp, gcno
    real(kind=dp_t) :: euler
    
    euler = 2.71828182845904523536E0

    rho_Hext = 0.0_dp_t
    Hmax = 0.0_dp_t
    
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(i,j,k,rho,xspec,T_in,T9,gpp,epp,xcno,gcno,ecno) &
    !$OMP SCHEDULE(DYNAMIC,1) 
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            rho = s(i,j,k,rho_comp)
            xspec = s(i,j,k,spec_comp:spec_comp+nspec-1) / rho
            if (drive_initial_convection) then
		T_in = tempbar_init(k)
	    else
	     	T_in = s(i,j,k,temp_comp)
	    endif
	    T9 = T_in * 10.0 ** (-9.0)
	   
	   ! approximate pp energy generation without screening according to Kippenhahn, Weigert, Weiss textbook
	    gpp = 1.0 + 3.82 * T9 + 1.51 * T9 ** TWO + 0.144 * T9 ** THREE - 0.0114 * T9 ** FOUR
	    epp = 2.57E4 * gpp * rho * xspec(1) * xspec(1) * T9 **( -TWO3RD) * euler ** (-3.381 * T9 ** (-THIRD))
	   ! approximate pp energy generation without screening according to maths.qmul.ac.uk
	   ! epp = 2.6 * xspec(1) * xspec(1) * rho * (T9 ** (4.5)) * (10.0 ** 3.5) 
	   
	   
	   ! approximate cno energy generation without screening according to Kippenhahn, Weigert, Weiss textbook
	    xcno = xspec(3)
	    gcno = 1.0 - 2.00 * T9 + 3.41 * T9 ** TWO - 2.43 * T9 ** THREE
	    ecno = 8.24E25 * gcno * xcno * xspec(1) * rho * T9 ** ( -TWO3RD) * euler ** (-15.231 * T9 ** (-THIRD) - (1.25 *T9 ) ** TWO )
	   
	   ! approximate cno energy generation without screening according to maths.qmul.ac.uk
	   ! xcno = xspec(3) + xspec(4)
	   ! ecno = 7.9 * xcno * xspec(1) * rho * (T9 ** (16.0)) * (10.0 ** 26.0)
	   

	    ! rho_Hext = sum of pp and cno energy generation times rho
	    rho_Hext(i,j,k) = rho * ( epp + ecno ) 
	    

	    enddo
       enddo
    enddo
    !$OMP END PARALLEL DO	  

  end subroutine get_rho_Hext_3d

end module heating_module
