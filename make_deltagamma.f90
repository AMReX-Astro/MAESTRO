module deltagamma_module

  use bl_types
  use bc_module
  use multifab_module
  use eos_module
  use network
  use variables

  implicit none

contains

  subroutine make_deltagamma (deltagamma,comp,s,rho0,p0,t0,rhoX0)

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: deltagamma
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: rho0(:), p0(:), t0(:), rhoX0(:,:)
    
    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: mp(:,:,:,:)
    real(kind=dp_t), pointer:: dp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i

    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sp => dataptr(s, i)
       dp => dataptr(deltagamma, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call makedeltagamma_2d(dp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng, rho0, p0, t0, rhoX0)
       case (3)
          call makedeltagamma_3d(dp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng, rho0, p0, t0, rhoX0)
       end select
    end do

  end subroutine make_deltagamma

  subroutine makedeltagamma_2d (deltagamma,s,lo,hi,ng,rho0,p0,t0,rhoX0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    rho0(lo(2):)
    real (kind = dp_t), intent(in   ) ::      p0(lo(2):)
    real (kind = dp_t), intent(in   ) ::      t0(lo(2):)
    real (kind = dp_t), intent(in   ) ::   rhoX0(lo(2):,:)

!     Local variables
    integer :: i, j

    integer :: input_flag
    real (kind = dp_t) :: gam10
      
    do_diag = .false.

    do j = lo(2), hi(2)

       ! first compute the base state gamma1
       den_row(1) = rho0(j)
       p_row(1) = p0(j)
       temp_row(1) = t0(j)   ! used as an initial guess

       xn_zone(:) = rhoX0(j,:)/rho0(j)

       input_flag = 4

       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, &
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                do_diag)

       gam10 = gam1_row(1)

       do i = lo(1), hi(1)

          ! compute the thermodynamics of the full state
          den_row(1) = s(i,j,rho_comp)
          h_row(1) = s(i,j,rhoh_comp) / s(i,j,rho_comp)
          temp_row(1) = t0(j)
            
          xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          input_flag = 2

          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)

          deltagamma(i,j) = gam1_row(1) - gam10

       enddo
    enddo

  end subroutine makedeltagamma_2d
    
  subroutine makedeltagamma_3d (deltagamma,s,lo,hi,ng,rho0,p0,t0,rhoX0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):,lo(3):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    rho0(lo(3):)
    real (kind = dp_t), intent(in   ) ::      p0(lo(3):)
    real (kind = dp_t), intent(in   ) ::      t0(lo(3):)
    real (kind = dp_t), intent(in   ) ::   rhoX0(lo(3):,:)


!     Local variables
    integer :: i, j, k
    integer :: input_flag
    real (kind = dp_t) :: gam10
      
    do_diag = .false.

    do k = lo(3), hi(3)

       ! first compute the base state gamma1
       den_row(1) = rho0(k)
       p_row(1) = p0(k)
       temp_row(1) = t0(k)   ! used as an initial guess

       xn_zone(:) = rhoX0(k,:)/rho0(k)

       input_flag = 4

       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, &
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                do_diag)

       gam10 = gam1_row(1)

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! compute the thermodynamics of the full state
             den_row(1) = s(i,j,k,rho_comp)
             h_row(1) = s(i,j,k,rhoh_comp) / s(i,j,k,rho_comp)
             temp_row(1) = t0(k)
            
             xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)
             
             input_flag = 2

             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, &
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             deltagamma(i,j,k) = gam1_row(1) - gam10

          enddo
       enddo
    enddo

  end subroutine makedeltagamma_3d

end module deltagamma_module
