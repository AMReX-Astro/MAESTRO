! create the coefficients needed for I_T (intra for the SDC
! predict_T_* enthalpy updates)


module make_intra_coeffs_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bl_constants_module

  implicit none

  private

  public :: make_intra_coeffs

contains 

  subroutine make_intra_coeffs(s,cp,xi)

    ! note: we explicitly fill the ghostcells by looping over them directly
    ! in the _1d, _2d, and _3d routines below.

    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: cp(:)
    type(multifab) , intent(inout) :: xi(:)

    ! local
    integer :: n,i,dm,nlevs
    integer :: ng_s,ng_cp,ng_xi
    integer :: lo(get_dim(s(1))),hi(get_dim(s(1)))

    real(kind=dp_t), pointer    :: sp(:,:,:,:)
    real(kind=dp_t), pointer    :: cpp(:,:,:,:),xip(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_intra_coeffs")

    dm = get_dim(s(1))
    nlevs = size(s)

    ng_s = nghost(s(1))
    ng_cp = nghost(cp(1))
    ng_xi = nghost(xi(1))

    do n=1,nlevs
       do i=1,nboxes(s(n))
          if (multifab_remote(s(n),i)) cycle
          sp       => dataptr(s(n),i)
          cpp  => dataptr(cp(n),i)
          xip  => dataptr(xi(n),i)

          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))

          select case (dm)
          case (1)
             call make_thermal_coeffs_1d(lo,hi,sp(:,1,1,:),ng_s,cpp(:,1,1,1),ng_cp, &
                                         xip(:,1,1,:),ng_xi)
          case (2)
             call make_thermal_coeffs_2d(lo,hi,sp(:,:,1,:),ng_s,cpp(:,:,1,1),ng_cp, &
                                         xip(:,:,1,:),ng_xi)
          case (3)
             call make_thermal_coeffs_3d(lo,hi,sp(:,:,:,:),ng_s,cpp(:,:,:,1),ng_cp, &
                                         xip(:,:,:,1),ng_xi)
          end select
       end do
    enddo

    call destroy(bpt)

  end subroutine make_intra_coeffs


  subroutine make_intra_coeffs_1d(lo,hi,s,ng_s,cp,ng_cp,xi,ng_xi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_cp,ng_xi
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::  cp(lo(1)-ng_cp:)
    real(kind=dp_t), intent(inout) ::  xi(lo(1)-ng_xi:,:)
    
    ! local
    integer :: i,comp    
    
    do i=lo(1)-1,hi(1)+1
          
       den_eos(1) = s(i,rho_comp)
       temp_eos(1) = s(i,temp_comp)
       xn_eos(1,:) = s(i,spec_comp:spec_comp+nspec-1)/den_eos(1)
       
       ! dens, temp, and xmass are inputs
       call eos(eos_input_rt, den_eos, temp_eos, &
                npts, &
                xn_eos, &
                p_eos, h_eos, e_eos, & 
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                .false.)
       
       cp(i) = cp_eos(1)

       do comp=1,nspec
          xi(i,comp) = dhdX_eos(1,comp)
       enddo

    enddo
    
  end subroutine make_intra_coeffs_1d
  

  subroutine make_intra_coeffs_2d(lo,hi,s,ng_s,cp,ng_cp,xi,ng_xi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_cp,ng_xi
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::  cp(lo(1)-ng_cp:,lo(2)-ng_cp:)
    real(kind=dp_t), intent(inout) ::  xi(lo(1)-ng_xi:,lo(2)-ng_xi:,:)
    
    ! local
    integer :: i,j,comp    
    
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          
          den_eos(1) = s(i,j,rho_comp)
          temp_eos(1) = s(i,j,temp_comp)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
          
          ! dens, temp, and xmass are inputs
          call eos(eos_input_rt, den_eos, temp_eos, &
                   npts, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, & 
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false.)
          
          cp(i,j) = cp_eos(1)
          
          do comp=1,nspec
             xi(i,j,comp) = dhdX_eos(1,comp)
          enddo

       enddo
    enddo
    
  end subroutine make_intra_coeffs_2d
  
  
  subroutine make_intra_coeffs_3d(lo,hi,s,ng_s,cp,ng_cp,xi,ng_xi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec
    
    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_cp,ng_xi
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::  cp(lo(1)-ng_cp:,lo(2)-ng_cp:,lo(3)-ng_cp:)
    real(kind=dp_t), intent(inout) ::  xi(lo(1)-ng_xi:,lo(2)-ng_xi:,lo(3)-ng_xi:,:)

    ! local
    integer :: i,j,k,comp
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,comp)
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             
             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = s(i,j,k,temp_comp)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
             
             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, & 
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false.)
             
             cp(i,j,k) = cp_eos(1)
             
             do comp=1,nspec
                xi(i,j,k,comp) = dhdX_eos(1,comp)
             enddo

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine make_intra_coeffs_3d

end module make_intra_coeffs_module
