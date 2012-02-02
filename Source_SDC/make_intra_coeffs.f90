! create the coefficients needed for I_T (intra for the SDC
! predict_T_* enthalpy updates).  We take an old and new state
! and return the time-centered quantities.


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

  subroutine make_intra_coeffs(sold,snew,cp,xi)

    ! note: we explicitly fill the ghostcells by looping over them directly
    ! in the _1d, _2d, and _3d routines below.

    type(multifab) , intent(in   ) :: sold(:), snew(:)
    type(multifab) , intent(inout) :: cp(:)
    type(multifab) , intent(inout) :: xi(:)

    ! local
    integer :: n,i,dm,nlevs
    integer :: ng_so,ng_sn,ng_cp,ng_xi
    integer :: lo(get_dim(sold(1))),hi(get_dim(sold(1)))

    real(kind=dp_t), pointer    :: sop(:,:,:,:),snp(:,:,:,:)
    real(kind=dp_t), pointer    :: cpp(:,:,:,:),xip(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_intra_coeffs")

    dm = get_dim(sold(1))
    nlevs = size(sold)

    ng_so = nghost(sold(1))
    ng_sn = nghost(snew(1))
    ng_cp = nghost(cp(1))
    ng_xi = nghost(xi(1))

    do n=1,nlevs
       do i=1,nboxes(sold(n))
          if (multifab_remote(sold(n),i)) cycle
          sop  => dataptr(sold(n),i)
          snp  => dataptr(snew(n),i)
          cpp  => dataptr(cp(n),i)
          xip  => dataptr(xi(n),i)

          lo = lwb(get_box(sold(n),i))
          hi = upb(get_box(sold(n),i))

          select case (dm)
          case (1)
             call make_intra_coeffs_1d(lo,hi,sop(:,1,1,:),ng_so,snp(:,1,1,:),ng_sn, &
                                       cpp(:,1,1,1),ng_cp,xip(:,1,1,:),ng_xi)
          case (2)
             call make_intra_coeffs_2d(lo,hi,sop(:,:,1,:),ng_so,snp(:,:,1,:),ng_sn, &
                                       cpp(:,:,1,1),ng_cp,xip(:,:,1,:),ng_xi)
          case (3)
             call make_intra_coeffs_3d(lo,hi,sop(:,:,:,:),ng_so,snp(:,:,:,:),ng_sn, &
                                       cpp(:,:,:,1),ng_cp,xip(:,:,:,:),ng_xi)
          end select
       end do
    enddo

    call destroy(bpt)

  end subroutine make_intra_coeffs


  subroutine make_intra_coeffs_1d(lo,hi,sold,ng_so,snew,ng_sn,cp,ng_cp,xi,ng_xi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec

    integer        , intent(in   ) :: lo(:),hi(:),ng_so,ng_sn,ng_cp,ng_xi
    real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_so:,:)
    real(kind=dp_t), intent(in   ) :: snew(lo(1)-ng_sn:,:)
    real(kind=dp_t), intent(inout) ::   cp(lo(1)-ng_cp:)
    real(kind=dp_t), intent(inout) ::   xi(lo(1)-ng_xi:,:)
    
    ! local
    integer :: i,comp    
    
    do i=lo(1)-1,hi(1)+1
          
       ! old state first
       den_eos  = sold(i,rho_comp)
       temp_eos = sold(i,temp_comp)
       xn_eos(:) = sold(i,spec_comp:spec_comp+nspec-1)/den_eos
       
       ! dens, temp, and xmass are inputs
       call eos(eos_input_rt, den_eos, temp_eos, &
                xn_eos, &
                p_eos, h_eos, e_eos, & 
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                .false.)
       
       cp(i) = cp_eos

       do comp=1,nspec
          xi(i,comp) = dhdX_eos(comp)
       enddo


       ! new state now -- average results
       den_eos  = snew(i,rho_comp)
       temp_eos = snew(i,temp_comp)
       xn_eos(:) = snew(i,spec_comp:spec_comp+nspec-1)/den_eos
       
       ! dens, temp, and xmass are inputs
       call eos(eos_input_rt, den_eos, temp_eos, &
                xn_eos, &
                p_eos, h_eos, e_eos, & 
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                .false.)
    
       ! average the current state with the old one
       cp(i) = HALF*(cp_eos + cp(i))

       do comp=1,nspec
          xi(i,comp) = HALF*(dhdX_eos(comp) + xi(i,comp))
       enddo

    enddo
    
  end subroutine make_intra_coeffs_1d
  

  subroutine make_intra_coeffs_2d(lo,hi,sold,ng_so,snew,ng_sn,cp,ng_cp,xi,ng_xi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec

    integer        , intent(in   ) :: lo(:),hi(:),ng_so,ng_sn,ng_cp,ng_xi
    real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real(kind=dp_t), intent(in   ) :: snew(lo(1)-ng_sn:,lo(2)-ng_sn:,:)
    real(kind=dp_t), intent(inout) ::   cp(lo(1)-ng_cp:,lo(2)-ng_cp:)
    real(kind=dp_t), intent(inout) ::   xi(lo(1)-ng_xi:,lo(2)-ng_xi:,:)
    
    ! local
    integer :: i,j,comp    
    
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
    
          ! old state first
          den_eos  = sold(i,j,rho_comp)
          temp_eos = sold(i,j,temp_comp)
          xn_eos(:) = sold(i,j,spec_comp:spec_comp+nspec-1)/den_eos
          
          ! dens, temp, and xmass are inputs
          call eos(eos_input_rt, den_eos, temp_eos, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, & 
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false.)
          
          cp(i,j) = cp_eos
          
          do comp=1,nspec
             xi(i,j,comp) = dhdX_eos(comp)
          enddo


          ! new state now -- average results
          den_eos  = snew(i,j,rho_comp)
          temp_eos = snew(i,j,temp_comp)
          xn_eos(:) = snew(i,j,spec_comp:spec_comp+nspec-1)/den_eos
          
          ! dens, temp, and xmass are inputs
          call eos(eos_input_rt, den_eos, temp_eos, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, & 
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false.)
          
          ! average the current state with the old one
          cp(i,j) = HALF*(cp_eos + cp(i,j))
          
          do comp=1,nspec
             xi(i,j,comp) = HALF*(dhdX_eos(comp) + xi(i,j,comp))
          enddo

       enddo
    enddo
    
  end subroutine make_intra_coeffs_2d
  
  
  subroutine make_intra_coeffs_3d(lo,hi,sold,ng_so,snew,ng_sn,cp,ng_cp,xi,ng_xi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec
    
    integer        , intent(in   ) :: lo(:),hi(:),ng_so,ng_sn,ng_cp,ng_xi
    real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real(kind=dp_t), intent(in   ) :: snew(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:,:)
    real(kind=dp_t), intent(inout) ::   cp(lo(1)-ng_cp:,lo(2)-ng_cp:,lo(3)-ng_cp:)
    real(kind=dp_t), intent(inout) ::   xi(lo(1)-ng_xi:,lo(2)-ng_xi:,lo(3)-ng_xi:,:)

    ! local
    integer :: i,j,k,comp
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,comp)
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             
             ! old state first
             den_eos  = sold(i,j,k,rho_comp)
             temp_eos = sold(i,j,k,temp_comp)
             xn_eos(:) = sold(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos
             
             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, den_eos, temp_eos, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, & 
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false.)
             
             cp(i,j,k) = cp_eos
             
             do comp=1,nspec
                xi(i,j,k,comp) = dhdX_eos(comp)
             enddo


             ! new state now -- average results
             den_eos  = snew(i,j,k,rho_comp)
             temp_eos = snew(i,j,k,temp_comp)
             xn_eos(:) = snew(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos
             
             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, den_eos, temp_eos, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, & 
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false.)
             
             cp(i,j,k) = HALF*(cp_eos + cp(i,j,k))
             
             do comp=1,nspec
                xi(i,j,k,comp) = HALF*(dhdX_eos(comp) + xi(i,j,k,comp))
             enddo

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine make_intra_coeffs_3d

end module make_intra_coeffs_module
