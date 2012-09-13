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
       do i=1,nfabs(sold(n))
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
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec

    integer        , intent(in   ) :: lo(:),hi(:),ng_so,ng_sn,ng_cp,ng_xi
    real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_so:,:)
    real(kind=dp_t), intent(in   ) :: snew(lo(1)-ng_sn:,:)
    real(kind=dp_t), intent(inout) ::   cp(lo(1)-ng_cp:)
    real(kind=dp_t), intent(inout) ::   xi(lo(1)-ng_xi:,:)
    
    ! local
    integer :: i,comp    
    type (eos_t) :: eos_state

    do i=lo(1)-1,hi(1)+1
          
       ! old state first
       eos_state%rho   = sold(i,rho_comp)
       eos_state%T     = sold(i,temp_comp)
       eos_state%xn(:) = sold(i,spec_comp:spec_comp+nspec-1)/eos_state%rho
       
       ! dens, temp, and xmass are inputs
       call eos(eos_input_rt, eos_state, .false.)
       
       cp(i) = eos_state%cp

       do comp=1,nspec
          xi(i,comp) = eos_state%dhdX(comp)
       enddo


       ! new state now -- average results
       eos_state%rho   = snew(i,rho_comp)
       eos_state%T     = snew(i,temp_comp)
       eos_state%xn(:) = snew(i,spec_comp:spec_comp+nspec-1)/eos_state%rho
       
       ! dens, temp, and xmass are inputs
       call eos(eos_input_rt, eos_state, .false.)
    
       ! average the current state with the old one
       cp(i) = HALF*(eos_state%cp + cp(i))

       do comp=1,nspec
          xi(i,comp) = HALF*(eos_state%dhdX(comp) + xi(i,comp))
       enddo

    enddo
    
  end subroutine make_intra_coeffs_1d
  

  subroutine make_intra_coeffs_2d(lo,hi,sold,ng_so,snew,ng_sn,cp,ng_cp,xi,ng_xi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec

    integer        , intent(in   ) :: lo(:),hi(:),ng_so,ng_sn,ng_cp,ng_xi
    real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real(kind=dp_t), intent(in   ) :: snew(lo(1)-ng_sn:,lo(2)-ng_sn:,:)
    real(kind=dp_t), intent(inout) ::   cp(lo(1)-ng_cp:,lo(2)-ng_cp:)
    real(kind=dp_t), intent(inout) ::   xi(lo(1)-ng_xi:,lo(2)-ng_xi:,:)
    
    ! local
    integer :: i,j,comp    
    type (eos_t) :: eos_state

    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
    
          ! old state first
          eos_state%rho   = sold(i,j,rho_comp)
          eos_state%T     = sold(i,j,temp_comp)
          eos_state%xn(:) = sold(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho
          
          ! dens, temp, and xmass are inputs
          call eos(eos_input_rt, eos_state, .false.)
          
          cp(i,j) = eos_state%cp
          
          do comp=1,nspec
             xi(i,j,comp) = eos_state%dhdX(comp)
          enddo


          ! new state now -- average results
          eos_state%rho   = snew(i,j,rho_comp)
          eos_state%T     = snew(i,j,temp_comp)
          eos_state%xn(:) = snew(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho
          
          ! dens, temp, and xmass are inputs
          call eos(eos_input_rt, eos_state, .false.)
          
          ! average the current state with the old one
          cp(i,j) = HALF*(eos_state%cp + cp(i,j))
          
          do comp=1,nspec
             xi(i,j,comp) = HALF*(eos_state%dhdX(comp) + xi(i,j,comp))
          enddo

       enddo
    enddo
    
  end subroutine make_intra_coeffs_2d
  
  
  subroutine make_intra_coeffs_3d(lo,hi,sold,ng_so,snew,ng_sn,cp,ng_cp,xi,ng_xi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    
    integer        , intent(in   ) :: lo(:),hi(:),ng_so,ng_sn,ng_cp,ng_xi
    real(kind=dp_t), intent(in   ) :: sold(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real(kind=dp_t), intent(in   ) :: snew(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:,:)
    real(kind=dp_t), intent(inout) ::   cp(lo(1)-ng_cp:,lo(2)-ng_cp:,lo(3)-ng_cp:)
    real(kind=dp_t), intent(inout) ::   xi(lo(1)-ng_xi:,lo(2)-ng_xi:,lo(3)-ng_xi:,:)

    ! local
    integer :: i,j,k,comp
    type (eos_t) :: eos_state
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,comp,eos_state)
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             
             ! old state first
             eos_state%rho   = sold(i,j,k,rho_comp)
             eos_state%T     = sold(i,j,k,temp_comp)
             eos_state%xn(:) = sold(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho
             
             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, eos_state, .false.)

             cp(i,j,k) = eos_state%cp
             
             do comp=1,nspec
                xi(i,j,k,comp) = eos_state%dhdX(comp)
             enddo


             ! new state now -- average results
             eos_state%rho   = snew(i,j,k,rho_comp)
             eos_state%T     = snew(i,j,k,temp_comp)
             eos_state%xn(:) = snew(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho
             
             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, eos_state, .false.)
             
             cp(i,j,k) = HALF*(eos_state%cp + cp(i,j,k))
             
             do comp=1,nspec
                xi(i,j,k,comp) = HALF*(eos_state%dhdX(comp) + xi(i,j,k,comp))
             enddo

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine make_intra_coeffs_3d

end module make_intra_coeffs_module
