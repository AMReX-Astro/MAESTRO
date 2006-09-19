module deltap_module

  use bl_types
  use bc_module
  use multifab_module
  use eos_module
  use network
  use variables
  
  implicit none

contains

  subroutine make_deltap (deltap,comp,s,p0,temp0)

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: deltap
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: p0(:),temp0(:)
    
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
       dp => dataptr(deltap, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call makedeltap_2d(dp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng, p0, temp0)
       case (3)
          call makedeltap_3d(dp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng, p0, temp0)
       end select
    end do

  end subroutine make_deltap

  subroutine makedeltap_2d (deltap,s,lo,hi,ng,p0,temp0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: deltap(lo(1):,lo(2):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::      p0(lo(2):)
    real (kind = dp_t), intent(in   ) ::   temp0(lo(2):)
    
!     Local variables
    integer :: i, j
    
    integer :: input_flag
      
    do_diag = .false.
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          den_row(1) = s(i,j,rho_comp)
          p_row(1) = p0(j)
          h_row(1) = s(i,j,rhoh_comp) / s(i,j,rho_comp)
          xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          ! (rho,h) --> T,p, etc
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

          deltap(i,j) = (p_row(1)-p0(j))/ p0(j)

       enddo
    enddo

  end subroutine makedeltap_2d

  subroutine makedeltap_3d (deltap,s,lo,hi,ng,p0,temp0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: deltap(lo(1):,lo(2):,lo(3):)  
    real (kind = dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::      p0(lo(3):)
    real (kind = dp_t), intent(in   ) ::   temp0(lo(3):)

!     Local variables
    integer :: i, j, k

    integer :: input_flag
    
    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_row(1) = s(i,j,k,rho_comp)
             p_row(1) = p0(k)
             h_row(1) = s(i,j,k,rhoh_comp) / s(i,j,k,rho_comp)
             xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

             ! (rho,h) --> T,p, etc
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

             deltap(i,j,k) = (p_row(1)-p0(k))/ p0(k)

          enddo
       enddo
    enddo
    
  end subroutine makedeltap_3d

end module deltap_module
