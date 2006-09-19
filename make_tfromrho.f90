module tfromrho_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_module
  use eos_module
  use network
  use variables

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_tfromrho (t,comp,s,t0,p0,time,dx)

    integer        , intent(inout) :: comp
    type(multifab) , intent(inout) :: t
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: t0(:),p0(:)
    real(kind=dp_t), intent(in   ) :: time,dx(:)

    real(kind=dp_t), pointer:: sp(:,:,:,:),tp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i

    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       tp => dataptr(t, i)
       sp => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call maketfromrho_2d(tp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng, t0, p0, time, dx)
       case (3)
          call maketfromrho_3d(tp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng, t0, p0, time, dx)
       end select
    end do

  end subroutine make_tfromrho

  subroutine maketfromrho_2d (t,s,lo,hi,ng,t0,p0,time,dx)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out) :: t(lo(1):,lo(2):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) :: t0(lo(2):)
    real (kind=dp_t), intent(in   ) :: p0(lo(2):)
    real (kind=dp_t), intent(in   ) :: time,dx(:)

    !     Local variables
    integer :: i, j

    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          den_row(1) = s(i,j,rho_comp)
          temp_row(1) = t0(j)
          p_row(1) = p0(j)
          xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          ! (rho,P) --> T,h
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
          
          t(i,j) = log(temp_row(1))/log(10.)
       enddo
    enddo

  end subroutine maketfromrho_2d

  subroutine maketfromrho_3d (t,s,lo,hi,ng,t0,p0,time,dx)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out) :: t(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) :: t0(lo(3):)
    real (kind=dp_t), intent(in   ) :: p0(lo(3):)
    real (kind=dp_t), intent(in   ) :: time,dx(:)

    !     Local variables
    integer :: i, j, k

    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_row(1) = s(i,j,k,rho_comp)
             temp_row(1) = t0(k)
             p_row(1) = p0(k)
             xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

             ! (rho,P) --> T,h
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

             t(i,j,k) = log(temp_row(1))/log(10.)

          enddo
       enddo
    enddo

  end subroutine maketfromrho_3d

end module tfromrho_module
