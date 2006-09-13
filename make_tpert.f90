module tpert_module

  use bl_types
  use bc_module
  use multifab_module
  use eos_module

  implicit none

contains

  subroutine make_tpert (T,comp,state,p0,temp0)

    integer        , intent(inout) :: comp
    type(multifab) , intent(inout) :: T
    type(multifab) , intent(in   ) :: state
    real(kind=dp_t), intent(in   ) ::    p0(:)
    real(kind=dp_t), intent(in   ) :: temp0(:)

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(state%dim),hi(state%dim),ng,dm
    integer :: i

    ng = state%ng
    dm = state%dim

    do i = 1, state%nboxes
       if ( multifab_remote(state, i) ) cycle
       sp => dataptr(state, i)
       tp => dataptr(T    , i)
       lo =  lwb(get_box(state, i))
       hi =  upb(get_box(state, i))
       select case (dm)
       case (2)
          call maketpert_2d(tp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng, p0, temp0)
       case (3)
          call maketpert_3d(tp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng, p0, temp0)
       end select
    end do

  end subroutine make_tpert

  subroutine maketpert_2d (T,state,lo,hi,ng,p0,temp0)

     implicit none
     integer, intent(in) :: lo(:), hi(:), ng
     real (kind = dp_t), intent(  out) ::     T(lo(1)   :,lo(2):)  
     real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,:)
     real (kind = dp_t), intent(in   ) ::    p0(lo(2):)
     real (kind = dp_t), intent(in   ) :: temp0(lo(2):)

     integer :: i, j

     do_diag = .false.

     do j = lo(2), hi(2)
       do i = lo(1), hi(1)
            
         den_row(1) = state(i,j,1)
        temp_row(1) = temp0(j)
           p_row(1) = p0(j)

         ! (rho, P) --> T
         input_flag = 4

         call eos(input_flag, den_row, temp_row, npts, nspec, &
              xmass, aion, zion, &
              p_row, h_row, e_row, &
              cv_row, cp_row, xne_row, eta_row, &
              pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
              s_row, do_diag)

         T(i,j) = temp_row(1) - temp0(j)
       enddo
     enddo

  end subroutine maketpert_2d

  subroutine maketpert_3d (T,state,lo,hi,ng,p0,temp0)

     implicit none
     integer, intent(in) :: lo(:), hi(:), ng
     real (kind = dp_t), intent(  out) ::     T(lo(1)   :,lo(2):   ,lo(3)   :  )  
     real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
     real (kind = dp_t), intent(in   ) ::    p0(lo(3):)
     real (kind = dp_t), intent(in   ) :: temp0(lo(3):)

     integer :: i, j, k

     do_diag = .false.

     do k = lo(3), hi(3)
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)
            
         den_row(1) = state(i,j,k,1)
        temp_row(1) = temp0(k)
           p_row(1) = p0(k)

         ! (rho, P) --> T
         input_flag = 4

         call eos(input_flag, den_row, temp_row, npts, nspec, &
              xmass, aion, zion, &
              p_row, h_row, e_row, &
              cv_row, cp_row, xne_row, eta_row, &
              pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
              s_row, do_diag)

         T(i,j,k) = temp_row(1) - temp0(k)
       enddo
       enddo
     enddo

  end subroutine maketpert_3d

end module tpert_module
