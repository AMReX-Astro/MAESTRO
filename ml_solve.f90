module ml_solve_module

   use bl_types
   use bl_constants_module
   use define_bc_module
   use multifab_module
   use boxarray_module
   use stencil_module
   use mg_module
   use list_box_module
   use ml_boxarray_module
   use itsol_module
   use sparse_solve_module
   use bl_mem_stat_module
   use bl_timer_module
   use box_util_module
   use bl_IO_module
   use fabio_module

   use ml_restriction_module
   use ml_prolongation_module
   use ml_interface_stencil_module
   use ml_util_module
   use bndry_reg_module
   use ml_cc_module
   use ml_nd_module
 
   implicit none

contains

   subroutine ml_cc_solve(mla,mgt,rh,full_soln,fine_flx,ref_ratio,do_diagnostics)

      type(ml_layout), intent(inout) :: mla
      type(mg_tower ), intent(inout) :: mgt(:)
      type(multifab ), intent(inout) :: rh(:)
      type(multifab ), intent(inout) :: full_soln(:)
      type(bndry_reg), intent(inout) :: fine_flx(2:)
      integer        , intent(in   ) :: ref_ratio(:,:)
      integer        , intent(in   ) :: do_diagnostics

      type(boxarray) :: bac

      type(lmultifab), pointer     :: fine_mask(:) => Null()
      integer                      :: i,dm,n,nlevs
      integer                      :: mglev
      real(dp_t)                   :: eps

      dm    = mla%dim
      nlevs = mla%nlevel

      eps = 1.d-12

      allocate(fine_mask(nlevs))
      do n = nlevs, 1, -1
        call lmultifab_build(fine_mask(n), mla%la(n), 1, 0)
        call setval(fine_mask(n), val = .true., all = .true.)
      end do
      do n = nlevs-1, 1, -1
        call copy(bac, get_boxarray(mla%la(n+1)))
        call boxarray_coarsen(bac, ref_ratio(n,:))
        call setval(fine_mask(n), .false., bac)
        call destroy(bac)
      end do


! ****************************************************************************

      call ml_cc(mla,mgt,rh,full_soln,fine_mask,ref_ratio,do_diagnostics,eps,.true.)

! ****************************************************************************

!   Put boundary conditions of soln in fine_flx to get correct grad(phi) at
!     crse-fine boundaries (after soln correctly interpolated in ml_cc)
    do n = 2,nlevs
       mglev = mgt(n)%nlevels
       do i = 1, dm
          call ml_fill_fine_fluxes(mgt(n)%ss(mglev), fine_flx(n)%bmf(i,0), &
                                   full_soln(n), mgt(n)%mm(mglev), -1, i)
          call ml_fill_fine_fluxes(mgt(n)%ss(mglev), fine_flx(n)%bmf(i,1), &
                                   full_soln(n), mgt(n)%mm(mglev),  1, i)
       end do
    end do

    do n = 1,nlevs
      call lmultifab_destroy(fine_mask(n))
    end do
    deallocate(fine_mask)

   end subroutine ml_cc_solve

   subroutine ml_nd_solve(mla,mgt,rh,full_soln,one_sided_ss,ref_ratio,do_diagnostics,eps_in)

       type(ml_layout), intent(inout) :: mla
       type(mg_tower) , intent(inout) :: mgt(:)
       type(multifab) , intent(inout) :: rh(:)
       type(multifab) , intent(inout) :: full_soln(:)
       type(multifab) , intent(in   ) :: one_sided_ss(2:)
       integer        , intent(in   ) :: ref_ratio(:,:)
       integer        , intent(in   ) :: do_diagnostics 
       real(dp_t)     , intent(in   ), optional :: eps_in

       type(lmultifab), pointer       :: fine_mask(:) => Null()

       integer                  :: nlevs
       integer                  :: n, dm

       integer, allocatable     :: lo(:), hi(:)
       logical, allocatable     :: nodal(:)

       real(dp_t)               :: eps

       if (present(eps_in)) then
         eps = eps_in
       else
         eps = 1.d-12
       end if

       nlevs = mla%nlevel
       dm = rh(nlevs)%dim
       allocate(nodal(dm), lo(dm), hi(dm))
       nodal = .true.

       allocate(fine_mask(nlevs))

!      We are only considering the dense stencils here (3 in 1d, 9 in 2d, 27 in 3d)

       do n = nlevs, 1, -1

          call lmultifab_build(fine_mask(n), mla%la(n), 1, 0, nodal)
          if ( n < nlevs ) then
             call create_nodal_mask(fine_mask(n), &
                                    mgt(n  )%mm(mgt(n  )%nlevels), &
                                    mgt(n+1)%mm(mgt(n+1)%nlevels), &
                                    ref_ratio(n,:))
          else
             call setval(fine_mask(n), val = .true., all = .true.)
          endif
       end do

       call ml_nd(mla,mgt,rh,full_soln,fine_mask,one_sided_ss,ref_ratio,do_diagnostics,eps)
     
       do n = 1,nlevs
          call lmultifab_destroy(fine_mask(n))
       end do

       deallocate(fine_mask)

   contains

        subroutine create_nodal_mask(mask,mm_crse,mm_fine,ref_ratio)

        type(lmultifab), intent(inout) :: mask
        type(imultifab), intent(in   ) :: mm_crse
        type(imultifab), intent(in   ) :: mm_fine
        integer        , intent(in   ) :: ref_ratio(:)
  
        type(box)        :: cbox,fbox
        logical, pointer :: mkp(:,:,:,:)
        integer, pointer :: cmp(:,:,:,:)
        integer, pointer :: fmp(:,:,:,:)

        integer :: loc(mask%dim),lof(mask%dim)
        integer :: i,j

        call setval(mask,.true.)

!       Note :          mm_fine is  in fine space
!       Note : mask and mm_crse are in crse space
  
        do j = 1,mask%nboxes

           cbox = get_ibox(mask,j)
           loc = lwb(cbox)

           do i = 1,mm_fine%nboxes

              fbox = get_ibox(mm_fine,i)
              lof = lwb(fbox)
              fbox = box_coarsen_v(fbox,ref_ratio)

              if (box_intersects(fbox,cbox)) then
                lo(:) = lwb(box_intersection(cbox,fbox))
                hi(:) = upb(box_intersection(cbox,fbox))

                mkp => dataptr(mask,j)
                cmp => dataptr(mm_crse,j)
                fmp => dataptr(mm_fine,i)
                select case (dm)
                case (2)
                   call create_nodal_mask_2d(mkp(:,:,1,1),cmp(:,:,1,1),loc,fmp(:,:,1,1),lof,lo,hi,ref_ratio)
                case (3)
                   call create_nodal_mask_3d(mkp(:,:,:,1),cmp(:,:,:,1),loc,fmp(:,:,:,1),lof,lo,hi,ref_ratio)
                end select
              end if
           end do
        end do

        end subroutine create_nodal_mask

        subroutine create_nodal_mask_2d(mask,mm_crse,loc,mm_fine,lof,lo,hi,ref_ratio)

             integer, intent(in   ) :: loc(:),lof(:)
             logical, intent(inout) ::    mask(loc(1):,loc(2):)
             integer, intent(inout) :: mm_crse(loc(1):,loc(2):)
             integer, intent(inout) :: mm_fine(lof(1):,lof(2):)
             integer, intent(in   ) :: lo(:),hi(:)
             integer, intent(in   ) :: ref_ratio(:)

             integer :: i,j

             do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (.not.  bc_dirichlet(mm_fine(i*ref_ratio(1),j*ref_ratio(2)),1,0) .or. &
                           bc_dirichlet(mm_crse(i             ,j             ),1,0 ) ) &
                    mask(i,j) = .false.
             end do
             end do


        end subroutine create_nodal_mask_2d

        subroutine create_nodal_mask_3d(mask,mm_crse,loc,mm_fine,lof,lo,hi,ref_ratio)

             integer, intent(in   ) :: loc(:),lof(:)
             logical, intent(inout) ::    mask(loc(1):,loc(2):,loc(3):)
             integer, intent(inout) :: mm_crse(loc(1):,loc(2):,loc(3):)
             integer, intent(inout) :: mm_fine(lof(1):,lof(2):,lof(3):)
             integer, intent(in   ) :: lo(:),hi(:)
             integer, intent(in   ) :: ref_ratio(:)

             integer :: i,j,k,i_fine,j_fine,k_fine

             do k = lo(3),hi(3)
             do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                i_fine = i*ref_ratio(1)
                j_fine = j*ref_ratio(2)
                k_fine = k*ref_ratio(3)
                if (.not. bc_dirichlet(mm_fine(i_fine,j_fine,k_fine),1,0) .or. &
                          bc_dirichlet(mm_crse(i     ,j     ,k     ),1,0) ) &
                    mask(i,j,k) = .false.
             end do
             end do
             end do

        end subroutine create_nodal_mask_3d

   end subroutine ml_nd_solve

! function ml_converged(res, sol, mask, bnorm, Anorm, eps) result(r)
!    logical :: r
!    type(multifab), intent(in) :: res(:), sol(:)
!    type(lmultifab), intent(in) :: mask(:)
!    real(dp_t), intent(in) :: Anorm, eps, bnorm
!    real(dp_t) :: ni_res, ni_sol
!    ni_res = ml_norm_inf(res, mask)
!    ni_sol = ml_norm_inf(sol, mask)
!    r =  ni_res <= eps*(Anorm*ni_sol + bnorm) .or. &
!         ni_res <= spacing(Anorm)
! end function ml_converged

! function ml_norm_inf(rr, mask) result(r)
!    real(dp_t)  :: r
!    type(multifab), intent(in) :: rr(:)
!    type(lmultifab), intent(in) :: mask(:)
!    integer :: n
!    r = 0
!    do n = 1, size(rr)
!       r = max(norm_inf(rr(n),mask(n)), r)
!       end do
! end function ml_norm_inf

! function ml_norm_l2(rr, mask) result(r)
!    real(dp_t)  :: r
!    type(multifab), intent(in) :: rr(:)
!    type(lmultifab), intent(in) :: mask(:)
!    integer :: n
!    r = 0
!    do n = 1, size(rr)
!       r = max(norm_l2(rr(n),mask(n)), r)
!    end do
! end function ml_norm_l2

end module ml_solve_module
