module fillpatch_module

  use layout_module
  use fab_module
  use bl_mem_stat_module
  use multifab_module
  use bc_module
  use setbc_module
  use interp_module

  implicit none

contains

  subroutine fillpatch(fine,crse,fine_domain,ng,ir,bc,icomp,bc_comp,nc)

    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    type(box     ), intent(in   ) :: fine_domain
    integer       , intent(in   ) :: ng
    integer       , intent(in   ) :: ir(:)
    integer       , intent(in   ) :: bc(:,:,:)
    integer       , intent(in   ) :: icomp,bc_comp,nc

    type(box) :: fbox,cstrip
    type(box) :: crse_domain
    integer   :: i,j,n,dm
    integer   :: iface,ishift
    integer   :: lo(4),hi(4),lo_f(3),lo_c(3),hi_f(3),hi_c(3),interior_lo(4)
    integer   :: fvcx_lo,fvcy_lo,fvcz_lo
    integer   :: cvcx_lo,cvcy_lo,cvcz_lo
    integer   :: cslope_lo(2),cslope_hi(2)
    logical   :: lim_slope, lin_limit

    real(kind=dp_t), allocatable :: cp(:,:,:,:)
    real(kind=dp_t), pointer     :: fp(:,:,:,:)
    real(kind=dp_t), pointer     :: fine_fab(:,:,:,:)
    real(kind=dp_t), allocatable :: fvcx(:),fvcy(:),fvcz(:)
    real(kind=dp_t), allocatable :: cvcx(:),cvcy(:),cvcz(:)
    integer                      :: local_bc(fine%dim,2,nc)

    integer         :: ng_of_crse
    real(kind=dp_t) :: dx(2)
    dx = ONE

    dm = fine%dim

    lo    = 1
    hi    = 1
    lo_c  = 1
    hi_c  = 1
    lo_f  = 1
    hi_f  = 1
    hi(4) = nc

    lim_slope = .true.
    lin_limit = .false.

    local_bc = INTERIOR

    crse_domain = box_coarsen_v(fine_domain,ir)

    do i = 1, fine%nboxes

      fbox   = get_ibox(fine,i)

      fp => dataptr(fine,i)

!       do iface = -1,1,2
!         do d = 1,dm
!            fstrip  = grow(fbox,ng)
!            ishift = iface * (ng + box_extent_d(fbox,d))
!            fstrip  = box_intersection(fstrip,box_shift_d(fstrip,ishift,d))
!            fstrip  = box_intersection(fstrip,fine_domain)

             local_bc = INTERIOR

!            if (.not. box_empty(fbox)) then

                cstrip = box_coarsen_v(fbox,ir)
                interior_lo(1:dm) = lwb(cstrip)

                if (cstrip%lo(1) == crse_domain%lo(1)) then
                   local_bc(1,1,1:nc) = bc(1,1,bc_comp:bc_comp+nc-1)
                end if
                if (cstrip%hi(1) == crse_domain%hi(1)) then
                   local_bc(1,2,1:nc) = bc(1,2,bc_comp:bc_comp+nc-1)
                end if
                if (dm > 1) then
                   if (cstrip%lo(2) == crse_domain%lo(2)) then
                      local_bc(2,1,1:nc) = bc(2,1,bc_comp:bc_comp+nc-1)
                   end if
                   if (cstrip%hi(2) == crse_domain%hi(2)) then
                      local_bc(2,2,1:nc) = bc(2,2,bc_comp:bc_comp+nc-1)
                   end if
                end if
                if (dm > 2) then
                   if (cstrip%lo(3) == crse_domain%lo(3)) then
                      local_bc(3,1,1:nc) = bc(3,1,bc_comp:bc_comp+nc-1)
                   end if
                   if (cstrip%hi(3) == crse_domain%hi(3)) then
                      local_bc(3,2,1:nc) = bc(3,2,bc_comp:bc_comp+nc-1)
                   end if
                end if

                lo_f(1:dm) = lwb(fbox)
                hi_f(1:dm) = upb(fbox)
                allocate(fp(lo_f(1):hi_f(1),lo_f(2):hi_f(2),lo(3):hi(3),lo(4):hi(4)))
                allocate(fvcx(lo_f(1):hi_f(1)+1))
                do j = lo_f(1),hi_f(1)+1
                   fvcx(j) = dble(j)
                end do

                if (dm > 1) then
                   allocate(fvcy(lo_f(2):hi_f(2)+1))
                   do j = lo_f(2),hi_f(2)+1
                      fvcy(j) = dble(j) 
                   end do
                end if
                if (dm > 2) then
                   allocate(fvcz(lo_f(3):hi_f(3)+1))
                   do j = lo_f(3),hi_f(3)+1
                      fvcy(j) = dble(j) 
                   end do
                end if

!               Note : these are deliberately done *before* cstrip is grown.
                cslope_lo(1:dm) = lwb(cstrip)
                cslope_hi(1:dm) = upb(cstrip)

                cstrip = grow(cstrip,1)

                lo_c(1:dm) = lwb(cstrip)
                hi_c(1:dm) = upb(cstrip)
                allocate(cp(lo_c(1):hi_c(1),lo_c(2):hi_c(2),lo_c(3):hi_c(3),lo(4):hi(4)))
                allocate(cvcx(lo_c(1):hi_c(1)+1))
                do j = lo_c(1),hi_c(1)+1
                   cvcx(j) = dble(j) * TWO
                end do
                if (dm > 1) then
                   allocate(cvcy(lo_c(2):hi_c(2)+1))
                   do j = lo_c(2),hi_c(2)+1
                      cvcy(j) = dble(j) * TWO
                   end do
                end if
                if (dm > 2) then
                   allocate(cvcz(lo_c(3):hi_c(3)+1))
                   do j = lo_c(3),hi_c(3)+1
                      cvcz(j) = dble(j) * TWO
                   end do
                end if
 
                call multifab_fab_copy(cp,lo_c,1,crse,icomp,nc)

!               Note: these calls assume that the crse array is one larger than
!                     the coarsened fine array in the lo- AND hi-directions in each dimension.
                select case (dm)
                case (2)

                   ng_of_crse = interior_lo(1) - lo_c(1)
                   do n = 1,nc
                      call setbc_2d(cp(:,:,1,n),interior_lo,ng_of_crse,local_bc(:,:,n),dx,bc_comp+n-1)
                   end do

                   fvcx_lo = lo_f(1)
                   fvcy_lo = lo_f(2)
                   cvcx_lo = lo_c(1)
                   cvcy_lo = lo_c(2)
                   call lin_cc_interp_2d(fp(:,:,1,:), lo_f, cp(:,:,1,:), lo_c, ir, local_bc, &
                                         fvcx, fvcx_lo, fvcy, fvcy_lo, &
                                         cvcx, cvcx_lo, cvcy, cvcy_lo, &
                                         cslope_lo, cslope_hi, lim_slope, lin_limit)
                case (3)
                   do n = 1,nc
                      call setbc_3d(cp(:,:,:,n),interior_lo,1,local_bc(:,:,n),dx,bc_comp+n-1)
                   end do
                   fvcx_lo = lo_f(1)
                   fvcy_lo = lo_f(2)
                   fvcz_lo = lo_f(3)
                   cvcx_lo = lo_c(1)
                   cvcy_lo = lo_c(2)
                   cvcz_lo = lo_c(3)
                   call lin_cc_interp_3d(fp(:,:,:,:), lo_f, cp(:,:,:,:), lo_c, ir, local_bc, &
                                         fvcx, fvcx_lo, fvcy, fvcy_lo, fvcz, fvcz_lo, &
                                         cvcx, cvcx_lo, cvcy, cvcy_lo, cvcz, cvcz_lo, &
                                         cslope_lo, cslope_hi, lim_slope, lin_limit)
                end select
   
                deallocate(cp)
                deallocate(cvcx)
                deallocate(fvcx)
                if (dm > 1) deallocate(cvcy)
                if (dm > 1) deallocate(fvcy)
                if (dm > 2) deallocate(cvcz)
                if (dm > 2) deallocate(fvcz)

                if (dm < 3) fbox%lo(3) = 1
                if (dm < 3) fbox%hi(3) = 1

                fine_fab => dataptr(fine,i)
                fine_fab(fbox%lo(1):fbox%hi(1),fbox%lo(2):fbox%hi(2), &
                         fbox%lo(3):fbox%hi(3),1:nc) = &
                      fp(fbox%lo(1):fbox%hi(1),fbox%lo(2):fbox%hi(2), &
                         fbox%lo(3):fbox%hi(3),1:nc)

!            end if
!        end do
!      end do
     end do

  end subroutine fillpatch

end module fillpatch_module
