module create_umac_grown_module

  use bl_constants_module
  use sparse_solve_module

  private
  public :: create_umac_grown, create_umac_grown_onesided

contains

  subroutine create_umac_grown(finelev,fine,crse)

    integer       , intent(in   ) :: finelev
    type(multifab), intent(inout) :: fine(:)
    type(multifab), intent(inout) :: crse(:)

    ! local
    integer        :: i,j,k,dm
    integer        :: c_lo(fine(1)%dim),c_hi(fine(1)%dim)
    integer        :: f_lo(fine(1)%dim),f_hi(fine(1)%dim)

    type(fgassoc)  :: fgasc
    type(boxarray) :: f_ba,c_ba,tba
    type(multifab) :: f_mf,c_mf,tcrse,tfine
    type(layout)   :: f_la,c_la,tla

    real(kind=dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "create_umac_grown")

    dm = fine(1)%dim

    ! fill_boundary on level 1
    ! the fill_boundary for levels 2 through nlev occur later in this function
    if (finelev .eq. 2) then
       do i=1,dm
          call multifab_fill_boundary(crse(i))
       end do
    end if

    ! Grab the cached boxarray of all ghost cells not covered by valid region.
    fgasc = layout_fgassoc(fine(1)%la, 1)

    call boxarray_build_copy(f_ba,fgasc%ba)
    call boxarray_build_copy(c_ba,fgasc%ba)

    do i=1,nboxes(f_ba)
       call set_box(c_ba,i,coarsen(get_box(f_ba,i),2))
       call set_box(f_ba,i,refine(get_box(c_ba,i),2))
    end do

    do i=1,dm

       call build(f_la,f_ba,get_pd(fine(i)%la),get_pmask(fine(i)%la))
       call build(c_la,c_ba,get_pd(crse(i)%la),get_pmask(crse(i)%la), &
                  explicit_mapping=get_proc(f_la))

       ! Create c_mf and f_mf on the same proc.
       call build(f_mf,f_la,1,0,fine(i)%nodal)
       call build(c_mf,c_la,1,0,crse(i)%nodal)

       ! Update c_mf with valid and ghost regions from crse.
       call boxarray_build_copy(tba, get_boxarray(crse(i)))
       do j=1,nboxes(tba)
          call set_box(tba,j,grow(get_box(crse(i),j),1))
       end do
       call build(tla, tba, get_pd(crse(i)%la), get_pmask(crse(i)%la), &
                  explicit_mapping=get_proc(crse(i)%la))
       call destroy(tba)
       call build(tcrse, tla, 1, 0, crse(i)%nodal)
       do j=1,nboxes(tcrse)
          if ( remote(tcrse,j) ) cycle
          fp => dataptr(tcrse,  j)
          cp => dataptr(crse(i),j)
          fp = cp
       end do

       call copy(c_mf, tcrse)

       call destroy(tcrse)
       call destroy(tla)

       ! Fill in some of the fine ghost cells from crse.
       do j=1,nboxes(f_mf)
          if ( remote(f_mf,j) ) cycle
          fp => dataptr(f_mf,j)
          cp => dataptr(c_mf,j)
          f_lo = lwb(get_box(f_mf,j))
          f_hi = upb(get_box(f_mf,j))
          c_lo = lwb(get_box(c_mf,j))
          c_hi = upb(get_box(c_mf,j))
          select case(dm)
          case (2)
             call pc_edge_interp_2d(i,f_lo,f_hi,c_lo,c_hi,fp(:,:,1,1),cp(:,:,1,1))
          case (3)
             call pc_edge_interp_3d(i,f_lo,f_hi,c_lo,c_hi,fp(:,:,:,1),cp(:,:,:,1))
          end select
       end do

       call copy(f_mf,fine(i))

       ! Fill in the rest of the fine ghost cells.
       do j=1,nboxes(f_mf)
          if ( remote(f_mf,j) ) cycle
          fp => dataptr(f_mf,j)
          cp => dataptr(c_mf,j)
          f_lo = lwb(get_box(f_mf,j))
          f_hi = upb(get_box(f_mf,j))
          c_lo = lwb(get_box(c_mf,j))
          c_hi = upb(get_box(c_mf,j))
          select case(dm)
          case (2)
             call lin_edge_interp_2d(i,f_lo,f_hi,c_lo,c_hi,fp(:,:,1,1))
          case (3)
             call lin_edge_interp_3d(i,f_lo,f_hi,c_lo,c_hi,fp(:,:,:,1))
          end select
       end do

       call destroy(c_mf)
       call destroy(c_la)

       ! Update ghost regions of fine where they overlap with f_mf.
       call boxarray_build_copy(tba, get_boxarray(fine(i)))
       do j = 1, nboxes(tba)
          call set_box(tba,j,grow(get_box(fine(i),j),1))
       end do
       call build(tla, tba, get_pd(fine(i)%la), get_pmask(fine(i)%la), &
                  explicit_mapping=get_proc(fine(i)%la))
       call destroy(tba)
       call build(tfine, tla, 1, 0, fine(i)%nodal)

       call copy(tfine, f_mf)

       call destroy(f_mf)
       call destroy(f_la)
 
       do j=1,nboxes(tfine)
          if ( remote(tfine,j) ) cycle
          call boxarray_box_diff(tba, get_ibox(tfine,j), get_ibox(fine(i),j))
          do k = 1, nboxes(tba)
             fp => dataptr(fine(i), j, get_box(tba,k))
             cp => dataptr(tfine,   j, get_box(tba,k))
             fp = cp
          end do
          call destroy(tba)
       end do

       call destroy(tfine)
       call destroy(tla)

       ! now fix up umac grown due to the low order interpolation we used
       do j=1,fine(i)%nboxes
          if ( multifab_remote(fine(i), j) ) cycle
          fp => dataptr(fine(i), j)
          f_lo = lwb(get_box(fine(i), j))
          f_hi = upb(get_box(fine(i), j))
          select case(dm)
          case (2)
             call correct_umac_grown_2d(fp(:,:,1,1),f_lo,f_hi,i)
          case (3)

          end select
       end do

       call multifab_fill_boundary(fine(i))
    end do

    call destroy(f_ba)
    call destroy(c_ba)
    call destroy(bpt)

  end subroutine create_umac_grown

  subroutine pc_edge_interp_2d(dir,f_lo,f_hi,c_lo,c_hi,fine,crse)

    integer,         intent(in   ) :: dir,f_lo(:),f_hi(:),c_lo(:),c_hi(:)
    real(kind=dp_t), intent(inout) :: fine(f_lo(1):,f_lo(2):)
    real(kind=dp_t), intent(inout) :: crse(c_lo(1):,c_lo(2):)

    ! local
    integer :: i,ii,j,jj

    if (dir .eq. 1) then

       do j=c_lo(2),c_hi(2)
          do i=c_lo(1),c_hi(1)+1
             do jj=0,1
                fine(2*i,2*j+jj) = crse(i,j)
             end do
          end do
       end do

    else

       do j=c_lo(2),c_hi(2)+1
          do i=c_lo(1),c_hi(1)
             do ii=0,1
                fine(2*i+ii,2*j) = crse(i,j)
             end do
          end do
       end do

    end if    

  end subroutine pc_edge_interp_2d

  subroutine pc_edge_interp_3d(dir,f_lo,f_hi,c_lo,c_hi,fine,crse)

    integer,         intent(in   ) :: dir,f_lo(:),f_hi(:),c_lo(:),c_hi(:)
    real(kind=dp_t), intent(inout) :: fine(f_lo(1):,f_lo(2):,f_lo(3):)
    real(kind=dp_t), intent(inout) :: crse(c_lo(1):,c_lo(2):,c_lo(3):)

    ! local
    integer :: i,ii,j,jj,k,kk

    if (dir .eq. 1) then

       do k=c_lo(3),c_hi(3)
          do j=c_lo(2),c_hi(2)
             do i=c_lo(1),c_hi(1)+1
                do kk=0,1
                   do jj=0,1
                      fine(2*i,2*j+jj,2*k+kk) = crse(i,j,k)
                   end do
                end do
             end do
          end do
       end do

    else if (dir .eq. 2) then

       do k=c_lo(3),c_hi(3)
          do j=c_lo(2),c_hi(2)+1
             do i=c_lo(1),c_hi(1)
                do kk=0,1
                   do ii=0,1
                      fine(2*i+ii,2*j,2*k+kk) = crse(i,j,k)
                   end do
                end do
             end do
          end do
       end do

    else

       do k=c_lo(3),c_hi(3)+1
          do j=c_lo(2),c_hi(2)
             do i=c_lo(1),c_hi(1)
                do jj=0,1
                   do ii=0,1
                      fine(2*i+ii,2*j+jj,2*k) = crse(i,j,k)
                   end do
                end do
             end do
          end do
       end do

    end if    

  end subroutine pc_edge_interp_3d

  subroutine lin_edge_interp_2d(dir,f_lo,f_hi,c_lo,c_hi,fine)

    integer,         intent(in   ) :: dir,f_lo(:),f_hi(:),c_lo(:),c_hi(:)
    real(kind=dp_t), intent(inout) :: fine(f_lo(1):,f_lo(2):)

    ! local
    integer :: i,ii,j,jj

    if (dir .eq. 1) then

       do j=c_lo(2),c_hi(2)
          do i=c_lo(1),c_hi(1)
             do jj=0,1
                fine(2*i+1,2*j+jj) = HALF*(fine(2*i,2*j+jj)+fine(2*i+2,2*j+jj))
             end do
          end do
       end do

    else

       do j=c_lo(2),c_hi(2)
          do i=c_lo(1),c_hi(1)
             do ii=0,1
                fine(2*i+ii,2*j+1) = HALF*(fine(2*i+ii,2*j)+fine(2*i+ii,2*j+2))
             end do
          end do
       end do

    end if    

  end subroutine lin_edge_interp_2d

  subroutine lin_edge_interp_3d(dir,f_lo,f_hi,c_lo,c_hi,fine)

    integer,         intent(in   ) :: dir,f_lo(:),f_hi(:),c_lo(:),c_hi(:)
    real(kind=dp_t), intent(inout) :: fine(f_lo(1):,f_lo(2):,f_lo(3):)

    ! local
    integer :: i,ii,j,jj,k,kk

    if (dir .eq. 1) then

       do k=c_lo(3),c_hi(3)
          do j=c_lo(2),c_hi(2)
             do i=c_lo(1),c_hi(1)
                do kk=0,1
                   do jj=0,1
                      fine(2*i+1,2*j+jj,2*k+kk) = &
                           HALF*(fine(2*i,2*j+jj,2*k+kk)+fine(2*i+2,2*j+jj,2*k+kk))
                   end do
                end do
             end do
          end do
       end do

    else if (dir .eq. 2) then

       do k=c_lo(3),c_hi(3)
          do j=c_lo(2),c_hi(2)
             do i=c_lo(1),c_hi(1)
                do kk=0,1
                   do ii=0,1
                      fine(2*i+ii,2*j+1,2*k+kk) = &
                           HALF*(fine(2*i+ii,2*j,2*k+kk)+fine(2*i+ii,2*j+2,2*k+kk))
                   end do
                end do
             end do
          end do
       end do

    else

       do k=c_lo(3),c_hi(3)
          do j=c_lo(2),c_hi(2)
             do i=c_lo(1),c_hi(1)
                do jj=0,1
                   do ii=0,1
                      fine(2*i+ii,2*j+jj,2*k+1) = &
                           HALF*(fine(2*i+ii,2*j+jj,2*k)+fine(2*i+ii,2*j+jj,2*k+2))
                   end do
                end do
             end do
          end do
       end do

    end if    

  end subroutine lin_edge_interp_3d

  subroutine correct_umac_grown_2d(vel,lo,hi,dir)
    
    integer        , intent(in   ) :: lo(:),hi(:),dir
    real(kind=dp_t), intent(inout) :: vel(lo(1)-1:,lo(2)-1:)

    ! local
    integer         :: i,j
    real(kind=dp_t) :: temp_lo,temp_hi

    if (dir .eq. 1) then

       do j=lo(2)-1,hi(2)+1
          vel(lo(1)-1,j) = TWO*vel(lo(1)-1,j) - vel(lo(1),j)
          vel(hi(1)+2,j) = TWO*vel(hi(1)+2,j) - vel(hi(1)+1,j)
       end do

       do j=lo(2)-1,hi(2)+1
          if (abs(mod(j,2)) .eq. 1) then
             temp_lo = vel(lo(1)-1,j)
             temp_hi = vel(hi(1)+2,j)
             vel(lo(1)-1,j) = (3.d0/4.d0)*vel(lo(1)-1,j) + FOURTH*vel(lo(1)-1,j+1)
             vel(hi(1)+2,j) = (3.d0/4.d0)*vel(hi(1)+2,j) + FOURTH*vel(hi(1)+2,j+1)
          else
             vel(lo(1)-1,j) = (3.d0/4.d0)*vel(lo(1)-1,j) + FOURTH*temp_lo
             vel(hi(1)+2,j) = (3.d0/4.d0)*vel(hi(1)+2,j) + FOURTH*temp_hi
          end if
       end do

       do i=lo(1),hi(1)+1
          vel(i,lo(2)-1) = (3.d0/4.d0)*vel(i,lo(2)-1) + EIGHTH*(vel(i,lo(2))+vel(i,lo(2)+1))
          vel(i,hi(2)+1) = (3.d0/4.d0)*vel(i,hi(2)+1) + EIGHTH*(vel(i,hi(2))+vel(i,hi(2)-1))
       end do

       do j=lo(2)-1,hi(2)+1
          vel(lo(1)-1,j) = HALF*(vel(lo(1)-1,j)+vel(lo(1),j))
          vel(hi(1)+2,j) = HALF*(vel(hi(1)+2,j)+vel(hi(1)+1,j))
       end do

    else if (dir .eq. 2) then

       do i=lo(1)-1,hi(1)+1
          vel(i,lo(2)-1) = TWO*vel(i,lo(2)-1) - vel(i,lo(2))
          vel(i,hi(2)+2) = TWO*vel(i,hi(2)+2) - vel(i,hi(2)+1)
       end do

       do i=lo(1)-1,hi(1)+1
          if (abs(mod(i,2)) .eq. 1) then
             temp_lo = vel(i,lo(2)-1)
             temp_hi = vel(i,hi(2)+2)
             vel(i,lo(2)-1) = (3.d0/4.d0)*vel(i,lo(2)-1) + FOURTH*vel(i+1,lo(2)-1)
             vel(i,hi(2)+2) = (3.d0/4.d0)*vel(i,hi(2)+2) + FOURTH*vel(i+1,hi(2)+2)
          else
             vel(i,lo(2)-1) = (3.d0/4.d0)*vel(i,lo(2)-1) + FOURTH*temp_lo
             vel(i,hi(2)+2) = (3.d0/4.d0)*vel(i,hi(2)+2) + FOURTH*temp_hi
          end if
       end do

       do j=lo(2),hi(2)+1
          vel(lo(1)-1,j) = (3.d0/4.d0)*vel(lo(1)-1,j) + EIGHTH*(vel(lo(1),j)+vel(lo(1)+1,j))
          vel(hi(1)+1,j) = (3.d0/4.d0)*vel(hi(1)+1,j) + EIGHTH*(vel(hi(1),j)+vel(hi(1)-1,j))
       end do

       do i=lo(1)-1,hi(1)+1
          vel(i,lo(2)-1) = HALF*(vel(i,lo(2)-1)+vel(i,lo(2)))
          vel(i,hi(2)+2) = HALF*(vel(i,hi(2)+2)+vel(i,hi(2)+1))

       end do

    end if

  end subroutine correct_umac_grown_2d

  subroutine create_umac_grown_onesided(nlevs,umac)

    integer       , intent(in   ) :: nlevs
    type(multifab), intent(inout) :: umac(:,:)

    integer :: i,n,dm
    integer :: lo(umac(1,1)%dim),hi(umac(1,1)%dim)

    real(kind=dp_t), pointer :: ump(:,:,:,:) 
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)

    dm = size(umac,dim=2)

    ! we only need to do this for fine levels
    do n=2,nlevs
       do i=1,umac(n,1)%nboxes
          if ( multifab_remote(umac(n,1), i) ) cycle
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)
          lo = lwb(get_box(umac(n,1), i))
          hi = upb(get_box(umac(n,1), i))
          select case (dm)
          case (2)
             call create_umac_grown_onesided_2d(ump(:,:,1,1),vmp(:,:,1,1),lo,hi)
          case (3)
             call create_umac_grown_onesided_3d(ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),lo,hi)
          end select
       end do
    end do

  end subroutine create_umac_grown_onesided

  subroutine create_umac_grown_onesided_2d(umac,vmac,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: umac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) :: vmac(lo(1)-1:,lo(2)-1:)

    integer i,j

    do j=lo(2),hi(2)
       umac(lo(1)-1,j) = TWO*umac(lo(1),j) - umac(lo(1)+1,j)
       umac(hi(1)+2,j) = TWO*umac(hi(1)+1,j) - umac(hi(1),j)
    end do
    
    do i=lo(1),hi(1)
       vmac(i,lo(2)-1) = TWO*vmac(i,lo(2)) - vmac(i,lo(2)+1)
       vmac(i,hi(2)+2) = TWO*vmac(i,hi(2)+1) - vmac(i,hi(2))
    end do

    ! use linear interpolation to fill fine level ghost cells needed for 
    ! transverse derivatives
    do i=lo(1),hi(1)+1
       umac(i,lo(2)-1) = TWO*umac(i,lo(2)) - umac(i,lo(2)+1)
       umac(i,hi(2)+1) = TWO*umac(i,hi(2)) - umac(i,hi(2)-1)
    end do
    
    do j=lo(2),hi(2)+1
       vmac(lo(1)-1,j) = TWO*vmac(lo(1),j) - vmac(lo(1)+1,j)
       vmac(hi(1)+1,j) = TWO*vmac(hi(1),j) - vmac(hi(1)-1,j)
    end do

  end subroutine create_umac_grown_onesided_2d

  subroutine create_umac_grown_onesided_3d(umac,vmac,wmac,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(inout) :: vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(inout) :: wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)

    integer i,j,k

    do j=lo(2),hi(2)
       do k=lo(3),hi(3)
          umac(lo(1)-1,j,k) = TWO*umac(lo(1),j,k) - umac(lo(1)+1,j,k)
          umac(hi(1)+2,j,k) = TWO*umac(hi(1)+1,j,k) - umac(hi(1),j,k)
       end do
    end do
    
    do i=lo(1),hi(1)
       do k=lo(3),hi(3)
          vmac(i,lo(2)-1,k) = TWO*vmac(i,lo(2),k) - vmac(i,lo(2)+1,k)
          vmac(i,hi(2)+2,k) = TWO*vmac(i,hi(2)+1,k) - vmac(i,hi(2),k)
       end do
    end do
    
    do i=lo(1),hi(1)
       do j=lo(2),hi(2)
          wmac(i,j,lo(3)-1) = TWO*wmac(i,j,lo(3)) - wmac(i,j,lo(3)+1)
          wmac(i,j,hi(3)+2) = TWO*wmac(i,j,hi(3)+1) - wmac(i,j,hi(3))
       end do
    end do
    
    ! use linear interpolation to fill fine level ghost cells needed for 
    ! transverse derivatives
    do i=lo(1),hi(1)+1
       do j=lo(2),hi(2)
          umac(i,j,lo(3)-1) = TWO*umac(i,j,lo(3)) - umac(i,j,lo(3)+1)
          umac(i,j,hi(3)+1) = TWO*umac(i,j,hi(3)) - umac(i,j,hi(3)-1)
       end do
    end do
    
    do i=lo(1),hi(1)+1
       do k=lo(3),hi(3)
          umac(i,lo(2)-1,k) = TWO*umac(i,lo(2),k) - umac(i,lo(2)+1,k)
          umac(i,hi(2)+1,k) = TWO*umac(i,hi(2),k) - umac(i,hi(2)-1,k)
       end do
    end do
    
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          vmac(i,j,lo(3)-1) = TWO*vmac(i,j,lo(3)) - vmac(i,j,lo(3)+1)
          vmac(i,j,hi(3)+1) = TWO*vmac(i,j,hi(3)) - vmac(i,j,hi(3)-1)
       end do
    end do
    
    do j=lo(2),hi(2)+1
       do k=lo(3),hi(3)
          vmac(lo(1)-1,j,k) = TWO*vmac(lo(1),j,k) - vmac(lo(1)+1,j,k)
          vmac(hi(1)+1,j,k) = TWO*vmac(hi(1),j,k) - vmac(hi(1)-1,j,k)
       end do
    end do

    do k=lo(3),hi(3)+1
       do i=lo(1),hi(1)
          wmac(i,lo(2)-1,k) = TWO*wmac(i,lo(2),k) - wmac(i,lo(2)+1,k)
          wmac(i,hi(2)+1,k) = TWO*wmac(i,hi(2),k) - wmac(i,hi(2)-1,k)
       end do
    end do

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          wmac(lo(1)-1,j,k) = TWO*wmac(lo(1),j,k) - wmac(lo(1)+1,j,k)
          wmac(hi(1)+1,j,k) = TWO*wmac(hi(1),j,k) - wmac(hi(1)-1,j,k)
       end do
    end do

  end subroutine create_umac_grown_onesided_3d

end module create_umac_grown_module
