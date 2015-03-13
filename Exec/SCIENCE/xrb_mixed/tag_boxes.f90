module tag_boxes_module

  use BoxLib
  use f2kcli
  use list_box_module
  use boxarray_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use box_util_module
  use bl_IO_module
  use cluster_module
  use ml_layout_module

  implicit none 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MUST set this to .true. if tagging uses ghost cells (e.g., tagging on gradient). !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical, save :: tagging_needs_ghost_cells = .false.

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    use variables, only: rho_comp, spec_comp
    use geometry, only: nr_fine, nr
    use probin_module, only: do_dens_tagging
    use network, only: network_species_index

    type( multifab)          , intent(in   ) :: mf
    type(lmultifab)          , intent(inout) :: tagboxes
    real(dp_t)               , intent(in   ) :: dx
    integer                  , intent(in   ) :: lev
    type( multifab), optional, intent(in   ) :: aux_tag_mf

    real(kind = dp_t), pointer :: sp(:,:,:,:),auxp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, j, lo(get_dim(mf)), ng_s, dm, ng_aux
    logical           ::      radialtag(0:nr_fine-1)
    logical           :: radialtag_proc(0:nr_fine-1)
    integer, parameter :: npad = 4

    integer :: ih1


    radialtag = .false.
    radialtag_proc = .false.

    dm = get_dim(mf)

    ih1 = network_species_index("hydrogen-1")

    ng_s = mf%ng
    if (present(aux_tag_mf)) ng_aux = aux_tag_mf%ng

    do i = 1, nfabs(mf)
       sp => dataptr(mf, i)
       if (present(aux_tag_mf)) auxp => dataptr(aux_tag_mf, i)
       lo =  lwb(get_box(tagboxes, i))
       select case (dm)
       case (2)

             if (present(aux_tag_mf) .and. (.not. do_dens_tagging)) then
                call update_radialtag_2d(radialtag_proc,&
                                         sp(:,:,1,rho_comp), ng_s, &
                                         auxp(:,:,1,1), ng_aux, lo)
             else
                call radialtag_2d(radialtag_proc, &
                                  sp(:,:,1,spec_comp-1+ih1),sp(:,:,1,rho_comp), &
                                  lo,ng_s,lev)
             endif
                  
       case  (3)

             if (present(aux_tag_mf) .and. (.not. do_dens_tagging)) then
                call update_radialtag_3d(radialtag_proc, &
                                         sp(:,:,:,rho_comp), ng_s, &
                                         auxp(:,:,:,1), ng_aux, lo)
             else
                call radialtag_3d(radialtag_proc, &
                                  sp(:,:,:,spec_comp-1+ih1),sp(:,:,:,rho_comp), &
                                  lo,ng_s,lev)
             endif

       end select
    end do


    ! gather radialtag
    call parallel_reduce(radialtag, radialtag_proc, MPI_LOR)

    ! apply some padding                                                                                                  
    do j = 1, npad

       ! pad the start of a tagged region                                                                                 
       do i = 1, nr(lev)-1
          if (radialtag(i) .and. .not. radialtag(i-1)) then
             ! found start of a tagged region                                                                             
             radialtag(i-1) = .true.
          endif
       enddo

       ! pad the end of a tagged region                                                                                   
       do i = nr(lev)-1, 1, -1
          if (radialtag(i) .and. .not. radialtag(i+1)) then
             ! found end of a tagged region                                                                               
             radialtag(i+1) = .true.
          endif
          enddo

    enddo


    do i = 1, nfabs(mf)
       tp => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       select case (dm)
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),radialtag,lo)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),radialtag,lo)
       end select
    end do

  end subroutine tag_boxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine radialtag_2d(radialtag,spec,rho,lo,ng,lev)

    use probin_module, only: tag_minval, tag_maxval, tag_xfac, &
                             base_cutoff_density, &
                             do_dens_tagging, lo_dens_tag, hi_dens_tag, &
                             dens_tag_lev_fac, max_levs

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) ::  spec(lo(1)-ng:,lo(2)-ng:)
    real(kind = dp_t), intent(in   ) :: rho(lo(1)-ng:,lo(2)-ng:)
    integer, optional, intent(in   ) :: lev

    ! local
    integer :: i,j,nx,ny
    real(kind=dp_t) :: Xspec
    integer :: llev

    llev = 1; if (present(lev)) llev = lev

    nx = size(spec,dim=1) - 2*ng
    ny = size(spec,dim=2) - 2*ng

    do j = lo(2),lo(2)+ny-1
       do i = lo(1),lo(1)+nx-1
          if (.not. do_dens_tagging) then
             Xspec = spec(i,j)/rho(i,j)
             if (Xspec .gt. tag_xfac*tag_minval .and. Xspec .lt. tag_maxval &
                  .and. rho(i,j) .gt. base_cutoff_density) then
                radialtag(j) = .true.
             endif
          else ! tag based on density

             if (rho(i,j) >= lo_dens_tag*dens_tag_lev_fac**(llev-max_levs) .and. &
                 rho(i,j) <= hi_dens_tag*dens_tag_lev_fac**(max_levs-llev)) then
                radialtag(j) = .true.
             endif

          endif
             
       end do
    enddo

  end subroutine radialtag_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine radialtag_3d(radialtag,spec,rho,lo,ng,lev)

    use probin_module, only: tag_minval, tag_maxval, tag_xfac, &
                             base_cutoff_density

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) ::  spec(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind = dp_t), intent(in   ) :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer, optional, intent(in   ) :: lev

    ! local
    integer :: i,j,k,nx,ny,nz
    real(kind=dp_t) :: Xspec

    nx = size(spec,dim=1) - 2*ng
    ny = size(spec,dim=2) - 2*ng
    nz = size(spec,dim=3) - 2*ng

    do k = lo(3),lo(3)+nz-1
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             Xspec = spec(i,j,k)/rho(i,j,k)
             if (Xspec .gt. tag_xfac*tag_minval .and. Xspec .lt. tag_maxval &
                  .and. rho(i,j,k) .gt. base_cutoff_density) then
                radialtag(k) = .true.
             endif
          end do
       enddo
    end do

  end subroutine radialtag_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_radialtag_2d(radialtag,rho,ng,aux,ng_aux,lo)

    use probin_module, only: tag_minval, base_cutoff_density

    integer          , intent(in   ) :: lo(:),ng,ng_aux
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) :: rho(lo(1)-ng:,lo(2)-ng:)
    real(kind = dp_t), intent(in   ) :: aux(lo(1)-ng_aux:,lo(2)-ng_aux:)

    ! local
    integer :: i, j, nx, ny
    real(kind = dp_t) :: Hdot

    nx = size(rho,dim=1) - 2*ng
    ny = size(rho,dim=2) - 2*ng

    do j = lo(2), lo(2)+ny-1
       do i = lo(1),lo(1)+nx-1
          Hdot = aux(i,j)/rho(i,j)
          if (Hdot .gt. tag_minval .and. rho(i,j) .gt. base_cutoff_density) &
               radialtag(j) = .true.
       enddo
    enddo

  end subroutine update_radialtag_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_radialtag_3d(radialtag,rho,ng,aux,ng_aux,lo)

    use probin_module, only: tag_minval, base_cutoff_density

    integer          , intent(in   ) :: lo(:),ng,ng_aux
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind = dp_t), intent(in   ) :: aux(lo(1)-ng_aux:,lo(2)-ng_aux:,lo(3)-ng_aux:)

    ! local
    integer :: i, j, k, nx, ny, nz
    real(kind = dp_t) :: Hdot

    nx = size(rho,dim=1) - 2*ng
    ny = size(rho,dim=2) - 2*ng
    nz = size(rho,dim=3) - 2*ng

    do k = lo(3), lo(3)+nz-1
       do j = lo(2), lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             Hdot = aux(i,j,k)/rho(i,j,k)
             if (Hdot .gt. tag_minval .and. &
                  rho(i,j,k) .gt. base_cutoff_density) &
                  radialtag(k) = .true.
          enddo
       enddo
    enddo

  end subroutine update_radialtag_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes_2d(tagbox,radialtag,lo)

    integer          , intent(in   ) :: lo(:)
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
    logical          , intent(in   ) :: radialtag(0:)

    integer :: j,ny

    ny = size(tagbox,dim=2)

    tagbox = .false.

    ! tag all boxes with radialtag = .true
    do j = lo(2),lo(2)+ny-1
       tagbox(:,j) = radialtag(j)
    enddo

  end subroutine tag_boxes_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes_3d(tagbox,radialtag,lo)

    integer          , intent(in   ) :: lo(:)
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    logical          , intent(in   ) :: radialtag(0:)

    integer :: k,nz

    nz = size(tagbox,dim=3)

    tagbox = .false.

    ! tag all boxes with radialtag = .true.
    do k = lo(3),lo(3)+nz-1
       tagbox(:,:,k) = radialtag(k)
    end do

  end subroutine tag_boxes_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tag_boxes_module
