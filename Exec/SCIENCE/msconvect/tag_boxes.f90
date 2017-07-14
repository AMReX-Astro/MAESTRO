!!!!!!!!!!!!!!!!!!!
!Tagging on the velocity magnitude. We allow for 4 indidual tagging velocities (tag_vel_1-4) for the first 4 tagging levels,
! if one of the values is negative then or we go deeper than 5 levels we use tag_vel_step to determine the next tagging velocities  
!!!!!!!!!!!!!!!!!!
module tag_boxes_module

  use multifab_module
  use bl_error_module
  use bl_constants_module
  
  implicit none 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MUST set this to .true. if tagging uses ghost cells (e.g., tagging on gradient). !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical, save :: tagging_needs_ghost_cells = .false.

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)
    
    use bl_constants_module
    use variables, only: rho_comp
    use parallel
    use multifab_module
    
    type( multifab)         , intent(in   ) :: mf
    type(lmultifab)         , intent(inout) :: tagboxes
    real(dp_t)              , intent(in   ) :: dx
    integer                 , intent(in   ) :: lev
    type(multifab), optional, intent(in   ) :: aux_tag_mf
    ! aux_tag_mf allows user to pass in additional multifabs for tagging logic
    
    ! local variables
    real(kind = dp_t), pointer :: mfp(:,:,:,:)
    real(kind = dp_t), pointer :: up(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i
    integer           :: lo(get_dim(mf)), hi(get_dim(mf))
    integer           :: tlo(4), mflo(4),uplo(4)
    type(mfiter)      :: mfi
    type(box)         :: bx
    real(kind = dp_t) :: maxvel, maxvel_grid, maxvel_proc
    real(kind = dp_t) :: maxdens, maxdens_grid, maxdens_proc
    real(kind = dp_t) :: dts_local, dts_global
    
    
    maxvel = ZERO
    
    if (present(aux_tag_mf)) then    
    
      call mfiter_build(mfi,tagboxes,.true.)
      
      maxvel_proc = ZERO
       
      do while(next_tile(mfi,i))
        bx = get_tilebox(mfi)
        lo =  lwb(bx)
        hi =  upb(bx)
        up => dataptr(aux_tag_mf,i)
       
        uplo = lbound(up)
        maxvel_grid = ZERO
       
        select case (get_dim(mf))
        case (2)
          call find_maxvel_2d(up(:,:,1,:),uplo(1:2),maxvel_grid,lo,hi)
        case (3)
          call find_maxvel_3d(up(:,:,:,:),uplo(1:3),maxvel_grid,lo,hi)
        end select 
        
        maxvel_proc = max(maxvel_grid,maxvel_proc)
      end do

      ! This sets maxvel to be the max of maxvel_proc over all processors.
      dts_local = maxvel_proc
      call parallel_reduce( dts_global,  dts_local, MPI_MAX)
      maxvel = dts_global
      
      
      do while(next_tile(mfi,i))
        bx = get_tilebox(mfi)
        lo =  lwb(bx)
        hi =  upb(bx)
        tp  => dataptr(tagboxes, i)
        up => dataptr(aux_tag_mf,i)

        uplo = lbound(up)
        tlo = lbound(tp)
         
       select case (get_dim(mf))
        case (2)
            call tag_boxes_2d(tp(:,:,1,1),tlo(1:2),up(:,:,1,:),uplo(1:2),maxvel,lo,hi,dx,lev)
        case  (3)
            call tag_boxes_3d(tp(:,:,:,1),tlo(1:3),up(:,:,:,:),uplo(1:3),maxvel,lo,hi,dx,lev)
        end select
      end do
    
    else !no velocity multifab given (initialisation) -> we do density tagging in a similar fashion to the velocity tagging. 
      call mfiter_build(mfi,tagboxes,.true.)

      maxdens_proc = ZERO
       
      do while(next_tile(mfi,i))
        bx = get_tilebox(mfi)
        lo =  lwb(bx)
        hi =  upb(bx)
        mfp => dataptr(mf,i)
        
        mflo = lbound(mfp)
        
        maxdens_grid = ZERO
       
        select case (get_dim(mf))
        case (2)
          call find_maxdens_2d(mfp(:,:,1,rho_comp),mflo(1:2),maxdens_grid,lo,hi)
        case (3)
          call find_maxdens_3d(mfp(:,:,:,rho_comp),mflo(1:3),maxdens_grid,lo,hi)
        end select 
        
        maxdens_proc = max(maxdens_grid,maxdens_proc)
        

      end do

      ! This sets maxvel to be the max of maxvel_proc over all processors.
      dts_local = maxdens_proc
      call parallel_reduce( dts_global,  dts_local, MPI_MAX)
      maxdens = dts_global
       
      do while(next_tile(mfi,i))
        bx = get_tilebox(mfi)
        lo =  lwb(bx)
        hi =  upb(bx)
        mfp => dataptr(mf, i)
        tp  => dataptr(tagboxes, i)
        

        mflo = lbound(mfp)
        tlo = lbound(tp)
         
       select case (get_dim(mf))
        case (2)
            call tag_boxes_init_2d(tp(:,:,1,1),tlo(1:2),mfp(:,:,1,rho_comp),mflo(1:2),maxdens,lo,hi,dx,lev)
        case  (3)
            call tag_boxes_init_3d(tp(:,:,:,1),tlo(1:3),mfp(:,:,:,rho_comp),mflo(1:3),maxdens,lo,hi,dx,lev)
        end select
      end do


    endif
    
  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,tlo,up,uplo,maxvel,lo,hi,dx,lev)
    
    use probin_module, only: tag_val_steep,tag_vel_extreme
    integer          , intent(in   ) :: lo(2),hi(2),tlo(2),uplo(2)
    logical          , intent(inout) :: tagbox( tlo(1):, tlo(2):)
    real(kind = dp_t), intent(in   ) ::     up(uplo(1):,uplo(2):,:)
    real(dp_t)       , intent(in   ) :: dx
    real(dp_t)       , intent(in   ) :: maxvel
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j
    integer :: tag_lev_switch
    real(kind = dp_t) :: tag_vel
    real(kind = dp_t)  :: magvel

    ! initially say that we do not want to tag any cells for refinement
    tagbox(lo(1):hi(1),lo(2):hi(2)) = .false.
    
    
    tag_vel = maxvel - maxvel/(dble(lev)*tag_val_steep)
    tag_vel = min(tag_vel,tag_vel_extreme)
    
    !$OMP PARALLEL DO PRIVATE(i,j,magvel)
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)
            
            magvel = sqrt(up(i,j,1)**2 + up(i,j,2)**2)
            if (magvel .ge. tag_vel) then
              tagbox(i,j) = .true.
            endif
            
      enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine tag_boxes_2d

  subroutine tag_boxes_3d(tagbox,tlo,up,uplo,maxvel,lo,hi,dx,lev)

    use probin_module, only: tag_val_steep,tag_vel_extreme
    
    integer          , intent(in   ) :: lo(3),hi(3),tlo(3),uplo(3)
    logical          , intent(inout) :: tagbox( tlo(1):, tlo(2):, tlo(3):)
    real(kind = dp_t), intent(in   ) ::     up(uplo(1):,uplo(2):,uplo(3):,:)
    real(dp_t)       , intent(in   ) :: dx
    real(dp_t)       , intent(in   ) :: maxvel
    integer          , intent(in   ) :: lev
   

    ! local variables
    integer :: i,j,k
    real(kind = dp_t)  :: tag_vel
    real(kind = dp_t)  :: magvel
     ! initially say that we do not want to tag any cells for refinement
    tagbox(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = .false.

    tag_vel = maxvel - maxvel/(dble(lev)*tag_val_steep)
    tag_vel = min(tag_vel,tag_vel_extreme)    
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,magvel)
    do k = lo(3),hi(3)
      do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            
              magvel = sqrt(up(i,j,k,1)**2 + up(i,j,k,2)**2 + up(i,j,k,3)**2)            
              if (magvel .ge. tag_vel) then
                tagbox(i,j,k) = .true.
              endif
            
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO    
  end subroutine tag_boxes_3d




  subroutine tag_boxes_init_2d(tagbox,tlo,mf_rho,mflo,maxdens,lo,hi,dx,lev)
    
    use probin_module, only: tag_val_steep
    
    integer          , intent(in   ) :: lo(2),hi(2),tlo(2), mflo(2)
    logical          , intent(inout) :: tagbox( tlo(1):, tlo(2):)
    real(kind = dp_t), intent(in   ) ::     mf_rho(mflo(1):,mflo(2):)
    real(dp_t)       , intent(in   ) :: dx
    real(dp_t)       , intent(in   ) :: maxdens
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j
    real(kind = dp_t) :: tag_dens

    ! initially say that we do not want to tag any cells for refinement
    tagbox(lo(1):hi(1),lo(2):hi(2)) = .false.
    
    
    tag_dens = maxdens - maxdens/(dble(lev)*tag_val_steep)
    
    
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)
            
            if (mf_rho(i,j) .ge. tag_dens) then
              tagbox(i,j) = .true.
            endif
            
      enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine tag_boxes_init_2d

  subroutine tag_boxes_init_3d(tagbox,tlo,mf_rho,mflo,maxdens,lo,hi,dx,lev)

    use probin_module, only: tag_val_steep
    
    integer          , intent(in   ) :: lo(3),hi(3),tlo(3),mflo(3)
    logical          , intent(inout) :: tagbox( tlo(1):, tlo(2):, tlo(3):)
    real(kind = dp_t), intent(in   ) ::     mf_rho(mflo(1):,mflo(2):,mflo(3):)
    real(dp_t)       , intent(in   ) :: dx
    real(dp_t)       , intent(in   ) :: maxdens
    integer          , intent(in   ) :: lev
   

    ! local variables
    integer :: i,j,k
    integer :: tag_lev_switch
    real(kind = dp_t)  :: tag_dens
     ! initially say that we do not want to tag any cells for refinement
    tagbox(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = .false.


    tag_dens = maxdens - maxdens/(dble(lev)*tag_val_steep)
        
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
      do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            
              if (mf_rho(i,j,k) .ge. tag_dens) then
                tagbox(i,j,k) = .true.
              endif
            
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO    
  end subroutine tag_boxes_init_3d





  subroutine find_maxvel_2d(up,uplo,maxvel,lo,hi)
    
    use bl_constants_module
    integer          , intent(in   ) :: lo(2),hi(2),uplo(2)
    real(kind = dp_t), intent(in   ) :: up(uplo(1):,uplo(2):,:)
    real(dp_t)       , intent(inout) :: maxvel
    
    ! local variables
    integer :: i,j
    real(kind = dp_t)  :: magvel
    
    
    
    maxvel = ZERO
    
    !$OMP PARALLEL DO PRIVATE(i,j,magvel) REDUCTION(MAX : maxvel)
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)
            
            magvel = sqrt(up(i,j,1)**2 + up(i,j,2)**2)            
            maxvel = max(maxvel, magvel)
            
      enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine find_maxvel_2d


  subroutine find_maxvel_3d(up,uplo,maxvel,lo,hi)
    
    use bl_constants_module
    integer          , intent(in   ) :: lo(3),hi(3),uplo(3)
    real(kind = dp_t), intent(in   ) :: up(uplo(1):,uplo(2):,uplo(3):,:)
    real(dp_t)       , intent(inout) :: maxvel
    
    ! local variables
    integer :: i,j,k
    real(kind = dp_t)  :: magvel
    
    
    maxvel = ZERO
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,magvel) REDUCTION(MAX : maxvel)
    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              
              magvel = sqrt(up(i,j,k,1)**2 + up(i,j,k,2)**2 + up(i,j,k,3)**2)          
              maxvel = max(maxvel, magvel)
              
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine find_maxvel_3d


  subroutine find_maxdens_2d(sn,slo,maxdens,lo,hi)
    
    use bl_constants_module
    integer          , intent(in   ) :: lo(2),hi(2),slo(2)
    real(kind = dp_t), intent(in   ) :: sn(slo(1):,slo(2):)
    real(dp_t)       , intent(out  ) :: maxdens
    
    ! local variables
    integer :: i,j
    
    
    
    maxdens = ZERO
    
    !$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(MAX : maxdens)
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)
            
            maxdens = max(maxdens, sn(i,j))
            
      enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine find_maxdens_2d


  subroutine find_maxdens_3d(sn,slo,maxdens,lo,hi)
    
    use bl_constants_module
    integer          , intent(in   ) :: lo(3),hi(3),slo(3)
    real(kind = dp_t), intent(in   ) :: sn(slo(1):,slo(2):,slo(3):)
    real(dp_t)       , intent(out  ) :: maxdens
    
    ! local variables
    integer :: i,j,k
    
    
    maxdens = ZERO
    
    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : maxdens)
    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              
              maxdens = max(maxdens, sn(i,j,k))
              
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine find_maxdens_3d



end module tag_boxes_module
