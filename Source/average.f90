module average_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use geometry
  use ml_layout_module
  
  implicit none
  
  ! if we set phibar to 0, then the evolution of the base state is effectively 
  ! turned off.  Here we set a parameter that allows us to do that.  We will 
  ! use the routines enable_base_evolution() and disable_base_evolution() to 
  ! control this.  By default, we have base state evolution enabled.
  logical, private, save :: evolve_base = .true.
  
contains
  
  subroutine enable_base_evolution()
    
    implicit none
    
    evolve_base = .true.
    
  end subroutine enable_base_evolution
  
  
  subroutine disable_base_evolution()
    
    implicit none
    
    evolve_base = .false.
    
  end subroutine disable_base_evolution
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average(mla,phi,phibar,dx,comp,ncomp)
    
    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: comp,ncomp
    type(multifab) , intent(inout) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar(:,0:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    
    ! local
    real(kind=dp_t), pointer     :: pp(:,:,:,:)
    integer                      :: lo(phi(1)%dim),hi(phi(1)%dim)
    integer                      :: i,k,n,nlevs,ng,dm,rr
    real(kind=dp_t), allocatable :: vol_grid(:,:),vol_proc(:,:)
    real(kind=dp_t), allocatable :: vol_tot(:,:),phibar_proc(:,:,:)
    real(kind=dp_t), allocatable :: source_buffer(:,:), target_buffer(:,:)
    
    dm = phi(1)%dim
    ng = phi(1)%ng
    nlevs = size(dx,dim=1)
    
    phibar(:,:,:) = ZERO    

    if (evolve_base) then
       
       if (spherical .eq. 1) then
          allocate(vol_grid(nlevs,0:nr(nlevs)-1))
       end if
       
       allocate(   vol_proc(nlevs,0:nr(nlevs)-1))
       allocate(    vol_tot(nlevs,0:nr(nlevs)-1))
       allocate(phibar_proc(nlevs,0:nr(nlevs)-1,ncomp))

       allocate(source_buffer(nlevs,ncomp))
       allocate(target_buffer(nlevs,ncomp))
       
       vol_proc(:,:) = ZERO
       vol_tot(:,:)  = ZERO
       
       phibar_proc(:,:,:) = ZERO
       

       if(spherical .eq. 0) then

          ! the first step is to compute the average as if the coarsest level were the
          ! only level in existence
          do i = 1, phi(1)%nboxes
             if ( multifab_remote(phi(1), i) ) cycle
             pp => dataptr(phi(1), i)
             lo =  lwb(get_box(phi(1), i))
             hi =  upb(get_box(phi(1), i))
             select case (dm)
             case (2)
                call average_coarsest_2d(pp(:,:,1,:),phibar_proc(1,:,:),lo,hi,ng, &
                                         comp,ncomp,dx(1,:))

                vol_proc(1,lo(2):hi(2)) = vol_proc(1,lo(2):hi(2)) &
                     + (hi(1)-lo(1)+1)*dx(1,1)*dx(1,2)
             case (3)
                   call average_coarsest_3d(pp(:,:,:,:),phibar_proc(1,:,:),lo,hi,ng, &
                                            comp,ncomp,dx(1,:))
                   
                   vol_proc(1,lo(3):hi(3)) = vol_proc(1,lo(3):hi(3)) &
                        + (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*dx(1,1)*dx(1,2)*dx(1,3)
             end select
          end do
          
          do k = 0,nr(1)-1
             call parallel_reduce(vol_tot(1,k),  vol_proc(1,k),      MPI_SUM)
             
             ! put all the components for the current k into a buffer arrar and reduce
             source_buffer(1,:) = phibar_proc(1,k,:)
             call parallel_reduce(target_buffer(1,:), source_buffer(1,:), MPI_SUM)
             phibar(1,k,:) = target_buffer(1,:)
             
             phibar(1,k,:) = phibar(1,k,:) / vol_tot(1,k)
          end do

          ! done computing the average at the coarsest level
          ! now we compute the average at the finer levels
          do n=2,nlevs

             rr = mla%mba%rr(n-1,dm)

             ! we begin by using piecewise linear interpolation to get the average at
             ! a fine level using the coarser data
             do k=0,nr(n)-1
                phibar(n,k,:) = phibar(n-1,int(k/rr),:)
             end do

          end do

       else if(spherical .eq. 1) then

          ! This does not work for multilevel yet.

          do n = 1, nlevs
             do i = 1, phi(n)%nboxes
                if ( multifab_remote(phi(n), i) ) cycle
                pp => dataptr(phi(n), i)
                lo =  lwb(get_box(phi(n), i))
                hi =  upb(get_box(phi(n), i))
                vol_grid(n,:) = ZERO
                call average_3d_sphr(n,pp(:,:,:,:),phibar_proc(n,:,:),lo,hi,ng,dx(n,:), &
                                     vol_grid(n,:),comp,ncomp)
                      
                vol_proc(n,:) = vol_proc(n,:) + vol_grid(n,:)
             end do
             
             do k = 0,nr(n)-1
                call parallel_reduce(vol_tot(n,k),  vol_proc(n,k),      MPI_SUM)
                
                ! put all the components for the current k into a buffer array and reduce
                source_buffer(n,:) = phibar_proc(n,k,:)
                call parallel_reduce(target_buffer(n,:), source_buffer(n,:), MPI_SUM)
                phibar(n,k,:) = target_buffer(n,:)
                
                if (vol_tot(n,k) .gt. ZERO) then
                   phibar(n,k,:) = phibar(n,k,:) / vol_tot(n,k)
                else
                   phibar(n,k,:) = ZERO
                end if
             end do

          enddo

          deallocate(vol_grid)
          
       end if

       deallocate(vol_proc,vol_tot,phibar_proc)
       deallocate(source_buffer, target_buffer)
       
    endif
    
  end subroutine average
  
  subroutine average_coarsest_2d(phi,phibar,lo,hi,ng,start_comp,ncomp,dx)
    
    integer         , intent(in   ) :: lo(:), hi(:), ng, start_comp, ncomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: phibar(0:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)
    
    ! Local variables
    integer          :: i, j, comp
    real (kind=dp_t) :: vol
    
    vol = dx(1)*dx(2)
    
    do comp = start_comp,start_comp+ncomp-1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             phibar(j,comp) = phibar(j,comp) + phi(i,j,comp)
          end do
          phibar(j,comp) = phibar(j,comp) * vol
       end do
    end do
    
  end subroutine average_coarsest_2d
  
  subroutine average_coarsest_3d(phi,phibar,lo,hi,ng,start_comp,ncomp,dx)
    
    integer         , intent(in   ) :: lo(:), hi(:), ng, start_comp, ncomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phibar(0:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)
    
    ! Local variables
    integer          :: i, j, k, comp
    real (kind=dp_t) :: vol
    
    vol = dx(1)*dx(2)*dx(3)
    
    do comp = start_comp,start_comp+ncomp-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                phibar(k,comp) = phibar(k,comp) + phi(i,j,k,comp)
             end do
          end do
          phibar(k,comp) = phibar(k,comp) * vol
       end do
    end do
    
  end subroutine average_coarsest_3d
  
  subroutine average_3d_sphr(n,phi,phibar,lo,hi,ng,dx,sum,start_comp,ncomp)
    
    integer         , intent(in   ) :: n
    integer         , intent(in   ) :: lo(:), hi(:), ng, start_comp, ncomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phibar(0:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)
    real (kind=dp_t), intent(inout) :: sum(0:)
    
    ! Local variables
    integer                       :: i, j, k, comp, index
    integer                       :: ii, jj, kk
    real (kind=dp_t)              :: x,y,z,radius,vol
    real (kind=dp_t)              :: xx, yy, zz
    real (kind=dp_t)              :: xmin, ymin, zmin
    integer :: nsub
    
    ! compute nsub such that we are always guaranteed to fill each of
    ! the base state radial bins
    nsub = int(dx(1)/dr(n)) + 1
    
    do k = lo(3),hi(3)
       zmin = dble(k)*dx(3) - center(3)
       
       do j = lo(2),hi(2)
          ymin = dble(j)*dx(2) - center(2)
          
          do i = lo(1),hi(1)
             xmin = dble(i)*dx(1) - center(1)
             
             do kk = 0, nsub-1
                zz = zmin + (dble(kk) + HALF)*dx(3)/nsub
                
                do jj = 0, nsub-1
                   yy = ymin + (dble(jj) + HALF)*dx(2)/nsub
                   
                   do ii = 0, nsub-1
                      xx = xmin + (dble(ii) + HALF)*dx(1)/nsub
                      
                      radius = sqrt(xx**2 + yy**2 + zz**2)
                      index = radius / dr(n) 
                      
                      if (index .lt. 0 .or. index .gt. nr(n)-1) then
                         print *,'RADIUS ',radius
                         print *,'BOGUS INDEX IN AVERAGE ',index
                         print *,'NOT IN RANGE 0 TO ',nr(n)-1
                         print *,'I J K ',i,j,k
                         print *,'X Y Z ',x,y,z
                         stop
                      end if
                      
                      vol = FOUR3RD*M_PI * dr(n) * &
                           (zl(n,index+1)**2 + zl(n,index+1)*zl(n,index) + zl(n,index)**2)
                      
                      do comp = start_comp,start_comp+ncomp-1
                         phibar(index,comp) = phibar(index,comp) + vol * phi(i,j,k,comp)
                      end do
                      
                      sum(index) = sum(index) + vol
                      
                   enddo
                enddo
             enddo
             
          end do
       end do
    end do
    
  end subroutine average_3d_sphr
  
end module average_module
