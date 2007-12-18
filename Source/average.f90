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
    
    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: comp,ncomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar(:,0:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    
    ! local
    real(kind=dp_t), pointer     :: pp(:,:,:,:)
    type(box)                    :: domain
    integer                      :: domlo(phi(1)%dim),domhi(phi(1)%dim)
    integer                      :: lo(phi(1)%dim),hi(phi(1)%dim)
    integer                      :: i,k,n,nlevs,ng,dm,rr,ncell
    real(kind=dp_t), allocatable :: vol_grid(:,:)
    real(kind=dp_t), allocatable :: vol_proc(:,:)
    real(kind=dp_t), allocatable :: vol_tot(:,:)
    real(kind=dp_t), allocatable :: phisum_proc(:,:,:)
    real(kind=dp_t), allocatable :: phisum(:,:,:)
    real(kind=dp_t), allocatable :: phipert_proc(:,:,:)
    real(kind=dp_t), allocatable :: phipert(:,:,:)
    real(kind=dp_t), allocatable :: source_buffer(:,:)
    real(kind=dp_t), allocatable :: target_buffer(:,:)


    dm = phi(1)%dim
    ng = phi(1)%ng
    nlevs = size(dx,dim=1)
    
    phibar(:,:,:) = ZERO    

    if (evolve_base) then
       
       if (spherical .eq. 1) then
          allocate(vol_grid(nlevs,0:nr(nlevs)-1))
       end if
       
       allocate(vol_proc(nlevs,0:nr(nlevs)-1))
       allocate( vol_tot(nlevs,0:nr(nlevs)-1))

       allocate( phisum_proc(nlevs,0:nr(nlevs)-1,ncomp))
       allocate(      phisum(nlevs,0:nr(nlevs)-1,ncomp))
       allocate(phipert_proc(nlevs,0:nr(nlevs)-1,ncomp))
       allocate(     phipert(nlevs,0:nr(nlevs)-1,ncomp))

       allocate(source_buffer(nlevs,ncomp))
       allocate(target_buffer(nlevs,ncomp))
       
       vol_proc(:,:) = ZERO
       vol_tot(:,:)  = ZERO
       
       phisum_proc(:,:,:) = ZERO
       phisum(:,:,:) = ZERO
       phipert_proc(:,:,:) = ZERO
       phipert(:,:,:) = ZERO

       if(spherical .eq. 0) then

          domain = layout_get_pd(phi(1)%la)
          domlo = lwb(domain)
          domhi = upb(domain)

          if(dm .eq. 2) then
             ncell = domhi(1)-domlo(1)+1
          else if(dm .eq. 3) then
             ncell = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
          end if

          ! the first step is to compute phibar assuming the coarsest level 
          ! is the only level in existence
          do i = 1, phi(1)%nboxes
             if ( multifab_remote(phi(1), i) ) cycle
             pp => dataptr(phi(1), i)
             lo =  lwb(get_box(phi(1), i))
             hi =  upb(get_box(phi(1), i))
             select case (dm)
             case (2)
                call sum_coarsest_2d(pp(:,:,1,:),phisum_proc(1,:,:),lo,hi,ng,comp,ncomp)
             case (3)
                call sum_coarsest_3d(pp(:,:,:,:),phisum_proc(1,:,:),lo,hi,ng,comp,ncomp)
             end select
          end do
          
          do k = 0, nr(1)-1
             ! gather phisum
             source_buffer(1,:) = phisum_proc(1,k,:)
             call parallel_reduce(target_buffer(1,:), source_buffer(1,:), MPI_SUM)
             phisum(1,k,:) = target_buffer(1,:)
             
             phibar(1,k,:) = phisum(1,k,:) / dble(ncell)
          end do

          ! now we compute the phibar at the finer levels
          do n = 2, nlevs

             rr = mla%mba%rr(n-1,dm)

             domain = layout_get_pd(phi(n)%la)
             domlo = lwb(domain)
             domhi = upb(domain)

             if(dm .eq. 2) then
                ncell = domhi(1)-domlo(1)+1
             else if(dm .eq. 3) then
                ncell = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
             end if

             ! compute phisum at next finer level
             ! begin by assuming piecewise linear interpolation
             do k = 0, nr(n)-1
                phisum(n,k,:) = phisum(n-1,k/rr,:)*rr**(dm-1)
             end do

             ! compute phipert_proc
             do i = 1, phi(n)%nboxes
                if ( multifab_remote(phi(n), i) ) cycle
                pp  => dataptr(phi(n), i)
                lo =  lwb(get_box(phi(n), i))
                hi =  upb(get_box(phi(n), i))
                select case (dm)
                case (2)
                   call compute_phipert_2d(pp(:,:,1,:),phipert_proc(n,:,:), &
                                           lo,hi,ng,comp,ncomp,rr)
                case (3)
                   call compute_phipert_3d(pp(:,:,:,:),phipert_proc(n,:,:), &
                                           lo,hi,ng,comp,ncomp,rr)
                end select
             end do

             do k = 0, nr(n)-1
                ! gather phipert
                source_buffer(n,:) = phipert_proc(n,k,:)
                call parallel_reduce(target_buffer(n,:), source_buffer(n,:), MPI_SUM)
                phipert(n,k,:) = target_buffer(n,:)

                ! update phisum and compute phibar
                phisum(n,k,:) = phisum(n,k,:) + phipert(n,k,:)
                phibar(n,k,:) = phisum(n,k,:) / dble(ncell)
             end do

          end do

       else if(spherical .eq. 1) then

          ! note: spherical does not work for multilevel yet.

          do n = 1, nlevs
             do i = 1, phi(n)%nboxes
                if ( multifab_remote(phi(n), i) ) cycle
                pp => dataptr(phi(n), i)
                lo =  lwb(get_box(phi(n), i))
                hi =  upb(get_box(phi(n), i))
                vol_grid(n,:) = ZERO
                call average_3d_sphr(n,pp(:,:,:,:),phisum_proc(n,:,:),lo,hi,ng,dx(n,:), &
                                     vol_grid(n,:),comp,ncomp)
                      
                vol_proc(n,:) = vol_proc(n,:) + vol_grid(n,:)
             end do
             
             do k = 0, nr(n)-1
                call parallel_reduce(vol_tot(n,k),  vol_proc(n,k),      MPI_SUM)
                
                ! put all the components for the current k into a buffer array and reduce
                source_buffer(n,:) = phisum_proc(n,k,:)
                call parallel_reduce(target_buffer(n,:), source_buffer(n,:), MPI_SUM)
                phisum(n,k,:) = target_buffer(n,:)
                
                if (vol_tot(n,k) .gt. ZERO) then
                   phibar(n,k,:) = phisum(n,k,:) / vol_tot(n,k)
                else
                   phibar(n,k,:) = ZERO
                end if
             end do

          enddo

          deallocate(vol_grid)
          
       end if

       deallocate(vol_proc,vol_tot)
       deallocate(phisum_proc,phisum)
       deallocate(phipert_proc,phipert)
       deallocate(source_buffer,target_buffer)
       
    endif
    
  end subroutine average
  
  subroutine sum_coarsest_2d(phi,phibar,lo,hi,ng,start_comp,ncomp)
    
    integer         , intent(in   ) :: lo(:), hi(:), ng, start_comp, ncomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: phibar(0:,:)
    
    ! Local variables
    integer          :: i, j, comp
    
    do comp = start_comp,start_comp+ncomp-1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             phibar(j,comp) = phibar(j,comp) + phi(i,j,comp)
          end do
       end do
    end do
    
  end subroutine sum_coarsest_2d
  
  subroutine sum_coarsest_3d(phi,phibar,lo,hi,ng,start_comp,ncomp)
    
    integer         , intent(in   ) :: lo(:), hi(:), ng, start_comp, ncomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phibar(0:,:)
    
    ! Local variables
    integer          :: i, j, k, comp
    
    do comp = start_comp,start_comp+ncomp-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                phibar(k,comp) = phibar(k,comp) + phi(i,j,k,comp)
             end do
          end do
       end do
    end do
    
  end subroutine sum_coarsest_3d

  subroutine compute_phipert_2d(phi,phipert,lo,hi,ng,start_comp,ncomp,rr)

    integer         , intent(in   ) :: lo(:), hi(:), ng, start_comp, ncomp, rr
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: phipert(0:,:)

    ! Local variables
    integer          :: i, j, icrse, jcrse, comp
    real (kind=dp_t) :: crseval

    do comp = start_comp,start_comp+ncomp-1

       ! loop over coarse cell index
       do jcrse = lo(2)/rr,hi(2)/rr
          do icrse = lo(1)/rr,hi(1)/rr

             crseval = ZERO

             ! compute coarse cell value by taking average of fine cells
             do j=0,rr-1
                do i=0,rr-1
                   crseval = crseval + phi(icrse*rr+i,jcrse*rr+j,comp)
                end do
             end do
             crseval = crseval / dble(rr**2)

             ! compute phipert
             do j=0,rr-1
                do i=0,rr-1
                   phipert(jcrse*rr+j,comp) = phipert(jcrse*rr+j,comp) &
                        + phi(icrse*rr+i,jcrse*rr+j,comp) - crseval
                end do
             end do

          end do
       end do

    end do

  end subroutine compute_phipert_2d

  subroutine compute_phipert_3d(phi,phipert,lo,hi,ng,start_comp,ncomp,rr)

    integer         , intent(in   ) :: lo(:), hi(:), ng, start_comp, ncomp, rr
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phipert(0:,:)

    ! Local variables
    integer          :: i, j, k, icrse, jcrse, kcrse, comp
    real (kind=dp_t) :: crseval
    
    do comp = start_comp,start_comp+ncomp-1

       ! loop over coarse cell index
       do kcrse = lo(3)/rr,hi(3)/rr
          do jcrse = lo(2)/rr,hi(2)/rr
             do icrse = lo(1)/rr,hi(1)/rr

                crseval = ZERO

                ! compute coarse cell value by taking average of fine cells
                do k=0,rr-1
                   do j=0,rr-1
                      do i=0,rr-1
                         crseval = crseval + phi(icrse*rr+i,jcrse*rr+j,kcrse*rr+k,comp)
                      end do
                   end do
                end do
                crseval = crseval / dble(rr**3)

                ! compute phipert
                do k=0,rr-1
                   do j=0,rr-1
                      do i=0,rr-1
                         phipert(kcrse*rr+k,comp) = phipert(kcrse*rr+k,comp) &
                              + phi(icrse*rr+i,jcrse*rr+j,kcrse*rr+k,comp) - crseval
                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

  end subroutine compute_phipert_3d

  
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
