module rk_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private
  public :: update_rk

contains

  subroutine update_rk(snew,sold,stemp,szero,step,mla,derivative_mode)

    use bl_prof_module
    use bl_constants_module
    
    integer, intent(in)  :: step
    type(multifab)    , intent(inout) :: snew(:)
    type(multifab)    , intent(inout) :: sold(:)
    type(multifab)    , intent(inout) :: stemp(:)
    type(multifab)    , intent(in   ) :: szero(:)    
    type(ml_layout)   , intent(inout) :: mla
    logical           , intent(in   ), optional :: derivative_mode
 
    ! local
    integer :: i,n
    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: ng_sn,ng_so,ng_st,ng_sz
    logical :: deriv
    real(kind=dp_t) :: fac, dtfac
    real(kind=dp_t), pointer:: snp(:,:,:,:)
    real(kind=dp_t), pointer:: sop(:,:,:,:)
    real(kind=dp_t), pointer:: stp(:,:,:,:)
    real(kind=dp_t), pointer:: szp(:,:,:,:)
        
    dm = mla%dim
    nlevs = mla%nlevel
    
    if (present(derivative_mode)) then
        deriv = derivative_mode
    else
        deriv = .false.
    endif
    
    
    select case(step)
        case(1)
                fac = SIXTH 
                dtfac = HALF
        case(2)
                fac = THIRD
                dtfac = HALF
        case(3)
                fac = THIRD
                dtfac = ONE
        case(4)
                fac = THIRD
                dtfac = ZERO
    end select    
    
    
    ng_sn = nghost(snew(1))
    ng_so = nghost(sold(1))
    ng_st = nghost(stemp(1))
    ng_sz = nghost(szero(1))    

    do n = 1, nlevs

       do i = 1, nfabs(snew(n))
          snp  => dataptr(snew(n),i)
          sop  => dataptr(sold(n),i)
          stp  => dataptr(stemp(n),i)
          szp  => dataptr(szero(n),i)

          lo = lwb(get_box(snew(n),i))
          hi = upb(get_box(snew(n),i))
          select case (dm)
          case (1)
             call rk_update_1d(snp(:,1,1,:),ng_sn,sop(:,1,1,:),ng_so, &
                                stp(:,1,1,:),ng_st,szp(:,1,1,:),ng_sz, &
                                fac,dtfac,lo,hi,deriv)
          case (2)
             call rk_update_2d(snp(:,:,1,:),ng_sn,sop(:,:,1,:),ng_so, &
                                stp(:,:,1,:),ng_st,szp(:,:,1,:),ng_sz, &
                                fac,dtfac,lo,hi,deriv)
          case (3)
             call rk_update_3d(snp(:,:,:,:),ng_sn,sop(:,:,:,:),ng_so, &
                                stp(:,:,:,:),ng_st,szp(:,:,:,:),ng_sz, &
                                fac,dtfac,lo,hi,deriv)
          end select
       end do

    enddo

  end subroutine update_rk

  subroutine rk_update_1d(snew,ng_sn,sold,ng_so,stemp,ng_st,szero,ng_sz, &
                         fac,dtfac,lo,hi,derivative)
    use probin_module, only: rk_timestep_fac 
    
    integer, intent(in) :: lo(:), hi(:), ng_sn, ng_so, ng_st, ng_sz  
    real (kind = dp_t), intent(inout) ::   snew(lo(1)-ng_sn:,:)  
    real (kind = dp_t), intent(inout) ::   sold(lo(1)-ng_so:,:)  
    real (kind = dp_t), intent(inout) ::   stemp(lo(1)-ng_st:,:)  
    real (kind = dp_t), intent(in   ) ::   szero(lo(1)-ng_sz:,:)
    real (kind = dp_t), intent(in   ) :: fac, dtfac
    logical           , intent(in   ) :: derivative
    integer :: i

    do i = lo(1)-ng_sn, hi(1)+ng_sn
            if (.not. derivative) then
                snew(i,:) = snew(i,:) - sold(i,:)
            endif
            stemp(i,:) = stemp(i,:)  + snew(i,:) * fac * rk_timestep_fac 
            sold(i,:) = szero(i,:) + snew(i,:) * dtfac * rk_timestep_fac 
    enddo

  end subroutine rk_update_1d


  subroutine rk_update_2d(snew,ng_sn,sold,ng_so,stemp,ng_st,szero,ng_sz, &
                         fac,dtfac,lo,hi,derivative)
    use probin_module, only: rk_timestep_fac

    integer, intent(in) :: lo(:), hi(:), ng_sn, ng_so, ng_st, ng_sz  
    real (kind = dp_t), intent(inout) ::   snew(lo(1)-ng_sn:,lo(2)-ng_sn:,:)  
    real (kind = dp_t), intent(inout) ::   sold(lo(1)-ng_so:,lo(2)-ng_so:,:)  
    real (kind = dp_t), intent(inout) ::   stemp(lo(1)-ng_st:,lo(2)-ng_st:,:)  
    real (kind = dp_t), intent(in   ) ::   szero(lo(1)-ng_sz:,lo(2)-ng_sz:,:)
    real (kind = dp_t), intent(in   ) :: fac, dtfac
    logical           , intent(in   ) :: derivative

    integer :: i,j
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j = lo(2)-ng_sn, hi(2)+ng_sn
            do i = lo(1)-ng_sn, hi(1)+ng_sn
                    if (.not. derivative) then
                            snew(i,j,:) = snew(i,j,:) - sold(i,j,:)
                    endif
                    stemp(i,j,:) = stemp(i,j,:)  + snew(i,j,:) * fac * rk_timestep_fac
                    sold(i,j,:) = szero(i,j,:) + snew(i,j,:) * dtfac * rk_timestep_fac
            enddo
    enddo
    !$OMP END PARALLEL DO    
  end subroutine rk_update_2d

  subroutine rk_update_3d(snew,ng_sn,sold,ng_so,stemp,ng_st,szero,ng_sz, &
                         fac,dtfac,lo,hi,derivative)
    use probin_module, only: rk_timestep_fac

    integer, intent(in) :: lo(:), hi(:), ng_sn, ng_so, ng_st, ng_sz  
    real (kind = dp_t), intent(inout) ::   snew(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:,:)  
    real (kind = dp_t), intent(inout) ::   sold(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)  
    real (kind = dp_t), intent(inout) ::   stemp(lo(1)-ng_st:,lo(2)-ng_st:,lo(3)-ng_st:,:)  
    real (kind = dp_t), intent(in   ) ::   szero(lo(1)-ng_sz:,lo(2)-ng_sz:,lo(3)-ng_sz:,:)
    real (kind = dp_t), intent(in   ) :: fac, dtfac
    logical           , intent(in   ) :: derivative

    integer :: i,j,k
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3)-ng_sn, hi(3)+ng_sn
            do j = lo(2)-ng_sn, hi(2)+ng_sn
                    do i = lo(1)-ng_sn, hi(1)+ng_sn
                            if (.not. derivative) then
                                    snew(i,j,k,:) = snew(i,j,k,:) - sold(i,j,k,:)
                            endif
                            stemp(i,j,k,:) = stemp(i,j,k,:)  + snew(i,j,k,:) &
                                                * fac   * rk_timestep_fac
                            sold(i,j,k,:) = szero(i,j,k,:) + snew(i,j,k,:) &
                                                * dtfac * rk_timestep_fac
                    enddo
            enddo
    enddo
    !$OMP END PARALLEL DO    
  end subroutine rk_update_3d
  
end module rk_module
