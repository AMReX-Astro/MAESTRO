!Incompressible shear jet velocity initialization.

!This module provides for the initialization of fluid velocity data
!
!The 'u' multifab holds this data, including u_x, u_y, and u_z (=0 for 2D)

module init_vel_module
  !Global modules
  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use fill_3d_module
  use eos_module
  use variables
  use network
  use geometry
  use ml_layout_module
  use ml_cc_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initveldata

contains

  subroutine initveldata(u,s0_init,p0_init,dx,bc,mla)
    !--- Variable declaration / initialization ---
    !Args
    type(multifab) , intent(inout) :: u(:)              !Velocity multifab
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)   !Base state 's' multifab
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)     !Base state pressure
    real(kind=dp_t), intent(in   ) :: dx(:,:)           !Grid size array
                                                        !(lev,dim)
    type(bc_level) , intent(in   ) :: bc(:)             !BCs at diff levels
    type(ml_layout), intent(inout) :: mla               !Primary ml_layout.
                                                        !Represents the
                                                        !entire computational
                                                        !grid.

    !Local
    real(kind=dp_t), pointer:: uop(:,:,:,:)             !Velocity pointer
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n,dm,nlevs
   
    !Initialize 
    dm = mla%dim
    nlevs = mla%nlevel
    ng = u(1)%ng        !Number of ghost cells at level 1

    !--- Execution ---
    !Initialize all boxes in the 'u' multifab for all levels
    do n=1,nlevs
       do i = 1, nfabs(u(n))

          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:))
          case (3) 
             if (spherical .eq. 1) then
                call bl_error("ERROR: spherical not implemented for the incompressible jet problem")
             else
                call bl_error("ERROR: 3d not implemented for the incompressible jet problem")
             end if
          end select
       end do
    enddo

    !Use initialized data to fill ghost cells
    if (nlevs .eq. 1) then
       ! Fill ghost cells for any interior boundaries between adjacent grids.
       ! This includes using any periodic boundary conditions at the physical
       ! (non-interior) boundary of the domain.
       ! Note: We can have multiple grids even with just 1 level because the
       !    computational domain is broken up by the MAESTRO algorithm to
       !    facilitate parallel computing.  This may happen even in the case of
       !    serial problems depending on how MAESTRO parameters are initialized.
       call multifab_fill_boundary(u(nlevs))

       ! Use (non-periodic) boundary conditions to fill ghost cells at the
       ! physical boundaries of the computational grid
       call multifab_physbc(u(nlevs),1,1,dm,bc(nlevs))
    else
       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
          ! Diagram (accurate?):
          !     -------------     --------------
          !     | c1  |f1|f2|     | c1  |avg(f)|
          !     |     |f3|f4|     |     |      |
          !     -------------  =  --------------
          !     | c3  | c4  |     | c3  | c4   | 
          !     |     |     |     |     |      |
          !     -------------     --------------
          call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,1,dm,fill_crse_input=.false.)
       enddo
    end if
  end subroutine initveldata

  subroutine initveldata_2d(u,lo,hi,ng,dx,s0_init,p0_init)
    !--- Variable declaration / initialization ---
    !Modules
    use probin_module,       only: prob_lo, prob_hi, yt, delx
    use bl_constants_module, only: M_PI
    
    !Args
    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind=dp_t)  , intent(in   ) :: s0_init(0:,:)
    real (kind=dp_t)  , intent(in   ) :: p0_init(0:)
    
    !Local variables
    integer          :: i, j, n
    real (kind=dp_t) :: x, y, ymid, L_x, yshr1, yshr2

    !--- Execution ---
    ymid = prob_lo(2) + HALF*(prob_hi(2)- prob_lo(2))
    L_x = (prob_hi(1) - prob_lo(1))
    yshr1 = prob_lo(2) + 0.25d0*(prob_hi(2)- prob_lo(2))
    yshr2 = prob_lo(2) + 0.75d0*(prob_hi(2)- prob_lo(2))

    ! initialize the velocity
    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j)+HALF) * dx(2)
       if ( parallel_IOProcessor() .and. n == 1) then
         print *, ' '
         print *, 'y = ', y
         print *, '----------'
       endif

       do i = lo(1), hi(1)
          x = prob_lo(1) + (dble(i)+HALF) * dx(1)

          !x velocity as a function of y
          if(y<ymid) then
             u(i,j,1) = tanh((y - yshr1)/yt) 
          else
             u(i,j,1) = tanh((yshr2 - y)/yt) 
          endif

          !y velocity as a function of x
          u(i,j,2) = delx*sin(2.0d0*M_PI*x/L_x)
       enddo
    enddo
  end subroutine initveldata_2d

  subroutine initveldata_3d(u,lo,hi,ng,dx,s0_init,p0_init)

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)
    
    u = ZERO

  end subroutine initveldata_3d

end module init_vel_module
