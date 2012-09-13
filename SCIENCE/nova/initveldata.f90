! Adapted from the initdata.f90 file in the xrb/ directory

module init_vel_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use eos_module
  use variables
  use network
  use geometry, only: spherical
  use ml_layout_module
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none

  real(dp_t), save :: pert_height

  private
  public :: initveldata

contains

  subroutine initveldata(u,s0_init,p0_init,dx,bc,mla)
    
    use probin_module, only: prob_lo, prob_hi, num_vortices

    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng,dm,nlevs
    integer :: i,n

    real(kind=dp_t) :: xloc_vortices(num_vortices)
    real(kind=dp_t) :: offset

    dm = mla%dim
    nlevs = mla%nlevel

    if (mod(num_vortices,2) == 1) then
       call bl_error("ERROR: num_vortices must be even")
    endif

    ! for now, this is calculated even if we don't use the velocity field in
    ! the initveldata_2d routine below
    offset = (prob_hi(1) - prob_lo(1)) / (num_vortices)

    do i = 1, num_vortices
       xloc_vortices(i) = (dble(i-1) + 0.5d0) * offset + prob_lo(1)
    enddo
    
    ng = nghost(u(1))

    do n=1,nlevs
       do i = 1, nfabs(u(n))
          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:), xloc_vortices)
          ! 3d doesn't currently have vortice information coded !
          case (3) 
             call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:))
          end select
       end do
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(u(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(u(nlevs),1,1,dm,bc(nlevs))
    else
    
       ! the loop over nlevs must count backwards to make sure the finer grids
       ! are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering 
          ! it
          call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))
          
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,1,dm, &
                                         fill_crse_input=.false.)

       enddo

    end if

  end subroutine initveldata




  subroutine initveldata_2d(u,lo,hi,ng,dx,s0_init,p0_init,xloc_vortices)

    use probin_module, only: prob_lo, apply_vel_field, velpert_scale, &
                             velpert_amplitude, velpert_height_loc

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: xloc_vortices(:)

    ! local variables
    real(kind=dp_t) :: xloc(2), upert(2)
    integer :: i, j, vortex
    real(kind=dp_t) :: xdist, ydist, r


    u = ZERO

    if (apply_vel_field) then

       do j = lo(2), hi(2)

          xloc(2) = prob_lo(2) + (dble(j)+HALF)*dx(2)

          ydist = xloc(2) - velpert_height_loc
          
          do i = lo(1), hi(1)

             upert = ZERO

             xloc(1) = prob_lo(1) + (dble(i)+HALF)*dx(1)

             ! loop over each vortex
             do vortex = 1, size(xloc_vortices, dim=1)

                xdist = xloc(1) - xloc_vortices(vortex)

                r = xdist**2 + ydist**2
                r = sqrt(r)

                ! e.g. Calder et al. ApJSS 143, 201-229 (2002)
                ! we set things up so that every other vortex has the same
                ! orientation
                upert(1) = upert(1) - (ydist/velpert_scale) * &
                     velpert_amplitude * exp( -r**2/(TWO*velpert_scale**2)) &
                     * (-ONE)**vortex

                upert(2) = upert(2) + (xdist/velpert_scale) * &
                     velpert_amplitude * exp(-r**2/(TWO*velpert_scale**2)) &
                     * (-ONE)**vortex

             enddo

             u(i,j,:) = u(i,j,:) + upert(:)

          enddo

       enddo
       
                
    endif

  end subroutine initveldata_2d




  subroutine initveldata_3d(u,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: apply_vel_field

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! initial the velocity
    u = ZERO

    ! we don't support the vortices in 3d yet...
    if (apply_vel_field) call bl_error("apply_vel_field not supported for 3d")

    
  end subroutine initveldata_3d

end module init_vel_module
