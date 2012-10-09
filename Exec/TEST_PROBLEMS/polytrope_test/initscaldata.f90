module init_scalar_module

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
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initscalardata, initscalardata_on_level

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla)

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n,dm,nlevs
    
    dm = mla%dim
    nlevs = mla%nlevel

    ng = s(1)%ng

    do n=1,nlevs
       do i = 1, nfabs(s(n))
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
                                    p0_init(n,:))
          case (3)
             if (spherical .eq. 1) then
                call bl_error("ERROR: spherical not implemented in initdata")
             else
                call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
                                       p0_init(n,:))
             end if
          end select
       end do
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(s(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),rho_comp,dm+rho_comp,nscal,bc(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),rho_comp,dm+rho_comp,nscal, &
                                         fill_crse_input=.false.)

       enddo

    end if

  end subroutine initscalardata

  subroutine initscalardata_on_level(n,s,s0_init,p0_init,dx,bc)

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    ! local
    integer                  :: ng,i,dm
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    real(kind=dp_t), pointer :: sop(:,:,:,:)

    ng = nghost(s)
    dm = get_dim(s)

    do i = 1, nfabs(s)
       sop => dataptr(s,i)
       lo =  lwb(get_box(s,i))
       hi =  upb(get_box(s,i))
       select case (dm)
       case (2)
          call initscalardata_2d(sop(:,:,1,:),lo,hi,ng,dx,s0_init,p0_init)
       case (3)
          if (spherical .eq. 1) then
             call bl_error("ERROR: spherical not implemented in initdata")
          else
             call bl_error("ERROR: 3d not implemented in initdata")
          end if
       end select
    end do

    call multifab_fill_boundary(s)

    call multifab_physbc(s,rho_comp,dm+rho_comp,nscal,bc)

  end subroutine initscalardata_on_level

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, prob_hi, perturb_model, &
                             perturb_type, yhigh, ylow, amp, kmode
    use geometry, only: nr_fine

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables
    integer         :: i,j,comp
    real(kind=dp_t) :: x,y
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)
    real(kind=dp_t) :: ymid, ywidth, Lx, Ly, fy, ys

    real(kind=dp_t) :: rhoX_1,rhoX_2

    ! initial the domain with the base state
    s = ZERO

    if (perturb_type .eq. 1) then

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             s(i,j,rho_comp) = s0_init(j,rho_comp)
             s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
             s(i,j,temp_comp) = s0_init(j,temp_comp)
             s(i,j,spec_comp:spec_comp+nspec-1) = s0_init(j,spec_comp:spec_comp+nspec-1)
             s(i,j,trac_comp:trac_comp+ntrac-1) = s0_init(j,trac_comp:trac_comp+ntrac-1)

          end do
       end do

       Lx = prob_hi(1) - prob_lo(1)

       do j = lo(2), hi(2)
          y = (j+HALF)*dx(2)+prob_lo(2)
          do i = lo(1), hi(1)
             x = (i+HALF)*dx(1)+prob_lo(1)

             ymid = HALF * (yhigh - ylow)
             ywidth = yhigh - ymid

             s(i,j,rho_comp) = s(i,j,rho_comp)*(1.d0+amp*sin(2.d0*M_PI*x/Lx)*exp(-(y-ymid)*(y-ymid)/ywidth/ywidth))

          end do
       end do

    else if (perturb_type .eq. 2) then

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             s(i,j,rho_comp)  = s0_init(j,rho_comp)
             s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
             s(i,j,temp_comp) = s0_init(j,temp_comp)
             s(i,j,spec_comp:spec_comp+nspec-1) = s0_init(j,spec_comp:spec_comp+nspec-1)
             s(i,j,trac_comp:trac_comp+ntrac-1) = s0_init(j,trac_comp:trac_comp+ntrac-1)
          enddo
       enddo

    else if (perturb_type .eq. 3) then

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             s(i,j,rho_comp) = s0_init(j,rho_comp)
             s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
             s(i,j,temp_comp) = s0_init(j,temp_comp)
             s(i,j,spec_comp:spec_comp+nspec-1) = s0_init(j,spec_comp:spec_comp+nspec-1)
             s(i,j,trac_comp:trac_comp+ntrac-1) = s0_init(j,trac_comp:trac_comp+ntrac-1)

          end do
       end do

       Lx = prob_hi(1) - prob_lo(1)
       Ly = prob_hi(2) - prob_lo(2)

       do j = lo(2), hi(2)
          y = (j+HALF)*dx(2)+prob_lo(2)
          do i = lo(1), hi(1)
             x = (i+HALF)*dx(1)+prob_lo(1)

             s(i,j,rho_comp) = s(i,j,rho_comp)*(1.d0+amp*sin(2.d0*M_PI*x/Lx*kmode)*sin(M_PI*y/Ly))

          end do
       end do

    else if (perturb_type .eq. 4) then

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             s(i,j,rho_comp) = s0_init(j,rho_comp)
             s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
             s(i,j,temp_comp) = s0_init(j,temp_comp)
             s(i,j,spec_comp:spec_comp+nspec-1) = s0_init(j,spec_comp:spec_comp+nspec-1)
             s(i,j,trac_comp:trac_comp+ntrac-1) = s0_init(j,trac_comp:trac_comp+ntrac-1)

          end do
       end do

       Lx = prob_hi(1) - prob_lo(1)
       Ly = prob_hi(2) - prob_lo(2)

       do j = lo(2), hi(2)
          y = (j+HALF)*dx(2)+prob_lo(2)
          do i = lo(1), hi(1)
             x = (i+HALF)*dx(1)+prob_lo(1)

             ys = (2.0d0*y - prob_hi(2) - prob_lo(2))/Ly

             fy = 1.0d0 - 3.0d0*ys**2 + 3.0d0*ys**4 - ys**6 + 2.5d0*ys - 2.5d0*ys**3

             s(i,j,rho_comp) = s(i,j,rho_comp)*(1.d0+amp*sin(2.d0*M_PI*x/Lx)*fy)

          end do
       end do

    end if

  end subroutine initscalardata_2d

  subroutine initscalardata_3d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, prob_hi, perturb_type, yhigh, ylow &
                                , amp, perturb_model, n_cellz, n_celly
    use mt19937_module    

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    integer :: i, j, k, iseed
    real (kind = dp_t) :: x, y, z
    real (kind = dp_t) :: pertheightarray(lo(1):hi(1),lo(2):hi(2))

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             s(i,j,k,:) = s0_init(k,:)
          end do
       end do
    end do

    if ( perturb_type .eq. 1 ) then

       do k = lo(3), hi(3)
          z = prob_lo(3) + (k + HALF)*dx(3)
          do j = lo(2), hi(2)
             y = prob_lo(2) + (j + HALF)*dx(2)
             do i = lo(1), hi(1)
                x = prob_lo(1) + (i + HALF)*dx(1)

                if ( (z .gt. ylow) .and. (z .lt. yhigh) ) then
                   s(i,j,k,rho_comp) = s(i,j,k,rho_comp)*amp
                end if

             end do
          end do
       end do

    end if


  end subroutine initscalardata_3d

end module init_scalar_module
