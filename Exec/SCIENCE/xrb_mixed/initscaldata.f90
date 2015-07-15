module init_scalar_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use variables
  use network
  use geometry, only: nr, spherical
  use ml_layout_module
  use ml_cc_restriction_module
  use multifab_fill_ghost_module

  implicit none

  real(dp_t), save :: pert_height

  private
  public :: initscalardata, initscalardata_on_level

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla)

    use probin_module, only: prob_lo, prob_hi, perturb_model, xrb_pert_height

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng,dm,nlevs
    integer :: i,n,r

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

             ! make sure our perturbation is within the domain
             if (perturb_model .and. &
                  (xrb_pert_height < prob_lo(2) .or. &
                   xrb_pert_height > prob_hi(2))) &
                   call bl_error("perturbation outside of the domain")

             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), &
                  s0_init(n,:,:), p0_init(n,:))
          case (3)

             ! make sure our perturbation is within the domain
             if (perturb_model .and. &
                  (xrb_pert_height < prob_lo(3) .or. &
                   xrb_pert_height > prob_hi(3))) &
                   call bl_error("perturbation outside of the domain")

             call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx(n,:), &
                  s0_init(n,:,:), p0_init(n,:))
          end select

       enddo
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(s(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),rho_comp,dm+rho_comp,nscal,bc(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids
       ! are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering 
          ! it
          call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),rho_comp,dm+rho_comp,nscal)

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
    integer                    :: ng,i,r
    integer                    :: lo(get_dim(s)),hi(get_dim(s)),dm
    real(kind=dp_t), pointer   :: sop(:,:,:,:)

    dm = get_dim(s)

    ng = nghost(s)

    do i = 1, nfabs(s)

       sop => dataptr(s,i)
       lo =  lwb(get_box(s,i))
       hi =  upb(get_box(s,i))
       select case (dm)
       case (2)
          call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx, s0_init, &
               p0_init)
       case (3)
          call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx, s0_init, &
               p0_init)
       end select

    end do

    call multifab_fill_boundary(s)

    call multifab_physbc(s,rho_comp,dm+rho_comp,nscal,bc)

  end subroutine initscalardata_on_level



  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, prob_hi, perturb_model, xrb_pert_height

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables
    integer         :: i,j
    real(kind=dp_t) :: x,y
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    real(kind=dp_t) :: xcen, ycen, dist

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          s(i,j,rho_comp)  = s0_init(j,rho_comp)
          s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
          s(i,j,temp_comp) = s0_init(j,temp_comp)
          s(i,j,spec_comp:spec_comp+nspec-1) = &
               s0_init(j,spec_comp:spec_comp+nspec-1)
          s(i,j,trac_comp:trac_comp+ntrac-1) = &
               s0_init(j,trac_comp:trac_comp+ntrac-1)

       enddo
    enddo
    
    ! add an optional perturbation
    if (perturb_model) then

       xcen = (prob_lo(1) + prob_hi(1)) / TWO
       ycen = xrb_pert_height

       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF) * dx(2)
          
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i)+HALF) * dx(1)
          
             dist = sqrt((x-xcen)**2 + (y-ycen)**2)

             call perturb(dist, p0_init(j), s0_init(j,:), &
                          dens_pert, rhoh_pert, rhoX_pert, temp_pert, &
                          trac_pert)

             s(i,j,rho_comp) = dens_pert
             s(i,j,rhoh_comp) = rhoh_pert
             s(i,j,temp_comp) = temp_pert
             s(i,j,spec_comp:spec_comp+nspec-1) = rhoX_pert(1:)
             s(i,j,trac_comp:trac_comp+ntrac-1) = trac_pert(:)

          enddo
       enddo
    endif
    
  end subroutine initscalardata_2d




  subroutine initscalardata_3d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, prob_hi, perturb_model, xrb_pert_height
    
    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! local
    integer         :: i,j,k
    real(kind=dp_t) :: x, y, z
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    real(kind=dp_t)            :: xcen, ycen, zcen, dist

    ! initial the domain with the base state
    s = ZERO
  
    if (spherical .eq. 1) then

       call bl_error('Error: initdata does not handle the spherical case')

    else 

       ! initialize the scalars
       do k = lo(3), hi(3)

          do j = lo(2), hi(2)

             do i = lo(1), hi(1)
                
                s(i,j,k,rho_comp)  = s0_init(k,rho_comp)
                s(i,j,k,rhoh_comp) = s0_init(k,rhoh_comp)
                s(i,j,k,temp_comp) = s0_init(k,temp_comp)
                s(i,j,k,spec_comp:spec_comp+nspec-1) = &
                     s0_init(k,spec_comp:spec_comp+nspec-1)
                s(i,j,k,trac_comp:trac_comp+ntrac-1) = &
                     s0_init(k,trac_comp:trac_comp+ntrac-1)

             enddo

          enddo

       enddo
       
       if (perturb_model) then

          xcen = (prob_lo(1) + prob_hi(1)) / TWO
          ycen = (prob_lo(2) + prob_hi(2)) / TWO
          zcen = xrb_pert_height

          ! add an optional perturbation
          do k = lo(3), hi(3)
             z = prob_lo(3) + (dble(k)+HALF) * dx(3)
             
             do j = lo(2), hi(2)
                y = prob_lo(2) + (dble(j)+HALF) * dx(2)
                
                do i = lo(1), hi(1)
                   x = prob_lo(1) + (dble(i)+HALF) * dx(1)

                   dist = sqrt((x-xcen)**2 + (y-ycen)**2 + (z-zcen)**2)
                   
                   call perturb(dist, p0_init(k), s0_init(k,:), &
                                dens_pert, rhoh_pert, rhoX_pert, temp_pert, &
                                trac_pert)

                   s(i,j,k,rho_comp) = dens_pert
                   s(i,j,k,rhoh_comp) = rhoh_pert
                   s(i,j,k,temp_comp) = temp_pert
                   s(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(:)
                   s(i,j,k,trac_comp:trac_comp+ntrac-1) = trac_pert(:)

                enddo

             enddo

          enddo

       endif

    end if
    
  end subroutine initscalardata_3d

  subroutine perturb(distance, p0_init, s0_init,  &
                     dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    use probin_module, ONLY: xrb_pert_size, xrb_pert_factor, xrb_pert_type
    use eos_module, only: eos_input_tp, eos_input_rp, eos
    use eos_type_module

    real(kind=dp_t), intent(in ) :: distance
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp,dens,t0,d0,rad_pert

    integer, parameter :: perturb_temp = 1, perturb_dens = 2
    integer :: eos_input_flag

    type (eos_t) :: eos_state

    rad_pert = -xrb_pert_size**2 / (FOUR*log(HALF))

    select case (xrb_pert_type)
    case(perturb_temp)

       t0 = s0_init(temp_comp)

       temp = t0 * (ONE + xrb_pert_factor * dexp(-distance**2 / rad_pert) )

       dens = s0_init(rho_comp)

       eos_input_flag = eos_input_tp

    case(perturb_dens)
          
       d0 = s0_init(rho_comp)
       
       dens = d0 * (ONE + xrb_pert_factor * dexp(-distance**2 / rad_pert) )
       
       temp = s0_init(temp_comp)

       eos_input_flag = eos_input_rp

    end select

    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = dens
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_flag, eos_state, .false.)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho * eos_state%h
    rhoX_pert(:) = dens_pert*eos_state%xn(:)

    temp_pert = eos_state%T
    
    trac_pert(:) = ZERO

  end subroutine perturb

end module init_scalar_module
