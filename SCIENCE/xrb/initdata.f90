module init_module

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
  public :: initscalardata, initveldata, scalar_diags

contains

  subroutine initscalardata(nlevs,s,s0,p0,dx,bc,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(s(1)%dim),hi(s(1)%dim),ng,dm
    integer :: i,n
    
    ng = s(1)%ng
    dm = s(1)%dim

    do n=1,nlevs

       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s0(n,:,:), p0(n,:))
          case (3)
             call initscalardata_3d(n,sop(:,:,:,:), lo, hi, ng, dx(n,:), s0(n,:,:), p0(n,:))
          end select
       end do
       
       call multifab_fill_boundary(s(n))
       call multifab_physbc(s(n),rho_comp,dm+rho_comp,nscal,dx(n,:),bc(n))

    enddo

    do n=nlevs,2,-1
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                       bc(n-1),bc(n),1,dm+rho_comp,nscal)
    enddo

  end subroutine initscalardata

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0,p0)

    use probin_module, only: prob_lo_x, prob_hi_x,    &
                             prob_lo_y,               &
                             perturb_model

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0(0:)

    ! Local variables
    integer         :: i,j,n
    real(kind=dp_t) :: x,y
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    real(kind=dp_t) :: xcen, ycen, dist
    real(kind=dp_t), parameter :: he4_pert = 1.e-2_dp_t
    integer         :: he4_comp, pert_index

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          s(i,j,rho_comp)  = s0(j,rho_comp)
          s(i,j,rhoh_comp) = s0(j,rhoh_comp)
          s(i,j,temp_comp) = s0(j,temp_comp)
          s(i,j,spec_comp:spec_comp+nspec-1) = &
               s0(j,spec_comp:spec_comp+nspec-1)
          s(i,j,trac_comp:trac_comp+ntrac-1) = &
               s0(j,trac_comp:trac_comp+ntrac-1)
       enddo
    enddo
    
    ! add an optional perturbation
    if (perturb_model) then

       xcen = (prob_lo_x + prob_hi_x) / TWO

       he4_comp = network_species_index('helium-4')

       do j = lo(2), hi(2)
          if (s0(j,spec_comp+he4_comp-1)/s0(j,rho_comp) .gt. he4_pert) then
             pert_index = j
             exit
          endif
       enddo

       ycen = prob_lo_y + (dble(pert_index)+HALF) * dx(2)

       do j = lo(2), hi(2)
          y = prob_lo_y + (dble(j)+HALF) * dx(2)
          
          do i = lo(1), hi(1)
             x = prob_lo_x + (dble(i)+HALF) * dx(1)
          
             dist = sqrt((x-xcen)**2 + (y-ycen)**2)

             call perturb(dist, p0(j), s0(j,:), &
                          dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

             s(i,j,rho_comp) = dens_pert
             s(i,j,rhoh_comp) = rhoh_pert
             s(i,j,temp_comp) = temp_pert
             s(i,j,spec_comp:spec_comp+nspec-1) = rhoX_pert(1:)
             s(i,j,trac_comp:trac_comp+ntrac-1) = trac_pert(:)
          enddo
       enddo
    endif
    
  end subroutine initscalardata_2d

  subroutine initscalardata_3d(n,s,lo,hi,ng,dx,s0,p0)

    use probin_module, only: prob_lo_x, prob_hi_x,    &
                             prob_lo_y, prob_hi_y,    &
                             prob_lo_z,               &
                             perturb_model
    
    integer           , intent(in   ) :: n,lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0(0:)

    !     Local variables
    integer         :: i,j,k,comp
    real(kind=dp_t) :: x, y, z
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    real(kind=dp_t)            :: xcen, ycen, zcen, dist
    real(kind=dp_t), parameter :: he4_pert = 1.e-2
    integer                    :: he4_comp, pert_index

    ! initial the domain with the base state
    s = ZERO
  
    if (spherical .eq. 1) then

       ! initialize the scalars
       call fill_3d_data(n,s(:,:,:,rho_comp), s0(:,rho_comp), lo,hi,dx,ng)
       call fill_3d_data(n,s(:,:,:,rhoh_comp),s0(:,rhoh_comp),lo,hi,dx,ng)
       call fill_3d_data(n,s(:,:,:,temp_comp),s0(:,temp_comp),lo,hi,dx,ng)

       do comp = spec_comp, spec_comp+nspec-1
          call fill_3d_data(n,s(:,:,:,comp),s0(:,comp),lo,hi,dx,ng)
       end do

       do comp = trac_comp, trac_comp+ntrac-1
          call fill_3d_data(n,s(:,:,:,comp),s0(:,comp),lo,hi,dx,ng)
       end do

    else 

       ! initialize the scalars
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                s(i,j,k,rho_comp)  = s0(k,rho_comp)
                s(i,j,k,rhoh_comp) = s0(k,rhoh_comp)
                s(i,j,k,temp_comp) = s0(k,temp_comp)
                s(i,j,k,spec_comp:spec_comp+nspec-1) = &
                     s0(k,spec_comp:spec_comp+nspec-1)
                s(i,j,k,trac_comp:trac_comp+ntrac-1) = &
                     s0(k,trac_comp:trac_comp+ntrac-1)
             enddo
          enddo
       enddo
       
       if (perturb_model) then

          xcen = (prob_lo_x + prob_hi_x) / TWO
          ycen = (prob_lo_y + prob_hi_y) / TWO

          he4_comp = network_species_index('helium-4')
          do k = lo(3), hi(3)
             if (s0(k,spec_comp+he4_comp-1)/s0(j,rho_comp) .gt. he4_pert) then
                pert_index = j
                exit
             endif
          enddo

          zcen = prob_lo_z + (dble(pert_index)+HALF)*dx(3)

          ! add an optional perturbation
          do k = lo(3), hi(3)
             z = prob_lo_z + (dble(k)+HALF) * dx(3)
             
             do j = lo(2), hi(2)
                y = prob_lo_y + (dble(j)+HALF) * dx(2)
                
                do i = lo(1), hi(1)
                   x = prob_lo_x + (dble(i)+HALF) * dx(1)

                   dist = sqrt((x-xcen)**2 + (y-ycen)**2 + (z-zcen)**2)
                   
                   call perturb(dist, p0(k), s0(k,:), &
                                dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

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

  subroutine initveldata(nlevs,u,s0,p0,dx,bc,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(u(1)%dim),hi(u(1)%dim),ng,dm
    integer :: i,n
    
    ng = u(1)%ng
    dm = u(1)%dim

    do n=1,nlevs

       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), &
                                 s0(n,:,:), p0(n,:))
          case (3) 
             call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                 s0(n,:,:), p0(n,:))
          end select
       end do

       call multifab_fill_boundary(u(n))
       call multifab_physbc(u(n),1,1,dm,dx(n,:),bc(n))
       
    enddo
    
    do n=nlevs,2,-1
       call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))
       call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                      bc(n-1),bc(n),1,1,dm)
    enddo

  end subroutine initveldata

  subroutine initveldata_2d(u,lo,hi,ng,dx,s0,p0)

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0(0:)

    ! Local variables

    ! initial the velocity
    u = ZERO

  end subroutine initveldata_2d

  subroutine initveldata_3d(u,lo,hi,ng,dx,s0,p0)

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0(0:)

    ! Local variables

    ! initial the velocity
    u = ZERO
    
  end subroutine initveldata_3d


  subroutine perturb(distance, p0, s0,  &
                     dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

    use probin_module, ONLY: pert_size, pert_factor, pert_type
    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: distance
    real(kind=dp_t), intent(in ) :: p0, s0(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp,dens,t0,d0,rad_pert

    integer, parameter :: perturb_temp = 1, perturb_dens = 2
    integer :: eos_input_flag

    rad_pert = pert_size**2 / (FOUR * log(TWO))

    select case (pert_type)
    case(perturb_temp)

       t0 = s0(temp_comp)

       temp = t0 * (ONE + pert_factor) * dexp(-distance**2 / rad_pert)
       temp = max(t0, temp)
       
       dens = s0(rho_comp)

       eos_input_flag = eos_input_tp

    case(perturb_dens)
          
       d0 = s0(rho_comp)
       
       dens = d0 * (ONE + pert_factor) * dexp(-distance**2 / rad_pert)
       dens = max(d0, dens)
       
       temp = s0(temp_comp)

       eos_input_flag = eos_input_rp

    end select

    temp_eos(1) = temp
    p_eos(1) = p0
    den_eos(1) = dens
    xn_eos(1,:) = s0(spec_comp:spec_comp+nspec-1)/s0(rho_comp)

    call eos(eos_input_flag, den_eos, temp_eos, &
             npts, nspec, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

    dens_pert = den_eos(1)
    rhoh_pert = den_eos(1)*h_eos(1)
    rhoX_pert(:) = dens_pert*xn_eos(1,:)

    temp_pert = temp_eos(1)
    
    trac_pert(:) = ZERO

  end subroutine perturb


  subroutine scalar_diags (istep,s,s0,p0,dx)

    integer        , intent(in   ) :: istep
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in)    :: s0(:,:)
    real(kind=dp_t), intent(in)    :: p0(:)
    real(kind=dp_t), intent(in)    :: dx(:)

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i,n
    
    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sop => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))

       select case (dm)
       case (2)
          call scalar_diags_2d(istep, sop(:,:,1,:), lo, hi, ng, dx, s0, p0)
       case (3)
!         call scalar_diags_3d(istep, sop(:,:,:,:), lo, hi, ng, dx, s0)
       end select
    end do

  end subroutine scalar_diags

  subroutine scalar_diags_2d (istep, s,lo,hi,ng,dx,s0,p0)

    use probin_module, only: grav_const

    integer, intent(in) :: istep, lo(:), hi(:), ng
    real (kind = dp_t), intent(in) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in) :: dx(:)
    real(kind=dp_t)   , intent(in) :: s0(0:,:)
    real(kind=dp_t)   , intent(in) :: p0(0:)

    ! Local variables
    integer :: i, j, n
    real(kind=dp_t) :: fac, mass, mass0
    real(kind=dp_t), allocatable ::  rhoavg(:)
    real(kind=dp_t), allocatable :: rhopert(:)
    real(kind=dp_t), allocatable ::    pavg(:)
    character(len=11) :: file_name
    character(len=10) :: file_name2
    character(len= 8) :: file_name3

    allocate(rhopert(lo(2):hi(2)))
    allocate( rhoavg(lo(2):hi(2)))
    allocate(   pavg(lo(2):hi(2)))

    write(unit=file_name ,fmt='("rhopert",i4.4)') istep
    write(unit=file_name2,fmt='("rhoavg",i4.4)') istep
    write(unit=file_name3,fmt='("pavg",i4.4)') istep
    open(90,file=file_name)
    open(91,file=file_name2)
    open(92,file=file_name3)

    fac = ONE / dble(hi(1)-lo(1)+1)
    mass  = ZERO
    mass0 = ZERO
    do j = lo(2), hi(2)
      rhoavg(j) = ZERO
      rhopert(j) = ZERO
      do i = lo(1), hi(1)
         rhopert(j) = rhopert(j) + (s(i,j,rho_comp) - s0(j,rho_comp))
         rhoavg(j) = rhoavg(j) +  s(i,j,rho_comp)
      enddo
      rhoavg(j)  = rhoavg(j) * fac
      rhopert(j)  = rhopert(j) * fac
      write(90,*) (dble(j)+HALF)*dx(2),rhopert(j)
      write(91,*) (dble(j)+HALF)*dx(2),rhoavg(j)
      mass  = mass  + rhoavg(j)
      mass0 = mass0 + s0(j,rho_comp)
    enddo

!   print *,'TOTAL MASS ',istep, mass, mass0

    pavg(hi(2)) = p0(hi(2))
    do j = hi(2)-1,lo(2),-1
      pavg(j) = pavg(j+1) + 0.5d0 * (rhoavg(j+1)+rhoavg(j))*abs(grav_const)*dx(2)
    enddo
    do j = lo(2),hi(2)
      write(92,*) (dble(j)+HALF)*dx(2),p0(j),pavg(j)
    enddo

    deallocate(rhoavg,rhopert,pavg)

  end subroutine scalar_diags_2d


end module init_module
