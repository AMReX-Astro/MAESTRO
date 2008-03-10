! If flag = .true. then convert the state array species from (rho X)
! to X.  If flag = .false. then convert X to (rho X).  Note, this only
! applies when we are coming in with the full state (not the
! perturbational state).  This does not touch the base state.

module convert_rhoX_to_X_module

  use multifab_module

  implicit none

  private

  public :: convert_rhoX_to_X, convert_rhoh_to_h
  
contains

  subroutine convert_rhoX_to_X(nlevs,s,dx,flag,mla,the_bc_level)

    use network, only: nspec
    use variables, only: spec_comp, foextrap_comp, rho_comp
    use ml_layout_module
    use define_bc_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module
    use multifab_physbc_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    logical        , intent(in   ) :: flag
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    integer :: lo(s(1)%dim),hi(s(1)%dim)
    integer :: dm,n,comp,bc_comp

    dm = s(1)%dim

    if (flag) then
       do n=1,nlevs
          do comp=spec_comp,spec_comp+nspec-1
             call multifab_div_div_c(s(n),comp,s(n),rho_comp,1)
          end do
       end do
    else
       do n=1,nlevs
          do comp=spec_comp,spec_comp+nspec-1
             call multifab_mult_mult_c(s(n),comp,s(n),rho_comp,1)
          end do
       end do
    end if

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s(nlevs),spec_comp,nspec)

       do comp = spec_comp,spec_comp+nspec-1

          if (flag) then
             bc_comp = foextrap_comp
          else
             bc_comp = dm+comp
          end if

          ! fill non-periodic domain boundary ghost cells
          call multifab_physbc(s(nlevs),comp,bc_comp,1,the_bc_level(nlevs))
       end do

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
          
          do comp = spec_comp,spec_comp+nspec-1
             if (flag) then
                bc_comp = foextrap_comp
             else
                bc_comp = dm+comp
             end if
             
             ! fill level n ghost cells using interpolation from level n-1 data
             ! note that multifab_fill_boundary and multifab_physbc are called for
             ! both levels n-1 and n
             call multifab_fill_ghost_cells(s(n),s(n-1), &
                                            s(n)%ng,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1),the_bc_level(n), &
                                            comp,bc_comp,1)
          end do

       end do

    end if
    
  end subroutine convert_rhoX_to_X

  subroutine convert_rhoh_to_h(nlevs,s,dx,flag,mla,the_bc_level)

    use variables, only: rho_comp, rhoh_comp, foextrap_comp
    use ml_layout_module
    use define_bc_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module
    use multifab_physbc_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    logical        , intent(in   ) :: flag
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    integer :: lo(s(1)%dim),hi(s(1)%dim)
    integer :: dm,n,bc_comp

    dm = s(1)%dim

    if (flag) then
       do n=1,nlevs
          call multifab_div_div_c(s(n),rhoh_comp,s(n),rho_comp,1)
       end do
    else
       do n=1,nlevs
          call multifab_mult_mult_c(s(n),rhoh_comp,s(n),rho_comp,1)
       end do
    end if
    
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s(nlevs),rhoh_comp,1)

       if (flag) then
          bc_comp = foextrap_comp
       else
          bc_comp = dm+rhoh_comp
       end if

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),rhoh_comp,bc_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
          
          if (flag) then
             bc_comp = foextrap_comp
          else
             bc_comp = dm+rhoh_comp
          end if
             
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1), &
                                         s(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         rhoh_comp,bc_comp,1)
       end do

    end if
    
  end subroutine convert_rhoh_to_h


  subroutine make_edge_spec_1d(s,sedgex,umac,force,dx,dt)

     use geometry, only: nr
     use probin_module, only: slope_order
     use bl_constants_module
     
     real(kind=dp_t), intent(in   ) ::      s(:)
     real(kind=dp_t), intent(inout) :: sedgex(:)
     real(kind=dp_t), intent(in   ) ::   umac(:)
     real(kind=dp_t), intent(in   ) ::  force(:)
     real(kind=dp_t), intent(in   ) :: dx,dt
     
     real(kind=dp_t), allocatable::  slopex(:)
     real(kind=dp_t), allocatable::  s_l(:),s_r(:)
     real(kind=dp_t), allocatable:: dxscr(:,:)
     real(kind=dp_t) :: dmin,dpls,ds,del,slim,sflag
     real(kind=dp_t) :: ubardth, dth, savg
     real(kind=dp_t) :: abs_eps, eps, umax, u
     
     integer :: i,is,ie,hi,lo
     integer        , parameter :: cen = 1, lim = 2, flag = 3, fromm = 4
     real(kind=dp_t), parameter :: fourthirds = 4.0_dp_t / 3.0_dp_t
     
     lo = 1
     hi = lo + size(force,dim=1) - 1
     
     allocate(s_l(lo-1:hi+2),s_r(lo-1:hi+2))
     allocate(slopex(lo:hi))
     allocate(dxscr(lo:hi,4))
     
     abs_eps = 1.0d-8
     
     dth = HALF*dt
     
     is = lo
     ie = hi
     
     umax = ZERO
     do i = is,ie+1
        umax = max(umax,abs(umac(i)))
     end do
     
     eps = abs_eps * umax

     if (slope_order .eq. 0) then

        slopex = ZERO

     else if (slope_order .eq. 2) then

        do i = is+1,ie-1
           del = half*(s(i+1) - s(i-1))
           dpls = two*(s(i+1) - s(i  ))
           dmin = two*(s(i  ) - s(i-1))
           slim = min(abs(dpls), abs(dmin))
           slim = merge(slim, zero, dpls*dmin.gt.ZERO)
           sflag = sign(one,del)
           slopex(i)= sflag*min(slim,abs(del))
        enddo
     
        slopex(is) = ZERO
        slopex(ie) = ZERO

     else if (slope_order .eq. 4) then
     
        do i = is+1,ie-1
           dxscr(i,cen) = half*(s(i+1)-s(i-1))
           dpls = two*(s(i+1)-s(i  ))
           dmin = two*(s(i  )-s(i-1))
           dxscr(i,lim)= min(abs(dmin),abs(dpls))
           dxscr(i,lim) = merge(dxscr(i,lim),zero,dpls*dmin.gt.ZERO)
           dxscr(i,flag) = sign(one,dxscr(i,cen))
           dxscr(i,fromm)= dxscr(i,flag)*min(dxscr(i,lim),abs(dxscr(i,cen)))
        enddo
     
        dxscr(is,fromm) = ZERO
        dxscr(ie,fromm) = ZERO
     
        do i = is+1,ie-1
           ds = fourthirds * dxscr(i,cen) - sixth * (dxscr(i+1,fromm) + dxscr(i-1,fromm))
           slopex(i) = dxscr(i,flag)*min(abs(ds),dxscr(i,lim))
        enddo
     
        slopex(is) = ZERO
        slopex(ie) = ZERO

     end if
        
     ! Compute edge values using slopes and forcing terms.
     do i = is,ie
        
        u = HALF * (umac(i) + umac(i+1))
        ubardth = dth*u/dx
        
        s_l(i+1)= s(i) + (HALF-ubardth)*slopex(i) + dth * force(i)
        s_r(i  )= s(i) - (HALF+ubardth)*slopex(i) + dth * force(i)
        
     enddo
     
     sedgex(is  ) = s_r(is  )
     sedgex(ie+1) = s_l(ie+1)
     
     do i = is+1, ie 
        sedgex(i)=merge(s_l(i),s_r(i),umac(i).gt.ZERO)
        savg = HALF*(s_r(i) + s_l(i))
        sedgex(i)=merge(savg,sedgex(i),abs(umac(i)) .lt. eps)
     enddo
     
  end subroutine make_edge_spec_1d
  
end module convert_rhoX_to_X_module
