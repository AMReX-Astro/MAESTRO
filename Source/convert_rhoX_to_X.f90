! If flag = .true. then convert the state array species from (rho X)
! to X.  If flag = .false. then convert X to (rho X).  Note, this only
! applies when we are coming in with the full state (not the
! perturbational state).  This does not touch the base state.

module convert_rhoX_to_X_module

  use multifab_module

  implicit none

  private

  public :: convert_rhoX_to_X, make_edge_rhoX_from_X
  
contains

  subroutine convert_rhoX_to_X(nlevs,s,dx,flag,mla,the_bc_level)

    use geometry, only: spherical
    use network, only: nspec
    use variables, only: spec_comp, foextrap_comp, nscal
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
    integer :: i,ng,dm,n,comp,bc_comp

    ng = s(1)%ng
    dm = s(1)%dim

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sp => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call convert_rhoX_to_X_2d(n,sp(:,:,1,:),lo,hi,ng,flag)
          case (3)
             call convert_rhoX_to_X_3d(n,sp(:,:,:,:),lo,hi,ng,flag)
          end select
       end do
    end do

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

  subroutine convert_rhoX_to_X_2d(n,s,lo,hi,ng,flag)

    use network, only: nspec
    use variables, only: spec_comp, rho_comp
    use geometry, only: nr

    integer        , intent(in   ) :: n
    integer        , intent(in   ) ::  lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,r,comp

    if (flag) then

       ! convert (rho X) -> X
       do comp = spec_comp, spec_comp+nspec-1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                s(i,j,comp) = s(i,j,comp)/s(i,j,rho_comp)
             end do
          end do
       end do

    else

       ! convert X -> (rho X)
       do comp = spec_comp, spec_comp+nspec-1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                s(i,j,comp) = s(i,j,rho_comp)*s(i,j,comp)
             end do
          end do
       end do

    end if

  end subroutine convert_rhoX_to_X_2d

  subroutine convert_rhoX_to_X_3d(n,s,lo,hi,ng,flag)

    use network, only: nspec
    use variables, only: spec_comp, rho_comp
    use geometry, only: nr

    integer        , intent(in   ) :: n
    integer        , intent(in   ) ::  lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,k,r,comp

    if (flag) then

       ! convert (rho X) -> X
       do comp = spec_comp, spec_comp+nspec-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp)/s(i,j,k,rho_comp)
                end do
             end do
          end do
       end do

    else

       ! convert X -> (rho X)
       do comp = spec_comp, spec_comp+nspec-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,rho_comp)*s(i,j,k,comp)
                end do
             end do
          end do
       end do

    end if

  end subroutine convert_rhoX_to_X_3d



  subroutine make_edge_rhoX_from_X(nlevs,u,sedge, &
                                   s0_predicted_edge, &
                                   the_bc_level,dx)


    ! here we take the edge states of X' and rho', and together with the base
    ! edges states (s0_predicted_edge, as returned from advect_base) rho0 and
    ! X_0, we construct the edge state for (rho X)'.  Note, the base state 
    ! quantity X_0 is really the favre average, (rho X)_0/rho_0.
    !
    ! here, (rho X)' = (rho_0 + rho')X' + (rho' X_0) on edges.

    use bl_prof_module
    use bl_constants_module
    use geometry
    use variables
    use probin_module, only: predict_X_at_edges
    use network, only: nspec
    use fill_3d_module, only: fill_3d_data
    use define_bc_module
    use multifab_physbc_module
    
    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    real(kind=dp_t), intent(in   ) :: s0_predicted_edge(:,0:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    
    ! local
    integer :: i,dm,n
    integer :: lo(u(1)%dim),hi(u(1)%dim)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_edge_rhoX_from_X")

    dm = u(1)%dim


    do n=1,nlevs

       do i=1,u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          sepx => dataptr(sedge(n,1), i)
          sepy => dataptr(sedge(n,2), i)
          lo = lwb(get_box(u(n),i))
          hi = upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call make_edge_rhoX_from_X_2d(sepx(:,:,1,:), sepy(:,:,1,:), &
                                           s0_predicted_edge(n,:,:), &
                                           lo, hi)

          case (3)
             sepz => dataptr(sedge(n,3),i)
             if (spherical .eq. 1) then

               ! for spherical, we need to create a routine that takes 
               ! s0_predicted_edge (the edge-centered, predicted 1/2 time 
               ! base state quantities, and put these onto a Cartesian grid, 
               ! on the edges.

               if (spherical .eq. 1) &
                    call bl_error("ERROR: spherical not yet implemented in make_edge_rhoX_from_X")

             else
               call make_edge_rhoX_from_X_3d_cart(sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                                  s0_predicted_edge(n,:,:), &
                                                  lo, hi)
             end if
          end select
       end do

    end do

    call destroy(bpt)
    
  end subroutine make_edge_rhoX_from_X

  subroutine make_edge_rhoX_from_X_2d(sx,sy, &
                                      s0_predicted_edge, &
                                      lo,hi)

    use bl_constants_module
    use variables,     only: rho_comp, spec_comp
    use network, only: nspec

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(in   ) :: s0_predicted_edge(0:,:)

    real(kind=dp_t) :: rho0_edge, X0_edge
    integer :: i, j, comp    

    ! x-interfaces
    do comp = spec_comp, spec_comp+nspec-1

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
                    
             ! s0_predicted_edge is on the y-edges.  Here, we simply average
             ! the two vertical edge states to find the edge state on the
             ! x-edges
             X0_edge = HALF*(s0_predicted_edge(j,  comp) + &
                             s0_predicted_edge(j+1,comp))

             rho0_edge = HALF*(s0_predicted_edge(j,  rho_comp) + &
                               s0_predicted_edge(j+1,rho_comp))

             sx(i,j,comp) = (rho0_edge + sx(i,j,rho_comp))*sx(i,j,comp) - &
                            rho0_edge*X0_edge
                    
          enddo
       enddo

    enddo


    ! y-interfaces
    do comp = spec_comp, spec_comp+nspec-1

       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
          
             ! s0_predicted_edge is on the y-edges, so we use them directly.
             sy(i,j,comp) = (s0_predicted_edge(j,rho_comp) + &
                             sy(i,j,rho_comp))*sy(i,j,comp) - &
                            (s0_predicted_edge(j,rho_comp)* &
                             s0_predicted_edge(j,comp))
             
          enddo
       enddo

    enddo
    
  end subroutine make_edge_rhoX_from_X_2d
  

  subroutine make_edge_rhoX_from_X_3d_cart(sx,sy,sz, &
                                           s0_predicted_edge, &
                                           lo,hi)

    use variables,     only: rho_comp, spec_comp
    use network, only: nspec
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sz(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(in   ) :: s0_predicted_edge(0:,:)

    real(kind=dp_t) :: rho0_edge, X0_edge    
    integer :: i, j, k, comp
    
    ! x edge
    do comp = spec_comp, spec_comp+nspec-1

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                ! s0_predicted_edge is on the z-edges.  Here, we simply average
                ! the two vertical edge states to find the edge state on the
                ! x-edges
                X0_edge = HALF*(s0_predicted_edge(k,  comp) + &
                                s0_predicted_edge(k+1,comp))

                rho0_edge = HALF*(s0_predicted_edge(k,  rho_comp) + &
                                  s0_predicted_edge(k+1,rho_comp))

                sx(i,j,k,comp) = (rho0_edge + &
                                  sx(i,j,k,rho_comp))*sx(i,j,k,comp) - &
                                 rho0_edge*X0_edge
                
                
             enddo
          enddo
       enddo

    enddo


    ! y edge
    do comp = spec_comp, spec_comp+nspec-1

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                ! s0_predicted_edge is on the z-edges.  Here, we simply average
                ! the two vertical edge states to find the edge state on the
                ! y-edges
                X0_edge = HALF*(s0_predicted_edge(k,  comp) + &
                                s0_predicted_edge(k+1,comp))

                rho0_edge = HALF*(s0_predicted_edge(k,  rho_comp) + &
                                  s0_predicted_edge(k+1,rho_comp))

                sy(i,j,k,comp) = (rho0_edge + &
                                  sy(i,j,k,rho_comp))*sy(i,j,k,comp) - &
                                 rho0_edge*X0_edge
                
             enddo
          enddo
       enddo

    enddo


    ! z edge
    do comp = spec_comp, spec_comp+nspec-1

       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! s0_predicted_edge is on the z-edges, so we use them directly.
                sz(i,j,k,comp) = (s0_predicted_edge(k,rho_comp) + &
                                  sz(i,j,k,rho_comp))*sz(i,j,k,comp) - &
                                 (s0_predicted_edge(k,rho_comp)* &
                                  s0_predicted_edge(k,comp))
             
             enddo
          enddo
       enddo
    
    enddo

  end subroutine make_edge_rhoX_from_X_3d_cart

end module convert_rhoX_to_X_module
