module update_scal_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_scal

contains

  subroutine update_scal(nlevs,which_step,nstart,nstop,sold,snew,umac,w0,w0_cart_vec,eta, &
                         sedge,sflux,scal_force,s0_old,s0_edge_old,s0_new,s0_edge_new, &
                         s0_old_cart,s0_new_cart,dx,dt,evolve_base_state,the_bc_level,mla)

    use bl_constants_module
    use geometry,  only: spherical
    use variables, only: spec_comp, rho_comp
    use network,   only: nspec
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module

    integer           , intent(in   ) :: nlevs, which_step, nstart, nstop
    type(multifab)    , intent(in   ) :: sold(:)
    type(multifab)    , intent(inout) :: snew(:)
    type(multifab)    , intent(in   ) :: umac(:,:)
    real(kind=dp_t)   , intent(in   ) :: w0(:,0:)
    type(multifab)    , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t)   , intent(inout) :: eta(:,0:,:)
    type(multifab)    , intent(in   ) :: sedge(:,:)
    type(multifab)    , intent(in   ) :: sflux(:,:)
    type(multifab)    , intent(in   ) :: scal_force(:)
    real(kind = dp_t) , intent(in   ) :: s0_old(:,0:,:), s0_edge_old(:,0:,:)
    real(kind = dp_t) , intent(in   ) :: s0_new(:,0:,:), s0_edge_new(:,0:,:)
    type(multifab)    , intent(in   ) :: s0_old_cart(:)
    type(multifab)    , intent(in   ) :: s0_new_cart(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:),dt
    logical           , intent(in   ) :: evolve_base_state
    type(bc_level)    , intent(in   ) :: the_bc_level(:)
    type(ml_layout)   , intent(inout) :: mla

    ! local
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: sfpx(:,:,:,:)
    real(kind=dp_t), pointer :: sfpy(:,:,:,:)
    real(kind=dp_t), pointer :: sfpz(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)
    real(kind=dp_t), pointer :: s0op(:,:,:,:)
    real(kind=dp_t), pointer :: s0np(:,:,:,:)

    type(box) :: domain

    integer :: domlo(sold(1)%dim),domhi(sold(1)%dim)
    integer :: lo(sold(1)%dim),hi(sold(1)%dim)
    integer :: i,ng,dm,n

    dm = sold(1)%dim
    ng = sold(1)%ng

    do n=1,nlevs

       domain = layout_get_pd(sold(n)%la)
       domlo = lwb(domain)
       domhi = upb(domain)

       do i = 1, sold(n)%nboxes
          if ( multifab_remote(sold(n),i) ) cycle
          sop => dataptr(sold(n),i)
          snp => dataptr(snew(n),i)
          ump => dataptr(umac(n,1),i)
          vmp => dataptr(umac(n,2),i)
          sepx => dataptr(sedge(n,1),i)
          sepy => dataptr(sedge(n,2),i)
          sfpx => dataptr(sflux(n,1),i)
          sfpy => dataptr(sflux(n,2),i)
          fp => dataptr(scal_force(n),i)
          lo =  lwb(get_box(sold(n),i))
          hi =  upb(get_box(sold(n),i))
          select case (dm)
          case (2)
             call update_scal_2d(which_step, nstart, nstop, &
                                 sop(:,:,1,:), snp(:,:,1,:), &
                                 ump(:,:,1,1), vmp(:,:,1,1), &
                                 w0(n,:), eta(n,:,:), &
                                 sepx(:,:,1,:), sepy(:,:,1,:), &
                                 sfpx(:,:,1,:), sfpy(:,:,1,:), &
                                 fp(:,:,1,:), &
                                 s0_old(n,:,:), s0_edge_old(n,:,:), &
                                 s0_new(n,:,:), s0_edge_new(n,:,:), &
                                 lo, hi, ng, dx(n,:), dt, evolve_base_state)
          case (3)
             wmp => dataptr(umac(n,3),i)
             sepz => dataptr(sedge(n,3),i)
             sfpz => dataptr(sflux(n,3),i)
             if (spherical .eq. 0) then
                call update_scal_3d_cart(which_step, nstart, nstop, &
                                         sop(:,:,:,:), snp(:,:,:,:), &
                                         ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                         w0(n,:), eta(n,:,:), &
                                         sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                         sfpx(:,:,:,:), sfpy(:,:,:,:), sfpz(:,:,:,:), &
                                         fp(:,:,:,:), &
                                         s0_old(n,:,:), s0_edge_old(n,:,:), &
                                         s0_new(n,:,:), s0_edge_new(n,:,:), &
                                         lo, hi, ng, dx(n,:), dt, evolve_base_state)
             else
                s0op => dataptr(s0_old_cart(n), i)
                s0np => dataptr(s0_new_cart(n), i)
                w0p => dataptr(w0_cart_vec(n),i)
                call update_scal_3d_sphr(which_step, nstart, nstop, &
                                         sop(:,:,:,:), snp(:,:,:,:), &
                                         ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                         w0p(:,:,:,:), &
                                         sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                         sfpx(:,:,:,:), sfpy(:,:,:,:), sfpz(:,:,:,:), &
                                         fp(:,:,:,:), &
                                         s0_old(n,:,:), s0_new(n,:,:), &
                                         s0op(:,:,:,:), s0np(:,:,:,:), &
                                         lo, hi, domlo, domhi, ng, dx(n,:), dt, &
                                         evolve_base_state)
             end if
          end select
       end do

       call multifab_fill_boundary_c(snew(n),nstart,nstop-nstart+1)
       call multifab_physbc(snew(n),nstart,dm+nstart,nstop-nstart+1,dx(n,:), &
                            the_bc_level(n))

       if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
          call multifab_fill_boundary_c(snew(n),rho_comp,1)
          call multifab_physbc(snew(n),rho_comp,dm+rho_comp,1,dx(n,:),the_bc_level(n))
       endif

    end do

    do n=nlevs,2,-1
       call ml_cc_restriction_c(snew(n-1),nstart,snew(n),nstart,mla%mba%rr(n-1,:), &
                                nstop-nstart+1)
       call multifab_fill_ghost_cells(snew(n),snew(n-1),ng,mla%mba%rr(n-1,:), &
                                      the_bc_level(n-1),the_bc_level(n  ), &
                                      nstart,dm+nstart,nstop-nstart+1)

       if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
          call ml_cc_restriction_c(snew(n-1),rho_comp,snew(n),rho_comp,mla%mba%rr(n-1,:),1)
          call multifab_fill_ghost_cells(snew(n),snew(n-1),ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         rho_comp,dm+rho_comp,1)
       endif
    enddo

  end subroutine update_scal

  subroutine update_scal_2d(which_step,nstart,nstop,sold,snew,umac,vmac,w0,eta, &
                            sedgex,sedgey,sfluxx,sfluxy,force,base_old,base_old_edge, &
                            base_new,base_new_edge,lo,hi,ng,dx,dt,evolve_base_state)

    use network, only: nspec
    use variables, only: spec_comp, rho_comp
    use bl_constants_module

    ! update each scalar in time.  Here, it is assumed that the edge
    ! states (sedgex and sedgey) are for the perturbational quantities.

    integer           , intent(in   ) :: which_step, nstart, nstop, lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)
    real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)
    real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)
    real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxx(lo(1)   :,lo(2)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxy(lo(1)   :,lo(2)   :,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   base_old(0:,:), base_old_edge(0:,:)
    real (kind = dp_t), intent(in   ) ::   base_new(0:,:), base_new_edge(0:,:)
    real (kind = dp_t), intent(in   ) :: w0(0:)
    real (kind = dp_t), intent(inout) :: eta(0:,:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: evolve_base_state

    integer :: i, j, comp, comp2
    real (kind = dp_t) :: delta_base,divterm
    real (kind = dp_t) :: delta,frac,sum,fac
    real (kind = dp_t), allocatable :: smin(:),smax(:)

    fac = ONE / dble(hi(1)-lo(1)+1)

    if (evolve_base_state) then
       if (which_step .eq. 1) then
          do comp = nstart, nstop
             eta(:,comp) = ZERO
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)
                   eta(j,comp) = eta(j,comp) + vmac(i,j)*sedgey(i,j,comp)
                end do
                eta(j,comp) = eta(j,comp) * fac
             end do
          end do
       end if
    end if

    do comp = nstart, nstop
       do j = lo(2), hi(2)

          delta_base = base_new(j,comp) - base_old(j,comp)

          do i = lo(1), hi(1)

             divterm = (sfluxx(i+1,j,comp) - sfluxx(i,j,comp))/dx(1) &
                  + (sfluxy(i,j+1,comp) - sfluxy(i,j,comp))/dx(2)

             snew(i,j,comp) = sold(i,j,comp) + delta_base &
                  - dt * divterm + dt * force(i,j,comp)

             if (evolve_base_state) then
                if (which_step .eq. 2) then
                   snew(i,j,comp) = snew(i,j,comp) + dt / dx(2) * (eta(j+1,comp)-eta(j,comp))
                end if
             end if

          enddo
       enddo
    enddo

    ! Define the update to rho as the sum of the updates to (rho X)_i
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

       allocate(smin(nstart:nstop),smax(nstart:nstop))
       smin(:) =  1.e20
       smax(:) = -1.e20

       snew(:,:,rho_comp) = sold(:,:,rho_comp)

       do comp = nstart, nstop
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                snew(i,j,rho_comp) = snew(i,j,rho_comp) + (snew(i,j,comp)-sold(i,j,comp))
                smin(comp) = min(smin(comp),snew(i,j,comp))
                smax(comp) = max(smax(comp),snew(i,j,comp))
             enddo
          enddo
       enddo

       if (which_step .eq. 1) then
          eta(:,rho_comp) = ZERO
          do comp = nstart, nstop
             do j = lo(2), hi(2)+1
                eta(j,rho_comp) = eta(j,rho_comp) + eta(j,comp)
             end do
          end do
       end if

    end if

    ! Do not allow the species to leave here negative.
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       do comp = nstart, nstop
          if (smin(comp) .lt. ZERO) then 
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (snew(i,j,comp) .lt. ZERO) then
                      delta = -snew(i,j,comp)
                      sum = ZERO 
                      do comp2 = nstart, nstop
                         if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                            sum = sum + snew(i,j,comp2)
                         end if
                      enddo
                      do comp2 = nstart, nstop
                         if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                            frac = snew(i,j,comp2) / sum
                            snew(i,j,comp2) = snew(i,j,comp2) - frac * delta
                         end if
                      enddo
                      snew(i,j,comp) = ZERO
                   end if
                enddo
             enddo
          end if
       enddo
       deallocate(smin,smax)
    end if

  end subroutine update_scal_2d

  subroutine update_scal_3d_cart(which_step,nstart,nstop,sold,snew,umac,vmac,wmac,w0, &
                                 eta,sedgex,sedgey,sedgez,sfluxx,sfluxy,sfluxz,force, &
                                 base_old,base_old_edge,base_new,base_new_edge,lo,hi, &
                                 ng,dx,dt,evolve_base_state)
    use network, only: nspec
    use variables, only: spec_comp, rho_comp
    use bl_constants_module

    integer           , intent(in   ) :: which_step, nstart, nstop, lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(inout) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(inout) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(inout) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxx(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxy(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxz(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   base_old(0:,:), base_old_edge(0:,:)
    real (kind = dp_t), intent(in   ) ::   base_new(0:,:), base_new_edge(0:,:)
    real (kind = dp_t), intent(in   ) :: w0(0:)
    real (kind = dp_t), intent(inout) :: eta(0:,:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: evolve_base_state

    integer :: i, j, k, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: delta,frac,sum,delta_base,fac
    real (kind = dp_t), allocatable :: smin(:),smax(:)

    fac = ONE / dble( (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1) )

    if (evolve_base_state) then
       if (which_step .eq. 1) then
          do comp = nstart, nstop
             eta(:,comp) = ZERO
             do k = lo(3), hi(3)+1
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      eta(k,comp) = eta(k,comp) + wmac(i,j,k)*sedgez(i,j,k,comp)
                   end do
                end do
                eta(k,comp) = eta(k,comp) * fac
             end do
          end do
       end if
    end if

    do comp = nstart, nstop

       do k = lo(3), hi(3)

          delta_base = base_new(k,comp) - base_old(k,comp)

          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                     + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                     + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)

                snew(i,j,k,comp) = sold(i,j,k,comp) + delta_base &
                     - dt * divterm + dt * force(i,j,k,comp)

                if (evolve_base_state) then
                   if (which_step .eq. 2) then
                      snew(i,j,k,comp) = snew(i,j,k,comp) &
                           + dt / dx(3) * (eta(k+1,comp)-eta(k,comp))
                   end if
                end if

             enddo
          enddo
       enddo
    end do

    ! Define the update to rho as the sum of the updates to (rho X)_i
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

       allocate(smin(nstart:nstop),smax(nstart:nstop))
       smin(:) =  1.e20
       smax(:) = -1.e20

       snew(:,:,:,rho_comp) = sold(:,:,:,rho_comp)

       do comp = nstart, nstop
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) &
                        + (snew(i,j,k,comp)-sold(i,j,k,comp))
                   smin(comp) = min(smin(comp),snew(i,j,k,comp))
                   smax(comp) = max(smax(comp),snew(i,j,k,comp))
                enddo
             enddo
          enddo
       enddo
    end if

    ! Do not allow the species to leave here negative.
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       do comp = nstart, nstop
          if (smin(comp) .lt. ZERO) then
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      if (snew(i,j,k,comp) .lt. ZERO) then
                         delta = -snew(i,j,k,comp)
                         sum = ZERO
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) &
                                 sum = sum + snew(i,j,k,comp2)
                         enddo
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                               frac = snew(i,j,k,comp2) / sum
                               snew(i,j,k,comp2) = snew(i,j,k,comp2) - frac * delta
                            end if
                         enddo
                         snew(i,j,k,comp) = ZERO
                      end if
                   enddo
                enddo
             enddo
          end if
       enddo
       deallocate(smin,smax)
    end if

  end subroutine update_scal_3d_cart

  subroutine update_scal_3d_sphr(which_step,nstart,nstop,sold,snew,umac,vmac,wmac, &
                                 w0_cart,sedgex,sedgey,sedgez,sfluxx,sfluxy,sfluxz, &
                                 force,base_old,base_new,base_old_cart,base_new_cart, &
                                 lo,hi,domlo,domhi,ng,dx,dt,evolve_base_state)
    use network, only: nspec
    use variables, only: spec_comp, rho_comp
    use bl_constants_module

    integer           , intent(in   ) :: which_step, nstart, nstop
    integer           , intent(in   ) :: lo(:), hi(:), domlo(:), domhi(:), ng
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(inout) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(inout) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(inout) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxx(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxy(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxz(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   base_old(0:,:)
    real (kind = dp_t), intent(in   ) ::   base_new(0:,:)
    real (kind = dp_t), intent(in   ) ::   base_old_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   base_new_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) :: w0_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: evolve_base_state

    integer :: i, j, k, comp, comp2
    real (kind = dp_t) :: divsu,divbaseu,mult,divterm
    real (kind = dp_t) :: delta,frac,sum
    real (kind = dp_t) :: bc_lox,bc_loy,bc_loz
    real (kind = dp_t) :: bc_hix,bc_hiy,bc_hiz
    real (kind = dp_t), allocatable :: smin(:),smax(:)

    ! is spherical

    do comp = nstart, nstop

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                        + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                        + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)

                   snew(i,j,k,comp) = sold(i,j,k,comp) + &
                        (base_new_cart(i,j,k,comp)-base_old_cart(i,j,k,comp)) &
                        - dt * divterm + dt * force(i,j,k,comp)

                enddo
             enddo
          enddo

    end do

    ! Define the update to rho as the sum of the updates to (rho X)_i
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

       allocate(smin(nstart:nstop),smax(nstart:nstop))
       smin(:) =  1.e20
       smax(:) = -1.e20

       snew(:,:,:,rho_comp) = sold(:,:,:,rho_comp)

       do comp = nstart, nstop
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) &
                        + (snew(i,j,k,comp)-sold(i,j,k,comp))
                   smin(comp) = min(smin(comp),snew(i,j,k,comp))
                   smax(comp) = max(smax(comp),snew(i,j,k,comp))
                enddo
             enddo
          enddo
       enddo
    end if

    ! Do not allow the species to leave here negative.
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       do comp = nstart, nstop
          if (smin(comp) .lt. ZERO) then
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      if (snew(i,j,k,comp) .lt. ZERO) then
                         delta = -snew(i,j,k,comp)
                         sum = ZERO
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) &
                                 sum = sum + snew(i,j,k,comp2)
                         enddo
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                               frac = snew(i,j,k,comp2) / sum
                               snew(i,j,k,comp2) = snew(i,j,k,comp2) - frac * delta
                            end if
                         enddo
                         snew(i,j,k,comp) = ZERO
                      end if
                   enddo
                enddo
             enddo
          end if
       enddo
       deallocate(smin,smax)
    end if

  end subroutine update_scal_3d_sphr

end module update_scal_module
