module update_scal_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_scal

contains

  subroutine update_scal(nlevs,nstart,nstop,sold,snew,umac,w0,w0_cart_vec, &
                         sedge,sflux,scal_force,s0_old,s0_edge_old,s0_new,s0_edge_new, &
                         s0_old_cart,s0_new_cart,dx,dt,the_bc_level,mla)

    use bl_prof_module
    use bl_constants_module
    use geometry,  only: spherical
    use variables, only: spec_comp, rho_comp
    use network,   only: nspec
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module

    integer           , intent(in   ) :: nlevs, nstart, nstop
    type(multifab)    , intent(in   ) :: sold(:)
    type(multifab)    , intent(inout) :: snew(:)
    type(multifab)    , intent(in   ) :: umac(:,:)
    real(kind=dp_t)   , intent(in   ) :: w0(:,0:)
    type(multifab)    , intent(in   ) :: w0_cart_vec(:)
    type(multifab)    , intent(in   ) :: sedge(:,:)
    type(multifab)    , intent(in   ) :: sflux(:,:)
    type(multifab)    , intent(in   ) :: scal_force(:)
    real(kind = dp_t) , intent(in   ) :: s0_old(:,0:,:), s0_edge_old(:,0:,:)
    real(kind = dp_t) , intent(in   ) :: s0_new(:,0:,:), s0_edge_new(:,0:,:)
    type(multifab)    , intent(in   ) :: s0_old_cart(:)
    type(multifab)    , intent(in   ) :: s0_new_cart(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:),dt
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
    integer :: i,ng_s,dm,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "update_scal")

    dm = sold(1)%dim
    ng_s = sold(1)%ng

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
             call update_scal_2d(nstart, nstop, &
                                 sop(:,:,1,:), snp(:,:,1,:), &
                                 ump(:,:,1,1), vmp(:,:,1,1), &
                                 w0(n,:), &
                                 sepx(:,:,1,:), sepy(:,:,1,:), &
                                 sfpx(:,:,1,:), sfpy(:,:,1,:), &
                                 fp(:,:,1,:), &
                                 s0_old(n,:,:), s0_edge_old(n,:,:), &
                                 s0_new(n,:,:), s0_edge_new(n,:,:), &
                                 lo, hi, ng_s, dx(n,:), dt)
          case (3)
             wmp => dataptr(umac(n,3),i)
             sepz => dataptr(sedge(n,3),i)
             sfpz => dataptr(sflux(n,3),i)
             if (spherical .eq. 0) then
                call update_scal_3d_cart(nstart, nstop, &
                                         sop(:,:,:,:), snp(:,:,:,:), &
                                         ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                         w0(n,:), &
                                         sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                         sfpx(:,:,:,:), sfpy(:,:,:,:), sfpz(:,:,:,:), &
                                         fp(:,:,:,:), &
                                         s0_old(n,:,:), s0_edge_old(n,:,:), &
                                         s0_new(n,:,:), s0_edge_new(n,:,:), &
                                         lo, hi, ng_s, dx(n,:), dt)
             else
                s0op => dataptr(s0_old_cart(n), i)
                s0np => dataptr(s0_new_cart(n), i)
                w0p => dataptr(w0_cart_vec(n),i)
                call update_scal_3d_sphr(nstart, nstop, &
                                         sop(:,:,:,:), snp(:,:,:,:), &
                                         ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                         w0p(:,:,:,:), &
                                         sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                         sfpx(:,:,:,:), sfpy(:,:,:,:), sfpz(:,:,:,:), &
                                         fp(:,:,:,:), &
                                         s0_old(n,:,:), s0_new(n,:,:), &
                                         s0op(:,:,:,:), s0np(:,:,:,:), &
                                         lo, hi, domlo, domhi, ng_s, dx(n,:), dt)
             end if
          end select
       end do

    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(snew(nlevs),nstart,nstop-nstart+1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(snew(nlevs),nstart,dm+nstart,nstop-nstart+1, &
                            the_bc_level(nlevs))

       ! do the same for density if we updated the species
       if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
          call multifab_fill_boundary_c(snew(nlevs),rho_comp,1)
          call multifab_physbc(snew(nlevs),rho_comp,dm+rho_comp,1,the_bc_level(nlevs))

       endif

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(snew(n-1),nstart,snew(n),nstart,mla%mba%rr(n-1,:), &
                                   nstop-nstart+1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(snew(n),snew(n-1),ng_s,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         nstart,dm+nstart,nstop-nstart+1)

          ! do the same for density if we updated the species
          if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
             call ml_cc_restriction_c(snew(n-1),rho_comp,snew(n),rho_comp, &
                                      mla%mba%rr(n-1,:),1)
             call multifab_fill_ghost_cells(snew(n),snew(n-1),ng_s,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1),the_bc_level(n), &
                                            rho_comp,dm+rho_comp,1)
          endif

       enddo

    end if

    call destroy(bpt)

  end subroutine update_scal

  subroutine update_scal_2d(nstart,nstop,sold,snew,umac,vmac,w0, &
                            sedgex,sedgey,sfluxx,sfluxy,force,base_old,base_old_edge, &
                            base_new,base_new_edge,lo,hi,ng_s,dx,dt)

    use network,       only: nspec
    use probin_module, only: enthalpy_pred_type
    use variables,     only: spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac
    use pred_parameters
    use bl_constants_module


    ! update each scalar in time.  Here, it is assumed that the edge
    ! states (sedgex and sedgey) are for the perturbational quantities.

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:), ng_s
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_s:,lo(2)-ng_s:,:)
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
    real (kind = dp_t), intent(in   ) :: dt,dx(:)

    integer            :: i, j, comp, comp2
    real (kind = dp_t) :: delta_base,divterm
    real (kind = dp_t) :: delta,frac,sum
    real (kind = dp_t) :: smin(nstart:nstop),smax(nstart:nstop)
    logical :: test

    do comp = nstart, nstop

       ! test = T means the edge states are NOT in perturbational form
       test = ( (comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1) ) &
         .or. ( (comp.eq.rhoh_comp).and. &
                     ( enthalpy_pred_type.eq.predict_h .or. &
                       enthalpy_pred_type.eq.predict_T_then_h ) ) &
         .or. ( (comp.ge.trac_comp).and.(comp.le.trac_comp+ntrac-1) )

       if (test) then

          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                divterm = (sfluxx(i+1,j,comp) - sfluxx(i,j,comp))/dx(1) &
                        + (sfluxy(i,j+1,comp) - sfluxy(i,j,comp))/dx(2)

                snew(i,j,comp) = sold(i,j,comp) + dt*(-divterm + force(i,j,comp))

             end do
          end do

       else

          do j=lo(2),hi(2)

             delta_base = base_new(j,comp) - base_old(j,comp)

             do i = lo(1), hi(1)

                divterm = (sfluxx(i+1,j,comp) - sfluxx(i,j,comp))/dx(1) &
                        + (sfluxy(i,j+1,comp) - sfluxy(i,j,comp))/dx(2)

                snew(i,j,comp) = sold(i,j,comp) &
                     + delta_base + dt*(-divterm + force(i,j,comp))
                
             end do
          end do

       end if
    enddo

    ! Define the update to rho as the sum of the updates to (rho X)_i
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

       smin(:) =  HUGE(smin)
       smax(:) = -HUGE(smax)

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
    end if

  end subroutine update_scal_2d

  subroutine update_scal_3d_cart(nstart,nstop,sold,snew,umac,vmac,wmac,w0, &
                                 sedgex,sedgey,sedgez,sfluxx,sfluxy,sfluxz,force, &
                                 base_old,base_old_edge,base_new,base_new_edge,lo,hi, &
                                 ng_s,dx,dt)
    use network,       only: nspec
    use probin_module, only: enthalpy_pred_type
    use variables,     only: spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac
    use pred_parameters
    use bl_constants_module


    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:), ng_s
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
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
    real (kind = dp_t), intent(in   ) :: dt,dx(:)

    integer            :: i, j, k, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: delta,frac,sum,delta_base
    real (kind = dp_t) :: smin(nstart:nstop),smax(nstart:nstop)
    logical            :: test

    do comp = nstart, nstop
    
       ! test = T means the edge states are NOT in perturbational form
       test = ( (comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1) ) &
         .or. ( (comp.eq.rhoh_comp).and. &
                     ( enthalpy_pred_type.eq.predict_h .or. &
                       enthalpy_pred_type.eq.predict_T_then_h ) ) &
         .or. ( (comp.ge.trac_comp).and.(comp.le.trac_comp+ntrac-1) )

       if (test) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
   
                   divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                           + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                           + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)
   
                   snew(i,j,k,comp) = sold(i,j,k,comp) &
                        + dt * (-divterm + force(i,j,k,comp))
   
                enddo
             enddo
          enddo

       else

          do k = lo(3), hi(3)

             delta_base = base_new(k,comp) - base_old(k,comp)

             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                           + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                           + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)
   
                   snew(i,j,k,comp) = sold(i,j,k,comp) + delta_base &
                        + dt * (-divterm + force(i,j,k,comp))
   
                enddo
             enddo
          enddo

       end if
    end do

    ! Define the update to rho as the sum of the updates to (rho X)_i
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

       smin(:) =  HUGE(smin)
       smax(:) = -HUGE(smax)

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
    end if

  end subroutine update_scal_3d_cart

  subroutine update_scal_3d_sphr(nstart,nstop,sold,snew,umac,vmac,wmac, &
                                 w0_cart,sedgex,sedgey,sedgez,sfluxx,sfluxy,sfluxz, &
                                 force,base_old,base_new,base_old_cart,base_new_cart, &
                                 lo,hi,domlo,domhi,ng_s,dx,dt)
    use network,       only: nspec
    use probin_module, only: enthalpy_pred_type
    use variables,     only: spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop
    integer           , intent(in   ) :: lo(:), hi(:), domlo(:), domhi(:), ng_s
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
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

    integer            :: i, j, k, comp, comp2
    real (kind = dp_t) :: divterm,delta_base
    real (kind = dp_t) :: delta,frac,sum
    real (kind = dp_t) :: smin(nstart:nstop),smax(nstart:nstop)
    logical            :: test

    ! is spherical

    do comp = nstart, nstop
    
       ! test = T means the edge states are NOT in perturbational form
       test = ( (comp.ge.spec_comp).and.(comp.le.spec_comp+nspec-1) ) &
         .or. ( (comp.eq.rhoh_comp).and. &
                     ( enthalpy_pred_type.eq.predict_h .or. &
                       enthalpy_pred_type.eq.predict_T_then_h ) ) &
         .or. ( (comp.ge.trac_comp).and.(comp.le.trac_comp+ntrac-1) )

       if (test) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                           + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                           + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)

                   snew(i,j,k,comp) = sold(i,j,k,comp) &
                        + dt * (-divterm + force(i,j,k,comp))
                enddo
             enddo
          enddo

       else

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   delta_base = base_new_cart(i,j,k,comp) - base_old_cart(i,j,k,comp)

                   divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                           + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                           + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)

                   snew(i,j,k,comp) = sold(i,j,k,comp) + delta_base &
                        + dt * (-divterm + force(i,j,k,comp))

                enddo
             enddo
          enddo

       end if
    end do

    ! Define the update to rho as the sum of the updates to (rho X)_i
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

       smin(:) =  HUGE(smin)
       smax(:) = -HUGE(smax)

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
    end if

  end subroutine update_scal_3d_sphr

end module update_scal_module
