module update_scal_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_scal

contains

  subroutine update_scal(nlevs,nstart,nstop,sold,snew,sflux,scal_force,rhoh0_old,rhoh0_new, &
                         rhoh0_old_cart,rhoh0_new_cart,p0,tempbar,dx,dt,the_bc_level,mla)

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
    type(multifab)    , intent(in   ) :: sflux(:,:)
    type(multifab)    , intent(in   ) :: scal_force(:)
    real(kind = dp_t) , intent(in   ) :: rhoh0_old(:,0:)
    real(kind = dp_t) , intent(in   ) :: rhoh0_new(:,0:)
    real(kind = dp_t) , intent(in   ) ::       p0(:,0:)
    real(kind = dp_t) , intent(in   ) ::  tempbar(:,0:)
    type(multifab)    , intent(in   ) :: rhoh0_old_cart(:)
    type(multifab)    , intent(in   ) :: rhoh0_new_cart(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:),dt
    type(bc_level)    , intent(in   ) :: the_bc_level(:)
    type(ml_layout)   , intent(inout) :: mla

    ! local
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: sfpx(:,:,:,:)
    real(kind=dp_t), pointer :: sfpy(:,:,:,:)
    real(kind=dp_t), pointer :: sfpz(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: rhoh0op(:,:,:,:)
    real(kind=dp_t), pointer :: rhoh0np(:,:,:,:)

    integer :: lo(sold(1)%dim),hi(sold(1)%dim)
    integer :: i,ng_s,dm,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "update_scal")

    dm = sold(1)%dim
    ng_s = sold(1)%ng

    do n=1,nlevs

       do i = 1, sold(n)%nboxes
          if ( multifab_remote(sold(n),i) ) cycle
          sop => dataptr(sold(n),i)
          snp => dataptr(snew(n),i)
          sfpx => dataptr(sflux(n,1),i)
          sfpy => dataptr(sflux(n,2),i)
          fp => dataptr(scal_force(n),i)
          lo =  lwb(get_box(sold(n),i))
          hi =  upb(get_box(sold(n),i))
          select case (dm)
          case (2)
             call update_scal_2d(nstart, nstop, &
                                 sop(:,:,1,:), snp(:,:,1,:), &
                                 sfpx(:,:,1,:), sfpy(:,:,1,:), &
                                 fp(:,:,1,:), &
                                 rhoh0_old(n,:), rhoh0_new(n,:), &
                                 p0(n,:), tempbar(n,:), &
                                 lo, hi, ng_s, dx(n,:), dt)
          case (3)
             sfpz => dataptr(sflux(n,3),i)
             if (spherical .eq. 0) then
                call update_scal_3d_cart(nstart, nstop, &
                                         sop(:,:,:,:), snp(:,:,:,:), &
                                         sfpx(:,:,:,:), sfpy(:,:,:,:), sfpz(:,:,:,:), &
                                         fp(:,:,:,:), &
                                         rhoh0_old(n,:), rhoh0_new(n,:), &
                                         p0(n,:), tempbar(n,:), &
                                         lo, hi, ng_s, dx(n,:), dt)
             else
                rhoh0op => dataptr(rhoh0_old_cart(n), i)
                rhoh0np => dataptr(rhoh0_new_cart(n), i)
                call update_scal_3d_sphr(nstart, nstop, &
                                         sop(:,:,:,:), snp(:,:,:,:), &
                                         sfpx(:,:,:,:), sfpy(:,:,:,:), sfpz(:,:,:,:), &
                                         fp(:,:,:,:), &
                                         rhoh0op(:,:,:,1), rhoh0np(:,:,:,1), &
                                         p0(n,:), tempbar(n,:), &
                                         lo, hi, ng_s, dx(n,:), dt)
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

  subroutine update_scal_2d(nstart,nstop,sold,snew,sfluxx,sfluxy,force, &
                            rhoh0_old,rhoh0_new,p0,tempbar,lo,hi,ng_s,dx,dt)

    use network,       only: nspec
    use eos_module
    use probin_module, only: enthalpy_pred_type, do_eos_h_above_cutoff, base_cutoff_density
    use variables,     only: spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:), ng_s
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxx(lo(1)   :,lo(2)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxy(lo(1)   :,lo(2)   :,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   rhoh0_old(0:)
    real (kind = dp_t), intent(in   ) ::   rhoh0_new(0:)
    real (kind = dp_t), intent(in   ) ::       p0(0:)
    real (kind = dp_t), intent(in   ) ::  tempbar(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)

    integer            :: i, j, comp, comp2
    real (kind = dp_t) :: delta_rhoh0,divterm
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

             delta_rhoh0 = rhoh0_new(j) - rhoh0_old(j)

             do i = lo(1), hi(1)

                divterm = (sfluxx(i+1,j,comp) - sfluxx(i,j,comp))/dx(1) &
                        + (sfluxy(i,j+1,comp) - sfluxy(i,j,comp))/dx(2)

                snew(i,j,comp) = sold(i,j,comp) &
                     + delta_rhoh0 + dt*(-divterm + force(i,j,comp))
                
             end do
          end do

       end if
    enddo

    if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then

       do j = lo(2), hi(2)
       do i = lo(1), hi(1)
 
          if (snew(i,j,rho_comp) < base_cutoff_density) then
            den_eos(1) = snew(i,j,rho_comp)
            temp_eos(1) = tempbar(j)
            p_eos(1) = p0(j)
            xn_eos(1,:) = snew(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

            ! (rho,P) --> T,h
            call eos(eos_input_rp, den_eos, temp_eos, &
                     npts, nspec, &
                     xn_eos, &
                     p_eos, h_eos, e_eos, &
                     cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                     dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                     dpdX_eos, dhdX_eos, &
                     gam1_eos, cs_eos, s_eos, &
                     dsdt_eos, dsdr_eos, &
                     do_diag)
   
            snew(i,j,rhoh_comp) = snew(i,j,rho_comp) * h_eos(1)

          end if
 
       enddo
       enddo

    end if

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

  subroutine update_scal_3d_cart(nstart,nstop,sold,snew,sfluxx,sfluxy,sfluxz,force, &
                                 rhoh0_old,rhoh0_new,p0,tempbar,lo,hi,ng_s,dx,dt)

    use network,       only: nspec
    use eos_module
    use probin_module, only: enthalpy_pred_type, do_eos_h_above_cutoff, base_cutoff_density
    use variables,     only: spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac
    use pred_parameters
    use bl_constants_module


    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:), ng_s
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxx(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxy(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxz(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   rhoh0_old(0:)
    real (kind = dp_t), intent(in   ) ::   rhoh0_new(0:)
    real (kind = dp_t), intent(in   ) ::       p0(0:)
    real (kind = dp_t), intent(in   ) ::  tempbar(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)

    integer            :: i, j, k, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: delta,frac,sum,delta_rhoh0
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

             delta_rhoh0 = rhoh0_new(k) - rhoh0_old(k)

             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                           + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                           + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)
   
                   snew(i,j,k,comp) = sold(i,j,k,comp) + delta_rhoh0 &
                        + dt * (-divterm + force(i,j,k,comp))
   
                enddo
             enddo
          enddo

       end if
    end do


    if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then

       do k = lo(3), hi(3)
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          if (snew(i,j,k,rho_comp) < base_cutoff_density) then
            den_eos(1) = snew(i,j,k,rho_comp)
            temp_eos(1) = tempbar(k)
            p_eos(1) = p0(k)
            xn_eos(1,:) = snew(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

            ! (rho,P) --> T,h
            call eos(eos_input_rp, den_eos, temp_eos, &
                     npts, nspec, &
                     xn_eos, &
                     p_eos, h_eos, e_eos, &
                     cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                     dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                     dpdX_eos, dhdX_eos, &
                     gam1_eos, cs_eos, s_eos, &
                     dsdt_eos, dsdr_eos, &
                     do_diag)

            snew(i,j,k,rhoh_comp) = snew(i,j,k,rho_comp) * h_eos(1)

          end if

       enddo
       enddo
       enddo

    end if


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

  subroutine update_scal_3d_sphr(nstart,nstop,sold,snew,sfluxx,sfluxy,sfluxz,force, &
                                 rhoh0_old_cart,rhoh0_new_cart,p0,tempbar,lo,hi,ng_s,dx,dt)

    use network,       only: nspec
    use eos_module
    use probin_module, only: enthalpy_pred_type, do_eos_h_above_cutoff, base_cutoff_density
    use variables,     only: spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop
    integer           , intent(in   ) :: lo(:), hi(:), ng_s
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxx(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxy(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sfluxz(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   rhoh0_old_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(in   ) ::   rhoh0_new_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(in   ) ::       p0(0:)
    real (kind = dp_t), intent(in   ) ::  tempbar(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)

    integer            :: i, j, k, comp, comp2
    real (kind = dp_t) :: divterm,delta_rhoh0
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

                   delta_rhoh0 = rhoh0_new_cart(i,j,k) - rhoh0_old_cart(i,j,k)

                   divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                           + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                           + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)

                   snew(i,j,k,comp) = sold(i,j,k,comp) + delta_rhoh0 &
                        + dt * (-divterm + force(i,j,k,comp))

                enddo
             enddo
          enddo

       end if

    end do

    if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then

       do k = lo(3), hi(3) 
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          if (snew(i,j,k,rho_comp) < base_cutoff_density) then
            den_eos(1) = snew(i,j,k,rho_comp)
            temp_eos(1) = tempbar(k)
            p_eos(1) = p0(k)
            xn_eos(1,:) = snew(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

            ! (rho,P) --> T,h
            call eos(eos_input_rp, den_eos, temp_eos, &
                     npts, nspec, &
                     xn_eos, &
                     p_eos, h_eos, e_eos, &
                     cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                     dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                     dpdX_eos, dhdX_eos, &
                     gam1_eos, cs_eos, s_eos, &
                     dsdt_eos, dsdr_eos, &
                     do_diag)

            snew(i,j,k,rhoh_comp) = snew(i,j,k,rho_comp) * h_eos(1)

          end if

       enddo
       enddo
       enddo

    end if


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
