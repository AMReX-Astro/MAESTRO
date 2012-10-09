module update_scal_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_scal

contains

  subroutine update_scal(mla,nstart,nstop,sold,snew,sflux,scal_force,p0,p0_new_cart, &
                         dx,dt,the_bc_level)

    use bl_prof_module
    use bl_constants_module
    use geometry,  only: spherical
    use variables, only: spec_comp, rho_comp
    use network,   only: nspec
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module

    type(ml_layout)   , intent(inout) :: mla
    integer           , intent(in   ) :: nstart, nstop
    type(multifab)    , intent(in   ) :: sold(:)
    type(multifab)    , intent(inout) :: snew(:)
    type(multifab)    , intent(in   ) :: sflux(:,:)
    type(multifab)    , intent(in   ) :: scal_force(:)
    real(kind = dp_t) , intent(in   ) :: p0(:,0:)
    type(multifab)    , intent(in   ) :: p0_new_cart(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:),dt
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: sfpx(:,:,:,:)
    real(kind=dp_t), pointer :: sfpy(:,:,:,:)
    real(kind=dp_t), pointer :: sfpz(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: p0np(:,:,:,:)

    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: i,n
    integer :: ng_so,ng_sn,ng_sf,ng_f,ng_p

    type(bl_prof_timer), save :: bpt

    call build(bpt, "update_scal")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_so = nghost(sold(1))
    ng_sn = nghost(snew(1))
    ng_sf = nghost(sflux(1,1))
    ng_f  = nghost(scal_force(1))
    ng_p  = nghost(p0_new_cart(1))

    do n=1,nlevs

       do i = 1, nboxes(sold(n))
          if ( multifab_remote(sold(n),i) ) cycle
          sop => dataptr(sold(n),i)
          snp => dataptr(snew(n),i)
          sfpx => dataptr(sflux(n,1),i)
          fp => dataptr(scal_force(n),i)
          lo =  lwb(get_box(sold(n),i))
          hi =  upb(get_box(sold(n),i))
          select case (dm)
          case (1)
             call update_scal_1d(nstart, nstop, &
                                  sop(:,1,1,:), ng_so, &
                                  snp(:,1,1,:), ng_sn, &
                                 sfpx(:,1,1,:), ng_sf, &
                                   fp(:,1,1,:), ng_f, &
                                 p0(n,:), lo, hi, dx(n,:), dt)
          case (2)
             sfpy => dataptr(sflux(n,2),i)
             call update_scal_2d(nstart, nstop, &
                                 sop(:,:,1,:), ng_so, snp(:,:,1,:), ng_sn, &
                                 sfpx(:,:,1,:), sfpy(:,:,1,:), ng_sf, &
                                 fp(:,:,1,:), ng_f, &
                                 p0(n,:), lo, hi, dx(n,:), dt)
          case (3)
             sfpy => dataptr(sflux(n,2),i)
             sfpz => dataptr(sflux(n,3),i)
             if (spherical .eq. 0) then
                call update_scal_3d_cart(nstart, nstop, &
                                         sop(:,:,:,:), ng_so, snp(:,:,:,:), ng_sn, &
                                         sfpx(:,:,:,:), sfpy(:,:,:,:), sfpz(:,:,:,:), &
                                         ng_sf, fp(:,:,:,:), ng_f, &
                                         p0(n,:), lo, hi, dx(n,:), dt)
             else
                p0np => dataptr(p0_new_cart(n), i)
                call update_scal_3d_sphr(nstart, nstop, &
                                         sop(:,:,:,:), ng_so, snp(:,:,:,:), ng_sn, &
                                         sfpx(:,:,:,:), sfpy(:,:,:,:), sfpz(:,:,:,:), &
                                         ng_sf, fp(:,:,:,:), ng_f, &
                                         p0np(:,:,:,1), ng_p, lo, hi, dx(n,:), dt)
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
          call multifab_fill_ghost_cells(snew(n),snew(n-1),ng_sn,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         nstart,dm+nstart,nstop-nstart+1, &
                                         fill_crse_input=.false.)

          ! do the same for density if we updated the species
          if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

             call ml_cc_restriction_c(snew(n-1),rho_comp,snew(n),rho_comp, &
                                      mla%mba%rr(n-1,:),1)

             call multifab_fill_ghost_cells(snew(n),snew(n-1),ng_sn,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1),the_bc_level(n), &
                                            rho_comp,dm+rho_comp,1,fill_crse_input=.false.)

          endif

       enddo

    end if

    call destroy(bpt)

  end subroutine update_scal

  subroutine update_scal_1d(nstart,nstop,sold,ng_so,snew,ng_sn,sfluxx,ng_sf, &
                            force,ng_f,p0,lo,hi,dx,dt)

    use network,       only: nspec
    use eos_module
    use probin_module, only: enthalpy_pred_type, do_eos_h_above_cutoff, &
         base_cutoff_density, prob_lo
    use variables,     only: spec_comp, rho_comp, rhoh_comp, trac_comp, &
         ntrac, temp_comp
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:)
    integer           , intent(in   ) :: ng_so, ng_sn, ng_sf, ng_f
    real (kind = dp_t), intent(in   ) ::   sold(lo(1)-ng_so:,:)
    real (kind = dp_t), intent(  out) ::   snew(lo(1)-ng_sn:,:)
    real (kind = dp_t), intent(in   ) :: sfluxx(lo(1)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  force(lo(1)-ng_f :,:)
    real (kind = dp_t), intent(in   ) ::        p0(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)

    logical            :: has_negative_species
    integer            :: i, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: delta,frac,sumX
    real (kind = dp_t) :: x
    
    do comp = nstart, nstop

       do i=lo(1),hi(1)
             
          divterm = (sfluxx(i+1,comp) - sfluxx(i,comp))/dx(1)

          snew(i,comp) = sold(i,comp) + dt*(-divterm + force(i,comp))
          
       end do

    enddo

    if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then
       
       do i = lo(1), hi(1)
             
          if (snew(i,rho_comp) .le. base_cutoff_density) then
             den_eos(1) = snew(i,rho_comp)
             temp_eos(1) = sold(i,temp_comp)
             p_eos(1) = p0(i)
             xn_eos(1,:) = snew(i,spec_comp:spec_comp+nspec-1)/den_eos(1)
             
             pt_index_eos(:) = (/i, -1, -1/)
             
             x = prob_lo(1) + (dble(i)+HALF) * dx(1)

             ! (rho,P) --> T,h
             call eos(x, eos_input_rp, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)
             
             snew(i,rhoh_comp) = snew(i,rho_comp) * h_eos(1)
             
          end if
          
       enddo
       
    end if
    
    ! update density
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       
       snew(:,rho_comp) = sold(:,rho_comp)
       
       do i = lo(1), hi(1)

          has_negative_species = .false.

          ! define the update to rho as the sum of the updates to (rho X)_i
          do comp = nstart, nstop
             snew(i,rho_comp) = snew(i,rho_comp) + (snew(i,comp)-sold(i,comp))
             if (snew(i,comp) .lt. ZERO) has_negative_species = .true.
          enddo

          ! enforce a density floor
          if (snew(i,rho_comp) .lt. 0.5d0*base_cutoff_density) then
             do comp = nstart, nstop
                snew(i,comp) = snew(i,comp) * 0.5d0*base_cutoff_density/snew(i,rho_comp)
             end do
             snew(i,rho_comp) = 0.5d0*base_cutoff_density
          end if

          ! do not allow the species to leave here negative.
          if (has_negative_species) then
             do comp = nstart, nstop
                if (snew(i,comp) .lt. ZERO) then
                   delta = -snew(i,comp)
                   sumX = ZERO 
                   do comp2 = nstart, nstop
                      if (comp2 .ne. comp .and. snew(i,comp2) .ge. ZERO) then
                         sumX = sumX + snew(i,comp2)
                      end if
                   enddo
                   do comp2 = nstart, nstop
                      if (comp2 .ne. comp .and. snew(i,comp2) .ge. ZERO) then
                         frac = snew(i,comp2) / sumX
                         snew(i,comp2) = snew(i,comp2) - frac * delta
                      end if
                   enddo
                   snew(i,comp) = ZERO
                end if
             end do
          end if

       enddo
       
    end if

  end subroutine update_scal_1d

  subroutine update_scal_2d(nstart,nstop,sold,ng_so,snew,ng_sn,sfluxx,sfluxy,ng_sf, &
                            force,ng_f,p0,lo,hi,dx,dt)

    use network,       only: nspec
    use eos_module
    use probin_module, only: enthalpy_pred_type, do_eos_h_above_cutoff, &
         base_cutoff_density, prob_lo
    use variables,     only: spec_comp, rho_comp, rhoh_comp, trac_comp, &
         ntrac, temp_comp
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:)
    integer           , intent(in   ) :: ng_so, ng_sn, ng_sf, ng_f
    real (kind = dp_t), intent(in   ) ::   sold(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real (kind = dp_t), intent(  out) ::   snew(lo(1)-ng_sn:,lo(2)-ng_sn:,:)
    real (kind = dp_t), intent(in   ) :: sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) :: sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real (kind = dp_t), intent(in   ) ::        p0(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)

    integer            :: i, j, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: y
    real (kind = dp_t) :: delta,frac,sumX
    logical            :: has_negative_species

    do comp = nstart, nstop

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             
             divterm = (sfluxx(i+1,j,comp) - sfluxx(i,j,comp))/dx(1) &
                     + (sfluxy(i,j+1,comp) - sfluxy(i,j,comp))/dx(2)

             snew(i,j,comp) = sold(i,j,comp) + dt*(-divterm + force(i,j,comp))
             
          end do
       end do

    enddo

    if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then
       
       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF) * dx(2)

          do i = lo(1), hi(1)

             if (snew(i,j,rho_comp) .le. base_cutoff_density) then
                den_eos(1) = snew(i,j,rho_comp)
                temp_eos(1) = sold(i,j,temp_comp)
                p_eos(1) = p0(j)
                xn_eos(1,:) = snew(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

                pt_index_eos(:) = (/i, j, -1/)
                
                ! (rho,P) --> T,h
                call eos(y, eos_input_rp, den_eos, temp_eos, &
                         npts, &
                         xn_eos, &
                         p_eos, h_eos, e_eos, &
                         cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                         dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                         dpdX_eos, dhdX_eos, &
                         gam1_eos, cs_eos, s_eos, &
                         dsdt_eos, dsdr_eos, &
                         .false., &
                         pt_index_eos)
                
                snew(i,j,rhoh_comp) = snew(i,j,rho_comp) * h_eos(1)
                
             end if
             
          enddo
       enddo
       
    end if
    
    ! update density
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       
       snew(:,:,rho_comp) = sold(:,:,rho_comp)
           
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             has_negative_species = .false.

             ! define the update to rho as the sum of the updates to (rho X)_i  
             do comp = nstart, nstop
                snew(i,j,rho_comp) = snew(i,j,rho_comp) + (snew(i,j,comp)-sold(i,j,comp))
                if (snew(i,j,comp) .lt. ZERO) has_negative_species = .true.
             enddo

             ! enforce a density floor
             if (snew(i,j,rho_comp) .lt. 0.5d0*base_cutoff_density) then
                do comp = nstart, nstop
                   snew(i,j,comp) = snew(i,j,comp) * &
                        0.5d0*base_cutoff_density/snew(i,j,rho_comp)
                end do
                snew(i,j,rho_comp) = 0.5d0*base_cutoff_density
             end if

             ! do not allow the species to leave here negative.
             if (has_negative_species) then
                do comp = nstart, nstop
                   if (snew(i,j,comp) .lt. ZERO) then
                      delta = -snew(i,j,comp)
                      sumX = ZERO 
                      do comp2 = nstart, nstop
                         if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                            sumX = sumX + snew(i,j,comp2)
                         end if
                      enddo
                      do comp2 = nstart, nstop
                         if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                            frac = snew(i,j,comp2) / sumX
                            snew(i,j,comp2) = snew(i,j,comp2) - frac * delta
                         end if
                      enddo
                      snew(i,j,comp) = ZERO
                   end if
                end do
             end if

          enddo
       enddo

    end if

  end subroutine update_scal_2d

  subroutine update_scal_3d_cart(nstart,nstop,sold,ng_so,snew,ng_sn,sfluxx,sfluxy,sfluxz, &
                                 ng_sf,force,ng_f,p0,lo,hi,dx,dt)

    use network,       only: nspec
    use eos_module
    use probin_module, only: enthalpy_pred_type, do_eos_h_above_cutoff, base_cutoff_density, prob_lo
    use variables,     only: spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac, temp_comp
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:)
    integer           , intent(in   ) :: ng_so, ng_sn, ng_sf, ng_f
    real (kind = dp_t), intent(in   ) ::   sold(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real (kind = dp_t), intent(  out) ::   snew(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:,:)
    real (kind = dp_t), intent(in   ) :: sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) :: sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) :: sfluxz(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real (kind = dp_t), intent(in   ) ::        p0(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)

    integer            :: i, j, k, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: z
    real (kind = dp_t) :: delta,frac,sumX
    logical            :: has_negative_species

    do comp = nstart, nstop

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                
                divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                        + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                        + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)
   
                snew(i,j,k,comp) = sold(i,j,k,comp) + dt * (-divterm + force(i,j,k,comp))
                
             enddo
          enddo
       enddo

    end do
    
    
    if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then

       do k = lo(3), hi(3)
          z = prob_lo(3) + (dble(k)+HALF) * dx(3)

          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (snew(i,j,k,rho_comp) .le. base_cutoff_density) then
                   den_eos(1) = snew(i,j,k,rho_comp)
                   temp_eos(1) = sold(i,j,k,temp_comp)
                   p_eos(1) = p0(k)
                   xn_eos(1,:) = snew(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

                   pt_index_eos(:) = (/i, j, k/)

                   ! (rho,P) --> T,h
                   call eos(z, eos_input_rp, den_eos, temp_eos, &
                            npts, &
                            xn_eos, &
                            p_eos, h_eos, e_eos, &
                            cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                            dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                            dpdX_eos, dhdX_eos, &
                            gam1_eos, cs_eos, s_eos, &
                            dsdt_eos, dsdr_eos, &
                            .false., &
                            pt_index_eos)

                   snew(i,j,k,rhoh_comp) = snew(i,j,k,rho_comp) * h_eos(1)

                end if

             enddo
          enddo
       enddo

    end if


    ! update density
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

       snew(:,:,:,rho_comp) = sold(:,:,:,rho_comp)
       
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                has_negative_species = .false.

                ! define the update to rho as the sum of the updates to (rho X)_i
                do comp = nstart, nstop
                   snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) &
                        + (snew(i,j,k,comp)-sold(i,j,k,comp))
                   if (snew(i,j,k,comp) .lt. ZERO) has_negative_species = .true.
                enddo

                ! enforce a density floor
                if (snew(i,j,k,rho_comp) .lt. 0.5d0*base_cutoff_density) then
                   do comp = nstart, nstop
                      snew(i,j,k,comp) = snew(i,j,k,comp) * &
                           0.5d0*base_cutoff_density/snew(i,j,k,rho_comp)
                   end do
                   snew(i,j,k,rho_comp) = 0.5d0*base_cutoff_density
                end if

                ! do not allow the species to leave here negative.
                if (has_negative_species) then
                   do comp = nstart, nstop
                      if (snew(i,j,k,comp) .lt. ZERO) then
                         delta = -snew(i,j,k,comp)
                         sumX = ZERO 
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                               sumX = sumX + snew(i,j,k,comp2)
                            end if
                         enddo
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                               frac = snew(i,j,k,comp2) / sumX
                               snew(i,j,k,comp2) = snew(i,j,k,comp2) - frac * delta
                            end if
                         enddo
                         snew(i,j,k,comp) = ZERO
                      end if
                   end do
                end if

             enddo
          enddo
       enddo

    end if

  end subroutine update_scal_3d_cart

  subroutine update_scal_3d_sphr(nstart,nstop,sold,ng_so,snew,ng_sn,sfluxx,sfluxy,sfluxz, &
                                 ng_sf,force,ng_f,p0_new_cart,ng_p,lo,hi,dx,dt)

    use network,       only: nspec
    use eos_module
    use probin_module, only: enthalpy_pred_type, do_eos_h_above_cutoff, base_cutoff_density, prob_lo
    use variables,     only: spec_comp, rho_comp, rhoh_comp, trac_comp, ntrac, temp_comp
    use pred_parameters
    use bl_constants_module
    use geometry,      only: center

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:)
    integer           , intent(in   ) :: ng_so, ng_sn, ng_sf, ng_f, ng_p
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxz(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real (kind = dp_t), intent(in   )::   p0_new_cart(lo(1)-ng_p :,lo(2)-ng_p :,lo(3)-ng_p :)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)

    integer            :: i, j, k, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: delta,frac,sumX
    logical            :: has_negative_species
    real (kind = dp_t) :: r0, x, y, z

    ! is spherical
    do comp = nstart, nstop

       !$OMP PARALLEL DO PRIVATE(i,j,k,divterm)       
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                        + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                        + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)

                snew(i,j,k,comp) = sold(i,j,k,comp) + dt * (-divterm + force(i,j,k,comp))

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    end do

    if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,r0)
       do k = lo(3), hi(3) 
          z = prob_lo(3) + (dble(k)+HALF) * dx(3)

          do j = lo(2), hi(2)
             y = prob_lo(2) + (dble(j)+HALF) * dx(2)
             
             do i = lo(1), hi(1)
                x = prob_lo(1) + (dble(i)+HALF) * dx(1)

                if (snew(i,j,k,rho_comp) .le. base_cutoff_density) then
                   den_eos(1) = snew(i,j,k,rho_comp)
                   temp_eos(1) = sold(i,j,k,temp_comp)
                   p_eos(1) = p0_new_cart(i,j,k)
                   xn_eos(1,:) = snew(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

                   pt_index_eos(:) = (/i, j, k/)
                   r0 = sqrt( (x-center(1))**2 + (y-center(2))**2 + &
                        (z-center(3))**2 ) 

                   ! (rho,P) --> T,h
                   call eos(r0, eos_input_rp, den_eos, temp_eos, &
                            npts, &
                            xn_eos, &
                            p_eos, h_eos, e_eos, &
                            cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                            dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                            dpdX_eos, dhdX_eos, &
                            gam1_eos, cs_eos, s_eos, &
                            dsdt_eos, dsdr_eos, &
                            .false., &
                            pt_index_eos)

                   snew(i,j,k,rhoh_comp) = snew(i,j,k,rho_comp) * h_eos(1)

                end if

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    end if

    ! update density
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

       snew(:,:,:,rho_comp) = sold(:,:,:,rho_comp)

       !$OMP PARALLEL DO PRIVATE(i,j,k,has_negative_species,comp,delta,sumX,comp2,frac)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                has_negative_species = .false.

                ! define the update to rho as the sum of the updates to (rho X)_i
                do comp = nstart, nstop
                   snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) &
                        + (snew(i,j,k,comp)-sold(i,j,k,comp))
                   if (snew(i,j,k,comp) .lt. ZERO) has_negative_species = .true.
                enddo

                ! enforce a density floor
                if (snew(i,j,k,rho_comp) .lt. 0.5d0*base_cutoff_density) then
                   do comp = nstart, nstop
                      snew(i,j,k,comp) = snew(i,j,k,comp) * &
                           0.5d0*base_cutoff_density/snew(i,j,k,rho_comp)
                   end do
                   snew(i,j,k,rho_comp) = 0.5d0*base_cutoff_density
                end if

                ! do not allow the species to leave here negative.
                if (has_negative_species) then
                   do comp = nstart, nstop
                      if (snew(i,j,k,comp) .lt. ZERO) then
                         delta = -snew(i,j,k,comp)
                         sumX = ZERO 
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                               sumX = sumX + snew(i,j,k,comp2)
                            end if
                         enddo
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                               frac = snew(i,j,k,comp2) / sumX
                               snew(i,j,k,comp2) = snew(i,j,k,comp2) - frac * delta
                            end if
                         enddo
                         snew(i,j,k,comp) = ZERO
                      end if
                   end do
                end if

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

  end subroutine update_scal_3d_sphr

end module update_scal_module
