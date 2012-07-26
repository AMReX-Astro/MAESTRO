module react_state_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: react_state, instantaneous_reaction_rates

contains

  subroutine react_state(mla,tempbar,sold,snew,rho_Hext,p0,dt,dx,sdc_source,the_bc_level)

    use probin_module, only: use_tfromp, do_heating, do_burning, drive_initial_convection
    use variables, only: temp_comp, rhoh_comp, rho_comp,nscal

    use multifab_fill_ghost_module
    use multifab_physbc_module, only : multifab_physbc
    use ml_restriction_module , only : ml_cc_restriction
    use heating_module        , only : get_rho_Hext 
    use rhoh_vs_t_module      , only : makeTfromRhoP, makeTfromRhoH
    use bl_constants_module   , only: ZERO

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(dp_t)     , intent(in   ) :: p0(:,0:)
    real(dp_t)     , intent(in   ) :: tempbar(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:,:)
    type(multifab) , intent(inout) :: sdc_source(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local
    type(bl_prof_timer), save :: bpt

    integer :: n,nlevs,dm

    call build(bpt, "react_state")

    if (drive_initial_convection) then
       call bl_error("react_state.f90:drive_initial_convection not supported")
    end if

    nlevs = mla%nlevel
    dm = mla%dim

    ! apply heating term
    if (do_heating) then

       call get_rho_Hext(mla,tempbar,sold,rho_Hext,the_bc_level,dx,dt)

       if (do_burning) then

          ! add to sdc_source for enthalpy
          do n=1,nlevs
             call multifab_plus_plus_c(sdc_source(n),rhoh_comp,rho_Hext(n),1,1,0)
          end do
          
       else

          ! if we aren't burning, then we should just copy the old state to the
          ! new and only update the rhoh component with the heating term
          do n = 1, nlevs
             call multifab_copy(snew(n),sold(n),nghost(sold(n)))
             call multifab_mult_mult_s(rho_Hext(n),dt)
             call multifab_plus_plus_c(snew(n),rhoh_comp,rho_Hext(n),1,1)
             call multifab_div_div_s(rho_Hext(n),dt)
          enddo

       endif

    else

       ! no heating, so we ZERO rho_Hext
       do n=1,nlevs
          call setval(rho_Hext(n),ZERO,all=.true.)
       enddo

    endif

    ! apply burning term
    if (do_burning) then

       call burner_loop(mla,sold,snew,dx,dt,the_bc_level,sdc_source,p0)

       ! pass temperature through for seeding the temperature update eos call
       do n=1,nlevs
          call multifab_copy_c(snew(n),temp_comp,sold(n),temp_comp,1, &
                               nghost(sold(n)))
       end do

    end if

    ! if we aren't doing any heating/burning, then just copy the old to the new
    if (.not. (do_heating .or. do_burning)) then
       do n = 1, nlevs
          call multifab_copy(snew(n),sold(n),nghost(sold(n)))
       enddo
    endif

    ! let's fill the boundaries
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(snew(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(snew(nlevs),rho_comp,dm+rho_comp,nscal,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(snew(n-1),snew(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(snew(n),snew(n-1),nghost(snew(n)), &
                                         mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         rho_comp,dm+rho_comp,nscal,fill_crse_input=.false.)
       enddo

    end if

    ! now update temperature
    if (use_tfromp) then
       call makeTfromRhoP(snew,p0,mla,the_bc_level,dx)
    else
       call makeTfromRhoH(snew,p0,mla,the_bc_level,dx)
    end if

    call destroy(bpt)

  end subroutine react_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop(mla,sold,snew,dx,dt,the_bc_level,sdc_source,p0)

    use bl_constants_module, only: ZERO
    use variables, only: foextrap_comp
    use network, only: nspec
    use geometry, only: spherical
    use fill_3d_module, only: put_1d_array_on_cart

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(in   ) :: sdc_source(:)
    real(dp_t)     , intent(in   ) :: p0(:,0:)

    ! Local
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer ::  rp(:,:,:,:)
    real(kind=dp_t), pointer ::  hnp(:,:,:,:)
    real(kind=dp_t), pointer ::  hep(:,:,:,:)
    real(kind=dp_t), pointer ::  p0p(:,:,:,:)
    real(kind=dp_t), pointer ::  sdcp(:,:,:,:)
    logical        , pointer ::   mp(:,:,:,:)

    type(multifab) :: p0_cart(mla%nlevel)

    integer :: lo(mla%dim),hi(mla%dim)
    integer :: ng_si,ng_so,ng_rw,ng_hn,ng_p0,ng_sd
    integer :: dm,nlevs
    integer :: i,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "burner_loop")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_si = nghost(sold(1))
    ng_so = nghost(snew(1))
    ng_sd = nghost(sdc_source(1))

    if (spherical == 1) then
       do n=1, nlevs
          call build(p0_cart(n),mla%la(n),1,0)
          call setval(p0_cart(n),ZERO,all=.true.)
       enddo

       call put_1d_array_on_cart(p0,p0_cart, &
                                 foextrap_comp,.false.,.false.,dx, &
                                 the_bc_level,mla)
    endif

    do n = 1, nlevs
       do i = 1, nboxes(sold(n))
          if ( multifab_remote(sold(n), i) ) cycle
          snp => dataptr(sold(n) , i)
          sop => dataptr(snew(n), i)
          sdcp => dataptr(sdc_source(n), i)
          lo =  lwb(get_box(sold(n), i))
          hi =  upb(get_box(sold(n), i))
          select case (dm)
          case (1)
             call burner_loop_1d(snp(:,1,1,:),ng_si,sop(:,1,1,:),ng_so, &
                                 sdcp(:,1,1,:),ng_sd,p0(n,:),dt,lo,hi)
          case (2)
             if (n .eq. nlevs) then
                call burner_loop_2d(snp(:,:,1,:),ng_si,sop(:,:,1,:),ng_so, &
                                    sdcp(:,:,1,:),ng_sd,p0(n,:),dt,lo,hi)
             else
                mp => dataptr(mla%mask(n), i)
                call burner_loop_2d(snp(:,:,1,:),ng_si,sop(:,:,1,:),ng_so, &
                                    sdcp(:,:,1,:),ng_sd,p0(n,:), &
                                    dt,lo,hi,mp(:,:,1,1))
             end if
          case (3)
             if (spherical == 1) then
                p0p => dataptr(p0_cart(n), i)
                ng_p0 = nghost(p0_cart(n))
                if (n .eq. nlevs) then
                   call burner_loop_3d_sph(snp(:,:,:,:),ng_si,sop(:,:,:,:),ng_so, &
                                           sdcp(:,:,:,:),ng_sd,p0p(:,:,:,1),ng_p0, &
                                           dt,lo,hi)
                else
                   mp => dataptr(mla%mask(n), i)
                   call burner_loop_3d_sph(snp(:,:,:,:),ng_si,sop(:,:,:,:),ng_so, &
                                           sdcp(:,:,:,:),ng_sd,p0p(:,:,:,1),ng_p0, &
                                           dt,lo,hi,mp(:,:,:,1))
                end if
             else
                if (n .eq. nlevs) then
                   call burner_loop_3d(snp(:,:,:,:),ng_si,sop(:,:,:,:),ng_so, &
                                       sdcp(:,:,:,:),ng_sd,p0(n,:),dt,lo,hi)
                else
                   mp => dataptr(mla%mask(n), i)
                   call burner_loop_3d(snp(:,:,:,:),ng_si,sop(:,:,:,:),ng_so, &
                                       sdcp(:,:,:,:),ng_sd,p0(n,:), &
                                       dt,lo,hi,mp(:,:,:,1))
                end if
             endif
          end select
       end do
    end do

    if (spherical .eq. 1) then
       do n = 1, nlevs
          call destroy(p0_cart(n))
       enddo
    endif

    call destroy(bpt)

  end subroutine burner_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_1d(sold,ng_si,snew,ng_so,sdc_source,ng_sd,p0,dt,lo,hi)

    use bl_constants_module
    use burner_module
    use variables, only: rho_comp, spec_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, network_species_index
    use probin_module, ONLY: burning_cutoff_density, burner_threshold_species, &
         burner_threshold_cutoff

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_sd
    real(kind=dp_t), intent(in   ) ::         sold(lo(1)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,:)
    real(kind=dp_t), intent(in   ) ::   sdc_source(lo(1)-ng_sd:,:)
    real(dp_t)     , intent(in   ) :: p0(0:)
    real(kind=dp_t), intent(in   ) :: dt

    !     Local variables
    integer            :: i
    real (kind = dp_t) :: rho_in,rho_out,rhoh_in,rhoh_out
    real (kind = dp_t) :: rhox_in(nspec)
    real (kind = dp_t) :: rhox_out(nspec)
    real (kind = dp_t) :: x_test
    integer, save      :: ispec_threshold
    logical, save      :: firstCall = .true.

    real (kind = dp_t) :: sdc_rhoX(nspec)
    real (kind = dp_t) :: sdc_rhoh
    real (kind = dp_t) :: p0_in

    if (firstCall) then
       ispec_threshold = network_species_index(burner_threshold_species)
       firstCall = .false.
    endif

    do i=lo(1),hi(1)
       
       sdc_rhoX(:) = sdc_source(i,spec_comp:spec_comp+nspec-1)
       sdc_rhoh = sdc_source(i,rhoh_comp)

       p0_in = p0(i)

       rho_in = sold(i,rho_comp)
       rhox_in(1:nspec) = sold(i,spec_comp:spec_comp+nspec-1)
       rhoh_in = sold(i,rhoh_comp)

       ! Fortran doesn't guarantee short-circuit evaluation of logicals so
       ! we need to test the value of ispec_threshold before using it 
       ! as an index in x_in
       if (ispec_threshold > 0) then
          x_test = rhox_in(ispec_threshold) / rho_in
       else
          x_test = ZERO
       endif

       ! if the threshold species is not in the network, then we burn
       ! normally.  if it is in the network, make sure the mass
       ! fraction is above the cutoff.
       if (rho_in > burning_cutoff_density .and.           &
            ( ispec_threshold < 0 .or.                  &
            (ispec_threshold > 0 .and.                  &
            x_test > burner_threshold_cutoff ))) then
          call burner(rhox_in, rhoh_in, dt, rho_out, rhox_out, rhoh_out, &
                      sdc_rhoX, sdc_rhoh, p0_in)
       else
          rho_out = rho_in + sum(sdc_rhoX(1:nspec))*dt
          rhox_out = rhox_in + sdc_rhoX*dt
          rhoh_out = rhoh_in + sdc_rhoh*dt
       endif

       ! update the density
       snew(i,rho_comp) = rho_out

       ! update the species
       snew(i,spec_comp:spec_comp+nspec-1) = rhox_out(1:nspec)

       ! update the enthalpy -- include the change due to external heating
       snew(i,rhoh_comp) = rhoh_out

       ! pass the tracers through
       snew(i,trac_comp:trac_comp+ntrac-1) = sold(i,trac_comp:trac_comp+ntrac-1)

    enddo

  end subroutine burner_loop_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_2d(sold,ng_si,snew,ng_so,sdc_source,ng_sd,p0,dt,lo,hi,mask)

    use bl_constants_module
    use burner_module
    use variables, only: rho_comp, spec_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, network_species_index
    use probin_module, ONLY: burning_cutoff_density, burner_threshold_species, &
         burner_threshold_cutoff

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_sd
    real(kind=dp_t), intent(in   ) ::         sold(lo(1)-ng_si:,lo(2)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real(kind=dp_t), intent(in   ) ::   sdc_source(lo(1)-ng_sd:,lo(2)-ng_sd:,:)
    real(dp_t)     , intent(in   ) :: p0(0:)
    real(kind=dp_t), intent(in   ) :: dt
    logical        , intent(in   ), optional :: mask(lo(1):,lo(2):)

    !     Local variables
    integer            :: i, j
    real (kind = dp_t) :: rho_in,rho_out,rhoh_in,rhoh_out
    real (kind = dp_t) :: rhox_in(nspec)
    real (kind = dp_t) :: rhox_out(nspec)
    real (kind = dp_t) :: x_test
    logical            :: cell_valid
    integer, save      :: ispec_threshold
    logical, save      :: firstCall = .true.

    real (kind = dp_t) :: sdc_rhoX(nspec)
    real (kind = dp_t) :: sdc_rhoh
    real (kind = dp_t) :: p0_in

    if (firstCall) then
       ispec_threshold = network_species_index(burner_threshold_species)
       firstCall = .false.
    endif

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          
          ! make sure the cell isn't covered by finer cells
          cell_valid = .true.
          if ( present(mask) ) then
             if ( (.not. mask(i,j)) ) cell_valid = .false.
          end if

          if (cell_valid) then

             sdc_rhoX(:) = sdc_source(i,j,spec_comp:spec_comp+nspec-1)
             sdc_rhoh = sdc_source(i,j,rhoh_comp)

             p0_in = p0(j)

             rho_in = sold(i,j,rho_comp)
             rhox_in(1:nspec) = sold(i,j,spec_comp:spec_comp+nspec-1)
             rhoh_in = sold(i,j,rhoh_comp)
          
             ! Fortran doesn't guarantee short-circuit evaluation of logicals so
             ! we need to test the value of ispec_threshold before using it 
             ! as an index in x_in
             if (ispec_threshold > 0) then
                x_test = rhox_in(ispec_threshold) / rho_in
             else
                x_test = ZERO
             endif

             ! if the threshold species is not in the network, then we burn
             ! normally.  if it is in the network, make sure the mass
             ! fraction is above the cutoff.
             if (rho_in > burning_cutoff_density .and.           &
                  ( ispec_threshold < 0 .or.                  &
                  (ispec_threshold > 0 .and.                  &
                  x_test > burner_threshold_cutoff ))) then
                call burner(rhox_in, rhoh_in, dt, rho_out, rhox_out, rhoh_out, &
                            sdc_rhoX, sdc_rhoh, p0_in)
             else
                rho_out = rho_in + sum(sdc_rhoX(1:nspec))*dt
                rhox_out = rhox_in + sdc_rhoX*dt
                rhoh_out = rhoh_in + sdc_rhoh*dt
             endif
             
             ! update the density
             snew(i,j,rho_comp) = rho_out
             
             ! update the species
             snew(i,j,spec_comp:spec_comp+nspec-1) = rhox_out(1:nspec)
             
             ! update the enthalpy -- include the change due to external heating
             snew(i,j,rhoh_comp) = rhoh_out
             
             ! pass the tracers through
             snew(i,j,trac_comp:trac_comp+ntrac-1) = sold(i,j,trac_comp:trac_comp+ntrac-1)   

          end if
          
       enddo
    enddo

  end subroutine burner_loop_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_3d(sold,ng_si,snew,ng_so,sdc_source,ng_sd,p0,dt,lo,hi,mask)

    use bl_constants_module
    use burner_module
    use variables, only: rho_comp, spec_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, network_species_index
    use probin_module, ONLY: burning_cutoff_density, burner_threshold_species, &
         burner_threshold_cutoff

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_sd
    real(kind=dp_t), intent(in   ) ::         sold(lo(1)-ng_si:,lo(2)-ng_si:,lo(3)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real(kind=dp_t), intent(in   ) ::   sdc_source(lo(1)-ng_sd:,lo(2)-ng_sd:,lo(3)-ng_sd:,:)
    real(dp_t)     , intent(in   ) :: p0(0:)
    real(kind=dp_t), intent(in   ) :: dt
    logical        , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: i, j, k
    real (kind = dp_t) :: rho_in,rho_out,rhoh_in,rhoh_out,ldt
    real (kind = dp_t) :: rhox_in(nspec)
    real (kind = dp_t) :: rhox_out(nspec)
    real (kind = dp_t) :: x_test
    logical            :: cell_valid
    integer, save      :: ispec_threshold
    logical, save      :: firstCall = .true.

    real (kind = dp_t) :: sdc_rhoX(nspec)
    real (kind = dp_t) :: sdc_rhoh
    real (kind = dp_t) :: p0_in

    if (firstCall) then
       ispec_threshold = network_species_index(burner_threshold_species)
       firstCall = .false.
    endif

    ldt = dt

    !$OMP PARALLEL DO PRIVATE(i,j,k,cell_valid,sdc_rhoX,sdc_rhoh,p0_in,rho_in,rhox_in) &
    !$OMP PRIVATE(rhoh_in,x_test,rho_out,rhox_out,rhoh_out) &
    !$OMP FIRSTPRIVATE(ldt) &
    !$OMP SCHEDULE(DYNAMIC,1)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if ( present(mask) ) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             end if

             if (cell_valid) then

                sdc_rhoX(:) = sdc_source(i,j,k,spec_comp:spec_comp+nspec-1)
                sdc_rhoh = sdc_source(i,j,k,rhoh_comp)

                p0_in = p0(k)

                rho_in = sold(i,j,k,rho_comp)
                rhox_in(1:nspec) = sold(i,j,k,spec_comp:spec_comp+nspec-1)
                rhoh_in = sold(i,j,k,rhoh_comp)

                ! Fortran doesn't guarantee short-circuit evaluation of logicals so
                ! we need to test the value of ispec_threshold before using it 
                ! as an index in x_in
                if (ispec_threshold > 0) then
                   x_test = rhox_in(ispec_threshold) / rho_in
                else
                   x_test = ZERO
                endif

                ! if the threshold species is not in the network, then we burn
                ! normally.  if it is in the network, make sure the mass
                ! fraction is above the cutoff.
                if (rho_in > burning_cutoff_density .and.           &
                     ( ispec_threshold < 0 .or.                  &
                     (ispec_threshold > 0 .and.                  &
                     x_test > burner_threshold_cutoff ))) then
                   call burner(rhox_in, rhoh_in, dt, rho_out, rhox_out, rhoh_out, &
                               sdc_rhoX, sdc_rhoh, p0_in)
                else
                   rho_out = rho_in + sum(sdc_rhoX(1:nspec))*dt
                   rhox_out = rhox_in + sdc_rhoX*dt
                   rhoh_out = rhoh_in + sdc_rhoh*dt
                endif

                ! update the density
                snew(i,j,k,rho_comp) = rho_out

                ! update the species
                snew(i,j,k,spec_comp:spec_comp+nspec-1) = rhox_out(1:nspec)

                ! update the enthalpy -- include the change due to external heating
                snew(i,j,k,rhoh_comp) = rhoh_out

                ! pass the tracers through
                snew(i,j,k,trac_comp:trac_comp+ntrac-1) = sold(i,j,k,trac_comp:trac_comp+ntrac-1)   

             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine burner_loop_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_3d_sph(sold,ng_si,snew,ng_so, &
                                sdc_source,ng_sd, &
                                p0_cart,ng_p0, &
                                dt,lo,hi,mask)

    use bl_constants_module
    use burner_module
    use variables, only: rho_comp, spec_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, network_species_index
    use probin_module, ONLY: burning_cutoff_density, burner_threshold_species, &
         burner_threshold_cutoff

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_sd,ng_p0
    real(kind=dp_t), intent(in   ) ::         sold(lo(1)-ng_si:,lo(2)-ng_si:,lo(3)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real(kind=dp_t), intent(in   ) ::   sdc_source(lo(1)-ng_sd:,lo(2)-ng_sd:,lo(3)-ng_sd:,:)
    real(kind=dp_t), intent(in   ) ::      p0_cart(lo(1)-ng_p0:,lo(2)-ng_p0:,lo(3)-ng_p0:)
    real(kind=dp_t), intent(in   ) :: dt
    logical        , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

!     Local variables
    integer            :: i, j, k
    real (kind = dp_t) :: rho_in,rho_out,rhoh_in,rhoh_out,ldt
    real (kind = dp_t) :: rhox_in(nspec)
    real (kind = dp_t) :: rhox_out(nspec)
    real (kind = dp_t) :: x_test
    logical            :: cell_valid
    integer, save      :: ispec_threshold
    logical, save      :: firstCall = .true.

    real (kind = dp_t) :: sdc_rhoX(nspec)
    real (kind = dp_t) :: sdc_rhoh
    real (kind = dp_t) :: p0_in

    if (firstCall) then
       ispec_threshold = network_species_index(burner_threshold_species)
       firstCall = .false.
    endif

    ldt = dt

    !$OMP PARALLEL DO PRIVATE(i,j,k,cell_valid,sdc_rhoX,sdc_rhoh,p0_in,rho_in,rhox_in) &
    !$OMP PRIVATE(rhoh_in,x_test,rho_out,rhox_out,rhoh_out) &
    !$OMP FIRSTPRIVATE(ldt) &
    !$OMP SCHEDULE(DYNAMIC,1)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if ( present(mask) ) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             end if

             if (cell_valid) then

                sdc_rhoX(:) = sdc_source(i,j,k,spec_comp:spec_comp+nspec-1)
                sdc_rhoh = sdc_source(i,j,k,rhoh_comp)

                p0_in = p0_cart(i,j,k)

                rho_in = sold(i,j,k,rho_comp)
                rhox_in(1:nspec) = sold(i,j,k,spec_comp:spec_comp+nspec-1)
                rhoh_in = sold(i,j,k,rhoh_comp)

                ! Fortran doesn't guarantee short-circuit evaluation of logicals so
                ! we need to test the value of ispec_threshold before using it 
                ! as an index in x_in
                if (ispec_threshold > 0) then
                   x_test = rhox_in(ispec_threshold) / rho_in
                else
                   x_test = ZERO
                endif

                ! if the threshold species is not in the network, then we burn
                ! normally.  if it is in the network, make sure the mass
                ! fraction is above the cutoff.
                if (rho_in > burning_cutoff_density .and.           &
                     ( ispec_threshold < 0 .or.                  &
                     (ispec_threshold > 0 .and.                  &
                     x_test > burner_threshold_cutoff ))) then
                   call burner(rhox_in, rhoh_in, dt, rho_out, rhox_out, rhoh_out, &
                               sdc_rhoX, sdc_rhoh, p0_in)
                else
                   rho_out = rho_in + sum(sdc_rhoX(1:nspec))*dt
                   rhox_out = rhox_in + sdc_rhoX*dt
                   rhoh_out = rhoh_in + sdc_rhoh*dt
                endif

                ! update the density
                snew(i,j,k,rho_comp) = rho_out

                ! update the species
                snew(i,j,k,spec_comp:spec_comp+nspec-1) = rhox_out(1:nspec)

                ! update the enthalpy -- include the change due to external heating
                snew(i,j,k,rhoh_comp) = rhoh_out

                ! pass the tracers through
                snew(i,j,k,trac_comp:trac_comp+ntrac-1) = sold(i,j,k,trac_comp:trac_comp+ntrac-1)   

             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine burner_loop_3d_sph

  subroutine instantaneous_reaction_rates(mla,s,rho_omegadot,rho_Hnuc)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_omegadot(:)
    type(multifab) , intent(inout) :: rho_Hnuc(:)

    ! local
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: wp(:,:,:,:)
    real(kind=dp_t), pointer :: hp(:,:,:,:)

    integer :: lo(mla%dim),hi(mla%dim)
    integer :: dm,nlevs,i,n
    integer :: ng_s, ng_w, ng_h

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s = s(1)%ng
    ng_w = rho_omegadot(1)%ng
    ng_h = rho_Hnuc(1)%ng

    do n=1,nlevs
       do i=1,nboxes(s(n))
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n), i)
          wp => dataptr(rho_omegadot(n), i)
          hp => dataptr(rho_Hnuc(n), i)
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
          select case(dm)
          case (1)
             call instantaneous_reaction_rates_1d(sp(:,1,1,:),ng_s,wp(:,1,1,:),ng_w, &
                                                  hp(:,1,1,1),ng_h,lo,hi)

          case (2)
             call instantaneous_reaction_rates_2d(sp(:,:,1,:),ng_s,wp(:,:,1,:),ng_w, &
                                                  hp(:,:,1,1),ng_h,lo,hi)
          case (3)
             call bl_error("instantaneous_reaction_rates_3d not written yet")
          end select
       end do
    end do

  end subroutine instantaneous_reaction_rates

  subroutine instantaneous_reaction_rates_1d(s,ng_s,rho_omegadot,ng_w,rho_Hnuc,ng_h,lo,hi)

    use network, only: nspec
    use variables, only: spec_comp, rhoh_comp

    integer,         intent(in   ) :: ng_s,ng_w,ng_h,lo(:),hi(:)
    real(kind=dp_t), intent(in   ) ::            s(lo(1)-ng_s:,:)
    real(kind=dp_t), intent(inout) :: rho_omegadot(lo(1)-ng_w:,:)
    real(kind=dp_t), intent(inout) ::     rho_Hnuc(lo(1)-ng_h:)

    ! local
    integer :: i

    real(kind=dp_t) :: y(nspec+1)
    real(kind=dp_t) :: ydot(nspec+1)
    real(kind=dp_t) :: rho_Hnuc_out

    real(kind=dp_t) :: time_dummy, rpar_dummy
    integer :: ipar_dummy

    do i=lo(1),hi(1)

       y(1:nspec) = s(i,spec_comp:spec_comp+nspec-1)
       y(nspec+1) = s(i,rhoh_comp)

       call f_rhs_instantaneous_reaction_rates(nspec+1,time_dummy,y,ydot,rho_Hnuc_out,&
                                               rpar_dummy,ipar_dummy)

       rho_omegadot(i,1:nspec) = ydot(1:nspec)
       rho_Hnuc(i) = rho_Hnuc_out

    end do

  end subroutine instantaneous_reaction_rates_1d

  subroutine instantaneous_reaction_rates_2d(s,ng_s,rho_omegadot,ng_w,rho_Hnuc,ng_h,lo,hi)

    use network, only: nspec
    use variables, only: spec_comp, rhoh_comp

    integer,         intent(in   ) :: ng_s,ng_w,ng_h,lo(:),hi(:)
    real(kind=dp_t), intent(in   ) ::            s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(inout) :: rho_omegadot(lo(1)-ng_w:,lo(2)-ng_w:,:)
    real(kind=dp_t), intent(inout) ::     rho_Hnuc(lo(1)-ng_h:,lo(2)-ng_h:)

    ! local
    integer :: i,j

    real(kind=dp_t) :: y(nspec+1)
    real(kind=dp_t) :: ydot(nspec+1)
    real(kind=dp_t) :: rho_Hnuc_out

    real(kind=dp_t) :: time_dummy, rpar_dummy
    integer :: ipar_dummy

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          y(1:nspec) = s(i,j,spec_comp:spec_comp+nspec-1)
          y(nspec+1) = s(i,j,rhoh_comp)

          call f_rhs_instantaneous_reaction_rates(nspec+1,time_dummy,y,ydot,rho_Hnuc_out, &
                                                  rpar_dummy,ipar_dummy)

          rho_omegadot(i,j,1:nspec) = ydot(1:nspec)
          rho_Hnuc(i,j) = rho_Hnuc_out

       end do
    end do

  end subroutine instantaneous_reaction_rates_2d

end module react_state_module
