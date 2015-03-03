module react_state_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module
  use bl_constants_module

  implicit none

  private

  public :: react_state

contains

  function slope(q, qm, qp) result (dq)
    
    real (kind=dp_t), intent(in) :: q, qm, qp
    real (kind=dp_t) :: dq

    real (kind=dp_t) :: test

    test = (qp - q)*(q - qm)
    if (test > ZERO) then
       dq = min(HALF*abs(qp - qm), min(TWO*abs(qp-q), TWO*abs(q-qm)))* sign(ONE, qp-qm)
    else
       dq = ZERO
    endif

  end function slope

  function ppm_fit(q, qm2, qm1, qp1, qp2) result (qcoeff)

    real (kind=dp_t), intent(in) :: q, qm2, qm1, qp1, qp2
    real (kind=dp_t) :: qcoeff(3)

    ! qcoeff(1) is the ppm ql, qcoeff(2) is the ppm qr, and
    ! qcoeff(3) is the ppm q6
    ! no limiting is done!!

    ! find ql and qr
    qcoeff(1) = (7.0_dp_t/12.0_dp_t)*(q + qm1) - (1.0_dp_t/12.0_dp_t)*(qm2 + qp1)
    qcoeff(2) = (7.0_dp_t/12.0_dp_t)*(q + qp1) - (1.0_dp_t/12.0_dp_t)*(qm1 + qp2)

    ! construct q6
    qcoeff(3) = 6.0_dp_t*(q - 0.5_dp_t*(qcoeff(1) + qcoeff(2)))

  end function ppm_fit

  function ppm_eval(xi, ql, qr, q6) result (qval)
    
    real (kind=dp_t), intent(in) :: xi, ql, qr, q6
    real (kind=dp_t) :: qval

    qval = ql + xi*(qr - ql + q6*(1.0_dp_t - xi))

  end function ppm_eval
    

  subroutine react_state(mla,tempbar_init,sold,snew,rho_omegadot,rho_Hnuc,rho_Hext,p0, &
                         dt,dx,the_bc_level)

    use probin_module, only: use_tfromp, do_heating, do_burning
    use variables, only: temp_comp, rhoh_comp, rho_comp,nscal

    use multifab_fill_ghost_module
    use multifab_physbc_module, only : multifab_physbc
    use ml_cc_restriction_module , only : ml_cc_restriction
    use heating_module        , only : get_rho_Hext 
    use rhoh_vs_t_module      , only : makeTfromRhoP, makeTfromRhoH
    use bl_constants_module   , only: ZERO

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: rho_omegadot(:)
    type(multifab) , intent(inout) :: rho_Hnuc(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(dp_t)     , intent(in   ) :: p0(:,0:)
    real(dp_t)     , intent(in   ) :: tempbar_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local
    type(bl_prof_timer), save :: bpt

    integer :: n,nlevs,dm

    call build(bpt, "react_state")

    nlevs = mla%nlevel
    dm = mla%dim

    ! apply heating term
    if(do_heating) then
       call get_rho_Hext(mla,tempbar_init,sold,rho_Hext,the_bc_level,dx,dt)

       ! if we aren't burning, then we should just copy the old state to the
       ! new and only update the rhoh component with the heating term
       if (.not. do_burning) then
          do n = 1, nlevs
             call multifab_copy(snew(n),sold(n), nghost(sold(n)))
             
             ! add in the heating term*dt
             call multifab_mult_mult_s(rho_Hext(n),dt)
             call multifab_plus_plus_c(snew(n),rhoh_comp,rho_Hext(n),1,1)
             call multifab_div_div_s(rho_Hext(n),dt)
          enddo
       endif

    else ! not burning, so we ZERO rho_Hext
       do n = 1, nlevs
          call setval(rho_Hext(n),ZERO,all=.true.)
       enddo

    endif

    ! apply burning term
    if (do_burning) then
       ! we pass in rho_Hext so that we can add it to rhoh incase we 
       ! applied heating
       call burner_loop(mla,tempbar_init,sold,snew,rho_omegadot,rho_Hnuc, &
                        rho_Hext,dx,dt,the_bc_level)

       ! pass temperature through for seeding the temperature update eos call
       do n=1,nlevs
          call multifab_copy_c(snew(n),temp_comp,sold(n),temp_comp,1, &
                               nghost(sold(n)))
       end do

    else ! not burning, so we ZERO rho_omegadot and rho_Hnuc
       do n = 1, nlevs
          call setval(rho_omegadot(n),ZERO,all=.true.)
          call setval(rho_Hnuc(n),ZERO,all=.true.)
       enddo

    endif

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
          call ml_cc_restriction(snew(n-1)        ,snew(n)        ,mla%mba%rr(n-1,:))
          call ml_cc_restriction(rho_omegadot(n-1),rho_omegadot(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(rho_Hext(n-1)    ,rho_Hext(n)    ,mla%mba%rr(n-1,:))
          call ml_cc_restriction(rho_Hnuc(n-1)    ,rho_Hnuc(n)    ,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(snew(n),snew(n-1),nghost(snew(n)),mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         rho_comp,dm+rho_comp,nscal,fill_crse_boundary_input=.false.)
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

  subroutine burner_loop(mla,tempbar_init,sold,snew,rho_omegadot,rho_Hnuc,rho_Hext,dx,dt,the_bc_level)

    use bl_constants_module, only: ZERO
    use variables, only: foextrap_comp
    use network, only: nspec
    use probin_module, only: drive_initial_convection, do_subgrid_burning
    use geometry, only: spherical
    use fill_3d_module, only: put_1d_array_on_cart

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: rho_omegadot(:)
    type(multifab) , intent(inout) :: rho_Hnuc(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer ::  rp(:,:,:,:)
    real(kind=dp_t), pointer ::  hnp(:,:,:,:)
    real(kind=dp_t), pointer ::  hep(:,:,:,:)
    real(kind=dp_t), pointer ::  tcp(:,:,:,:)
    logical        , pointer ::   mp(:,:,:,:)

    type(multifab) :: tempbar_init_cart(mla%nlevel)

    integer :: lo(mla%dim),hi(mla%dim)
    integer :: ng_si,ng_so,ng_rw,ng_he,ng_hn,ng_tc
    integer :: dm,nlevs
    integer :: i,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "burner_loop")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_si = nghost(sold(1))
    ng_so = nghost(snew(1))
    ng_rw = nghost(rho_omegadot(1))
    ng_hn = nghost(rho_Hnuc(1))
    ng_he = nghost(rho_Hext(1))

    ! put tempbar_init on Cart
    if (spherical == 1) then
       do n=1, nlevs
          ! tempbar_init_cart will hold the initial tempbar on a Cartesian
          ! grid to be used if drive_initial_convection is true
          call build(tempbar_init_cart(n),mla%la(n),1,0)
          call setval(tempbar_init_cart(n), ZERO, all=.true.)
       enddo

       if (drive_initial_convection) then
          ! fill all components
          call put_1d_array_on_cart(tempbar_init,tempbar_init_cart, &
                                    foextrap_comp,.false.,.false.,dx, &
                                    the_bc_level,mla)
       endif
    endif

    do n = 1, nlevs
       do i = 1, nfabs(sold(n))
          print *, 'lev, ifab: ', n, i
          snp => dataptr(sold(n) , i)
          sop => dataptr(snew(n), i)
          rp => dataptr(rho_omegadot(n), i)
          hnp => dataptr(rho_Hnuc(n), i)
          hep => dataptr(rho_Hext(n), i)
          lo =  lwb(get_box(sold(n), i))
          hi =  upb(get_box(sold(n), i))
          select case (dm)
          case (1)
             call burner_loop_1d(tempbar_init(n,:), &
                                 snp(:,1,1,:),ng_si,sop(:,1,1,:),ng_so, &
                                 rp(:,1,1,:),ng_rw, &
                                 hnp(:,1,1,1),ng_hn,hep(:,1,1,1),ng_he, &
                                 dt,lo,hi)
          case (2)

             if (do_subgrid_burning) then
                if (n .eq. nlevs) then
                   call burner_loop_2d_sub(tempbar_init(n,:), &
                                           snp(:,:,1,:),ng_si,sop(:,:,1,:),ng_so, &
                                           rp(:,:,1,:),ng_rw, &
                                           hnp(:,:,1,1),ng_hn,hep(:,:,1,1),ng_he, &
                                           dx(n,:),dt,lo,hi)
                else
                   mp => dataptr(mla%mask(n), i)
                   call burner_loop_2d_sub(tempbar_init(n,:), &
                                           snp(:,:,1,:),ng_si,sop(:,:,1,:),ng_so, &
                                           rp(:,:,1,:),ng_rw, &
                                           hnp(:,:,1,1),ng_hn,hep(:,:,1,1),ng_he, &
                                           dx(n,:),dt,lo,hi,mp(:,:,1,1))
                end if

             else

                if (n .eq. nlevs) then
                   call burner_loop_2d(tempbar_init(n,:), &
                                       snp(:,:,1,:),ng_si,sop(:,:,1,:),ng_so, &
                                       rp(:,:,1,:),ng_rw, &
                                       hnp(:,:,1,1),ng_hn,hep(:,:,1,1),ng_he, &
                                       dt,lo,hi)
                else
                   mp => dataptr(mla%mask(n), i)
                   call burner_loop_2d(tempbar_init(n,:), &
                                       snp(:,:,1,:),ng_si,sop(:,:,1,:),ng_so, &
                                       rp(:,:,1,:),ng_rw, &
                                       hnp(:,:,1,1),ng_hn,hep(:,:,1,1),ng_he, &
                                       dt,lo,hi,mp(:,:,1,1))
                end if

             endif

          case (3)
             if (spherical == 1) then
                tcp => dataptr(tempbar_init_cart(n), i)
                ng_tc = nghost(tempbar_init_cart(1))
                if (n .eq. nlevs) then
                   call burner_loop_3d_sph(tcp(:,:,:,1),ng_tc, &
                                           snp(:,:,:,:),ng_si,sop(:,:,:,:),ng_so, &
                                           rp(:,:,:,:),ng_rw, &
                                           hnp(:,:,:,1),ng_hn,hep(:,:,:,1),ng_he, &
                                           dt,lo,hi)
                else
                   mp => dataptr(mla%mask(n), i)
                   call burner_loop_3d_sph(tcp(:,:,:,1),ng_tc, &
                                           snp(:,:,:,:),ng_si,sop(:,:,:,:),ng_so, &
                                           rp(:,:,:,:),ng_rw, &
                                           hnp(:,:,:,1),ng_hn,hep(:,:,:,1),ng_he, &
                                           dt,lo,hi,mp(:,:,:,1))
                end if
             else
                print *, 'Calling 3D vec loop...'
                call burner_loop_3d_vec(tempbar_init(n,:), &
                                    snp(:,:,:,:),ng_si,sop(:,:,:,:),ng_so, &
                                    rp(:,:,:,:),ng_rw, &
                                    hnp(:,:,:,1),ng_hn,hep(:,:,:,1),ng_he, &
                                    dt,lo,hi)
                !if (n .eq. nlevs) then
                !   call burner_loop_3d(tempbar_init(n,:), &
                !                       snp(:,:,:,:),ng_si,sop(:,:,:,:),ng_so, &
                !                       rp(:,:,:,:),ng_rw, &
                !                       hnp(:,:,:,1),ng_hn,hep(:,:,:,1),ng_he, &
                !                       dt,lo,hi)
                !else
                !   mp => dataptr(mla%mask(n), i)
                !   call burner_loop_3d(tempbar_init(n,:), &
                !                       snp(:,:,:,:),ng_si,sop(:,:,:,:),ng_so, &
                !                       rp(:,:,:,:),ng_rw, &
                !                       hnp(:,:,:,1),ng_hn,hep(:,:,:,1),ng_he, &
                !                       dt,lo,hi,mp(:,:,:,1))
                !end if
             endif
          end select
       end do
    end do

    call destroy(bpt)

    if (spherical == 1) then
       do n = 1, nlevs
          call destroy(tempbar_init_cart(n))
       enddo
    endif

  end subroutine burner_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_1d(tempbar_init,sold,ng_si,snew,ng_so,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,dt,lo,hi)

    use bl_constants_module
    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, network_species_index
    use probin_module, ONLY: burning_cutoff_density, burner_threshold_species, &
         burner_threshold_cutoff, drive_initial_convection

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_rw,ng_he,ng_hn
    real(kind=dp_t), intent(in   ) ::        sold (lo(1)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,:)
    real(kind=dp_t), intent(  out) :: rho_omegadot(lo(1)-ng_rw:,:)
    real(kind=dp_t), intent(  out) ::     rho_Hnuc(lo(1)-ng_hn:)
    real(kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(0:)
    real(kind=dp_t), intent(in   ) :: dt

    !     Local variables
    integer            :: i,n
    real (kind = dp_t) :: rho,T_in
    real (kind = dp_t) :: x_in(nspec)
    real (kind = dp_t) :: x_out(nspec)
    real (kind = dp_t) :: rhowdot(nspec)
    real (kind = dp_t) :: rhoH
    real (kind = dp_t) :: x_test
    integer, save      :: ispec_threshold
    logical, save      :: firstCall = .true.

    real (kind = dp_t) :: sumX
    real (kind = dp_t), parameter :: x_err = 1.d-10

    if (firstCall) then
       ispec_threshold = network_species_index(burner_threshold_species)
       firstCall = .false.
    endif

    do i = lo(1), hi(1)
          
          rho = sold(i,rho_comp)
          x_in(1:nspec) = sold(i,spec_comp:spec_comp+nspec-1) / rho

          if (drive_initial_convection) then
             T_in = tempbar_init(i)
          else
             T_in = sold(i,temp_comp)
          endif

          ! Fortran doesn't guarantee short-circuit evaluation of logicals so
          ! we need to test the value of ispec_threshold before using it 
          ! as an index in x_in
          if (ispec_threshold > 0) then
             x_test = x_in(ispec_threshold)
          else
             x_test = ZERO
          endif

          ! if the threshold species is not in the network, then we burn
          ! normally.  if it is in the network, make sure the mass
          ! fraction is above the cutoff.
          if (rho > burning_cutoff_density .and.                      &
              ( ispec_threshold < 0 .or.                              &
               (ispec_threshold > 0 .and.                             &
                x_test > burner_threshold_cutoff                      &
               )                                                      &
              )                                                       &
             ) then
             call burner(rho, T_in, x_in, dt, x_out, rhowdot, rhoH)
          else
             x_out = x_in
             rhowdot = 0.d0
             rhoH = 0.d0
          endif

          ! check if sum{X_k} = 1
          sumX = ZERO
          do n = 1, nspec
             sumX = sumX + x_out(n)
          enddo
          if (abs(sumX - ONE) > x_err) then
             call bl_error("ERROR: abundances do not sum to 1", abs(sumX-ONE))
          endif

          ! pass the density through
          snew(i,rho_comp) = sold(i,rho_comp)

          ! update the species
          snew(i,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho

          ! store the energy generation and species creation quantities
          rho_omegadot(i,1:nspec) = rhowdot(1:nspec)
          rho_Hnuc(i) = rhoH

          ! update the enthalpy -- include the change due to external heating
          snew(i,rhoh_comp) = sold(i,rhoh_comp) + dt*rho_Hnuc(i) + dt*rho_Hext(i)

          ! pass the tracers through
          snew(i,trac_comp:trac_comp+ntrac-1) = sold(i,trac_comp:trac_comp+ntrac-1)   
          
    enddo

  end subroutine burner_loop_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_2d(tempbar_init,sold,ng_si,snew,ng_so,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,dt,lo,hi,mask)

    use bl_constants_module
    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, network_species_index
    use probin_module, ONLY: burning_cutoff_density, burner_threshold_species, &
         burner_threshold_cutoff, drive_initial_convection

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_rw,ng_he,ng_hn
    real(kind=dp_t), intent(in   ) ::        sold (lo(1)-ng_si:,lo(2)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real(kind=dp_t), intent(  out) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,:)
    real(kind=dp_t), intent(  out) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:)
    real(kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(0:)
    real(kind=dp_t), intent(in   ) :: dt
    logical        , intent(in   ), optional :: mask(lo(1):,lo(2):)

    !     Local variables
    integer            :: i, j, n
    real (kind = dp_t) :: rho,T_in
    real (kind = dp_t) :: x_in(nspec)
    real (kind = dp_t) :: x_out(nspec)
    real (kind = dp_t) :: rhowdot(nspec)
    real (kind = dp_t) :: rhoH
    real (kind = dp_t) :: x_test
    logical            :: cell_valid
    integer, save      :: ispec_threshold
    logical, save      :: firstCall = .true.

    real (kind = dp_t) :: sumX
    real (kind = dp_t), parameter :: x_err = 1.d-10

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

             rho = sold(i,j,rho_comp)
             x_in(1:nspec) = sold(i,j,spec_comp:spec_comp+nspec-1) / rho
          
             if (drive_initial_convection) then
                T_in = tempbar_init(j)
             else
                T_in = sold(i,j,temp_comp)
             endif

             ! Fortran doesn't guarantee short-circuit evaluation of logicals so
             ! we need to test the value of ispec_threshold before using it 
             ! as an index in x_in
             if (ispec_threshold > 0) then
                x_test = x_in(ispec_threshold)
             else
                x_test = ZERO
             endif

             ! if the threshold species is not in the network, then we burn
             ! normally.  if it is in the network, make sure the mass
             ! fraction is above the cutoff.
             if (rho > burning_cutoff_density .and.           &
                  ( ispec_threshold < 0 .or.                  &
                  (ispec_threshold > 0 .and.                  &
                  x_test > burner_threshold_cutoff ))) then
                call burner(rho, T_in, x_in, dt, x_out, rhowdot, rhoH)
             else
                x_out = x_in
                rhowdot = 0.d0
                rhoH = 0.d0
             endif
             
             ! check if sum{X_k} = 1
             sumX = ZERO
             do n = 1, nspec
                sumX = sumX + x_out(n)
             enddo
             if (abs(sumX - ONE) > x_err) then
                call bl_error("ERROR: abundances do not sum to 1", abs(sumX-ONE))
             endif

             ! pass the density through
             snew(i,j,rho_comp) = sold(i,j,rho_comp)
             
             ! update the species
             snew(i,j,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho
             
             ! store the energy generation and species creation quantities
             rho_omegadot(i,j,1:nspec) = rhowdot(1:nspec)
             rho_Hnuc(i,j) = rhoH
             
             ! update the enthalpy -- include the change due to external heating
             snew(i,j,rhoh_comp) = sold(i,j,rhoh_comp) + dt*rho_Hnuc(i,j) + dt*rho_Hext(i,j)
             
             ! pass the tracers through
             snew(i,j,trac_comp:trac_comp+ntrac-1) = sold(i,j,trac_comp:trac_comp+ntrac-1)   

          end if
          
       enddo
    enddo

  end subroutine burner_loop_2d


  subroutine burner_loop_2d_sub(tempbar_init,sold,ng_si,snew,ng_so,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,dx,dt,lo,hi,mask)

    ! a special version of the 2d loop that subsamples in the vertical direction and averages
    ! the result -- this gets a better value for the average enuc in a cell.

    use bl_constants_module
    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, network_species_index
    use probin_module, ONLY: burning_cutoff_density, burner_threshold_species, &
         burner_threshold_cutoff, drive_initial_convection

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_rw,ng_he,ng_hn
    real(kind=dp_t), intent(in   ) ::        sold (lo(1)-ng_si:,lo(2)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real(kind=dp_t), intent(  out) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,:)
    real(kind=dp_t), intent(  out) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:)
    real(kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:), dt
    logical        , intent(in   ), optional :: mask(lo(1):,lo(2):)

    !     Local variables
    integer, parameter :: nsub = 4
    integer :: jj

    integer            :: i, j, n
    real (kind = dp_t) :: rho,T_in
    real (kind = dp_t) :: x_in(nspec)
    real (kind = dp_t) :: x_out(nspec)
    real (kind = dp_t) :: rhowdot(nspec)
    real (kind = dp_t) :: rhoH
    real (kind = dp_t) :: x_test
    logical            :: cell_valid
    integer, save      :: ispec_threshold
    logical, save      :: firstCall = .true.

    real (kind = dp_t) :: sumX
    real (kind = dp_t), parameter :: x_err = 1.d-10

    real (kind=dp_t) :: slope_rho, slope_T, slope_X
    real (kind=dp_t) :: Tcoeff(3)
    real (kind=dp_t) :: slope_X_temp
    real (kind=dp_t) :: x_out_temp(nspec), rhowdot_temp(nspec), rhoH_temp

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

             ! density
             slope_rho = slope(sold(i,j,rho_comp), sold(i,j-1,rho_comp), sold(i,j+1,rho_comp))

             ! temp -- do dT and T_0 separately
             if (j > 1 .and. j < size(tempbar_init,dim=1)-2) then
                slope_T = slope(sold(i,j  ,temp_comp)-tempbar_init(j), &
                                sold(i,j-1,temp_comp)-tempbar_init(j-1), &
                                sold(i,j+1,temp_comp)-tempbar_init(j+1))

                Tcoeff = ppm_fit(tempbar_init(j), &
                                 tempbar_init(j-2), tempbar_init(j-1), &
                                 tempbar_init(j+1), tempbar_init(j+2))
             else
                slope_T = slope(sold(i,j  ,temp_comp)-tempbar_init(j), &
                                sold(i,j-1,temp_comp)-tempbar_init(j), &
                                sold(i,j+1,temp_comp)-tempbar_init(j))

                Tcoeff = [tempbar_init(j), tempbar_init(j), 0.0d0]
             endif

             ! X -- we really need to do a group limit here -- for now do 0
             slope_X = ZERO

             ! subcycle over the zones
             x_out = ZERO
             rhowdot = ZERO
             rhoH = ZERO

             do jj = 0, nsub-1
                             
                rho = sold(i,j,rho_comp) + dble(jj - nsub/2 + HALF)*slope_rho/dx(2)

                x_in(1:nspec) = sold(i,j,spec_comp:spec_comp+nspec-1) / sold(i,j,rho_comp) + &
                     dble(jj - nsub/2 + HALF)*slope_X/dx(2)

                ! T is the sum of dT interpolated + T0 reconstructed via PPM
                T_in = (sold(i,j,temp_comp)-tempbar_init(j)) + dble(jj - nsub/2 + HALF)*slope_T/dx(2) + &
                     ppm_eval(dble(jj+HALF)/nsub, Tcoeff(1), Tcoeff(2), Tcoeff(3))


                if (i == 0) then
                   print *, 'j, T_ij, T0, T_in', j, sold(i,j,temp_comp), tempbar_init(j), T_in
                endif
                sumX = ZERO
                do n = 1, nspec
                   sumX = sumX + x_in(n)
                enddo
                if (abs(sumX - ONE) > x_err) then
                   print *, x_in
                   print *, slope_X
                   call bl_error("ERROR: before burn, abundances do not sum to 1", abs(sumX-ONE))
                endif

                ! Fortran doesn't guarantee short-circuit evaluation of logicals so
                ! we need to test the value of ispec_threshold before using it 
                ! as an index in x_in
                if (ispec_threshold > 0) then
                   x_test = x_in(ispec_threshold)
                else
                   x_test = ZERO
                endif

                ! if the threshold species is not in the network, then we burn
                ! normally.  if it is in the network, make sure the mass
                ! fraction is above the cutoff.
                if (rho > burning_cutoff_density .and.           &
                     ( ispec_threshold < 0 .or.                  &
                     (ispec_threshold > 0 .and.                  &
                     x_test > burner_threshold_cutoff ))) then
                   call burner(rho, T_in, x_in, dt, x_out_temp, rhowdot_temp, rhoH_temp)
                else
                   x_out_temp = x_in
                   rhowdot_temp = 0.d0
                   rhoH_temp = 0.d0
                endif
             
                ! check if sum{X_k} = 1
                sumX = ZERO
                do n = 1, nspec
                   sumX = sumX + x_out_temp(n)
                enddo
                if (abs(sumX - ONE) > x_err) then
                   print *, sold(i,j,spec_comp:spec_comp+nspec-1)
                   call bl_error("ERROR: abundances do not sum to 1", abs(sumX-ONE))
                endif

                x_out = x_out + x_out_temp
                rhowdot = rhowdot + rhowdot_temp
                rhoH = rhoH + rhoH_temp

             enddo

             ! normalize
             x_out = x_out/nsub
             rhowdot = rhowdot/nsub
             rhoH = rhoH/nsub


             ! pass the density through
             snew(i,j,rho_comp) = sold(i,j,rho_comp)
             
             ! update the species
             snew(i,j,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * sold(i,j,rho_comp)
             
             ! store the energy generation and species creation quantities
             rho_omegadot(i,j,1:nspec) = rhowdot(1:nspec)
             rho_Hnuc(i,j) = rhoH
             
             ! update the enthalpy -- include the change due to external heating
             snew(i,j,rhoh_comp) = sold(i,j,rhoh_comp) + dt*rho_Hnuc(i,j) + dt*rho_Hext(i,j)
             
             ! pass the tracers through
             snew(i,j,trac_comp:trac_comp+ntrac-1) = sold(i,j,trac_comp:trac_comp+ntrac-1)   

          end if
          
       enddo
    enddo

  end subroutine burner_loop_2d_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_3d(tempbar_init,sold,ng_si,snew,ng_so,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,dt,lo,hi,mask)

    use bl_constants_module
    use parallel
    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, network_species_index
    use probin_module, ONLY: burning_cutoff_density, burner_threshold_species, &
         burner_threshold_cutoff, drive_initial_convection

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_rw,ng_he,ng_hn
    real(kind=dp_t), intent(in   ) ::         sold(lo(1)-ng_si:,lo(2)-ng_si:,lo(3)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real(kind=dp_t), intent(  out) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,lo(3)-ng_rw:,:)
    real(kind=dp_t), intent(  out) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:,lo(3)-ng_hn:)
    real(kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:,lo(3)-ng_he:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(0:)
    real(kind=dp_t), intent(in   ) :: dt
    logical        , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: i, j, k, n
    real (kind = dp_t) :: rho,T_in,ldt

    real (kind = dp_t) :: x_in(nspec)
    real (kind = dp_t) :: x_out(nspec)
    real (kind = dp_t) :: rhowdot(nspec)
    real (kind = dp_t) :: rhoH
    real (kind = dp_t) :: x_test
    logical            :: cell_valid
    integer, save      :: ispec_threshold
    logical, save      :: firstCall = .true.

    real (kind = dp_t) :: sumX
    real(kind=dp_t) :: r1, r2
    real (kind = dp_t), parameter :: x_err = 1.d-10

    if (firstCall) then
       ispec_threshold = network_species_index(burner_threshold_species)
       firstCall = .false.
    endif

    ldt = dt
    r1 = parallel_wtime()

    !$OMP PARALLEL DO PRIVATE(i,j,k,cell_valid,rho,x_in,T_in,x_test,x_out,rhowdot,rhoH,sumX,n) FIRSTPRIVATE(ldt) &
    !$OMP SCHEDULE(DYNAMIC,1)
    do k = lo(3), hi(3)
       print *, 'k ', k
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             cell_valid = .true.
             if ( present(mask) ) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             end if

             if (cell_valid) then

                rho = sold(i,j,k,rho_comp)
                x_in = sold(i,j,k,spec_comp:spec_comp+nspec-1) / rho

                if (drive_initial_convection) then
                   T_in = tempbar_init(k)
                else
                   T_in = sold(i,j,k,temp_comp)
                endif
                
                ! Fortran doesn't guarantee short-circuit evaluation of logicals 
                ! so we need to test the value of ispec_threshold before using it 
                ! as an index in x_in
                if (ispec_threshold > 0) then
                   x_test = x_in(ispec_threshold)
                else
                   x_test = ZERO
                endif
                
                ! if the threshold species is not in the network, then we burn
                ! normally.  if it is in the network, make sure the mass
                ! fraction is above the cutoff.
                if (rho > burning_cutoff_density .and.                &
                     ( ispec_threshold < 0 .or.                       &
                     (ispec_threshold > 0 .and.                       &
                     x_test > burner_threshold_cutoff))) then
                   !print *, 'burn ', i, j, k
                   call burner(rho, T_in, x_in, ldt, x_out, rhowdot, rhoH)
                else
                   x_out = x_in
                   rhowdot = 0.d0
                   rhoH = 0.d0
                endif
                
                ! check if sum{X_k} = 1
                sumX = ZERO
                do n = 1, nspec
                   sumX = sumX + x_out(n)
                enddo
                if (abs(sumX - ONE) > x_err) then
                   call bl_error("ERROR: abundances do not sum to 1", abs(sumX-ONE))
                endif

                ! pass the density through
                snew(i,j,k,rho_comp) = sold(i,j,k,rho_comp)
                
                ! update the species
                snew(i,j,k,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho
                
                ! store the energy generation and species create quantities
                rho_omegadot(i,j,k,1:nspec) = rhowdot(1:nspec)
                rho_Hnuc(i,j,k) = rhoH
                
                ! update the enthalpy -- include the change due to external heating
                snew(i,j,k,rhoh_comp) = sold(i,j,k,rhoh_comp) &
                     + ldt*rho_Hnuc(i,j,k) + ldt*rho_Hext(i,j,k)
                
                ! pass the tracers through
                snew(i,j,k,trac_comp:trac_comp+ntrac-1) = &
                     sold(i,j,k,trac_comp:trac_comp+ntrac-1)

             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    r2 = parallel_wtime() - r1
    print *, 'Burner loop runtime:       ', r2
    print *, '  total wrap runtime:      ', wrap_total
    print *, '    total build runtime:   ', bdf_build_total
    print *, '    total advance runtime: ', bdf_advance_total
    print *, '      total update runtime: ', bdf_update_total
    print *, '      total predict runtime: ', bdf_predict_total
    print *, '      total solve runtime: ', bdf_solve_total
    print *, '      total check runtime: ', bdf_check_total
    print *, '      total correct runtime: ', bdf_correct_total
    print *, '      total adjust runtime: ', bdf_adjust_total
    wrap_total = 0.0_dp_t
    bdf_build_total = 0.0_dp_t
    bdf_advance_total = 0.0_dp_t
    bdf_update_total =  0.0_dp_t
    bdf_predict_total = 0.0_dp_t
    bdf_solve_total = 0.0_dp_t 
    bdf_check_total = 0.0_dp_t 
    bdf_correct_total = 0.0_dp_t
    bdf_adjust_total = 0.0_dp_t 



  end subroutine burner_loop_3d

  !This loop makes use of a vectorized burner.  So the loop is trivial,
  !creating a vector of data, and the heavy numerical lifting is deferred to the
  !burner/integrator.  Useful for burners that target GPUs.
  subroutine burner_loop_3d_vec(tempbar_init,sold,ng_si,snew,ng_so,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,dt,lo,hi)
                           ! rho_Hext,ng_he,dt,lo,hi,mask)

    use bl_constants_module
    use parallel
    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, network_species_index
    use probin_module, ONLY: burning_cutoff_density, burner_threshold_species, &
         burner_threshold_cutoff, drive_initial_convection

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_rw,ng_he,ng_hn
    real(kind=dp_t), intent(in   ) ::         sold(lo(1)-ng_si:,lo(2)-ng_si:,lo(3)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real(kind=dp_t), intent(  out) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,lo(3)-ng_rw:,:)
    real(kind=dp_t), intent(  out) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:,lo(3)-ng_hn:)
    real(kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:,lo(3)-ng_he:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(0:)
    real(kind=dp_t), intent(in   ) :: dt
    !logical        , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: vec_size, x_size, y_size, z_size
    integer            :: i, j, k, n, ii
    real (kind = dp_t) :: ldt
    real (kind = dp_t), allocatable :: rho(:)
    real (kind = dp_t), allocatable :: T_in(:)
    real (kind = dp_t), allocatable :: x_in(:,:)
    real (kind = dp_t), allocatable :: x_out(:,:)
    real (kind = dp_t), allocatable :: rhowdot(:,:)
    real (kind = dp_t), allocatable :: rhoH(:)
    real (kind = dp_t) :: x_test
    !logical            :: cell_valid
    integer, save      :: ispec_threshold
    logical, save      :: firstCall = .true.

    real (kind = dp_t) :: sumX
    real(kind=dp_t) :: r1, r2
    real (kind = dp_t), parameter :: x_err = 1.d-10

    if (firstCall) then
       ispec_threshold = network_species_index(burner_threshold_species)
       firstCall = .false.
    endif
   
    ldt = dt

    ! Build input and output vectors
    x_size = (hi(1) - lo(1) + 1)
    y_size = (hi(2) - lo(2) + 1)
    z_size = (hi(3) - lo(3) + 1)
    vec_size = x_size * y_size * z_size
    allocate(    rho(vec_size))
    allocate(   T_in(vec_size))
    allocate(   x_in(nspec,vec_size))
    allocate(  x_out(nspec,vec_size))
    allocate(rhowdot(nspec,vec_size))
    allocate(   rhoH(vec_size))
    
    r1 = parallel_wtime()
    !!$OMP PARALLEL DO PRIVATE(i,j,k,cell_valid,rho,x_in,T_in,x_test,x_out,rhowdot,rhoH,sumX,n) FIRSTPRIVATE(ldt) &
    !!$OMP PARALLEL DO PRIVATE(i,j,k,rho,x_in,T_in,x_test,x_out,rhowdot,rhoH,sumX,n) FIRSTPRIVATE(ldt) &
    !!$OMP SCHEDULE(DYNAMIC,1)
    !TODO: Worth it to OMP-ify this loop?
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ii = i - lo(1) + y_size * (j - lo(2) + z_size * (k-lo(3))) + 1
             rho(ii) = sold(i,j,k,rho_comp)
             x_in(:,ii) = sold(i,j,k,spec_comp:spec_comp+nspec-1) / rho(ii)

             if (drive_initial_convection) then
                T_in(ii) = tempbar_init(k)
             else
                T_in(ii) = sold(i,j,k,temp_comp)
             endif
          end do 
       end do 
    end do 
   
    ! Call vectorized burner 
    call burner_vec(rho, T_in, x_in, ldt, x_out, rhowdot, rhoH)

    ! Process output
    !TODO: Worth it to OMP-ify this loop?
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ii = i - lo(1) + y_size * (j - lo(2) + z_size * (k-lo(3))) + 1

             ! Fortran doesn't guarantee short-circuit evaluation of logicals 
             ! so we need to test the value of ispec_threshold before using it 
             ! as an index in x_in
             if (ispec_threshold > 0) then
                x_test = x_in(ispec_threshold,ii)
             else
                x_test = ZERO
             endif

             ! if the threshold species is in the network, make sure the mass
             ! fraction is above the cutoff.  Otherwise don't use the burner's
             ! output.  
             !    TODO: As of now we burn even if we're not going to
             !    use the output assuming that the vectorized burning will be
             !    a net win even if we do some unneeded burner.  Should make
             !    sure this assumption holds and may want develop a more clever
             !    vectorized burning strategy that doesn't do any wasted burning
             if (.not. (rho(ii) > burning_cutoff_density .and.                &
                  ( ispec_threshold < 0 .or.                       &
                  (ispec_threshold > 0 .and.                       &
                  x_test > burner_threshold_cutoff)))) then
                x_out(:,ii) = x_in(:,ii)
                rhowdot(:,ii) = 0.d0
                rhoH(ii) = 0.d0
             endif
             
             ! check if sum{X_k} = 1
             sumX = ZERO
             do n = 1, nspec
                sumX = sumX + x_out(n,ii)
             enddo
             if (abs(sumX - ONE) > x_err) then
                call bl_error("ERROR: abundances do not sum to 1", abs(sumX-ONE))
             endif

             ! pass the density through
             snew(i,j,k,rho_comp) = sold(i,j,k,rho_comp)
             
             ! update the species
             snew(i,j,k,spec_comp:spec_comp+nspec-1) = x_out(1:nspec,ii) * rho(ii)
             
             ! store the energy generation and species create quantities
             rho_omegadot(i,j,k,1:nspec) = rhowdot(1:nspec,ii)
             rho_Hnuc(i,j,k) = rhoH(ii)
             
             ! update the enthalpy -- include the change due to external heating
             snew(i,j,k,rhoh_comp) = sold(i,j,k,rhoh_comp) &
                  + ldt*rho_Hnuc(i,j,k) + ldt*rho_Hext(i,j,k)
             
             ! pass the tracers through
             snew(i,j,k,trac_comp:trac_comp+ntrac-1) = &
                  sold(i,j,k,trac_comp:trac_comp+ntrac-1)

          enddo
       enddo
    enddo
    
    r2 = parallel_wtime() - r1
    print *, 'Burner loop runtime:       ', r2
    print *, '  total wrap runtime:      ', wrap_total
    print *, '    total build runtime:   ', bdf_build_total
    print *, '    total advance runtime: ', bdf_advance_total
    print *, '      total update runtime: ', bdf_update_total
    print *, '      total predict runtime: ', bdf_predict_total
    print *, '      total solve runtime: ', bdf_solve_total
    print *, '      total check runtime: ', bdf_check_total
    print *, '      total correct runtime: ', bdf_correct_total
    print *, '      total adjust runtime: ', bdf_adjust_total
    wrap_total = 0.0_dp_t
    bdf_build_total = 0.0_dp_t
    bdf_advance_total = 0.0_dp_t
    bdf_update_total =  0.0_dp_t
    bdf_predict_total = 0.0_dp_t
    bdf_solve_total = 0.0_dp_t 
    bdf_check_total = 0.0_dp_t 
    bdf_correct_total = 0.0_dp_t
    bdf_adjust_total = 0.0_dp_t 

  end subroutine burner_loop_3d_vec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_3d_sph(tempbar_init_cart,ng_tc, &
                                sold,ng_si,snew,ng_so, &
                                rho_omegadot,ng_rw, &
                                rho_Hnuc,ng_hn,rho_Hext,ng_he, &
                                dt,lo,hi,mask)

    use bl_constants_module
    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec, network_species_index
    use probin_module, ONLY: burning_cutoff_density, burner_threshold_species, &
         burner_threshold_cutoff, drive_initial_convection

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_rw,ng_he,ng_hn,ng_tc
    real(kind=dp_t), intent(in   ) ::         sold(lo(1)-ng_si:,lo(2)-ng_si:,lo(3)-ng_si:,:)
    real(kind=dp_t), intent(in   ) :: tempbar_init_cart(lo(1)-ng_tc:,lo(2)-ng_tc:,lo(3)-ng_tc:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real(kind=dp_t), intent(  out) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,lo(3)-ng_rw:,:)
    real(kind=dp_t), intent(  out) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:,lo(3)-ng_hn:)
    real(kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:,lo(3)-ng_he:)
    real(kind=dp_t), intent(in   ) :: dt
    logical        , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: i, j, k, n
    real (kind = dp_t) :: rho,T_in,ldt

    real (kind = dp_t) :: x_in(nspec)
    real (kind = dp_t) :: x_out(nspec)
    real (kind = dp_t) :: rhowdot(nspec)
    real (kind = dp_t) :: rhoH
    real (kind = dp_t) :: x_test
    logical            :: cell_valid
    integer, save      :: ispec_threshold
    logical, save      :: firstCall = .true.

    real (kind = dp_t) :: sumX
    real (kind = dp_t), parameter :: x_err = 1.d-10

    if (firstCall) then
       ispec_threshold = network_species_index(burner_threshold_species)
       firstCall = .false.
    endif

    ldt = dt

    !$OMP PARALLEL DO PRIVATE(i,j,k,cell_valid,rho,x_in,T_in,x_test,x_out,rhowdot,rhoH,sumX,n) &
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

                rho = sold(i,j,k,rho_comp)
                x_in = sold(i,j,k,spec_comp:spec_comp+nspec-1) / rho

                if (drive_initial_convection) then
                   T_in = tempbar_init_cart(i,j,k)
                else
                   T_in = sold(i,j,k,temp_comp)
                endif
             
                ! Fortran doesn't guarantee short-circuit evaluation of logicals 
                ! so we need to test the value of ispec_threshold before using it 
                ! as an index in x_in
                if (ispec_threshold > 0) then
                   x_test = x_in(ispec_threshold)
                else
                   x_test = ZERO
                endif

                ! if the threshold species is not in the network, then we burn
                ! normally.  if it is in the network, make sure the mass
                ! fraction is above the cutoff.
                if (rho > burning_cutoff_density .and.                      &
                     ( ispec_threshold < 0 .or.                              &
                     (ispec_threshold > 0 .and.                             &
                     x_test > burner_threshold_cutoff)                     &
                     )                                                       &
                     ) then
                   call burner(rho, T_in, x_in, ldt, x_out, rhowdot, rhoH)
                else
                   x_out = x_in
                   rhowdot = 0.d0
                   rhoH = 0.d0
                endif
             
                ! check if sum{X_k} = 1
                sumX = ZERO
                do n = 1, nspec
                   sumX = sumX + x_out(n)
                enddo
                if (abs(sumX - ONE) > x_err) then
                   call bl_error("ERROR: abundances do not sum to 1", abs(sumX-ONE))
                endif

                ! pass the density through
                snew(i,j,k,rho_comp) = sold(i,j,k,rho_comp)
                
                ! update the species
                snew(i,j,k,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho
                
                ! store the energy generation and species create quantities
                rho_omegadot(i,j,k,1:nspec) = rhowdot(1:nspec)
                rho_Hnuc(i,j,k) = rhoH
                
                ! update the enthalpy -- include the change due to external heating
                snew(i,j,k,rhoh_comp) = sold(i,j,k,rhoh_comp) &
                     + ldt*rho_Hnuc(i,j,k) + ldt*rho_Hext(i,j,k)
                
                ! pass the tracers through
                snew(i,j,k,trac_comp:trac_comp+ntrac-1) = &
                     sold(i,j,k,trac_comp:trac_comp+ntrac-1)

             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine burner_loop_3d_sph

end module react_state_module
