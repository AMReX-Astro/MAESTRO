module react_state_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: react_state

contains

  subroutine react_state(mla,tempbar_init,sold,snew,rho_omegadot,rho_Hnuc,rho_Hext,p0, &
                         dt,dx,the_bc_level)

    use probin_module, only: use_tfromp, do_heating, do_burning
    use variables, only: temp_comp

    use multifab_fill_ghost_module
    use multifab_physbc_module, only: multifab_physbc
    use ml_restriction_module , only: ml_cc_restriction
    use heating_module        , only: get_rho_Hext 
    use rhoh_vs_t_module      , only: makeTfromRhoP, makeTfromRhoH
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
    if (do_heating) then
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


    ! fill the boundaries
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
                                         rho_comp,dm+rho_comp,nscal,fill_crse_input=.false.)

       enddo

    endif


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
    use variables, only: rho_comp, rhoh_comp, spec_comp, temp_comp, &
                         nscal, ntrac, trac_comp, foextrap_comp
    use network, only: nspec
    use probin_module, only: do_average_burn

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


    integer :: lo(dm),hi(dm),ng_si,ng_so,ng_rw,ng_he,ng_hn
    integer :: i,n,r,comp

    type(bl_prof_timer), save :: bpt

    ! for do_average_burn
    real(kind=dp_t), allocatable :: rho_avg(:,:), T_avg(:,:), X_in_avg(:,:,:), &
         X_out_avg(:,:,:), rhowdot_avg(:,:,:), rhoH_avg(:,:)


    call build(bpt, "burner_loop")

    if (do_average_burn) then
       allocate(    rho_avg(nlevs,0:nr_fine-1))
       allocate(      T_avg(nlevs,0:nr_fine-1))
       allocate(   X_in_avg(nlevs,0:nr_fine-1,nspec))
       allocate(  X_out_avg(nlevs,0:nr_fine-1,nspec))
       allocate(rhowdot_avg(nlevs,0:nr_fine-1,nspec))
       allocate(   rhoH_avg(nlevs,0:nr_fine-1))
    endif


    ng_si = sold(1)%ng
    ng_so = snew(1)%ng
    ng_rw = rho_omegadot(1)%ng
    ng_hn = rho_Hnuc(1)%ng
    ng_he = rho_Hext(1)%ng

    if (do_average_burn) then
       ! laterally average the density, temperature and mass fractions
       ! and compute the energy release and mass fraction change in
       ! 1-d, as a function of height.
       call average(mla,sold,rho_avg,dx,rho_comp)
       call average(mla,sold,T_avg,dx,temp_comp)

       ! create X from (rho X) -- we will use snew to do the division, so we don't 
       ! modify sold
       do n = 1, nlevs
          do comp = spec_comp,spec_comp+nspec-1
             call multifab_copy_c(snew(n),comp,sold(n),comp,1,0)
             call multifab_div_div_c(snew(n),comp,sold(n),rho_comp,1)
          enddo
       enddo

       ! average the mass fractions
       do comp = spec_comp,spec_comp+nspec-1
          call average(mla,snew,X_in_avg(:,:,comp-spec_comp+1),dx,comp)          
       enddo

       ! we don't bother converting back, since the species in snew will
       ! be overwritten after the burn


       ! now call the burner on each average quantity as a function of height
       select case (dm)
       case (2)
          do n=1,nlevs_radial
             do r = 0, nr(n)-1

                ! not every point in a multilevel run will be defined
                if (rho_avg(n,r) > ZERO .and. T_avg(n,r) > ZERO) then
                   call burner(rho_avg(n,r), T_avg(n,r), X_in_avg(n,r,:), dt, &
                        X_out_avg(n,r,:), rhowdot_avg(n,r,:), rhoH_avg(n,r))
                endif
                
             enddo
          enddo

       case(3)
          call bl_error("ERROR: do_average_burn not implemented in 3-d")
       end select
    endif

    do n = 1, nlevs
       do i = 1, sold(n)%nboxes
          if ( multifab_remote(sold(n), i) ) cycle
          snp => dataptr(sold(n) , i)
          sop => dataptr(snew(n), i)
          rp => dataptr(rho_omegadot(n), i)
          hnp => dataptr(rho_Hnuc(n), i)
          hep => dataptr(rho_Hext(n), i)
          lo =  lwb(get_box(sold(n), i))
          hi =  upb(get_box(sold(n), i))
          select case (dm)
          case (2)
             if (do_average_burn) then
                call burner_loop_2d_avg(snp(:,:,1,:),ng_si,sop(:,:,1,:),ng_so,rp(:,:,1,:),ng_rw, &
                                        hnp(:,:,1,1),ng_hn,hep(:,:,1,1),ng_he,dt,lo,hi, &
                                        rho_avg(n,:), T_avg(n,:), X_in_avg(n,:,:), &
                                        X_out_avg(n,:,:), rhowdot_avg(n,:,:), rhoH_avg(n,:))
             else
                call burner_loop_2d(snp(:,:,1,:),ng_si,sop(:,:,1,:),ng_so,rp(:,:,1,:),ng_rw, &
                                    hnp(:,:,1,1),ng_hn,hep(:,:,1,1),ng_he,dt,lo,hi)
             endif
          case (3)
             call burner_loop_3d(snp(:,:,:,:),ng_si,sop(:,:,:,:),ng_so,rp(:,:,:,:),ng_rw, &
                                 hnp(:,:,:,1),ng_hn,hep(:,:,:,1),ng_he,dt,lo,hi)
          end select
       end do
    end do

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

          ! density
          call multifab_fill_ghost_cells(snew(n),snew(n-1),ng_so,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         rho_comp,dm+rho_comp,1,fill_crse_input=.false.)

          ! enthalpy
          call multifab_fill_ghost_cells(snew(n),snew(n-1),ng_so,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         rhoh_comp,dm+rhoh_comp,1,fill_crse_input=.false.)

          ! species
          call multifab_fill_ghost_cells(snew(n),snew(n-1),ng_so,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         spec_comp,dm+spec_comp,nspec,fill_crse_input=.false.)

          ! tracers
          call multifab_fill_ghost_cells(snew(n),snew(n-1),ng_so,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         trac_comp,dm+trac_comp,ntrac,fill_crse_input=.false.)

       enddo

    end if

    call destroy(bpt)

  end subroutine burner_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_2d(sold,ng_si,snew,ng_so,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,dt,lo,hi)

    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec
    use probin_module, ONLY: do_burning, burning_cutoff_density

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_rw,ng_he,ng_hn
    real(kind=dp_t), intent(in   ) ::        sold (lo(1)-ng_si:,lo(2)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real(kind=dp_t), intent(  out) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,:)
    real(kind=dp_t), intent(  out) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:)
    real(kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:)
    real(kind=dp_t), intent(in   ) :: dt

    !     Local variables
    integer            :: i, j
    real (kind = dp_t) :: rho,T_in
    real (kind = dp_t) :: x_in(nspec)
    real (kind = dp_t) :: x_out(nspec)
    real (kind = dp_t) :: rhowdot(nspec)
    real (kind = dp_t) :: rhoH

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          
          rho = sold(i,j,rho_comp)
          x_in(1:nspec) = sold(i,j,spec_comp:spec_comp+nspec-1) / rho
          T_in = sold(i,j,temp_comp)

          if (do_burning .and. rho > burning_cutoff_density) then
             call burner(rho, T_in, x_in, dt, x_out, rhowdot, rhoH)
          else
             x_out = x_in
             rhowdot = 0.d0
             rhoH = 0.d0
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
          
       enddo
    enddo

  end subroutine burner_loop_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_2d_avg(sold,ng_si,snew,ng_so,rho_omegadot,ng_rw, &
                                rho_Hnuc,ng_hn,rho_Hext,ng_he,dt,lo,hi, &
                                rho_avg,T_avg,X_in_avg, &
                                X_out_avg,rhowdot_avg,rhoH_avg)

    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec
    use probin_module, ONLY: do_burning, burning_cutoff_density, &
         do_average_burn, transverse_tol
    
    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_rw,ng_he,ng_hn
    real(kind=dp_t), intent(in   ) ::        sold (lo(1)-ng_si:,lo(2)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real(kind=dp_t), intent(  out) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,:)
    real(kind=dp_t), intent(  out) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:)
    real(kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:)
    real(kind=dp_t), intent(in   ) :: dt
    real(kind=dp_t), intent(in   ) ::     rho_avg(0:)
    real(kind=dp_t), intent(in   ) ::       T_avg(0:)
    real(kind=dp_t), intent(in   ) ::    X_in_avg(0:,:)
    real(kind=dp_t), intent(in   ) ::   X_out_avg(0:,:)
    real(kind=dp_t), intent(in   ) :: rhowdot_avg(0:,:)
    real(kind=dp_t), intent(in   ) ::    rhoH_avg(0:)

    !     Local variables
    integer          :: i, j
    real (kind=dp_t) :: rho, T_in
    real (kind=dp_t) :: x_out(nspec)
    real (kind=dp_t) :: rhowdot(nspec)
    real (kind=dp_t) :: rhoH

    real (kind=dp_t) :: err_rho, err_T

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          
          rho = sold(i,j,rho_comp)
          T_in = sold(i,j,temp_comp)

          if (do_burning .and. rho > burning_cutoff_density) then
             ! check to make sure that the average was a good representation
             ! of this zone

             if (do_average_burn) then
                err_rho = abs(rho  - rho_avg(j))/rho
                err_T = abs(T_in - T_avg(j))/T_in

                if ( err_rho > transverse_tol .or. err_T > transverse_tol) then
                   print *, 'ERROR: transverse density error: ', err_rho
                   print *, 'ERROR: transverse temperature error: ', err_T
                   call bl_error("ERROR: transverse variations too large")
                endif
             endif

             x_out = X_out_avg(j,:)
             rhowdot = rhowdot_avg(j,:)
             rhoH = rhoH_avg(j)
          else
             x_out = sold(i,j,spec_comp:spec_comp+nspec-1) / rho
             rhowdot = 0.d0
             rhoH = 0.d0
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
          
       enddo
    enddo

  end subroutine burner_loop_2d_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop_3d(sold,ng_si,snew,ng_so,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,dt,lo,hi)

    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec
    use probin_module, ONLY: do_burning, burning_cutoff_density

    integer        , intent(in   ) :: lo(:),hi(:),ng_si,ng_so,ng_rw,ng_he,ng_hn
    real(kind=dp_t), intent(in   ) ::         sold(lo(1)-ng_si:,lo(2)-ng_si:,lo(3)-ng_si:,:)
    real(kind=dp_t), intent(  out) ::         snew(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real(kind=dp_t), intent(  out) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,lo(3)-ng_rw:,:)
    real(kind=dp_t), intent(  out) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:,lo(3)-ng_hn:)
    real(kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:,lo(3)-ng_he:)
    real(kind=dp_t), intent(in   ) :: dt

    !     Local variables
    integer            :: i, j, k
    real (kind = dp_t) :: rho,T_in

    real (kind = dp_t) :: x_in(nspec)
    real (kind = dp_t) :: x_out(nspec)
    real (kind = dp_t) :: rhowdot(nspec)
    real (kind = dp_t) :: rhoH

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rho = sold(i,j,k,rho_comp)
             x_in = sold(i,j,k,spec_comp:spec_comp+nspec-1) / rho
             T_in = sold(i,j,k,temp_comp)
             
             if (do_burning .and. rho > burning_cutoff_density) then
                call burner(rho, T_in, x_in, dt, x_out, rhowdot, rhoH)
             else
                x_out = x_in
                rhowdot = 0.d0
                rhoH = 0.d0
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
                  + dt*rho_Hnuc(i,j,k) + dt*rho_Hext(i,j,k)

             ! pass the tracers through
             snew(i,j,k,trac_comp:trac_comp+ntrac-1) = &
                  sold(i,j,k,trac_comp:trac_comp+ntrac-1)
             
          enddo
       enddo
    enddo

  end subroutine burner_loop_3d

end module react_state_module
