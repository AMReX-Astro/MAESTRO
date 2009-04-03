module react_state_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: react_state

contains

  subroutine react_state(mla,sold,snew,rho_omegadot,rho_Hnuc,rho_Hext,p0, &
                         dt,dx,the_bc_level,time)

    use heating_module
    use probin_module, only: use_tfromp
    use rhoh_vs_t_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: rho_omegadot(:)
    type(multifab) , intent(inout) :: rho_Hnuc(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(dp_t)     , intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:,:),time
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local
    type(bl_prof_timer), save :: bpt

    call build(bpt, "react_state")

    ! get heating term
    call get_rho_Hext(mla,sold,rho_Hext,dx,time,the_bc_level)

    ! do the burning
    call burner_loop(mla,sold,snew,rho_omegadot,rho_Hnuc,rho_Hext,dt,the_bc_level)

    ! now update temperature
    if (use_tfromp) then
       call makeTfromRhoP(snew,p0,sold,mla,the_bc_level,dx)
    else
       call makeTfromRhoH(snew,sold,mla,the_bc_level)
    end if

    call destroy(bpt)

  end subroutine react_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine burner_loop(mla,sold,snew,rho_omegadot,rho_Hnuc,rho_Hext,dt,the_bc_level)

    use geometry, only: dm, nlevs
    use variables, only: rho_comp, nscal
    use multifab_fill_ghost_module
    use ml_restriction_module
    use multifab_physbc_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: rho_omegadot(:)
    type(multifab) , intent(inout) :: rho_Hnuc(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local
    real(kind=dp_t), pointer:: snp(:,:,:,:)
    real(kind=dp_t), pointer:: sop(:,:,:,:)
    real(kind=dp_t), pointer::  rp(:,:,:,:)
    real(kind=dp_t), pointer::  hnp(:,:,:,:)
    real(kind=dp_t), pointer::  hep(:,:,:,:)

    integer :: lo(dm),hi(dm),ng_si,ng_so,ng_rw,ng_he,ng_hn
    integer :: i,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "burner_loop")

    ng_si = sold(1)%ng
    ng_so = snew(1)%ng
    ng_rw = rho_omegadot(1)%ng
    ng_hn = rho_Hnuc(1)%ng
    ng_he = rho_Hext(1)%ng

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
             call burner_loop_2d(snp(:,:,1,:),ng_si,sop(:,:,1,:),ng_so,rp(:,:,1,:),ng_rw, &
                                 hnp(:,:,1,1),ng_hn,hep(:,:,1,1),ng_he,dt,lo,hi)
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
          call multifab_fill_ghost_cells(snew(n),snew(n-1),ng_so,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         rho_comp,dm+rho_comp,nscal,fill_crse_input=.false.)
       enddo

    end if

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

    logical :: update_temp

    update_temp = .false.

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

    logical :: update_temp

    update_temp = .false.

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
