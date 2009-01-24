module react_state_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: react_state

contains

  subroutine react_state(mla,s_in,s_out,rho_omegadot,rho_Hnuc,rho_Hext,tempbar,p0, &
                         dt,dx,the_bc_level,time)

    use variables, only: rho_comp, nscal, foextrap_comp, temp_comp
    use network, only: nspec
    use bl_prof_module
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use heating_module
    use geometry, only: dm, nlevs
    use probin_module, only: use_tfromp
    use rhoh_vs_t_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s_in(:)
    type(multifab) , intent(inout) :: s_out(:)
    type(multifab) , intent(inout) :: rho_omegadot(:)
    type(multifab) , intent(inout) :: rho_Hnuc(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(dp_t)     , intent(in   ) :: tempbar(:,0:)
    real(dp_t)     , intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:,:),time
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local
    real(kind=dp_t), pointer:: sinp(:,:,:,:)
    real(kind=dp_t), pointer:: sotp(:,:,:,:)
    real(kind=dp_t), pointer::   rp(:,:,:,:)
    real(kind=dp_t), pointer::   hnp(:,:,:,:)
    real(kind=dp_t), pointer::   hep(:,:,:,:)

    integer :: lo(dm),hi(dm),ng_si,ng_so,ng_rw,ng_he,ng_hn
    integer :: i,n,ispec

    type(bl_prof_timer), save :: bpt

    call build(bpt, "react_state")

    ng_si = s_in(1)%ng
    ng_so = s_out(1)%ng
    ng_rw = rho_omegadot(1)%ng
    ng_hn = rho_Hnuc(1)%ng
    ng_he = rho_Hext(1)%ng

    call get_rho_Hext(mla,s_in,rho_Hext,dx,time)

    do n = 1, nlevs
       do i = 1, s_in(n)%nboxes
          if ( multifab_remote(s_in(n), i) ) cycle
          sinp => dataptr(s_in(n) , i)
          sotp => dataptr(s_out(n), i)
          rp => dataptr(rho_omegadot(n), i)
          hnp => dataptr(rho_Hnuc(n), i)
          hep => dataptr(rho_Hext(n), i)
          lo =  lwb(get_box(s_in(n), i))
          hi =  upb(get_box(s_in(n), i))
          select case (dm)
          case (2)
             call react_state_2d(sinp(:,:,1,:),ng_si,sotp(:,:,1,:),ng_so, &
                                 rp(:,:,1,:),ng_rw,hnp(:,:,1,1),ng_hn, &
                                 hep(:,:,1,1),ng_he,dt,lo,hi)
          case (3)
             call react_state_3d(sinp(:,:,:,:),ng_si,sotp(:,:,:,:),ng_so, &
                                 rp(:,:,:,:),ng_rw,hnp(:,:,:,1),ng_hn, &
                                 hep(:,:,:,1),ng_he,dt,lo,hi)
          end select
       end do
    end do

!    ! now update temperature
    if (use_tfromp) then
       ! option to update with eos
!       call makeTfromRhoP(s_out,p0,tempbar,mla,the_bc_level,dx)
       ! option to pass temperature through
       do n=1,nlevs
          call multifab_copy_c(s_out(n),temp_comp,s_in(n),temp_comp,1,3)
       end do
    else
       call makeTfromRhoH(s_out,tempbar,mla,the_bc_level,dx)
    end if

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(s_out(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s_out(nlevs),rho_comp,dm+rho_comp,nscal,the_bc_level(nlevs))

       ! since rho_omegadot, rho_Hnuc, and rho_Hext are going to be averaged later, we 
       ! need to also fill the ghostcells on those -- we'll use extrapolation
       call multifab_fill_boundary(rho_omegadot(nlevs))
       call multifab_fill_boundary(rho_Hnuc(nlevs))
       call multifab_fill_boundary(rho_Hext(nlevs))

       do ispec = 1, nspec
          call multifab_physbc(rho_omegadot(nlevs),ispec,foextrap_comp,1,the_bc_level(nlevs))
       enddo
       call multifab_physbc(rho_Hnuc(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))
       call multifab_physbc(rho_Hext(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(s_out(n-1)       ,s_out(n)       ,mla%mba%rr(n-1,:))
          call ml_cc_restriction(rho_omegadot(n-1),rho_omegadot(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(rho_Hext(n-1)    ,rho_Hext(n)    ,mla%mba%rr(n-1,:))
          call ml_cc_restriction(rho_Hnuc(n-1)    ,rho_Hnuc(n)    ,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s_out(n),s_out(n-1), &
                                         ng_so,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         rho_comp,dm+rho_comp,nscal)

          do ispec = 1, nspec
             call multifab_fill_ghost_cells(rho_omegadot(n),rho_omegadot(n-1), &
                                            ng_rw,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1), the_bc_level(n), &
                                            ispec,foextrap_comp,1)
          enddo

          call multifab_fill_ghost_cells(rho_Hnuc(n),rho_Hnuc(n-1), &
                                         ng_hn,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         1,foextrap_comp,1)

          call multifab_fill_ghost_cells(rho_Hext(n),rho_Hext(n-1), &
                                         ng_he,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         1,foextrap_comp,1)
       enddo

    end if

    call destroy(bpt)

  end subroutine react_state

  subroutine react_state_2d(s_in,ng_si,s_out,ng_so, &
                            rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,dt,lo,hi)

    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec
    use probin_module, ONLY: do_burning, burning_cutoff_density, enthalpy_pred_type, &
         use_tfromp
    use pred_parameters
    use bl_constants_module, only: zero
    use eos_module

    integer, intent(in) :: lo(:), hi(:), ng_si, ng_so, ng_rw, ng_he, ng_hn
    real (kind = dp_t), intent(in   ) ::        s_in (lo(1)-ng_si:,lo(2)-ng_si:,:)
    real (kind = dp_t), intent(  out) ::        s_out(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real (kind = dp_t), intent(  out) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,:)
    real (kind = dp_t), intent(  out) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:)
    real (kind = dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:)
    real (kind = dp_t), intent(in   ) :: dt

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
          
          rho = s_in(i,j,rho_comp)
          x_in(1:nspec) = s_in(i,j,spec_comp:spec_comp+nspec-1) / rho
          T_in = s_in(i,j,temp_comp)

          if (do_burning .and. rho > burning_cutoff_density) then
             call burner(rho, T_in, x_in, dt, x_out, rhowdot, rhoH)
          else
             x_out = x_in
             rhowdot = ZERO
             rhoH = ZERO
          endif

          ! pass the density through
          s_out(i,j,rho_comp) = s_in(i,j,rho_comp)

          ! update the species
          s_out(i,j,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho

          ! store the energy generation and species creation quantities
          rho_omegadot(i,j,1:nspec) = rhowdot(1:nspec)
          rho_Hnuc(i,j) = rhoH

          ! update the enthalpy -- include the change due to external heating
          s_out(i,j,rhoh_comp) = s_in(i,j,rhoh_comp) + dt*rho_Hnuc(i,j) + dt*rho_Hext(i,j)

          ! pass the tracers through
          s_out(i,j,trac_comp:trac_comp+ntrac-1) = s_in(i,j,trac_comp:trac_comp+ntrac-1)   
          
       enddo
    enddo

  end subroutine react_state_2d

  subroutine react_state_3d(s_in,ng_si,s_out,ng_so, &
                            rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,dt,lo,hi)

    use burner_module
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp, trac_comp, ntrac
    use network, only: nspec
    use probin_module, ONLY: do_burning, burning_cutoff_density, enthalpy_pred_type, &
         use_tfromp
    use pred_parameters
    use bl_constants_module, only: zero
    use eos_module

    integer, intent(in)             :: lo(:), hi(:), ng_si, ng_so, ng_rw, ng_he, ng_hn
    real (kind = dp_t),intent(in   )::         s_in(lo(1)-ng_si:,lo(2)-ng_si:,lo(3)-ng_si:,:)
    real (kind = dp_t),intent(  out)::        s_out(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real (kind = dp_t),intent(  out):: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,lo(3)-ng_rw:,:)
    real (kind = dp_t),intent(  out)::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:,lo(3)-ng_hn:)
    real (kind = dp_t),intent(in   )::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:,lo(3)-ng_he:)
    real (kind = dp_t),intent(in   ):: dt

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

             rho = s_in(i,j,k,rho_comp)
             x_in = s_in(i,j,k,spec_comp:spec_comp+nspec-1) / rho
             T_in = s_in(i,j,k,temp_comp)
             
             if (do_burning .and. rho > burning_cutoff_density) then
                call burner(rho, T_in, x_in, dt, x_out, rhowdot, rhoH)
             else
                x_out = x_in
                rhowdot = ZERO
                rhoH = ZERO
             endif
             
             ! pass the density through
             s_out(i,j,k,rho_comp) = s_in(i,j,k,rho_comp)

             ! update the species
             s_out(i,j,k,spec_comp:spec_comp+nspec-1) = x_out(1:nspec) * rho
             
             ! store the energy generation and species create quantities
             rho_omegadot(i,j,k,1:nspec) = rhowdot(1:nspec)
             rho_Hnuc(i,j,k) = rhoH

             ! update the enthalpy -- include the change due to external heating
             s_out(i,j,k,rhoh_comp) = s_in(i,j,k,rhoh_comp) &
                  + dt*rho_Hnuc(i,j,k) + dt*rho_Hext(i,j,k)

             ! pass the tracers through
             s_out(i,j,k,trac_comp:trac_comp+ntrac-1) = &
                  s_in(i,j,k,trac_comp:trac_comp+ntrac-1)
             
          enddo
       enddo
    enddo

  end subroutine react_state_3d

end module react_state_module
