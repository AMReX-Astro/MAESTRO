! Compute the source term to the divergence constraint, S.

module make_S_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: make_S

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_S(Source,delta_gamma1_term,delta_gamma1,state,u,rho_omegadot,rho_Hnuc, &
                     rho_Hext,thermal,p0,gamma1bar,delta_gamma1_termbar,psi,dx,mla, &
                     the_bc_level)

    use bl_constants_module
    use bl_prof_module
    use probin_module, only: use_delta_gamma1_term
    use ml_layout_module
    use average_module
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use variables, only: foextrap_comp
    use geometry, only: spherical

    type(multifab) , intent(inout) :: Source(:)
    type(multifab) , intent(inout) :: delta_gamma1_term(:)
    type(multifab) , intent(inout) :: delta_gamma1(:)
    type(multifab) , intent(in   ) :: state(:)
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(in   ) :: rho_omegadot(:)
    type(multifab) , intent(in   ) :: rho_Hnuc(:)
    type(multifab) , intent(in   ) :: rho_Hext(:)
    type(multifab) , intent(in   ) :: thermal(:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(inout) :: delta_gamma1_termbar(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    real(kind=dp_t), pointer:: srcp(:,:,:,:),dgtp(:,:,:,:),sp(:,:,:,:),up(:,:,:,:)
    real(kind=dp_t), pointer:: tp(:,:,:,:),dgp(:,:,:,:)
    real(kind=dp_t), pointer:: omegap(:,:,:,:), hep(:,:,:,:), hnp(:,:,:,:)

    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: i,n,ng_sr,ng_dt,ng_dg,ng_s,ng_u,ng_rw,ng_he,ng_hn,ng_th

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_S")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_sr = nghost(Source(1))
    ng_dt = nghost(delta_gamma1_term(1))
    ng_dg = nghost(delta_gamma1(1))
    ng_s  = nghost(state(1))
    ng_u  = nghost(u(1))
    ng_rw = nghost(rho_omegadot(1))
    ng_he = nghost(rho_Hext(1))
    ng_hn = nghost(rho_Hnuc(1))
    ng_th = nghost(thermal(1))

    do n = 1, nlevs
       do i = 1, nboxes(state(n))
          if ( multifab_remote(state(n), i) ) cycle
          srcp => dataptr(Source(n), i)
          dgtp     => dataptr(delta_gamma1_term(n), i)
          dgp    => dataptr(delta_gamma1(n), i)
          sp     => dataptr(state(n), i)
          up     => dataptr(u(n), i)
          omegap => dataptr(rho_omegadot(n), i)
          hep    => dataptr(rho_Hext(n), i)
          hnp    => dataptr(rho_Hnuc(n), i)
          tp     => dataptr(thermal(n), i)
          lo = lwb(get_box(state(n), i))
          hi = upb(get_box(state(n), i))
          select case (dm)
          case (1)
             call make_S_1d(n,lo, hi, srcp(:,1,1,1), ng_sr, dgtp(:,1,1,1), ng_dt, &
                            dgp(:,1,1,1), ng_dg, sp(:,1,1,:), ng_s, up(:,1,1,1), ng_u, &
                            omegap(:,1,1,:), ng_rw, hnp(:,1,1,1), ng_hn, &
                            hep(:,1,1,1), ng_he, &
                            tp(:,1,1,1), ng_th, p0(n,:), gamma1bar(n,:), dx(n,:))
          case (2)
             call make_S_2d(n,lo, hi, srcp(:,:,1,1), ng_sr, dgtp(:,:,1,1), ng_dt, &
                            dgp(:,:,1,1), ng_dg, sp(:,:,1,:), ng_s, up(:,:,1,:), ng_u, &
                            omegap(:,:,1,:), ng_rw, hnp(:,:,1,1), ng_hn, &
                            hep(:,:,1,1), ng_he, &
                            tp(:,:,1,1), ng_th, p0(n,:), gamma1bar(n,:), dx(n,:))
          case (3)
             if (spherical .eq. 1) then
                call make_S_3d_sphr(lo, hi, srcp(:,:,:,1), ng_sr, dgtp(:,:,:,1), ng_dt, &
                                    dgp(:,:,:,1), ng_dg, sp(:,:,:,:), ng_s, &
                                    omegap(:,:,:,:), ng_rw, hnp(:,:,:,1), ng_hn, &
                                    hep(:,:,:,1), ng_he, &
                                    tp(:,:,:,1), ng_th)
             else
                call make_S_3d(n,lo, hi, srcp(:,:,:,1), ng_sr, dgtp(:,:,:,1), ng_dt, &
                               dgp(:,:,:,1), ng_dg, sp(:,:,:,:), ng_s, up(:,:,:,:), ng_u, &
                               omegap(:,:,:,:), ng_rw, hnp(:,:,:,1), ng_hn, &
                               hep(:,:,:,1), ng_he, &
                               tp(:,:,:,1), ng_th, p0(n,:), gamma1bar(n,:), dx(n,:))
             end if
          end select
       end do
    enddo

    ! fill the ghostcells for delta_gamma1_term and Source
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(Source(nlevs))
       call multifab_fill_boundary(delta_gamma1_term(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(Source(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))
       call multifab_physbc(delta_gamma1_term(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(Source(n-1),Source(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(delta_gamma1_term(n-1),delta_gamma1_term(n), &
                                 mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(Source(n),Source(n-1), &
                                         ng_sr,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         1,foextrap_comp,1,fill_crse_input=.false.)
          
          call multifab_fill_ghost_cells(delta_gamma1_term(n),delta_gamma1_term(n-1), &
                                         ng_dt,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), the_bc_level(n), &
                                         1,foextrap_comp,1,fill_crse_input=.false.)
       enddo

    end if


    if (use_delta_gamma1_term) then

       call average(mla,delta_gamma1_term,delta_gamma1_termbar,dx,1)

       do n = 1, nlevs
          do i = 1, nboxes(state(n))
             if ( multifab_remote(state(n), i) ) cycle
             dgtp   => dataptr(delta_gamma1_term(n), i)
             dgp    => dataptr(delta_gamma1(n), i)
             select case (dm)
             case (1)
                call correct_delta_gamma1_term_1d(lo,hi,dgtp(:,1,1,1),ng_dt, &
                                                  dgp(:,1,1,1),ng_dg, &
                                                  gamma1bar(n,:),psi(n,:), &
                                                  delta_gamma1_termbar(n,:),p0(n,:))
             case (2)
                call correct_delta_gamma1_term_2d(lo,hi,dgtp(:,:,1,1),ng_dt, &
                                                  dgp(:,:,1,1),ng_dg, &
                                                  gamma1bar(n,:),psi(n,:), &
                                                  delta_gamma1_termbar(n,:),p0(n,:))
             case (3)
                if (spherical .eq. 1) then
                   call bl_error("ERROR: correct_delta_gamma1_term not implemented for spherical in make_S")
                else
                   call correct_delta_gamma1_term_3d(lo,hi,dgtp(:,:,:,1),ng_dt, &
                                                     dgp(:,:,:,1),ng_dg, &
                                                     gamma1bar(n,:),psi(n,:), &
                                                     delta_gamma1_termbar(n,:),p0(n,:))
                end if
             end select
          end do
       enddo

       ! fill ghostcells for delta_gamma1_term again, since it was just updated
       if (nlevs .eq. 1) then

          ! fill ghost cells for two adjacent grids at the same level
          ! this includes periodic domain boundary ghost cells
          call multifab_fill_boundary(delta_gamma1_term(nlevs))

          ! fill non-periodic domain boundary ghost cells
          call multifab_physbc(delta_gamma1_term(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

       else

          ! the loop over nlevs must count backwards to make sure the finer grids are done first
          do n=nlevs,2,-1

             ! set level n-1 data to be the average of the level n data covering it
             call ml_cc_restriction(delta_gamma1_term(n-1),delta_gamma1_term(n), &
                                    mla%mba%rr(n-1,:))

             call multifab_fill_ghost_cells(delta_gamma1_term(n),delta_gamma1_term(n-1), &
                                            ng_dt,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1), the_bc_level(n), &
                                            1,foextrap_comp,1,fill_crse_input=.false.)
          enddo

       end if

    end if

    call destroy(bpt)

  end subroutine make_S

  subroutine make_S_1d(n,lo,hi,Source,ng_sr,delta_gamma1_term,ng_dt,delta_gamma1,ng_dg, &
                       s,ng_s,u,ng_u,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                       rho_Hext,ng_he,thermal,ng_th, &
                       p0,gamma1bar,dx)

    use bl_constants_module
    use eos_module
    use network, only: nspec
    use variables, only: rho_comp, temp_comp, spec_comp
    use probin_module, only: use_delta_gamma1_term
    use geometry, only: anelastic_cutoff_coord, nr

    integer         , intent(in   ) :: n,lo(:),hi(:)
    integer         , intent(in   ) :: ng_sr,ng_dt,ng_dg,ng_s,ng_u,ng_rw,ng_he,ng_hn,ng_th
    real (kind=dp_t), intent(  out) ::            Source(lo(1)-ng_sr:)
    real (kind=dp_t), intent(  out) :: delta_gamma1_term(lo(1)-ng_dt:)
    real (kind=dp_t), intent(  out) ::      delta_gamma1(lo(1)-ng_dg:)
    real (kind=dp_t), intent(in   ) ::                 s(lo(1)-ng_s :,:)
    real (kind=dp_t), intent(in   ) ::                 u(lo(1)-ng_u :)
    real (kind=dp_t), intent(in   ) ::      rho_omegadot(lo(1)-ng_rw:,:)
    real (kind=dp_t), intent(in   ) ::          rho_Hnuc(lo(1)-ng_hn:)
    real (kind=dp_t), intent(in   ) ::          rho_Hext(lo(1)-ng_he:)
    real (kind=dp_t), intent(in   ) ::           thermal(lo(1)-ng_th:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer         :: i, comp
    real(kind=dp_t) :: sigma, xi_term, pres_term, gradp0

    Source = zero

    do i = lo(1), hi(1)

          den_eos = s(i,rho_comp)
          temp_eos = s(i,temp_comp)
          xn_eos(:) = s(i,spec_comp:spec_comp+nspec-1)/den_eos

          pt_index_eos(:) = (/i, -1, -1/)

          ! dens, temp, and xmass are inputs
          call eos(eos_input_rt, den_eos, temp_eos, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, & 
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false., &
                   pt_index_eos)

          sigma = dpdt_eos / (den_eos * cp_eos * dpdr_eos)

          xi_term = ZERO
          pres_term = ZERO
          do comp = 1, nspec
             xi_term = xi_term - &
                  dhdX_eos(comp)*rho_omegadot(i,comp)/den_eos 

             pres_term = pres_term + &
                  dpdX_eos(comp)*rho_omegadot(i,comp)/den_eos
          enddo

          Source(i) = (sigma/den_eos) * &
               ( rho_Hext(i) + rho_Hnuc(i) + thermal(i) ) &
               + sigma*xi_term &
               + pres_term/(den_eos*dpdr_eos)

          if (use_delta_gamma1_term .and. i < anelastic_cutoff_coord(n)) then
             if (i .eq. 0) then
                gradp0 = (p0(i+1) - p0(i))/dx(1)
             else if (i .eq. nr(n)-1) then
                gradp0 = (p0(i) - p0(i-1))/dx(1)
             else
                gradp0 = HALF*(p0(i+1) - p0(i-1))/dx(1)
             endif

             delta_gamma1(i) = gam1_eos - gamma1bar(i)

             delta_gamma1_term(i) = &
                  (gam1_eos - gamma1bar(i))*u(i)* &
                  gradp0/(gamma1bar(i)*gamma1bar(i)*p0(i))
          else
             delta_gamma1_term(i) = ZERO
             delta_gamma1(i) = ZERO
          endif

    enddo

  end subroutine make_S_1d

  subroutine make_S_2d(n,lo,hi,Source,ng_sr,delta_gamma1_term,ng_dt,delta_gamma1,ng_dg, &
                       s,ng_s,u,ng_u,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                       rho_Hext,ng_he,thermal,ng_th, &
                       p0,gamma1bar,dx)

    use bl_constants_module
    use eos_module
    use network, only: nspec
    use variables, only: rho_comp, temp_comp, spec_comp
    use probin_module, only: use_delta_gamma1_term
    use geometry, only: anelastic_cutoff_coord, nr

    integer         , intent(in   ) :: n,lo(:),hi(:)
    integer         , intent(in   ) :: ng_sr,ng_dt,ng_dg,ng_s,ng_u,ng_rw,ng_he,ng_hn,ng_th
    real (kind=dp_t), intent(  out) ::            Source(lo(1)-ng_sr:,lo(2)-ng_sr:)
    real (kind=dp_t), intent(  out) :: delta_gamma1_term(lo(1)-ng_dt:,lo(2)-ng_dt:)
    real (kind=dp_t), intent(  out) ::      delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:)
    real (kind=dp_t), intent(in   ) ::                 s(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real (kind=dp_t), intent(in   ) ::                 u(lo(1)-ng_u :,lo(2)-ng_u :,:)
    real (kind=dp_t), intent(in   ) ::      rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,:)
    real (kind=dp_t), intent(in   ) ::          rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:)
    real (kind=dp_t), intent(in   ) ::          rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:)
    real (kind=dp_t), intent(in   ) ::           thermal(lo(1)-ng_th:,lo(2)-ng_th:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer         :: i, j, comp
    real(kind=dp_t) :: sigma, xi_term, pres_term, gradp0

    Source = zero

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          den_eos = s(i,j,rho_comp)
          temp_eos = s(i,j,temp_comp)
          xn_eos(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos

          pt_index_eos(:) = (/i, j, -1/)

          ! dens, temp, and xmass are inputs
          call eos(eos_input_rt, den_eos, temp_eos, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, & 
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false., &
                   pt_index_eos)

          sigma = dpdt_eos / (den_eos * cp_eos * dpdr_eos)

          xi_term = ZERO
          pres_term = ZERO
          do comp = 1, nspec
             xi_term = xi_term - &
                  dhdX_eos(comp)*rho_omegadot(i,j,comp)/den_eos 

             pres_term = pres_term + &
                  dpdX_eos(comp)*rho_omegadot(i,j,comp)/den_eos
          enddo

          Source(i,j) = (sigma/den_eos) * &
               ( rho_Hext(i,j) + rho_Hnuc(i,j) + thermal(i,j) ) &
               + sigma*xi_term &
               + pres_term/(den_eos*dpdr_eos)

          if (use_delta_gamma1_term .and. j < anelastic_cutoff_coord(n)) then
             if (j .eq. 0) then
                gradp0 = (p0(j+1) - p0(j))/dx(2)
             else if (j .eq. nr(n)-1) then
                gradp0 = (p0(j) - p0(j-1))/dx(2)
             else
                gradp0 = HALF*(p0(j+1) - p0(j-1))/dx(2)
             endif

             delta_gamma1(i,j) = gam1_eos - gamma1bar(j)

             delta_gamma1_term(i,j) = &
                  (gam1_eos - gamma1bar(j))*u(i,j,2)* &
                  gradp0/(gamma1bar(j)*gamma1bar(j)*p0(j))
          else
             delta_gamma1_term(i,j) = ZERO
             delta_gamma1(i,j) = ZERO
          endif

       enddo
    enddo

  end subroutine make_S_2d

  subroutine make_S_3d(n,lo,hi,Source,ng_sr,delta_gamma1_term,ng_dt,delta_gamma1,ng_dg, &
                       s,ng_s,u,ng_u,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                       rho_Hext,ng_he,thermal,ng_th,p0,gamma1bar,dx)

    use bl_constants_module
    use eos_module
    use network, only: nspec
    use variables, only: rho_comp, temp_comp, spec_comp
    use probin_module, only: use_delta_gamma1_term
    use geometry, only: anelastic_cutoff_coord, nr

    integer         , intent(in   ) :: n,lo(:),hi(:)
    integer         , intent(in   ) :: ng_sr,ng_dt,ng_dg,ng_s,ng_u,ng_rw,ng_he,ng_hn,ng_th
    real (kind=dp_t), intent(  out) ::          Source(lo(1)-ng_sr:,lo(2)-ng_sr:,lo(3)-ng_sr:)
    real (kind=dp_t), intent(  out) :: delta_gamma1_term(lo(1)-ng_dt:,lo(2)-ng_dt:,lo(3)-ng_dt:)
    real (kind=dp_t), intent(  out) :: delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:) 
    real (kind=dp_t), intent(in   ) ::            s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real (kind=dp_t), intent(in   ) ::            u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)
    real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,lo(3)-ng_rw:,:)
    real (kind=dp_t), intent(in   ) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:,lo(3)-ng_hn:)
    real (kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:,lo(3)-ng_he:)
    real (kind=dp_t), intent(in   ) ::      thermal(lo(1)-ng_th:,lo(2)-ng_th:,lo(3)-ng_th:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer         :: i, j, k, comp
    real(kind=dp_t) :: sigma, xi_term, pres_term, gradp0

    Source = zero

    !$OMP PARALLEL DO PRIVATE(i,j,k,comp,sigma,xi_term,pres_term,gradp0)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos = s(i,j,k,rho_comp)
             temp_eos = s(i,j,k,temp_comp)
             xn_eos(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos

             pt_index_eos(:) = (/i, j, k/)

             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, den_eos, temp_eos, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, & 
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)

             sigma = dpdt_eos / (den_eos * cp_eos * dpdr_eos)

             xi_term = ZERO
             pres_term = ZERO
             do comp = 1, nspec
                xi_term = xi_term - &
                     dhdX_eos(comp)*rho_omegadot(i,j,k,comp)/den_eos 

                pres_term = pres_term + &
                     dpdX_eos(comp)*rho_omegadot(i,j,k,comp)/den_eos
             enddo

             Source(i,j,k) = (sigma/den_eos) * &
                  ( rho_Hext(i,j,k) + rho_Hnuc(i,j,k) + thermal(i,j,k) ) &
                  + sigma*xi_term &
                  + pres_term/(den_eos*dpdr_eos)

             if (use_delta_gamma1_term .and. k < anelastic_cutoff_coord(n)) then
                if (k .eq. 0) then
                   gradp0 = (p0(k+1) - p0(k))/dx(3)
                else if (k .eq. nr(n)-1) then
                   gradp0 = (p0(k) - p0(k-1))/dx(3)
                else
                   gradp0 = HALF*(p0(k+1) - p0(k-1))/dx(3)
                endif
                
                delta_gamma1(i,j,k) = gam1_eos - gamma1bar(k)
                
                delta_gamma1_term(i,j,k) = (gam1_eos - gamma1bar(k))*u(i,j,k,3)* &
                     gradp0/(gamma1bar(k)*gamma1bar(k)*p0(k))
             else
                delta_gamma1_term(i,j,k) = ZERO
                delta_gamma1(i,j,k) = ZERO
             endif
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_S_3d

  subroutine make_S_3d_sphr(lo,hi,Source,ng_sr,dg1_term,ng_dt,delta_gamma1, &
                            ng_dg,s,ng_s,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,thermal,ng_th)

    use bl_constants_module
    use eos_module
    use network, only: nspec
    use variables, only: rho_comp, temp_comp, spec_comp
    use probin_module, only: use_delta_gamma1_term
    use geometry, only: anelastic_cutoff_coord, nr

    integer         , intent(in   ) :: lo(:),hi(:)
    integer         , intent(in   ) :: ng_sr,ng_dt,ng_dg,ng_s,ng_rw,ng_he,ng_hn,ng_th
    real (kind=dp_t), intent(  out) ::       Source(lo(1)-ng_sr:,lo(2)-ng_sr:,lo(3)-ng_sr:)
    real (kind=dp_t), intent(  out) ::     dg1_term(lo(1)-ng_dt:,lo(2)-ng_dt:,lo(3)-ng_dt:)
    real (kind=dp_t), intent(  out) :: delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:) 
    real (kind=dp_t), intent(in   ) ::            s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,lo(3)-ng_rw:,:)
    real (kind=dp_t), intent(in   ) ::     rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:,lo(3)-ng_hn:)
    real (kind=dp_t), intent(in   ) ::     rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:,lo(3)-ng_he:)
    real (kind=dp_t), intent(in   ) ::      thermal(lo(1)-ng_th:,lo(2)-ng_th:,lo(3)-ng_th:)

    !     Local variables
    integer         :: i, j, k, comp
    real(kind=dp_t) :: sigma, xi_term, pres_term

    Source = zero

    !$OMP PARALLEL DO PRIVATE(i,j,k,comp,sigma,xi_term,pres_term)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos = s(i,j,k,rho_comp)
             temp_eos = s(i,j,k,temp_comp)
             xn_eos(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos

             pt_index_eos(:) = (/i, j, k/)

             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, den_eos, temp_eos, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, & 
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)

             sigma = dpdt_eos / (den_eos * cp_eos * dpdr_eos)

             xi_term = ZERO
             pres_term = ZERO
             do comp = 1, nspec
                xi_term = xi_term - &
                     dhdX_eos(comp)*rho_omegadot(i,j,k,comp)/den_eos 

                pres_term = pres_term + &
                     dpdX_eos(comp)*rho_omegadot(i,j,k,comp)/den_eos
             enddo

             Source(i,j,k) = (sigma/den_eos) * &
                  ( rho_Hext(i,j,k) + rho_Hnuc(i,j,k) + thermal(i,j,k) ) &
                  + sigma*xi_term &
                  + pres_term/(den_eos*dpdr_eos)


             if (use_delta_gamma1_term) then
                call bl_error("ERROR: use_delta_gamma1_term not implemented for spherical in make_S")
             else
                dg1_term(i,j,k) = ZERO
                delta_gamma1(i,j,k) = ZERO
             end if

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_S_3d_sphr

  subroutine correct_delta_gamma1_term_1d(lo,hi,delta_gamma1_term,ng_dt,delta_gamma1,ng_dg, &
                                          gamma1bar,psi,delta_gamma1_termbar,p0)

    integer         , intent(in   ) :: lo(:), hi(:), ng_dt, ng_dg
    real (kind=dp_t), intent(inout) :: delta_gamma1_term(lo(1)-ng_dt:)
    real (kind=dp_t), intent(in   ) ::      delta_gamma1(lo(1)-ng_dg:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: psi(0:)
    real (kind=dp_t), intent(in   ) :: delta_gamma1_termbar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    integer :: i

    do i = lo(1), hi(1)

       delta_gamma1_term(i) = delta_gamma1_term(i) - delta_gamma1_termbar(i) &
            + delta_gamma1(i)*psi(i)/(gamma1bar(i)**2*p0(i))

    end do

  end subroutine correct_delta_gamma1_term_1d

  subroutine correct_delta_gamma1_term_2d(lo,hi,delta_gamma1_term,ng_dt,delta_gamma1,ng_dg, &
                                          gamma1bar,psi,delta_gamma1_termbar,p0)

    integer         , intent(in   ) :: lo(:), hi(:), ng_dt, ng_dg
    real (kind=dp_t), intent(inout) :: delta_gamma1_term(lo(1)-ng_dt:,lo(2)-ng_dt:)
    real (kind=dp_t), intent(in   ) ::      delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: psi(0:)
    real (kind=dp_t), intent(in   ) :: delta_gamma1_termbar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          delta_gamma1_term(i,j) = delta_gamma1_term(i,j) - delta_gamma1_termbar(j) &
               + delta_gamma1(i,j)*psi(j)/(gamma1bar(j)**2*p0(j))

       end do
    end do

  end subroutine correct_delta_gamma1_term_2d

  subroutine correct_delta_gamma1_term_3d(lo,hi,delta_gamma1_term,ng_dt,delta_gamma1,ng_dg, &
                                          gamma1bar,psi,delta_gamma1_termbar,p0)

    integer         , intent(in   ) :: lo(:), hi(:), ng_dt, ng_dg
    real (kind=dp_t), intent(inout) :: delta_gamma1_term(lo(1)-ng_dt:,lo(2)-ng_dt:,lo(3)-ng_dt:)
    real (kind=dp_t), intent(in   ) ::    delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: psi(0:)
    real (kind=dp_t), intent(in   ) :: delta_gamma1_termbar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             delta_gamma1_term(i,j,k) = delta_gamma1_term(i,j,k) - delta_gamma1_termbar(k) &
                  + delta_gamma1(i,j,k)*psi(k)/(gamma1bar(k)**2*p0(k))

          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine correct_delta_gamma1_term_3d

end module make_S_module
