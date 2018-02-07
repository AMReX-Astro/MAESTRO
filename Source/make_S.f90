! Compute the source term to the divergence constraint, S.
!
! We also construct the delta_gamma1_term here which is used to
! correct for the approximation that Gamma_1 -> \bar{Gamma_1} in the
! constraint.  First we compute the main part of the term and then
! we average it.  Finally we correct the full version of the term to
! account for the psi term that comes from the base state.  

module make_S_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: make_S

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_S(Source,delta_gamma1_term,delta_gamma1, &
                    state,u, &
                    normal, &
                    rho_omegadot,rho_Hnuc,rho_Hext,thermal, &
                    p0,gamma1bar,delta_gamma1_termbar,psi, &
                    dx,mla,the_bc_level)

    use bl_constants_module
    use bl_prof_module
    use probin_module, only: use_delta_gamma1_term
    use ml_layout_module
    use average_module
    use ml_restrict_fill_module
    use variables, only: foextrap_comp
    use geometry, only: spherical, polar, nr_fine, dr
    use fill_3d_module, only : put_1d_array_on_cart

    type(multifab) , intent(inout) :: Source(:)
    type(multifab) , intent(inout) :: delta_gamma1_term(:)
    type(multifab) , intent(inout) :: delta_gamma1(:)
    type(multifab) , intent(in   ) :: state(:)
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(in   ) :: normal(:)
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
    real(kind=dp_t), pointer:: np(:,:,:,:)
    real(kind=dp_t), pointer:: gp0p(:,:,:,:), p0p(:,:,:,:), g1p(:,:,:,:), psp(:,:,:,:)

    integer :: lo(MAX_SPACEDIM),hi(MAX_SPACEDIM),dm,nlevs
    integer :: i,n, r
    integer :: ng_sr,ng_dt,ng_dg,ng_s,ng_u,ng_rw,ng_he,ng_hn,ng_th
    integer :: ng_gp, ng_p0, ng_g1, ng_n, ng_ps

    real(kind=dp_t) :: gradp0_rad(1,0:nr_fine-1)
    type(multifab) :: gradp0_cart(mla%nlevel)
    type(multifab) :: p0_cart(mla%nlevel)
    type(multifab) :: psi_cart(mla%nlevel)
    type(multifab) :: gamma1bar_cart(mla%nlevel)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_S")

    dm = mla%dim
    nlevs = mla%nlevel

    if (spherical == 1 .or. polar == 1) then
       do n = 1, nlevs
          call build(gradp0_cart(n), mla%la(n), 1, 0)
          call setval(gradp0_cart(n), ZERO, all=.true.)

          call build(p0_cart(n), mla%la(n), 1, 0)
          call setval(p0_cart(n), ZERO, all=.true.)

          call build(psi_cart(n), mla%la(n), 1, 0)
          call setval(psi_cart(n), ZERO, all=.true.)

          call build(gamma1bar_cart(n), mla%la(n), 1, 0)
          call setval(gamma1bar_cart(n), ZERO, all=.true.)
       enddo

       if (use_delta_gamma1_term) then
          ! compute gradp0 and put it on a cart
          do r = 0, nr_fine-1
             if (r == 0) then
                gradp0_rad(1,r) = (p0(1,r+1) - p0(1,r))/dr(1)
             else if (r == nr_fine-1) then
                gradp0_rad(1,r) = (p0(1,r) - p0(1,r-1))/dr(1)
             else
                gradp0_rad(1,r) = HALF*(p0(1,r+1) - p0(1,r-1))/dr(1)
             endif
          enddo

          call put_1d_array_on_cart(gradp0_rad,gradp0_cart,foextrap_comp, &
                                    .false.,.false.,dx,the_bc_level,mla)

          call put_1d_array_on_cart(p0,p0_cart,foextrap_comp, &
                                    .false.,.false.,dx,the_bc_level,mla)

          call put_1d_array_on_cart(psi,psi_cart,foextrap_comp, &
                                    .false.,.false.,dx,the_bc_level,mla)

          call put_1d_array_on_cart(gamma1bar,gamma1bar_cart,foextrap_comp, &
                                    .false.,.false.,dx,the_bc_level,mla)
       endif
    endif


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

       lo(:) = 1; hi(:) = 1
       
       do i = 1, nfabs(state(n))
          srcp => dataptr(Source(n), i)
          dgtp     => dataptr(delta_gamma1_term(n), i)
          dgp    => dataptr(delta_gamma1(n), i)
          sp     => dataptr(state(n), i)
          up     => dataptr(u(n), i)
          omegap => dataptr(rho_omegadot(n), i)
          hep    => dataptr(rho_Hext(n), i)
          hnp    => dataptr(rho_Hnuc(n), i)
          tp     => dataptr(thermal(n), i)
          lo(1:dm) = lwb(get_box(state(n), i))
          hi(1:dm) = upb(get_box(state(n), i))

          if (spherical == 1 .or. polar == 1) then

             gp0p => dataptr(gradp0_cart(n), i)
             ng_gp = nghost(gradp0_cart(1))
             
             p0p => dataptr(p0_cart(n), i)
             ng_p0 = nghost(p0_cart(1))
             
             g1p => dataptr(gamma1bar_cart(n), i)
             ng_g1 = nghost(gamma1bar_cart(1))

             np => dataptr(normal(n), i)
             ng_n = nghost(normal(1))
             
             select case(dm)
             case(2)
                call make_S_2d_polar(lo, hi, srcp(:,:,1,1), ng_sr, dgtp(:,:,1,1), ng_dt, &
                                 dgp(:,:,1,1), ng_dg, &
                                 sp(:,:,1,:), ng_s, up(:,:,1,:), ng_u, &
                                 omegap(:,:,1,:), ng_rw, hnp(:,:,1,1), ng_hn, &
                                 hep(:,:,1,1), ng_he, &
                                 tp(:,:,1,1), ng_th, &
                                 gp0p(:,:,1,1), ng_gp, &
                                 p0p(:,:,1,1), ng_p0, &
                                 g1p(:,:,1,1), ng_g1, &
                                 np(:,:,1,:), ng_n)
 
             case(3)
                call make_S_3d_sphr(lo, hi, srcp(:,:,:,1), ng_sr, dgtp(:,:,:,1), ng_dt, &
                                 dgp(:,:,:,1), ng_dg, &
                                 sp(:,:,:,:), ng_s, up(:,:,:,:), ng_u, &
                                 omegap(:,:,:,:), ng_rw, hnp(:,:,:,1), ng_hn, &
                                 hep(:,:,:,1), ng_he, &
                                 tp(:,:,:,1), ng_th, &
                                 gp0p(:,:,:,1), ng_gp, &
                                 p0p(:,:,:,1), ng_p0, &
                                 g1p(:,:,:,1), ng_g1, &
                                 np(:,:,:,:), ng_n)
             end select
          else
             call make_S_cart(dm, n, lo, hi, &
                              srcp(:,:,:,1), lbound(srcp), &
                              dgtp(:,:,:,1), lbound(dgtp), &
                              dgp(:,:,:,1), lbound(dgp), &
                              sp, lbound(sp), &
                              up, lbound(up), &
                              omegap, lbound(omegap), &
                              hnp(:,:,:,1), lbound(hnp), &
                              hep(:,:,:,1), lbound(hep), &
                              tp(:,:,:,1), lbound(tp), &
                              p0(n,:), gamma1bar(n,:), dx(n,:))
          end if

       end do
    enddo

    ! restrict data (has no ghost cells)
    call ml_restrict_and_fill(nlevs,Source,mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=Source(1)%ng)

    ! restrict data (has no ghost cells)
    call ml_restrict_and_fill(nlevs,delta_gamma1_term,mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=delta_gamma1_term(1)%ng)

    if (use_delta_gamma1_term) then

       call average(mla,delta_gamma1_term,delta_gamma1_termbar,dx,1)

       do n = 1, nlevs
          do i = 1, nfabs(state(n))
             dgtp   => dataptr(delta_gamma1_term(n), i)
             dgp    => dataptr(delta_gamma1(n), i)

             lo(1:dm) = lwb(get_box(state(n), i))
             hi(1:dm) = upb(get_box(state(n), i))

             select case (dm)
             case (1)
                call correct_delta_gamma1_term_1d(lo,hi,dgtp(:,1,1,1),ng_dt, &
                                                  dgp(:,1,1,1),ng_dg, &
                                                  gamma1bar(n,:),psi(n,:),p0(n,:))
             case (2)
                if (polar == 1) then
                    g1p => dataptr(gamma1bar_cart(n), i)
                    ng_g1 = nghost(gamma1bar_cart(1))

                    p0p => dataptr(p0_cart(n), i)
                    ng_p0 = nghost(p0_cart(1))

                    psp => dataptr(psi_cart(n), i)
                    ng_ps = nghost(psi_cart(1))

                    call correct_delta_gamma1_term_2d_polar(lo,hi,dgtp(:,:,1,1),ng_dt, &
                                                          dgp(:,:,1,1),ng_dg, &
                                                          g1p(:,:,1,1),ng_g1, &
                                                          psp(:,:,1,1),ng_ps, &
                                                          p0p(:,:,1,1),ng_p0)
               
                else
                    call correct_delta_gamma1_term_2d(lo,hi,dgtp(:,:,1,1),ng_dt, &
                                                  dgp(:,:,1,1),ng_dg, &
                                                  gamma1bar(n,:),psi(n,:),p0(n,:))
                end if
             case (3)
                if (spherical == 1) then
                   g1p => dataptr(gamma1bar_cart(n), i)
                   ng_g1 = nghost(gamma1bar_cart(1))

                   p0p => dataptr(p0_cart(n), i)
                   ng_p0 = nghost(p0_cart(1))

                   psp => dataptr(psi_cart(n), i)
                   ng_ps = nghost(psi_cart(1))

                   call correct_delta_gamma1_term_3d_sphr(lo,hi,dgtp(:,:,:,1),ng_dt, &
                                                          dgp(:,:,:,1),ng_dg, &
                                                          g1p(:,:,:,1),ng_g1, &
                                                          psp(:,:,:,1),ng_ps, &
                                                          p0p(:,:,:,1),ng_p0)
                else
                   call correct_delta_gamma1_term_3d(lo,hi,dgtp(:,:,:,1),ng_dt, &
                                                     dgp(:,:,:,1),ng_dg, &
                                                     gamma1bar(n,:),psi(n,:),p0(n,:))
                end if
             end select
          end do
       enddo

    end if

    if (spherical == 1 .or. polar == 1) then
       do n = 1, nlevs
          call destroy(gradp0_cart(n))
          call destroy(p0_cart(n))
          call destroy(psi_cart(n))
          call destroy(gamma1bar_cart(n))
       enddo
    endif

    call destroy(bpt)

  end subroutine make_S

  subroutine make_S_cart(dm, n, lo, hi, &
                         Source, srlo, &
                         delta_gamma1_term, dtlo, &
                         delta_gamma1, dglo, &
                         s, slo, &
                         u, ulo, &
                         rho_omegadot, rwlo, &
                         rho_Hnuc, hnlo, &
                         rho_Hext, helo, &
                         thermal, thlo, &
                         p0,gamma1bar,dx)

    use bl_constants_module
    use eos_module, only: eos, eos_input_rt
    use eos_type_module
    use network, only: nspec
    use variables, only: rho_comp, temp_comp, spec_comp
    use probin_module, only: use_delta_gamma1_term
    use geometry, only: anelastic_cutoff_coord, nr

    integer         , intent(in   ) :: dm, n, lo(:),hi(:)
    integer         , intent(in   ) :: srlo(4), dtlo(4), dglo(4)
    integer         , intent(in   ) :: slo(4), ulo(4), rwlo(4)
    integer         , intent(in   ) :: hnlo(4), helo(4), thlo(4)

    real (kind=dp_t), intent(  out) ::            Source(srlo(1): ,srlo(2): ,srlo(3): )
    real (kind=dp_t), intent(  out) :: delta_gamma1_term(dtlo(1): ,dtlo(2): ,dtlo(3): )
    real (kind=dp_t), intent(  out) ::      delta_gamma1(dglo(1): ,dglo(2): ,dglo(3): )
    real (kind=dp_t), intent(in   ) ::                 s( slo(1): , slo(2): , slo(3): , slo(4): )
    real (kind=dp_t), intent(in   ) ::                 u( ulo(1): , ulo(2): , ulo(3): , ulo(4): )
    real (kind=dp_t), intent(in   ) ::      rho_omegadot(rwlo(1): ,rwlo(2): ,rwlo(3): ,rwlo(4): )
    real (kind=dp_t), intent(in   ) ::          rho_Hnuc(hnlo(1): ,hnlo(2): ,hnlo(3): )
    real (kind=dp_t), intent(in   ) ::          rho_Hext(helo(1): ,helo(2): ,helo(3): )
    real (kind=dp_t), intent(in   ) ::           thermal(thlo(1): ,thlo(2): ,thlo(3): )
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer         :: i, j, k, r, comp
    real(kind=dp_t) :: sigma, xi_term, pres_term, gradp0

    integer :: pt_index(MAX_SPACEDIM)
    type(eos_t) :: eos_state

    Source = zero
    r = -1
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,r,comp,sigma,xi_term,pres_term,gradp0,eos_state,pt_index)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = s(i,j,k,rho_comp)
             eos_state%T     = s(i,j,k,temp_comp)
             eos_state%xn(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, eos_state, pt_index)

             sigma = eos_state%dpdt / &
                  (eos_state%rho * eos_state%cp * eos_state%dpdr)

             xi_term = ZERO
             pres_term = ZERO
             do comp = 1, nspec
                xi_term = xi_term - &
                     eos_state%dhdX(comp)*rho_omegadot(i,j,k,comp)/eos_state%rho 

                pres_term = pres_term + &
                     eos_state%dpdX(comp)*rho_omegadot(i,j,k,comp)/eos_state%rho
             enddo

             Source(i,j,k) = (sigma/eos_state%rho) * &
                  ( rho_Hext(i,j,k) + rho_Hnuc(i,j,k) + thermal(i,j,k) ) &
                  + sigma*xi_term &
                  + pres_term/(eos_state%rho*eos_state%dpdr)

             select case (dm)
             case (1)
                r = i
             case (2)
                r = j
             case (3)
                r = k
             end select
             
             if (use_delta_gamma1_term .and. r < anelastic_cutoff_coord(n)) then
                if (r .eq. 0) then
                   gradp0 = (p0(r+1) - p0(r))/dx(dm)
                else if (r .eq. nr(n)-1) then
                   gradp0 = (p0(r) - p0(r-1))/dx(dm)
                else
                   gradp0 = HALF*(p0(r+1) - p0(r-1))/dx(dm)
                endif
                
                delta_gamma1(i,j,k) = eos_state%gam1 - gamma1bar(r)
                
                delta_gamma1_term(i,j,k) = (eos_state%gam1 - gamma1bar(r))*u(i,j,k,dm)* &
                     gradp0/(gamma1bar(r)*gamma1bar(r)*p0(r))
             else
                delta_gamma1_term(i,j,k) = ZERO
                delta_gamma1(i,j,k) = ZERO
             endif
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_S_cart

  
  subroutine make_S_2d_polar(lo,hi,Source,ng_sr,dg1_term,ng_dt,delta_gamma1, &
                            ng_dg,s,ng_s,u,ng_u,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,thermal,ng_th, &
                            gradp0_cart,ng_gp,p0_cart,ng_p0, &
                            gamma1bar_cart,ng_g1,normal,ng_n)

    use bl_constants_module
    use eos_module, only: eos, eos_input_rt
    use eos_type_module
    use network, only: nspec
    use variables, only: rho_comp, temp_comp, spec_comp
    use probin_module, only: use_delta_gamma1_term

    integer         , intent(in   ) :: lo(:),hi(:)
    integer         , intent(in   ) :: ng_sr,ng_dt,ng_dg,ng_s,ng_u, &
                                       ng_rw,ng_he, &
                                       ng_hn,ng_th,ng_gp,ng_p0,ng_g1,ng_n

    real (kind=dp_t), intent(  out) ::         Source(lo(1)-ng_sr:,lo(2)-ng_sr:)
    real (kind=dp_t), intent(  out) ::       dg1_term(lo(1)-ng_dt:,lo(2)-ng_dt:)
    real (kind=dp_t), intent(  out) ::   delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:) 
    real (kind=dp_t), intent(in   ) ::              s(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real (kind=dp_t), intent(in   ) ::              u(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real (kind=dp_t), intent(in   ) ::   rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,:)
    real (kind=dp_t), intent(in   ) ::       rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:)
    real (kind=dp_t), intent(in   ) ::       rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:)
    real (kind=dp_t), intent(in   ) ::        thermal(lo(1)-ng_th:,lo(2)-ng_th:)
    real (kind=dp_t), intent(in   ) ::    gradp0_cart(lo(1)-ng_gp:,lo(2)-ng_gp:)
    real (kind=dp_t), intent(in   ) ::        p0_cart(lo(1)-ng_p0:,lo(2)-ng_p0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar_cart(lo(1)-ng_g1:,lo(2)-ng_g1:)
    real (kind=dp_t), intent(in   ) ::         normal(lo(1)-ng_n:,lo(2)-ng_n:,:)

    !     Local variables
    integer         :: i, j, comp
    real(kind=dp_t) :: sigma, xi_term, pres_term

    integer :: pt_index(MAX_SPACEDIM)
    type(eos_t) :: eos_state

    real(kind=dp_t) :: Ut_dot_er

    Source = zero

    !$OMP PARALLEL DO PRIVATE(i,j,comp,sigma,xi_term,pres_term,eos_state,pt_index,Ut_dot_er)
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)

            eos_state%rho   = s(i,j,rho_comp)
            eos_state%T     = s(i,j,temp_comp)
            eos_state%xn(:) = s(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho

            pt_index(:) = (/i, j, -1/)

            ! dens, temp, and xmass are inputs
            call eos(eos_input_rt, eos_state, pt_index)

            sigma = eos_state%dpdt / &
                (eos_state%rho * eos_state%cp * eos_state%dpdr)

            xi_term = ZERO
            pres_term = ZERO
            do comp = 1, nspec
            xi_term = xi_term - &
                    eos_state%dhdX(comp)*rho_omegadot(i,j,comp)/eos_state%rho 

            pres_term = pres_term + &
                    eos_state%dpdX(comp)*rho_omegadot(i,j,comp)/eos_state%rho
            enddo

            Source(i,j) = (sigma/eos_state%rho) * &
                ( rho_Hext(i,j) + rho_Hnuc(i,j) + thermal(i,j) ) &
                + sigma*xi_term &
                + pres_term/(eos_state%rho*eos_state%dpdr)


            if (use_delta_gamma1_term) then
            delta_gamma1(i,j) = eos_state%gam1 - gamma1bar_cart(i,j)
            
            Ut_dot_er = &
                    u(i,j,1)*normal(i,j,1) + &
                    u(i,j,2)*normal(i,j,2) 
                    
            dg1_term(i,j) = delta_gamma1(i,j)*Ut_dot_er* &
                    gradp0_cart(i,j)/ &
                        (gamma1bar_cart(i,j)**2*p0_cart(i,j))

            else
            dg1_term(i,j) = ZERO
            delta_gamma1(i,j) = ZERO
            end if

        enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_S_2d_polar
  
  
  subroutine make_S_3d_sphr(lo,hi,Source,ng_sr,dg1_term,ng_dt,delta_gamma1, &
                            ng_dg,s,ng_s,u,ng_u,rho_omegadot,ng_rw,rho_Hnuc,ng_hn, &
                            rho_Hext,ng_he,thermal,ng_th, &
                            gradp0_cart,ng_gp,p0_cart,ng_p0, &
                            gamma1bar_cart,ng_g1,normal,ng_n)

    use bl_constants_module
    use eos_module, only: eos, eos_input_rt
    use eos_type_module
    use network, only: nspec
    use variables, only: rho_comp, temp_comp, spec_comp
    use probin_module, only: use_delta_gamma1_term

    integer         , intent(in   ) :: lo(:),hi(:)
    integer         , intent(in   ) :: ng_sr,ng_dt,ng_dg,ng_s,ng_u, &
                                       ng_rw,ng_he, &
                                       ng_hn,ng_th,ng_gp,ng_p0,ng_g1,ng_n

    real (kind=dp_t), intent(  out) ::         Source(lo(1)-ng_sr:,lo(2)-ng_sr:,lo(3)-ng_sr:)
    real (kind=dp_t), intent(  out) ::       dg1_term(lo(1)-ng_dt:,lo(2)-ng_dt:,lo(3)-ng_dt:)
    real (kind=dp_t), intent(  out) ::   delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:) 
    real (kind=dp_t), intent(in   ) ::              s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real (kind=dp_t), intent(in   ) ::              u(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real (kind=dp_t), intent(in   ) ::   rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,lo(3)-ng_rw:,:)
    real (kind=dp_t), intent(in   ) ::       rho_Hnuc(lo(1)-ng_hn:,lo(2)-ng_hn:,lo(3)-ng_hn:)
    real (kind=dp_t), intent(in   ) ::       rho_Hext(lo(1)-ng_he:,lo(2)-ng_he:,lo(3)-ng_he:)
    real (kind=dp_t), intent(in   ) ::        thermal(lo(1)-ng_th:,lo(2)-ng_th:,lo(3)-ng_th:)
    real (kind=dp_t), intent(in   ) ::    gradp0_cart(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:)
    real (kind=dp_t), intent(in   ) ::        p0_cart(lo(1)-ng_p0:,lo(2)-ng_p0:,lo(3)-ng_p0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar_cart(lo(1)-ng_g1:,lo(2)-ng_g1:,lo(3)-ng_g1:)
    real (kind=dp_t), intent(in   ) ::         normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)

    !     Local variables
    integer         :: i, j, k, comp
    real(kind=dp_t) :: sigma, xi_term, pres_term

    integer :: pt_index(MAX_SPACEDIM)
    type(eos_t) :: eos_state

    real(kind=dp_t) :: Ut_dot_er

    Source = zero

    !$OMP PARALLEL DO PRIVATE(i,j,k,comp,sigma,xi_term,pres_term,eos_state,pt_index,Ut_dot_er)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = s(i,j,k,rho_comp)
             eos_state%T     = s(i,j,k,temp_comp)
             eos_state%xn(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, eos_state, pt_index)

             sigma = eos_state%dpdt / &
                  (eos_state%rho * eos_state%cp * eos_state%dpdr)

             xi_term = ZERO
             pres_term = ZERO
             do comp = 1, nspec
                xi_term = xi_term - &
                     eos_state%dhdX(comp)*rho_omegadot(i,j,k,comp)/eos_state%rho 

                pres_term = pres_term + &
                     eos_state%dpdX(comp)*rho_omegadot(i,j,k,comp)/eos_state%rho
             enddo

             Source(i,j,k) = (sigma/eos_state%rho) * &
                  ( rho_Hext(i,j,k) + rho_Hnuc(i,j,k) + thermal(i,j,k) ) &
                  + sigma*xi_term &
                  + pres_term/(eos_state%rho*eos_state%dpdr)


             if (use_delta_gamma1_term) then
                delta_gamma1(i,j,k) = eos_state%gam1 - gamma1bar_cart(i,j,k)
                
                Ut_dot_er = &
                     u(i,j,k,1)*normal(i,j,k,1) + &
                     u(i,j,k,2)*normal(i,j,k,2) + &
                     u(i,j,k,3)*normal(i,j,k,3)                

                dg1_term(i,j,k) = delta_gamma1(i,j,k)*Ut_dot_er* &
                     gradp0_cart(i,j,k)/ &
                          (gamma1bar_cart(i,j,k)**2*p0_cart(i,j,k))

             else
                dg1_term(i,j,k) = ZERO
                delta_gamma1(i,j,k) = ZERO
             end if

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_S_3d_sphr

  subroutine correct_delta_gamma1_term_1d(lo,hi,delta_gamma1_term,ng_dt, &
                                          delta_gamma1,ng_dg, &
                                          gamma1bar,psi,p0)

    integer         , intent(in   ) :: lo(:), hi(:), ng_dt, ng_dg
    real (kind=dp_t), intent(inout) :: delta_gamma1_term(lo(1)-ng_dt:)
    real (kind=dp_t), intent(in   ) ::      delta_gamma1(lo(1)-ng_dg:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: psi(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    integer :: i

    do i = lo(1), hi(1)

       delta_gamma1_term(i) = delta_gamma1_term(i) &
            + delta_gamma1(i)*psi(i)/(gamma1bar(i)**2*p0(i))

    end do

  end subroutine correct_delta_gamma1_term_1d

  subroutine correct_delta_gamma1_term_2d(lo,hi,delta_gamma1_term,ng_dt, &
                                          delta_gamma1,ng_dg, &
                                          gamma1bar,psi,p0)

    integer         , intent(in   ) :: lo(:), hi(:), ng_dt, ng_dg
    real (kind=dp_t), intent(inout) :: delta_gamma1_term(lo(1)-ng_dt:,lo(2)-ng_dt:)
    real (kind=dp_t), intent(in   ) ::      delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: psi(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          delta_gamma1_term(i,j) = delta_gamma1_term(i,j) &
               + delta_gamma1(i,j)*psi(j)/(gamma1bar(j)**2*p0(j))

       end do
    end do

  end subroutine correct_delta_gamma1_term_2d

  subroutine correct_delta_gamma1_term_2d_polar(lo,hi,delta_gamma1_term,ng_dt, &
                                               delta_gamma1,ng_dg, &
                                               gamma1bar_cart,ng_g1, &
                                               psi_cart,ng_ps, &
                                               p0_cart,ng_p0)

    integer         , intent(in   ) :: lo(:), hi(:), ng_dt, ng_dg, ng_g1, ng_ps, ng_p0
    real (kind=dp_t), intent(inout) :: delta_gamma1_term(lo(1)-ng_dt:,lo(2)-ng_dt:)
    real (kind=dp_t), intent(in   ) ::    delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:)
    real (kind=dp_t), intent(in   ) ::  gamma1bar_cart(lo(1)-ng_g1:,lo(2)-ng_g1:)
    real (kind=dp_t), intent(in   ) ::         p0_cart(lo(1)-ng_p0:,lo(2)-ng_p0:)
    real (kind=dp_t), intent(in   ) ::        psi_cart(lo(1)-ng_ps:,lo(2)-ng_ps:)

    integer :: i, j

    !$OMP PARALLEL DO PRIVATE(i,j)
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)

            delta_gamma1_term(i,j) = delta_gamma1_term(i,j) &
                + delta_gamma1(i,j)*psi_cart(i,j)/(gamma1bar_cart(i,j)**2*p0_cart(i,j))

        end do
    end do
    !$OMP END PARALLEL DO

  end subroutine correct_delta_gamma1_term_2d_polar  
  
  
  subroutine correct_delta_gamma1_term_3d(lo,hi,delta_gamma1_term,ng_dt, &
                                          delta_gamma1,ng_dg, &
                                          gamma1bar,psi,p0)

    integer         , intent(in   ) :: lo(:), hi(:), ng_dt, ng_dg
    real (kind=dp_t), intent(inout) :: delta_gamma1_term(lo(1)-ng_dt:,lo(2)-ng_dt:,lo(3)-ng_dt:)
    real (kind=dp_t), intent(in   ) ::    delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: psi(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             delta_gamma1_term(i,j,k) = delta_gamma1_term(i,j,k) &
                  + delta_gamma1(i,j,k)*psi(k)/(gamma1bar(k)**2*p0(k))

          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine correct_delta_gamma1_term_3d

  subroutine correct_delta_gamma1_term_3d_sphr(lo,hi,delta_gamma1_term,ng_dt, &
                                               delta_gamma1,ng_dg, &
                                               gamma1bar_cart,ng_g1, &
                                               psi_cart,ng_ps, &
                                               p0_cart,ng_p0)

    integer         , intent(in   ) :: lo(:), hi(:), ng_dt, ng_dg, ng_g1, ng_ps, ng_p0
    real (kind=dp_t), intent(inout) :: delta_gamma1_term(lo(1)-ng_dt:,lo(2)-ng_dt:,lo(3)-ng_dt:)
    real (kind=dp_t), intent(in   ) ::    delta_gamma1(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:)
    real (kind=dp_t), intent(in   ) ::  gamma1bar_cart(lo(1)-ng_g1:,lo(2)-ng_g1:,lo(3)-ng_g1:)
    real (kind=dp_t), intent(in   ) ::         p0_cart(lo(1)-ng_p0:,lo(2)-ng_p0:,lo(3)-ng_p0:)
    real (kind=dp_t), intent(in   ) ::        psi_cart(lo(1)-ng_ps:,lo(2)-ng_ps:,lo(3)-ng_ps:)

    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             delta_gamma1_term(i,j,k) = delta_gamma1_term(i,j,k) &
                  + delta_gamma1(i,j,k)*psi_cart(i,j,k)/(gamma1bar_cart(i,j,k)**2*p0_cart(i,j,k))

          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine correct_delta_gamma1_term_3d_sphr

end module make_S_module
