! Compute the source term to the divergence constraint, S.

module make_S_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: make_S

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_S(nlevs,Source,delta_gamma1_term,delta_gamma1,state,u,rho_omegadot, &
                    rho_Hext,thermal,p0,gamma1bar,delta_gamma1_termbar,psi,dx,mla)

    use bl_constants_module
    use bl_prof_module
    use probin_module, only: use_delta_gamma1_term
    use ml_layout_module
    use average_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: Source(:)
    type(multifab) , intent(inout) :: delta_gamma1_term(:)
    type(multifab) , intent(inout) :: delta_gamma1(:)
    type(multifab) , intent(in   ) :: state(:)
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(in   ) :: rho_omegadot(:)
    type(multifab) , intent(in   ) :: rho_Hext(:)
    type(multifab) , intent(in   ) :: thermal(:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(inout) :: delta_gamma1_termbar(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(inout) :: mla
    
    real(kind=dp_t), pointer:: srcp(:,:,:,:),dgtp(:,:,:,:),sp(:,:,:,:),up(:,:,:,:)
    real(kind=dp_t), pointer:: tp(:,:,:,:),dgp(:,:,:,:)
    real(kind=dp_t), pointer:: omegap(:,:,:,:), hp(:,:,:,:)

    integer :: lo(state(1)%dim),hi(state(1)%dim),ng,dm
    integer :: i,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_S")

    ng = state(1)%ng
    dm = state(1)%dim

    do n = 1, nlevs
       do i = 1, state(n)%nboxes
          if ( multifab_remote(state(n), i) ) cycle
          srcp => dataptr(Source(n), i)
          dgtp     => dataptr(delta_gamma1_term(n), i)
          dgp    => dataptr(delta_gamma1(n), i)
          sp     => dataptr(state(n), i)
          up     => dataptr(u(n), i)
          omegap => dataptr(rho_omegadot(n), i)
          hp     => dataptr(rho_Hext(n), i)
          tp     => dataptr(thermal(n), i)
          lo = lwb(get_box(state(n), i))
          hi = upb(get_box(state(n), i))
          select case (dm)
          case (2)
             call make_S_2d(n,lo, hi, srcp(:,:,1,1), dgtp(:,:,1,1), dgp(:,:,1,1), &
                            sp(:,:,1,:), up(:,:,1,:), omegap(:,:,1,:), hp(:,:,1,1), &
                            tp(:,:,1,1), ng, p0(n,:), gamma1bar(n,:), dx(n,:))
          case (3)
             call make_S_3d(n,lo, hi, srcp(:,:,:,1), dgtp(:,:,:,1), dgp(:,:,:,1), &
                            sp(:,:,:,:), up(:,:,:,:), omegap(:,:,:,:), hp(:,:,:,1), &
                            tp(:,:,:,1), ng, p0(n,:), gamma1bar(n,:), dx(n,:))
          end select
       end do
    enddo

    if (use_delta_gamma1_term) then

       call average(mla,delta_gamma1_term,delta_gamma1_termbar,dx,1)
       
       do n = 1, nlevs
          do i = 1, state(n)%nboxes
             if ( multifab_remote(state(n), i) ) cycle
             dgtp   => dataptr(delta_gamma1_term(n), i)
             dgp    => dataptr(delta_gamma1(n), i)
             select case (dm)
             case (2)
                call correct_delta_gamma1_term_2d(lo,hi,dgtp(:,:,1,1),dgp(:,:,1,1), &
                                                  gamma1bar(n,:),psi(n,:), &
                                                  delta_gamma1_termbar(n,:),p0(n,:))
             case (3)
                call correct_delta_gamma1_term_3d(lo,hi,dgtp(:,:,:,1),dgp(:,:,:,1), &
                                                  gamma1bar(n,:),psi(n,:), &
                                                  delta_gamma1_termbar(n,:),p0(n,:))
             end select
          end do
       enddo
       
    end if

    call destroy(bpt)

   end subroutine make_S


   subroutine make_S_2d(n,lo,hi,Source,delta_gamma1_term,delta_gamma1,s,u, &
                        rho_omegadot,rho_Hext,thermal,ng,p0,gamma1bar,dx)

      use bl_constants_module
      use eos_module
      use variables, only: rho_comp, temp_comp, spec_comp
      use probin_module, only: use_delta_gamma1_term
      use geometry, only: anelastic_cutoff_coord, r_start_coord, r_end_coord

      integer         , intent(in   ) :: n,lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) :: Source(lo(1):,lo(2):)
      real (kind=dp_t), intent(  out) :: delta_gamma1_term(lo(1):,lo(2):)
      real (kind=dp_t), intent(  out) :: delta_gamma1(lo(1):,lo(2):)
      real (kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(in   ) :: u(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1):,lo(2):,:)
      real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1):,lo(2):)
      real (kind=dp_t), intent(in   ) :: thermal(lo(1)-1:,lo(2)-1:)
      real (kind=dp_t), intent(in   ) :: p0(0:)
      real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
      real (kind=dp_t), intent(in   ) :: dx(:)

!     Local variables
      integer         :: i, j, comp
      real(kind=dp_t) :: sigma, react_term, pres_term, gradp0

      Source = zero

      do_diag = .false.

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           den_eos(1) = s(i,j,rho_comp)
           temp_eos(1) = s(i,j,temp_comp)
           xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
           
           ! dens, temp, and xmass are inputs
           call eos(eos_input_rt, den_eos, temp_eos, &
                    npts, nspec, &
                    xn_eos, &
                    p_eos, h_eos, e_eos, & 
                    cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                    dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                    dpdX_eos, dhdX_eos, &
                    gam1_eos, cs_eos, s_eos, &
                    dsdt_eos, dsdr_eos, &
                    do_diag)

           sigma = dpdt_eos(1) / (den_eos(1) * cp_eos(1) * dpdr_eos(1))

           react_term = ZERO
           pres_term = ZERO
           do comp = 1, nspec
              react_term = react_term - &
                   (dhdX_eos(1,comp) + ebin(comp))*rho_omegadot(i,j,comp)/den_eos(1)

              pres_term = pres_term + &
                   dpdX_eos(1,comp)*rho_omegadot(i,j,comp)/den_eos(1)
           enddo

           Source(i,j) = (sigma/den_eos(1)) * ( rho_Hext(i,j) + thermal(i,j) ) &
                        + sigma*react_term &
                        + pres_term/(den_eos(1)*dpdr_eos(1))

           if (use_delta_gamma1_term .and. j < anelastic_cutoff_coord(n)) then
              if (j .eq. r_start_coord(n)) then
                 gradp0 = (p0(j+1) - p0(j))/dx(2)
              else if (j .eq. r_end_coord(n)) then
                 gradp0 = (p0(j) - p0(j-1))/dx(2)
              else
                 gradp0 = HALF*(p0(j+1) - p0(j-1))/dx(2)
              endif

              delta_gamma1(i,j) = gam1_eos(1) - gamma1bar(j)

              delta_gamma1_term(i,j) = &
                   (gam1_eos(1) - gamma1bar(j))*u(i,j,2)* &
                   gradp0/(gamma1bar(j)*gamma1bar(j)*p0(j))
           else
              delta_gamma1_term(i,j) = ZERO
              delta_gamma1(i,j) = ZERO
           endif

        enddo
      enddo
 
   end subroutine make_S_2d

   subroutine make_S_3d(n,lo,hi,Source,delta_gamma1_term,delta_gamma1,s,u, &
                        rho_omegadot,rho_Hext,thermal,ng,p0,gamma1bar,dx)

      use bl_constants_module
      use eos_module
      use geometry, only: spherical
      use variables, only: rho_comp, temp_comp, spec_comp
      use probin_module, only: use_delta_gamma1_term
      use geometry, only: anelastic_cutoff_coord, r_start_coord, r_end_coord
     
      integer         , intent(in   ) :: n,lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) :: Source(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(  out) :: delta_gamma1_term(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(  out) :: delta_gamma1(lo(1):,lo(2):,lo(3):) 
      real (kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(in   ) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1):,lo(2):,lo(3):,:)
      real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1):,lo(2):,lo(3):)
      real (kind=dp_t), intent(in   ) :: thermal(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real (kind=dp_t), intent(in   ) :: p0(0:)
      real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
      real (kind=dp_t), intent(in   ) :: dx(:)

!     Local variables
      integer         :: i, j, k, comp
      real(kind=dp_t) :: sigma, react_term, pres_term, gradp0

      Source = zero

      do_diag = .false.

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              den_eos(1) = s(i,j,k,rho_comp)
              temp_eos(1) = s(i,j,k,temp_comp)
              xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

              ! dens, temp, and xmass are inputs
              call eos(eos_input_rt, den_eos, temp_eos, &
                       npts, nspec, &
                       xn_eos, &
                       p_eos, h_eos, e_eos, & 
                       cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                       dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                       dpdX_eos, dhdX_eos, &
                       gam1_eos, cs_eos, s_eos, &
                       dsdt_eos, dsdr_eos, &
                       do_diag)

              sigma = dpdt_eos(1) / (den_eos(1) * cp_eos(1) * dpdr_eos(1))

              react_term = ZERO
              pres_term = ZERO
              do comp = 1, nspec
                 react_term = react_term - &
                      (dhdX_eos(1,comp) + ebin(comp))*rho_omegadot(i,j,k,comp)/den_eos(1)

                 pres_term = pres_term + &
                      dpdX_eos(1,comp)*rho_omegadot(i,j,k,comp)/den_eos(1)
              enddo

              Source(i,j,k) = (sigma/den_eos(1)) * ( rho_Hext(i,j,k) + thermal(i,j,k) ) &
                           + sigma*react_term &
                           + pres_term/(den_eos(1)*dpdr_eos(1))


              if (use_delta_gamma1_term .and. spherical .eq. 1) then
                 call bl_error("ERROR: use_delta_gamma1_term not implemented for spherical in make_S")
              end if

              if (use_delta_gamma1_term .and. k < anelastic_cutoff_coord(n)) then
                 if (k .eq. r_start_coord(n)) then
                    gradp0 = (p0(k+1) - p0(k))/dx(3)
                 else if (k .eq. r_end_coord(n)) then
                    gradp0 = (p0(k) - p0(k-1))/dx(3)
                 else
                    gradp0 = HALF*(p0(k+1) - p0(k-1))/dx(3)
                 endif
                 
                 delta_gamma1(i,j,k) = gam1_eos(1) - gamma1bar(k)
                 
                 delta_gamma1_term(i,j,k) = &
                      (gam1_eos(1) - gamma1bar(k))*u(i,j,k,3)* &
                      gradp0/(gamma1bar(k)*gamma1bar(k)*p0(k))
              else
                 delta_gamma1_term(i,j,k) = ZERO
                 delta_gamma1(i,j,k) = ZERO
              endif

           enddo
        enddo
      enddo
 
   end subroutine make_S_3d

   subroutine correct_delta_gamma1_term_2d(lo,hi,delta_gamma1_term,delta_gamma1, &
                                           gamma1bar,psi,delta_gamma1_termbar,p0)

     integer         , intent(in   ) :: lo(:), hi(:)
     real (kind=dp_t), intent(inout) :: delta_gamma1_term(lo(1):,lo(2):)
     real (kind=dp_t), intent(in   ) :: delta_gamma1(lo(1):,lo(2):)
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

   subroutine correct_delta_gamma1_term_3d(lo,hi,delta_gamma1_term,delta_gamma1, &
                                           gamma1bar,psi,delta_gamma1_termbar,p0)

     integer         , intent(in   ) :: lo(:), hi(:)
     real (kind=dp_t), intent(inout) :: delta_gamma1_term(lo(1):,lo(2):,lo(3):)
     real (kind=dp_t), intent(in   ) :: delta_gamma1(lo(1):,lo(2):,lo(3):)
     real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
     real (kind=dp_t), intent(in   ) :: psi(0:)
     real (kind=dp_t), intent(in   ) :: delta_gamma1_termbar(0:)
     real (kind=dp_t), intent(in   ) :: p0(0:)
     
     integer :: i, j, k
     
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
           
              delta_gamma1_term(i,j,k) = delta_gamma1_term(i,j,k) - delta_gamma1_termbar(k) &
                   + delta_gamma1(i,j,k)*psi(k)/(gamma1bar(k)**2*p0(k))
           
           end do
        end do
     end do
     
   end subroutine correct_delta_gamma1_term_3d

end module make_S_module
