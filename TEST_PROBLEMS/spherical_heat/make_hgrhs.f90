! Create the righthand side to the elliptic equation that is solved in
! the final project step, \beta * (S - \bar{S}).  This quantity is 
! node-centered.  The computation is done is two steps -- first we 
! compute it on the cell-centers and then we average to the nodes.
!
! Note, we include the delta_gamma1_term here, to (possibly) account for
! the effect of replacing \Gamma_1 by {\Gamma_1}_0 in the constraint
! equation (see paper III).

module hgrhs_module

  use bl_types
  use multifab_module
  
  implicit none

  private

  public :: make_hgrhs, correct_hgrhs

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine make_hgrhs(the_bc_tower,mla,hgrhs,Source,delta_gamma1_term,Sbar, &
                        div_coeff,dx)

    use define_bc_module
    use ml_layout_module
    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical, dm, nlevs
    use fill_3d_module
    use variables, only: foextrap_comp
    use ml_restriction_module
    use multifab_fill_ghost_module
    use multifab_physbc_module
    
    type(bc_tower),  intent(in   ) :: the_bc_tower
    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: hgrhs(:)
    type(multifab) , intent(in   ) :: Source(:)
    type(multifab) , intent(in   ) :: delta_gamma1_term(:)
    real(kind=dp_t), intent(in   ) :: Sbar(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    
    type(multifab) :: rhs_cc(nlevs)
    type(multifab) :: Sbar_cart(nlevs)
    type(multifab) :: div_coeff_cart(nlevs)

    real(kind=dp_t), pointer:: hp(:,:,:,:),gp(:,:,:,:),rp(:,:,:,:)
    real(kind=dp_t), pointer:: dp(:,:,:,:),sp(:,:,:,:),sbp(:,:,:,:)
    integer :: lo(dm),hi(dm),i,n
    integer :: ng_rh,ng_sr,ng_dg,ng_dc,ng_sb,ng_hg

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_hgrhs")
    
    if(spherical .eq. 1) then
       do n = 1, nlevs
          call multifab_build(Sbar_cart(n),Source(n)%la,1,0)
          call multifab_build(div_coeff_cart(n),Source(n)%la,1,0)
          call setval(Sbar_cart(n),ZERO,all=.true.)
          call setval(div_coeff_cart(n),ZERO,all=.true.)
       end do
    end if
    
    if (spherical .eq. 1) then
       call put_1d_array_on_cart(div_coeff,div_coeff_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla)
       call put_1d_array_on_cart(Sbar,Sbar_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla)
    end if

    do n = 1, nlevs
       call multifab_build(rhs_cc(n),Source(n)%la,1,1)
       call setval(rhs_cc(n),ZERO,all=.true.)
    end do

    ng_rh = rhs_cc(1)%ng
    ng_sr = Source(1)%ng
    ng_dg = delta_gamma1_term(1)%ng
    ng_dc = div_coeff_cart(1)%ng
    ng_sb = Sbar_cart(1)%ng
    ng_hg = hgrhs(1)%ng

    do n = 1, nlevs
       do i = 1, nfabs(Source(n))
          rp => dataptr(rhs_cc(n), i)
          sp => dataptr(Source(n), i)
          gp => dataptr(delta_gamma1_term(n), i)
          lo =  lwb(get_box(Source(n), i))
          hi =  upb(get_box(Source(n), i))
          select case (dm)
          case (1)
             call make_rhscc_1d(lo,hi,rp(:,1,1,1),ng_rh,sp(:,1,1,1),ng_sr, &
                                gp(:,1,1,1),ng_dg,Sbar(n,:),div_coeff(n,:))
          case (2)
             call make_rhscc_2d(lo,hi,rp(:,:,1,1),ng_rh,sp(:,:,1,1),ng_sr, &
                                gp(:,:,1,1),ng_dg,Sbar(n,:),div_coeff(n,:))
          case (3)
             if (spherical .eq. 1) then
                dp => dataptr(div_coeff_cart(n), i)
                sbp => dataptr(Sbar_cart(n), i)
                call make_rhscc_3d_sphr(lo,hi,rp(:,:,:,1),ng_rh,sp(:,:,:,1),ng_sr, &
                                        gp(:,:,:,1),ng_dg,sbp(:,:,:,1),ng_sb, &
                                        dp(:,:,:,1),ng_dc)
             else
                call make_rhscc_3d_cart(lo,hi,rp(:,:,:,1),ng_rh,sp(:,:,:,1),ng_sr, &
                                        gp(:,:,:,1),ng_dg,Sbar(n,:),div_coeff(n,:))
             end if
          end select
       end do
    end do

    if (nlevs .eq. 1) then
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(rhs_cc(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rhs_cc(nlevs),1,foextrap_comp,1, &
                            the_bc_tower%bc_tower_array(nlevs))
    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(rhs_cc(n-1),rhs_cc(n),mla%mba%rr(n-1,:))
 
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(rhs_cc(n),rhs_cc(n-1),1,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n),1,foextrap_comp,1, &
                                         fill_crse_input=.false.)
       end do

    end if
       
    do n=1,nlevs
       call setval(hgrhs(n),ZERO,all=.true.)
       do i = 1, nfabs(Source(n))
          hp => dataptr(hgrhs(n), i)
          rp => dataptr(rhs_cc(n), i)
          lo =  lwb(get_box(Source(n), i))
          hi =  upb(get_box(Source(n), i))
          select case (dm)
          case (1)
             call make_hgrhs_1d(lo,hi,hp(:,1,1,1),ng_hg,rp(:,1,1,1),ng_rh)
          case (2)
             call make_hgrhs_2d(lo,hi,hp(:,:,1,1),ng_hg,rp(:,:,1,1),ng_rh)
          case (3)
             call make_hgrhs_3d(lo,hi,hp(:,:,:,1),ng_hg,rp(:,:,:,1),ng_rh)
          end select
       end do
    end do ! end loop over levels
    
    do n = 1, nlevs
       call destroy(rhs_cc(n))
       if (spherical .eq. 1) then
          call destroy(Sbar_cart(n))
          call destroy(div_coeff_cart(n))
       end if
    end do

    call destroy(bpt)
    
  end subroutine make_hgrhs
  
  subroutine make_rhscc_1d(lo,hi,rhs_cc,ng_rh,Source,ng_sr,delta_gamma1_term,ng_dg,Sbar, &
                           div_coeff)

    integer         , intent(in   ) :: lo(:), hi(:), ng_rh, ng_sr, ng_dg
    real (kind=dp_t), intent(  out) ::            rhs_cc(lo(1)-ng_rh:)
    real (kind=dp_t), intent(in   ) ::            Source(lo(1)-ng_sr:)
    real (kind=dp_t), intent(in   ) :: delta_gamma1_term(lo(1)-ng_dg:)
    real (kind=dp_t), intent(in   ) ::      Sbar(0:)
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    
    ! Local variables
    integer :: i
    
    do i = lo(1),hi(1)
       rhs_cc(i) = div_coeff(i) * (Source(i) - Sbar(i) + delta_gamma1_term(i))
    end do
    
  end subroutine make_rhscc_1d
  
  subroutine make_rhscc_2d(lo,hi,rhs_cc,ng_rh,Source,ng_sr,delta_gamma1_term,ng_dg,Sbar, &
                           div_coeff)

    integer         , intent(in   ) :: lo(:), hi(:), ng_rh, ng_sr, ng_dg
    real (kind=dp_t), intent(  out) ::            rhs_cc(lo(1)-ng_rh:,lo(2)-ng_rh:)
    real (kind=dp_t), intent(in   ) ::            Source(lo(1)-ng_sr:,lo(2)-ng_sr:)
    real (kind=dp_t), intent(in   ) :: delta_gamma1_term(lo(1)-ng_dg:,lo(2)-ng_dg:)  
    real (kind=dp_t), intent(in   ) ::      Sbar(0:)
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    
    ! Local variables
    integer :: i, j
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          rhs_cc(i,j) = div_coeff(j) * (Source(i,j) - Sbar(j) + delta_gamma1_term(i,j))
       end do
    end do
    
  end subroutine make_rhscc_2d
  
  subroutine make_rhscc_3d_cart(lo,hi,rhs_cc,ng_rh,Source,ng_sr,delta_gamma1_term,ng_dg, &
                                Sbar,div_coeff)

    integer         , intent(in   ) :: lo(:), hi(:), ng_rh, ng_sr, ng_dg
    real (kind=dp_t), intent(  out) ::         rhs_cc(lo(1)-ng_rh:,lo(2)-ng_rh:,lo(3)-ng_rh:)
    real (kind=dp_t), intent(in   ) ::         Source(lo(1)-ng_sr:,lo(2)-ng_sr:,lo(3)-ng_sr:)
    real (kind=dp_t), intent(in) :: delta_gamma1_term(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:)
    real (kind=dp_t), intent(in   ) ::      Sbar(0:)
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    
    ! Local variables
    integer :: i,j,k
    
!$omp parallel do private(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhs_cc(i,j,k) = div_coeff(k) * (Source(i,j,k) - Sbar(k) + &
                  delta_gamma1_term(i,j,k))
          end do
       end do
    end do
!$omp end parallel do
    
  end subroutine make_rhscc_3d_cart
   
  subroutine make_rhscc_3d_sphr(lo,hi,rhs_cc,ng_rh,Source,ng_sr,delta_gamma1_term,ng_dg, &
                                Sbar_cart,ng_sb,div_coeff_cart,ng_dc)

    integer         , intent(in   ) :: lo(:), hi(:), ng_rh, ng_sr, ng_dg, ng_sb, ng_dc
    real (kind=dp_t), intent(  out) ::         rhs_cc(lo(1)-ng_rh:,lo(2)-ng_rh:,lo(3)-ng_rh:)
    real (kind=dp_t), intent(in   ) ::         Source(lo(1)-ng_sr:,lo(2)-ng_sr:,lo(3)-ng_sr:)
    real (kind=dp_t), intent(in) :: delta_gamma1_term(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:)
    real (kind=dp_t), intent(in   ) ::      Sbar_cart(lo(1)-ng_sb:,lo(2)-ng_sb:,lo(3)-ng_sb:)
    real (kind=dp_t), intent(in   ) :: div_coeff_cart(lo(1)-ng_dc:,lo(2)-ng_dc:,lo(3)-ng_dc:)
    
    ! Local variables
    integer :: i,j,k
    
!$omp parallel do private(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhs_cc(i,j,k) = div_coeff_cart(i,j,k) * (Source(i,j,k) - Sbar_cart(i,j,k) + &
                  delta_gamma1_term(i,j,k))
             
          end do
       end do
    end do
!$omp end parallel do
    
  end subroutine make_rhscc_3d_sphr
  
  subroutine make_hgrhs_1d(lo,hi,rhs,ng_hg,rhs_cc,ng_rh)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), ng_hg, ng_rh
    real (kind=dp_t), intent(  out) ::    rhs(lo(1)-ng_hg:)
    real (kind=dp_t), intent(in   ) :: rhs_cc(lo(1)-ng_rh:)
    
    ! Local variables
    integer :: i
    
    do i = lo(1), hi(1)+1
       rhs(i) = HALF * ( rhs_cc(i) + rhs_cc(i-1) )
    end do
    
  end subroutine make_hgrhs_1d
  
  subroutine make_hgrhs_2d(lo,hi,rhs,ng_hg,rhs_cc,ng_rh)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), ng_hg, ng_rh
    real (kind=dp_t), intent(  out) ::    rhs(lo(1)-ng_hg:,lo(2)-ng_hg:)  
    real (kind=dp_t), intent(in   ) :: rhs_cc(lo(1)-ng_rh:,lo(2)-ng_rh:)
    
    ! Local variables
    integer :: i, j
    
    do j = lo(2),hi(2)+1
       do i = lo(1), hi(1)+1
          rhs(i,j) = FOURTH * ( rhs_cc(i,j  ) + rhs_cc(i-1,j  ) &
                              + rhs_cc(i,j-1) + rhs_cc(i-1,j-1) )
       end do
    end do
    
  end subroutine make_hgrhs_2d
  
  subroutine make_hgrhs_3d(lo,hi,rhs,ng_hg,rhs_cc,ng_rh)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), ng_hg, ng_rh
    real (kind=dp_t), intent(  out) ::    rhs(lo(1)-ng_hg:,lo(2)-ng_hg:,lo(3)-ng_hg:)  
    real (kind=dp_t), intent(in   ) :: rhs_cc(lo(1)-ng_rh:,lo(2)-ng_rh:,lo(3)-ng_rh:)  
    
    ! Local variables
    integer :: i, j,k
    
!$omp parallel do private(i,j,k)
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)+1
!             rhs(i,j,k) = EIGHTH * ( rhs_cc(i,j  ,k-1) + rhs_cc(i-1,j  ,k-1) &
!                  +rhs_cc(i,j-1,k-1) + rhs_cc(i-1,j-1,k-1) &
!                  +rhs_cc(i,j  ,k  ) + rhs_cc(i-1,j  ,k  ) &
!                  +rhs_cc(i,j-1,k  ) + rhs_cc(i-1,j-1,k  ) )
! hack to force U-tilde to be zero
             rhs(i,j,k) = ZERO
          end do
       end do
    end do
!$omp end parallel do
    
  end subroutine make_hgrhs_3d

  subroutine correct_hgrhs(the_bc_tower,mla,rho0,hgrhs,div_coeff,dx,dt,gamma1bar,p0, &
                           delta_p_term)

    use define_bc_module
    use ml_layout_module
    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical, dm, nlevs
    use fill_3d_module
    use variables, only: foextrap_comp, rho_comp
    use ml_restriction_module
    use multifab_fill_ghost_module
    use multifab_physbc_module
    use probin_module, only: nodal
    
    type(bc_tower),  intent(in   ) :: the_bc_tower
    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: hgrhs(:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    type(multifab) , intent(in   ) :: delta_p_term(:)
    
    type(multifab) :: correction_cc(nlevs)
    type(multifab) :: correction_nodal(nlevs)

    type(multifab) :: gamma1bar_cart(nlevs)
    type(multifab) :: p0_cart(nlevs)
    type(multifab) :: div_coeff_cart(nlevs)
    type(multifab) :: rho0_cart(nlevs)

    real(kind=dp_t), pointer :: ptp(:,:,:,:), ccp(:,:,:,:), cnp(:,:,:,:)
    real(kind=dp_t), pointer :: gbp(:,:,:,:), p0p(:,:,:,:), dcp(:,:,:,:)
    real(kind=dp_t), pointer :: r0p(:,:,:,:)
    integer :: lo(dm),hi(dm),i,n
    integer :: ng_cc, ng_dp, ng_gb, ng_p0, ng_dc, ng_r0, ng_cn

    type(bl_prof_timer), save :: bpt

    call build(bpt, "correct_hgrhs")
    
    if(spherical .eq. 1) then
       do n = 1, nlevs
          call multifab_build(gamma1bar_cart(n),delta_p_term(n)%la,1,0)
          call multifab_build(p0_cart(n),delta_p_term(n)%la,1,0)
          call multifab_build(div_coeff_cart(n),delta_p_term(n)%la,1,0)
          call multifab_build(rho0_cart(n),delta_p_term(n)%la,1,0)
          call setval(gamma1bar_cart(n),ZERO,all=.true.)
          call setval(p0_cart(n),ZERO,all=.true.)
          call setval(div_coeff_cart(n),ZERO,all=.true.)
          call setval(rho0_cart(n),ZERO,all=.true.)
       end do
    end if
    
    if (spherical .eq. 1) then
       call put_1d_array_on_cart(gamma1bar,gamma1bar_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla)
       call put_1d_array_on_cart(p0,p0_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla)
       call put_1d_array_on_cart(div_coeff,div_coeff_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla)
       call put_1d_array_on_cart(rho0,rho0_cart,dm+rho_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla)

    end if

    do n = 1, nlevs
       call multifab_build(correction_cc(n),delta_p_term(n)%la,1,1)
       call setval(correction_cc(n),ZERO,all=.true.)
       call multifab_build(correction_nodal(n),delta_p_term(n)%la,1,0,nodal)
       call setval(correction_nodal(n),ZERO,all=.true.)
    end do

    ng_cc = correction_cc(1)%ng
    ng_dp = delta_p_term(1)%ng
    ng_gb = gamma1bar_cart(1)%ng
    ng_p0 = p0_cart(1)%ng
    ng_dc = div_coeff_cart(1)%ng
    ng_r0 = rho0_cart(1)%ng
    ng_cn = correction_nodal(1)%ng

    do n = 1, nlevs
       do i = 1, nfabs(delta_p_term(n))
          ccp => dataptr(correction_cc(n), i)
          ptp => dataptr(delta_p_term(n), i)
          lo =  lwb(get_box(delta_p_term(n), i))
          hi =  upb(get_box(delta_p_term(n), i))
          select case (dm)
          case (1)
             call create_correction_cc_1d(n,lo,hi,ccp(:,1,1,1),ng_cc,ptp(:,1,1,1),ng_dp, &
                                          div_coeff(n,:),gamma1bar(n,:),p0(n,:),dt)
          case (2)
             call create_correction_cc_2d(n,lo,hi,ccp(:,:,1,1),ng_cc,ptp(:,:,1,1),ng_dp, &
                                          div_coeff(n,:),gamma1bar(n,:),p0(n,:),dt)
          case (3)
             if (spherical .eq. 1) then
                gbp  => dataptr(gamma1bar_cart(n),i)
                p0p  => dataptr(p0_cart(n),i)
                dcp  => dataptr(div_coeff_cart(n),i)
                r0p  => dataptr(rho0_cart(n),i)
                call create_correction_cc_3d_sphr(lo,hi,ccp(:,:,:,1),ng_cc,ptp(:,:,:,1), &
                                                  ng_dp,dcp(:,:,:,1),ng_dc,gbp(:,:,:,1), &
                                                  ng_gb,p0p(:,:,:,1),ng_p0,r0p(:,:,:,1), &
                                                  ng_r0,dt)

             else
                call create_correction_cc_3d_cart(n,lo,hi,ccp(:,:,:,1),ng_cc,ptp(:,:,:,1), &
                                                  ng_dp,div_coeff(n,:),gamma1bar(n,:), &
                                                  p0(n,:),dt)
             end if
          end select
       end do
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(correction_cc(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(correction_cc(nlevs),1,foextrap_comp,1, &
                            the_bc_tower%bc_tower_array(nlevs))
    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(correction_cc(n-1),correction_cc(n),mla%mba%rr(n-1,:))
 
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(correction_cc(n),correction_cc(n-1),1, &
                                         mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n),1,foextrap_comp,1, &
                                         fill_crse_input=.false.)
       end do

    end if
       
    do n=1,nlevs
       call setval(correction_nodal(n),ZERO,all=.true.)
       do i = 1, nfabs(delta_p_term(n))
          cnp => dataptr(correction_nodal(n), i)
          ccp => dataptr(correction_cc(n), i)
          lo =  lwb(get_box(delta_p_term(n), i))
          hi =  upb(get_box(delta_p_term(n), i))
          select case (dm)
          case (1)
             call create_correction_nodal_1d(lo,hi,cnp(:,1,1,1),ng_cn,ccp(:,1,1,1),ng_cc)
          case (2)
             call create_correction_nodal_2d(lo,hi,cnp(:,:,1,1),ng_cn,ccp(:,:,1,1),ng_cc)
          case (3)
             call create_correction_nodal_3d(lo,hi,cnp(:,:,:,1),ng_cn,ccp(:,:,:,1),ng_cc)
          end select
       end do
    end do ! end loop over levels
    
    ! add correction term
    do n=1,nlevs
       call multifab_plus_plus_c(hgrhs(n),1,correction_nodal(n),1,1)
    end do

    do n = 1, nlevs
       call destroy(correction_cc(n))
       call destroy(correction_nodal(n))
       if (spherical .eq. 1) then
          call destroy(gamma1bar_cart(n))
          call destroy(p0_cart(n))
          call destroy(div_coeff_cart(n))
          call destroy(rho0_cart(n))
       end if
    end do

    call destroy(bpt)
    
  end subroutine correct_hgrhs
  
  subroutine create_correction_cc_1d(n,lo,hi,correction_cc,ng_cc,delta_p_term,ng_dp, &
                                     div_coeff,gamma1bar,p0,dt)

    use probin_module, only: dpdt_factor
    use geometry, only: base_cutoff_density_coord

    integer         , intent(in   ) :: n, lo(:), hi(:), ng_cc, ng_dp
    real (kind=dp_t), intent(  out) :: correction_cc(lo(1)-ng_cc:)
    real (kind=dp_t), intent(in   ) ::  delta_p_term(lo(1)-ng_dp:)
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)    
    real (kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: i
    real(kind=dp_t) :: correction_factor
    
    do i = lo(1),hi(1)
       if(i .lt. base_cutoff_density_coord(n)) then
          correction_factor = div_coeff(i)*(dpdt_factor/(gamma1bar(i)*p0(i))) / dt
       else
          correction_factor = 0.0d0
       end if
       correction_cc(i) = correction_factor*delta_p_term(i)
    end do
    
  end subroutine create_correction_cc_1d
  
  subroutine create_correction_cc_2d(n,lo,hi,correction_cc,ng_cc,delta_p_term,ng_dp, &
                                     div_coeff,gamma1bar,p0,dt)

    use probin_module, only: dpdt_factor
    use geometry, only: base_cutoff_density_coord

    integer         , intent(in   ) :: n, lo(:), hi(:), ng_cc, ng_dp
    real (kind=dp_t), intent(  out) :: correction_cc(lo(1)-ng_cc:,lo(2)-ng_cc:)
    real (kind=dp_t), intent(in   ) ::  delta_p_term(lo(1)-ng_dp:,lo(2)-ng_dp:)
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)    
    real (kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: i, j
    real(kind=dp_t) :: correction_factor
    
    do j = lo(2),hi(2)
       if(j .lt. base_cutoff_density_coord(n)) then
          correction_factor = div_coeff(j)*(dpdt_factor/(gamma1bar(j)*p0(j))) / dt
       else
          correction_factor = 0.0d0
       end if
       do i = lo(1),hi(1)
          correction_cc(i,j) = correction_factor*delta_p_term(i,j)
       end do
    end do
    
  end subroutine create_correction_cc_2d

  subroutine create_correction_cc_3d_cart(n,lo,hi,correction_cc,ng_cc,delta_p_term,ng_dp, &
                                          div_coeff,gamma1bar,p0,dt)

    use probin_module, only: dpdt_factor
    use geometry, only: base_cutoff_density_coord

    integer         , intent(in   ) :: n, lo(:), hi(:),ng_cc,ng_dp
    real (kind=dp_t), intent(  out) :: correction_cc(lo(1)-ng_cc:,lo(2)-ng_cc:,lo(3)-ng_cc:)
    real (kind=dp_t), intent(in   ) ::  delta_p_term(lo(1)-ng_dp:,lo(2)-ng_dp:,lo(3)-ng_dp:)
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)    
    real (kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: i, j, k
    real(kind=dp_t) :: correction_factor
    
!$omp parallel do private(i,j,k,correction_factor)
    do k = lo(3),hi(3)
       if(k .lt. base_cutoff_density_coord(n)) then
          correction_factor = div_coeff(k)*(dpdt_factor/(gamma1bar(k)*p0(k))) / dt
       else
          correction_factor = 0.0d0
       end if
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             correction_cc(i,j,k) = correction_factor*delta_p_term(i,j,k)
          end do
       end do
    end do
!$omp end parallel do
    
  end subroutine create_correction_cc_3d_cart

  subroutine create_correction_cc_3d_sphr(lo,hi,correction_cc,ng_cc,delta_p_term,ng_dp, &
                                          div_coeff_cart,ng_dc,gamma1bar_cart,ng_gb, &
                                          p0_cart,ng_p0,rho0_cart,ng_r0,dt)

    use probin_module, only: dpdt_factor, base_cutoff_density

    integer         , intent(in   ) :: lo(:),hi(:),ng_cc,ng_dp,ng_dc,ng_gb,ng_p0,ng_r0
    real (kind=dp_t), intent(  out) ::  correction_cc(lo(1)-ng_cc:,lo(2)-ng_cc:,lo(3)-ng_cc:)
    real (kind=dp_t), intent(in   ) ::   delta_p_term(lo(1)-ng_dp:,lo(2)-ng_dp:,lo(3)-ng_dp:)
    real (kind=dp_t), intent(in   ) :: div_coeff_cart(lo(1)-ng_dc:,lo(2)-ng_dc:,lo(3)-ng_dc:)
    real (kind=dp_t), intent(in   ) :: gamma1bar_cart(lo(1)-ng_gb:,lo(2)-ng_gb:,lo(3)-ng_gb:)
    real (kind=dp_t), intent(in   ) ::        p0_cart(lo(1)-ng_p0:,lo(2)-ng_p0:,lo(3)-ng_p0:)
    real (kind=dp_t), intent(in   ) ::      rho0_cart(lo(1)-ng_r0:,lo(2)-ng_r0:,lo(3)-ng_r0:)
    real (kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: i, j, k
    real(kind=dp_t) :: correction_factor
    
!$omp parallel do private(i,j,k,correction_factor)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if(rho0_cart(i,j,k) .gt. base_cutoff_density) then
                correction_factor = div_coeff_cart(i,j,k) * &
                     (dpdt_factor/(gamma1bar_cart(i,j,k)*p0_cart(i,j,k))) / dt
             else
                correction_factor = 0.0d0
             end if
             correction_cc(i,j,k) = correction_factor*delta_p_term(i,j,k)
          end do
       end do
    end do
!$omp end parallel do
    
  end subroutine create_correction_cc_3d_sphr
  
  subroutine create_correction_nodal_1d(lo,hi,correction_nodal,ng_cn,correction_cc,ng_cc)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), ng_cn, ng_cc
    real (kind=dp_t), intent(  out) :: correction_nodal(lo(1)-ng_cn:)
    real (kind=dp_t), intent(in   ) ::    correction_cc(lo(1)-ng_cc:)
    
    ! Local variables
    integer :: i
    
    do i = lo(1), hi(1)+1
       correction_nodal(i) = HALF * ( correction_cc(i) + correction_cc(i-1) )
    end do
    
  end subroutine create_correction_nodal_1d

  
  subroutine create_correction_nodal_2d(lo,hi,correction_nodal,ng_cn,correction_cc,ng_cc)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), ng_cn, ng_cc
    real (kind=dp_t), intent(  out) :: correction_nodal(lo(1)-ng_cn:,lo(2)-ng_cn:)  
    real (kind=dp_t), intent(in   ) ::    correction_cc(lo(1)-ng_cc:,lo(2)-ng_cc:)
    
    ! Local variables
    integer :: i, j
    
    do j = lo(2),hi(2)+1
       do i = lo(1), hi(1)+1
          correction_nodal(i,j) = FOURTH * &
               (correction_cc(i,j  ) + correction_cc(i-1,j  ) &
               +correction_cc(i,j-1) + correction_cc(i-1,j-1) )
       end do
    end do
    
  end subroutine create_correction_nodal_2d

  subroutine create_correction_nodal_3d(lo,hi,correction_nodal,ng_cn,correction_cc,ng_cc)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:), ng_cn, ng_cc
    real (kind=dp_t), intent(out) :: correction_nodal(lo(1)-ng_cn:,lo(2)-ng_cn:,lo(3)-ng_cn:)
    real (kind=dp_t), intent(in ) ::    correction_cc(lo(1)-ng_cc:,lo(2)-ng_cc:,lo(3)-ng_cc:)
    
    ! Local variables
    integer :: i, j,k
    
!$omp parallel do private(i,j,k)
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)+1
             correction_nodal(i,j,k) = EIGHTH * &
                  (correction_cc(i,j  ,k-1) + correction_cc(i-1,j  ,k-1) &
                  +correction_cc(i,j-1,k-1) + correction_cc(i-1,j-1,k-1) &
                  +correction_cc(i,j  ,k  ) + correction_cc(i-1,j  ,k  ) &
                  +correction_cc(i,j-1,k  ) + correction_cc(i-1,j-1,k  ) )
          end do
       end do
    end do
!$omp end parallel do
    
  end subroutine create_correction_nodal_3d
  
end module hgrhs_module
