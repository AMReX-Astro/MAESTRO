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
  
  subroutine make_hgrhs(nlevs,the_bc_tower,mla,hgrhs,Source,delta_gamma1_term,Sbar, &
                        div_coeff,dx)

    use define_bc_module
    use ml_layout_module
    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical
    use fill_3d_module
    use variables, only: foextrap_comp
    use ml_restriction_module
    use multifab_fill_ghost_module
    use multifab_physbc_module
    
    integer        , intent(in   ) :: nlevs
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
    integer :: lo(Source(1)%dim),hi(Source(1)%dim)
    integer :: i,dm,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_hgrhs")
    
    dm = Source(1)%dim

    if(spherical .eq. 1) then
       do n = 1, nlevs
          call multifab_build(Sbar_cart(n),Source(n)%la,1,0)
          call multifab_build(div_coeff_cart(n),Source(n)%la,1,0)
          call setval(Sbar_cart(n),ZERO,all=.true.)
          call setval(div_coeff_cart(n),ZERO,all=.true.)
       end do
    end if
    
    if (spherical .eq. 1) then
       call put_1d_array_on_cart(nlevs,div_coeff,div_coeff_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla,1)
       call put_1d_array_on_cart(nlevs,Sbar,Sbar_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla,1)
    end if

    do n = 1, nlevs
       call multifab_build(rhs_cc(n),Source(n)%la,1,1)
       call setval(rhs_cc(n),ZERO,all=.true.)
    end do

    do n = 1, nlevs
       do i = 1, Source(n)%nboxes
          if ( multifab_remote(Source(n), i) ) cycle
          rp => dataptr(rhs_cc(n), i)
          sp => dataptr(Source(n), i)
          gp => dataptr(delta_gamma1_term(n), i)
          lo =  lwb(get_box(Source(n), i))
          hi =  upb(get_box(Source(n), i))
          select case (dm)
          case (2)
             call make_rhscc_2d(lo,hi,rp(:,:,1,1),sp(:,:,1,1),gp(:,:,1,1),Sbar(n,:), &
                                div_coeff(n,:))
          case (3)
             if (spherical .eq. 1) then
                dp => dataptr(div_coeff_cart(n), i)
                sbp => dataptr(Sbar_cart(n), i)
                call make_rhscc_3d_sphr(lo,hi,rp(:,:,:,1),sp(:,:,:,1),gp(:,:,:,1), &
                                        sbp(:,:,:,1),dp(:,:,:,1))
             else
                call make_rhscc_3d_cart(lo,hi,rp(:,:,:,1),sp(:,:,:,1),gp(:,:,:,1), &
                                        Sbar(n,:),div_coeff(n,:))
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
                                         the_bc_tower%bc_tower_array(n),1,foextrap_comp,1)
       end do

    end if
       
    do n=1,nlevs
       call setval(hgrhs(n),ZERO,all=.true.)
       do i = 1, Source(n)%nboxes
          if ( multifab_remote(Source(n), i) ) cycle
          hp => dataptr(hgrhs(n), i)
          rp => dataptr(rhs_cc(n), i)
          lo =  lwb(get_box(Source(n), i))
          hi =  upb(get_box(Source(n), i))
          select case (dm)
          case (2)
             call make_hgrhs_2d(lo,hi,hp(:,:,1,1),rp(:,:,1,1))
          case (3)
             call make_hgrhs_3d(lo,hi,hp(:,:,:,1),rp(:,:,:,1))
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
  
  subroutine make_rhscc_2d(lo,hi,rhs_cc,Source,delta_gamma1_term,Sbar,div_coeff)

    integer         , intent(in   ) :: lo(:), hi(:)
    real (kind=dp_t), intent(  out) ::            rhs_cc(lo(1)-1:,lo(2)-1:)
    real (kind=dp_t), intent(in   ) ::            Source(lo(1)  :,lo(2)  :)
    real (kind=dp_t), intent(in   ) :: delta_gamma1_term(lo(1)  :,lo(2)  :)  
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
  
  subroutine make_rhscc_3d_cart(lo,hi,rhs_cc,Source,delta_gamma1_term,Sbar,div_coeff)

    integer         , intent(in   ) :: lo(:), hi(:)
    real (kind=dp_t), intent(  out) ::            rhs_cc(lo(1)-1:,lo(2)-1:,lo(3)-1:)  
    real (kind=dp_t), intent(in   ) ::            Source(lo(1)  :,lo(2)  :,lo(3)  :)  
    real (kind=dp_t), intent(in   ) :: delta_gamma1_term(lo(1)  :,lo(2)  :,lo(3)  :)  
    real (kind=dp_t), intent(in   ) ::      Sbar(0:)
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    
    ! Local variables
    integer :: i,j,k
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhs_cc(i,j,k) = div_coeff(k) * (Source(i,j,k) - Sbar(k) + &
                  delta_gamma1_term(i,j,k))
          end do
       end do
    end do
    
  end subroutine make_rhscc_3d_cart
   
  subroutine make_rhscc_3d_sphr(lo,hi,rhs_cc,Source,delta_gamma1_term,Sbar_cart, &
                                div_coeff_cart)

    integer         , intent(in   ) :: lo(:), hi(:)
    real (kind=dp_t), intent(  out) ::            rhs_cc(lo(1)-1:,lo(2)-1:,lo(3)-1:)  
    real (kind=dp_t), intent(in   ) ::            Source(lo(1)  :,lo(2)  :,lo(3)  :)  
    real (kind=dp_t), intent(in   ) :: delta_gamma1_term(lo(1)  :,lo(2)  :,lo(3)  :)  
    real (kind=dp_t), intent(in   ) ::         Sbar_cart(lo(1)  :,lo(2)  :,lo(3)  :)  
    real (kind=dp_t), intent(in   ) ::    div_coeff_cart(lo(1)  :,lo(2)  :,lo(3)  :)  
    
    ! Local variables
    integer :: i, j,k
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhs_cc(i,j,k) = div_coeff_cart(i,j,k) * (Source(i,j,k) - Sbar_cart(i,j,k) + &
                  delta_gamma1_term(i,j,k))
             
          end do
       end do
    end do
    
  end subroutine make_rhscc_3d_sphr
  
  subroutine make_hgrhs_2d(lo,hi,rhs,rhs_cc)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:)
    real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):)  
    real (kind=dp_t), intent(in   ) :: rhs_cc(lo(1)-1:,lo(2)-1:)
    
    ! Local variables
    integer :: i, j
    
    do j = lo(2),hi(2)+1
       do i = lo(1), hi(1)+1
          rhs(i,j) = FOURTH * ( rhs_cc(i,j  ) + rhs_cc(i-1,j  ) &
               + rhs_cc(i,j-1) + rhs_cc(i-1,j-1) )
       end do
    end do
    
  end subroutine make_hgrhs_2d
  
  subroutine make_hgrhs_3d(lo,hi,rhs,rhs_cc)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:)
    real (kind=dp_t), intent(  out) ::    rhs(lo(1)  :,lo(2)  :,lo(3)  :)  
    real (kind=dp_t), intent(in   ) :: rhs_cc(lo(1)-1:,lo(2)-1:,lo(3)-1:)  
    
    ! Local variables
    integer :: i, j,k
    
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)+1
             rhs(i,j,k) = EIGHTH * ( rhs_cc(i,j  ,k-1) + rhs_cc(i-1,j  ,k-1) &
                  +rhs_cc(i,j-1,k-1) + rhs_cc(i-1,j-1,k-1) &
                  +rhs_cc(i,j  ,k  ) + rhs_cc(i-1,j  ,k  ) &
                  +rhs_cc(i,j-1,k  ) + rhs_cc(i-1,j-1,k  ) )
          end do
       end do
    end do
    
  end subroutine make_hgrhs_3d

  subroutine correct_hgrhs(nlevs,the_bc_tower,mla,hgrhs,div_coeff,dx,dt,gamma1bar,p0, &
                           ptherm,pthermbar)

    use define_bc_module
    use ml_layout_module
    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical
    use fill_3d_module
    use variables, only: foextrap_comp
    use ml_restriction_module
    use multifab_fill_ghost_module
    use multifab_physbc_module
    
    integer        , intent(in   ) :: nlevs
    type(bc_tower),  intent(in   ) :: the_bc_tower
    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: hgrhs(:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    type(multifab) , intent(in   ) :: ptherm(:)
    real(kind=dp_t), intent(in   ) :: pthermbar(:,0:)
    
    type(multifab) :: correction_cc(nlevs)
    type(multifab) :: correction_nodal(nlevs)

    type(multifab) :: gamma1bar_cart(nlevs)
    type(multifab) :: p0_cart(nlevs)
    type(multifab) :: pthermbar_cart(nlevs)

    real(kind=dp_t), pointer :: rhp(:,:,:,:), ptp(:,:,:,:), ccp(:,:,:,:), cnp(:,:,:,:)
    real(kind=dp_t), pointer :: gbp(:,:,:,:), p0p(:,:,:,:)
    integer :: lo(ptherm(1)%dim),hi(ptherm(1)%dim)
    integer :: i,dm,n
    logical :: nodal(ptherm(1)%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "correct_hgrhs")
    
    dm = ptherm(1)%dim
    nodal = .true.

    if(spherical .eq. 1) then
       do n = 1, nlevs
          call multifab_build(gamma1bar_cart(n),ptherm(n)%la,1,0)
          call multifab_build(p0_cart(n),ptherm(n)%la,1,0)
          call multifab_build(pthermbar_cart(n),ptherm(n)%la,1,0)
          call setval(gamma1bar_cart(n),ZERO,all=.true.)
          call setval(p0_cart(n),ZERO,all=.true.)
          call setval(pthermbar_cart(n),ZERO,all=.true.)
       end do
    end if
    
    if (spherical .eq. 1) then
       call put_1d_array_on_cart(nlevs,gamma1bar,gamma1bar_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla,1)
       call put_1d_array_on_cart(nlevs,p0,p0_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla,1)
       call put_1d_array_on_cart(nlevs,pthermbar,pthermbar_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla,1)
    end if

    do n = 1, nlevs
       call multifab_build(correction_cc(n),ptherm(n)%la,1,1)
       call setval(correction_cc(n),ZERO,all=.true.)
       call multifab_build(correction_nodal(n),ptherm(n)%la,1,0,nodal)
       call setval(correction_nodal(n),ZERO,all=.true.)
    end do

    do n = 1, nlevs
       do i = 1, ptherm(n)%nboxes
          if ( multifab_remote(ptherm(n), i) ) cycle
          ccp => dataptr(correction_cc(n), i)
          ptp => dataptr(ptherm(n), i)
          lo =  lwb(get_box(ptherm(n), i))
          hi =  upb(get_box(ptherm(n), i))
          select case (dm)
          case (2)
             call create_correction_cc_2d(lo,hi,ccp(:,:,1,1),ptp(:,:,1,1),div_coeff(n,:), &
                                          gamma1bar(n,:),p0(n,:),pthermbar(n,:),dt)
          case (3)
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
                                         the_bc_tower%bc_tower_array(n),1,foextrap_comp,1)
       end do

    end if
       
    do n=1,nlevs
       call setval(correction_nodal(n),ZERO,all=.true.)
       do i = 1, ptherm(n)%nboxes
          if ( multifab_remote(ptherm(n), i) ) cycle
          cnp => dataptr(correction_nodal(n), i)
          ccp => dataptr(correction_cc(n), i)
          lo =  lwb(get_box(ptherm(n), i))
          hi =  upb(get_box(ptherm(n), i))
          select case (dm)
          case (2)
             call create_correction_nodal_2d(lo,hi,cnp(:,:,1,1),ccp(:,:,1,1))
          case (3)
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
          call destroy(pthermbar_cart(n))
       end if
    end do

    call destroy(bpt)
    
  end subroutine correct_hgrhs
  
  subroutine create_correction_cc_2d(lo,hi,correction_cc,ptherm,div_coeff,gamma1bar, &
                                     p0,pthermbar,dt)

    use probin_module, only: dpdt_factor

    integer         , intent(in   ) :: lo(:), hi(:)
    real (kind=dp_t), intent(  out) :: correction_cc(lo(1)-1:,lo(2)-1:)
    real (kind=dp_t), intent(in   ) :: ptherm(lo(1):,lo(2):)
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)    
    real (kind=dp_t), intent(in   ) :: pthermbar(0:)
    real (kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: i, j
    real(kind=dp_t) :: temp
    
    do j = lo(2),hi(2)
       temp = div_coeff(j)*(dpdt_factor/(gamma1bar(j)*p0(j))) / dt
       do i = lo(1),hi(1)
          correction_cc(i,j) = temp*(ptherm(i,j)-pthermbar(j))
       end do
    end do
    
  end subroutine create_correction_cc_2d
  
  subroutine create_correction_nodal_2d(lo,hi,correction_nodal,correction_cc)

    use bl_constants_module

    integer         , intent(in   ) :: lo(:), hi(:)
    real (kind=dp_t), intent(  out) :: correction_nodal(lo(1):,lo(2):)  
    real (kind=dp_t), intent(in   ) :: correction_cc(lo(1)-1:,lo(2)-1:)
    
    ! Local variables
    integer :: i, j
    
    do j = lo(2),hi(2)+1
       do i = lo(1), hi(1)+1
          correction_nodal(i,j) = FOURTH * ( correction_cc(i,j  ) + correction_cc(i-1,j  ) &
               + correction_cc(i,j-1) + correction_cc(i-1,j-1) )
       end do
    end do
    
  end subroutine create_correction_nodal_2d
  
end module hgrhs_module
