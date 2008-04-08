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

  public :: make_hgrhs

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
       call fill_3d_data_c(nlevs,dx,the_bc_tower%bc_tower_array,mla, &
                           div_coeff_cart,div_coeff,foextrap_comp)
       call fill_3d_data_c(nlevs,dx,the_bc_tower%bc_tower_array,mla, &
                           Sbar_cart,Sbar,foextrap_comp)
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
  
end module hgrhs_module
