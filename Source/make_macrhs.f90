! Create the righthand side to the elliptic equation that is solved in 
! the MAC project step, \beta * (S - \bar{S}).  For the MAC projection, 
! this quantity is cell-centered.
!
! Note, we include the delta_gamma1_term here, to (possibly) account for
! the effect of replacing \Gamma_1 by {\Gamma_1}_0 in the constraint
! equation (see paper III).

module macrhs_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: make_macrhs

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_macrhs(nlevs,macrhs,rho0,Source,delta_gamma1_term,Sbar,div_coeff,dx, &
                         gamma1bar_old,gamma1bar_new,p0_old,p0_new,ptherm_old,ptherm_new, &
                         pthermbar_old,pthermbar_new,dt)

    use bl_prof_module
    use bl_constants_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: macrhs(:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    type(multifab) , intent(in   ) :: Source(:)
    type(multifab) , intent(in   ) :: delta_gamma1_term(:)
    real(kind=dp_t), intent(in   ) :: Sbar(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    real(kind=dp_t), intent(in   ) :: gamma1bar_old(:,0:), gamma1bar_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    type(multifab) , intent(in   ) :: ptherm_old(:), ptherm_new(:)
    real(kind=dp_t), intent(in   ) :: pthermbar_old(:,0:), pthermbar_new(:,0:)

    real(kind=dp_t), pointer:: mp(:,:,:,:),sp(:,:,:,:),gp(:,:,:,:),pop(:,:,:,:),pnp(:,:,:,:)

    integer :: lo(Source(1)%dim),hi(Source(1)%dim)
    integer :: i,dm,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_macrhs")

    dm = Source(1)%dim

    do n = 1, nlevs

       do i = 1, Source(n)%nboxes
          if ( multifab_remote(Source(n), i) ) cycle
          mp => dataptr(macrhs(n), i)
          sp => dataptr(Source(n), i)
          gp => dataptr(delta_gamma1_term(n), i)
          pop => dataptr(ptherm_old(n), i)
          pnp => dataptr(ptherm_new(n), i)
          lo =  lwb(get_box(Source(n), i))
          hi =  upb(get_box(Source(n), i))
          select case (dm)
          case (2)
             call make_macrhs_2d(lo,hi,rho0(n,:),mp(:,:,1,1),sp(:,:,1,1),gp(:,:,1,1), &
                                 Sbar(n,:), &
                                 div_coeff(n,:),gamma1bar_old(n,:),gamma1bar_new(n,:), &
                                 p0_old(n,:),p0_new(n,:),pop(:,:,1,1),pnp(:,:,1,1), &
                                 pthermbar_old(n,:),pthermbar_new(n,:),dt)
          case (3)
             call make_macrhs_3d(n,lo,hi,mp(:,:,:,1),sp(:,:,:,1),gp(:,:,:,1),Sbar(n,:), &
                                 div_coeff(n,:),dx(n,:), &
                                 gamma1bar_old(n,:),gamma1bar_new(n,:), &
                                 p0_old(n,:),p0_new(n,:),pop(:,:,:,1),pnp(:,:,:,1), &
                                 pthermbar_old(n,:),pthermbar_new(n,:),dt)
          end select
       end do

    enddo

    call destroy(bpt)

  end subroutine make_macrhs

  subroutine make_macrhs_2d(lo,hi,rho0,rhs,Source,delta_gamma1_term,Sbar,div_coeff, &
                            gamma1bar_old,gamma1bar_new,p0_old,p0_new, &
                            ptherm_old,ptherm_new,pthermbar_old,pthermbar_new,dt)

    use probin_module, only: dpdt_factor, base_cutoff_density

    integer         , intent(in   ) :: lo(:), hi(:)
    real (kind=dp_t), intent(in   ) :: rho0(0:) 
    real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):)  
    real (kind=dp_t), intent(in   ) :: Source(lo(1):,lo(2):)  
    real (kind=dp_t), intent(in   ) :: delta_gamma1_term(lo(1):,lo(2):)  
    real (kind=dp_t), intent(in   ) :: Sbar(0:)  
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar_old(0:),gamma1bar_new(0:)
    real (kind=dp_t), intent(in   ) :: p0_old(0:),p0_new(0:)
    real (kind=dp_t), intent(in   ) :: ptherm_old(lo(1):,lo(2):),ptherm_new(lo(1):,lo(2):)
    real (kind=dp_t), intent(in   ) :: pthermbar_old(0:),pthermbar_new(0:)
    real (kind=dp_t), intent(in   ) :: dt

    !     Local variables
    integer :: i, j
    real(kind=dp_t) :: gamma1bar_p0_avg, ptherm_diff

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          rhs(i,j) = div_coeff(j) * (Source(i,j) - Sbar(j) + delta_gamma1_term(i,j))
       end do
    end do

    if (dpdt_factor .gt. 0.0d0) then
       do j = lo(2),hi(2)
          gamma1bar_p0_avg = 0.25d0*(gamma1bar_old(j)+gamma1bar_new(j))*(p0_old(j)+p0_new(j))
          do i = lo(1),hi(1)
             ptherm_diff = &
                  0.5d0*(ptherm_old(i,j)+ptherm_new(i,j)-pthermbar_old(j)-pthermbar_new(j))
             if (rho0(j) .gt. base_cutoff_density) then
                rhs(i,j) = rhs(i,j) + div_coeff(j) * &
                     (dpdt_factor / gamma1bar_p0_avg) * (ptherm_diff / dt)
             end if
          end do
       end do
    end if

  end subroutine make_macrhs_2d

  subroutine make_macrhs_3d(n,lo,hi,rhs,Source,delta_gamma1_term,Sbar,div_coeff,dx, &
                            gamma1bar_old,gamma1bar_new,p0_old,p0_new, &
                            ptherm_old,ptherm_new,pthermbar_old,pthermbar_new,dt)

    use geometry, only: spherical
    use probin_module, only: dpdt_factor
    use fill_3d_module

    integer         , intent(in   ) :: n,lo(:), hi(:)
    real (kind=dp_t), intent(  out) ::               rhs(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) ::            Source(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) :: delta_gamma1_term(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) ::      Sbar(0:)  
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)  
    real (kind=dp_t), intent(in   ) :: dx(:)
    real (kind=dp_t), intent(in   ) :: gamma1bar_old(0:),gamma1bar_new(0:)
    real (kind=dp_t), intent(in   ) :: p0_old(0:),p0_new(0:)
    real (kind=dp_t), intent(in   ) :: ptherm_old(lo(1):,lo(2):,lo(3):)
    real (kind=dp_t), intent(in   ) :: ptherm_new(lo(1):,lo(2):,lo(3):)
    real (kind=dp_t), intent(in   ) :: pthermbar_old(0:),pthermbar_new(0:)
    real (kind=dp_t), intent(in   ) :: dt

    !     Local variables
    integer :: i, j, k
    real(kind=dp_t) :: gamma1bar_p0_avg, ptherm_diff
    real(kind=dp_t), allocatable :: div_cart(:,:,:,:),Sbar_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: gamma1bar_old_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: gamma1bar_new_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: p0_old_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: p0_new_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: pthermbar_old_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: pthermbar_new_cart(:,:,:,:)

    if (spherical .eq. 1) then

       allocate(div_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1,div_coeff,div_cart, &
                                         lo,hi,dx,0)

       allocate(Sbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1,Sbar,Sbar_cart, &
                                         lo,hi,dx,0)

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhs(i,j,k) = div_cart(i,j,k,1) * (Source(i,j,k) - Sbar_cart(i,j,k,1) + &
                     delta_gamma1_term(i,j,k))
             end do
          end do
       end do

       deallocate(Sbar_cart,div_cart)

       if (dpdt_factor .ge. 0.0d0) then

          allocate(gamma1bar_old_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1, &
                                            gamma1bar_old,gamma1bar_old_cart,lo,hi,dx,0)

          allocate(gamma1bar_new_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1, &
                                            gamma1bar_new,gamma1bar_new_cart,lo,hi,dx,0)

          allocate(p0_old_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1, &
                                            p0_old,p0_old_cart,lo,hi,dx,0)

          allocate(p0_new_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1, &
                                            p0_new,p0_new_cart,lo,hi,dx,0)

          allocate(pthermbar_old_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1, &
                                            pthermbar_old,pthermbar_old_cart,lo,hi,dx,0)

          allocate(pthermbar_new_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1, &
                                            pthermbar_new,pthermbar_new_cart,lo,hi,dx,0)


          do k = lo(3),hi(3)
             do j = lo(2),hi(2)                
                do i = lo(1),hi(1)

                   gamma1bar_p0_avg = 0.25d0 * &
                        (gamma1bar_old_cart(i,j,k,1) + gamma1bar_new_cart(i,j,k,1)) * &
                        (p0_old_cart(i,j,k,1) + p0_new_cart(i,j,k,1))

                   ptherm_diff = 0.5d0*(ptherm_old(i,j,k) + ptherm_new(i,j,k) &
                        - pthermbar_old_cart(i,j,k,1) - pthermbar_new_cart(i,j,k,1))
                   rhs(i,j,k) = rhs(i,j,k) + div_cart(i,j,k,1) * &
                        (dpdt_factor / gamma1bar_p0_avg) * (ptherm_diff / dt)
                end do
             end do
          end do

          deallocate(gamma1bar_old_cart,gamma1bar_new_cart,p0_old_cart,p0_new_cart)
          deallocate(pthermbar_old_cart,pthermbar_new_cart)

       end if

    else

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhs(i,j,k) = div_coeff(k) * &
                     (Source(i,j,k) - Sbar(k) + delta_gamma1_term(i,j,k))
             end do
          end do
       end do

       if (dpdt_factor .gt. 0.0d0) then
          do k = lo(3),hi(3)
             gamma1bar_p0_avg = 0.25d0 * (gamma1bar_old(k) + gamma1bar_new(k)) * &
                  (p0_old(k) + p0_new(k))
             do j = lo(2),hi(2)                
                do i = lo(1),hi(1)
                   ptherm_diff = 0.5d0*(ptherm_old(i,j,k) + ptherm_new(i,j,k) &
                        - pthermbar_old(k) - pthermbar_new(k))
                   rhs(i,j,k) = rhs(i,j,k) + div_coeff(k) * &
                        (dpdt_factor / gamma1bar_p0_avg) * (ptherm_diff / dt)
                end do
             end do
          end do
       end if
       
    end if

  end subroutine make_macrhs_3d

end module macrhs_module
