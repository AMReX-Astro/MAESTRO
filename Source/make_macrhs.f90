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

  subroutine make_macrhs(macrhs,rho0,Source,delta_gamma1_term,Sbar,div_coeff,dx, &
                         gamma1bar_old,gamma1bar_new,p0_old,p0_new,delta_p_term,dt)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: dm, nlevs

    type(multifab) , intent(inout) :: macrhs(:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    type(multifab) , intent(in   ) :: Source(:)
    type(multifab) , intent(in   ) :: delta_gamma1_term(:)
    real(kind=dp_t), intent(in   ) :: Sbar(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    real(kind=dp_t), intent(in   ) :: gamma1bar_old(:,0:), gamma1bar_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    type(multifab) , intent(in   ) :: delta_p_term(:)

    real(kind=dp_t), pointer:: mp(:,:,:,:),sp(:,:,:,:),gp(:,:,:,:),pop(:,:,:,:)

    integer :: lo(dm),hi(dm)
    integer :: i,n,ng_rh,ng_sr,ng_dg,ng_dp

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_macrhs")

    ng_rh = macrhs(1)%ng
    ng_sr = Source(1)%ng
    ng_dg = delta_gamma1_term(1)%ng
    ng_dp = delta_p_term(1)%ng

    do n = 1, nlevs

       do i = 1, Source(n)%nboxes
          if ( multifab_remote(Source(n), i) ) cycle
          mp => dataptr(macrhs(n), i)
          sp => dataptr(Source(n), i)
          gp => dataptr(delta_gamma1_term(n), i)
          pop => dataptr(delta_p_term(n), i)
          lo =  lwb(get_box(Source(n), i))
          hi =  upb(get_box(Source(n), i))
          select case (dm)
          case (2)
             call make_macrhs_2d(n,lo,hi,rho0(n,:),mp(:,:,1,1),ng_rh,sp(:,:,1,1),ng_sr, &
                                 gp(:,:,1,1),ng_dg,Sbar(n,:),div_coeff(n,:), &
                                 gamma1bar_old(n,:),gamma1bar_new(n,:),p0_old(n,:), &
                                 p0_new(n,:),pop(:,:,1,1),ng_dp,dt)
          case (3)
             call make_macrhs_3d(n,lo,hi,rho0(n,:),mp(:,:,:,1),ng_rh,sp(:,:,:,1),ng_sr, &
                                 gp(:,:,:,1),ng_dg,Sbar(n,:),div_coeff(n,:),dx(n,:), &
                                 gamma1bar_old(n,:),gamma1bar_new(n,:), &
                                 p0_old(n,:),p0_new(n,:),pop(:,:,:,1),ng_dp,dt)
          end select
       end do

    enddo

    call destroy(bpt)

  end subroutine make_macrhs

  subroutine make_macrhs_2d(n,lo,hi,rho0,rhs,ng_rh,Source,ng_sr,delta_gamma1_term,ng_dg, &
                            Sbar,div_coeff,gamma1bar_old,gamma1bar_new,p0_old,p0_new, &
                            delta_p_term,ng_dp,dt)

    use probin_module, only: dpdt_factor
    use geometry, only: base_cutoff_density_coord

    integer         , intent(in   ) :: n, lo(:), hi(:), ng_rh, ng_sr, ng_dg, ng_dp
    real (kind=dp_t), intent(in   ) :: rho0(0:) 
    real (kind=dp_t), intent(  out) ::               rhs(lo(1)-ng_rh:,lo(2)-ng_rh:)  
    real (kind=dp_t), intent(in   ) ::            Source(lo(1)-ng_sr:,lo(2)-ng_sr:)  
    real (kind=dp_t), intent(in   ) :: delta_gamma1_term(lo(1)-ng_dg:,lo(2)-ng_dg:)  
    real (kind=dp_t), intent(in   ) :: Sbar(0:)  
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar_old(0:),gamma1bar_new(0:)
    real (kind=dp_t), intent(in   ) :: p0_old(0:),p0_new(0:)
    real (kind=dp_t), intent(in   ) ::      delta_p_term(lo(1)-ng_dp:,lo(2)-ng_dp:)
    real (kind=dp_t), intent(in   ) :: dt

    !     Local variables
    integer :: i, j
    real(kind=dp_t) :: gamma1bar_p0_avg

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          rhs(i,j) = div_coeff(j) * (Source(i,j) - Sbar(j) + delta_gamma1_term(i,j))
       end do
    end do

    if (dpdt_factor .gt. 0.0d0) then
       do j = lo(2),hi(2)
          if (j .lt. base_cutoff_density_coord(n)) then
             gamma1bar_p0_avg = 0.25d0 * &
                  (gamma1bar_old(j)+gamma1bar_new(j))*(p0_old(j)+p0_new(j))
             do i = lo(1),hi(1)
                rhs(i,j) = rhs(i,j) + div_coeff(j) * &
                     (dpdt_factor / gamma1bar_p0_avg) * (delta_p_term(i,j) / dt)
             end do
          end if
       end do
    end if

  end subroutine make_macrhs_2d

  subroutine make_macrhs_3d(n,lo,hi,rho0,rhs,ng_rh,Source,ng_sr,delta_gamma1_term,ng_dg, &
                            Sbar,div_coeff,dx,gamma1bar_old,gamma1bar_new,p0_old,p0_new, &
                            delta_p_term,ng_dp,dt)

    use geometry, only: spherical, base_cutoff_density_coord
    use probin_module, only: dpdt_factor, base_cutoff_density
    use fill_3d_module

    integer         , intent(in   ) :: n,lo(:),hi(:),ng_rh,ng_sr,ng_dg,ng_dp
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(  out) ::            rhs(lo(1)-ng_rh:,lo(2)-ng_rh:,lo(3)-ng_rh:)
    real (kind=dp_t), intent(in   ) ::         Source(lo(1)-ng_sr:,lo(2)-ng_sr:,lo(3)-ng_sr:)
    real (kind=dp_t), intent(in) :: delta_gamma1_term(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:)
    real (kind=dp_t), intent(in   ) ::      Sbar(0:)  
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)  
    real (kind=dp_t), intent(in   ) :: dx(:)
    real (kind=dp_t), intent(in   ) :: gamma1bar_old(0:),gamma1bar_new(0:)
    real (kind=dp_t), intent(in   ) :: p0_old(0:),p0_new(0:)
    real (kind=dp_t), intent(in   ) ::   delta_p_term(lo(1)-ng_dp:,lo(2)-ng_dp:,lo(3)-ng_dp:)
    real (kind=dp_t), intent(in   ) :: dt

    !     Local variables
    integer :: i, j, k
    real(kind=dp_t) :: gamma1bar_p0_avg
    real(kind=dp_t), allocatable ::           div_cart(:,:,:,:)
    real(kind=dp_t), allocatable ::          Sbar_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: gamma1bar_old_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: gamma1bar_new_cart(:,:,:,:)
    real(kind=dp_t), allocatable ::        p0_old_cart(:,:,:,:)
    real(kind=dp_t), allocatable ::        p0_new_cart(:,:,:,:)
    real(kind=dp_t), allocatable ::          rho0_cart(:,:,:,:)

    if (spherical .eq. 1) then

       allocate(div_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,div_coeff,div_cart, &
                                         lo,hi,dx,0,0)

       allocate(Sbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,Sbar,Sbar_cart, &
                                         lo,hi,dx,0,0)

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhs(i,j,k) = div_cart(i,j,k,1) * (Source(i,j,k) - Sbar_cart(i,j,k,1) + &
                     delta_gamma1_term(i,j,k))
             end do
          end do
       end do

       deallocate(Sbar_cart)

       if (dpdt_factor .ge. 0.0d0) then

          allocate(gamma1bar_old_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false., &
                                            gamma1bar_old,gamma1bar_old_cart,lo,hi,dx,0,0)

          allocate(gamma1bar_new_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false., &
                                            gamma1bar_new,gamma1bar_new_cart,lo,hi,dx,0,0)

          allocate(p0_old_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false., &
                                            p0_old,p0_old_cart,lo,hi,dx,0,0)

          allocate(p0_new_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false., &
                                            p0_new,p0_new_cart,lo,hi,dx,0,0)

          allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
          call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,rho0,rho0_cart,lo,hi,dx,0,0)


          do k = lo(3),hi(3)
             do j = lo(2),hi(2)                
                do i = lo(1),hi(1)
                   if(rho0_cart(i,j,k,1) .gt. base_cutoff_density) then
                      gamma1bar_p0_avg = 0.25d0 * &
                           (gamma1bar_old_cart(i,j,k,1) + gamma1bar_new_cart(i,j,k,1)) * &
                           (p0_old_cart(i,j,k,1) + p0_new_cart(i,j,k,1))
                      rhs(i,j,k) = rhs(i,j,k) + div_cart(i,j,k,1) * &
                           (dpdt_factor / gamma1bar_p0_avg) * (delta_p_term(i,j,k) / dt)
                   end if
                end do
             end do
          end do

          deallocate(gamma1bar_old_cart,gamma1bar_new_cart,p0_old_cart,p0_new_cart,rho0_cart)

       end if

       deallocate(div_cart)

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
             if (k .lt. base_cutoff_density_coord(n)) then
                gamma1bar_p0_avg = 0.25d0 * (gamma1bar_old(k) + gamma1bar_new(k)) * &
                     (p0_old(k) + p0_new(k))
                do j = lo(2),hi(2)                
                   do i = lo(1),hi(1)
                      rhs(i,j,k) = rhs(i,j,k) + div_coeff(k) * &
                           (dpdt_factor / gamma1bar_p0_avg) * (delta_p_term(i,j,k) / dt)
                   end do
                end do
             end if
          end do
       end if
       
    end if

  end subroutine make_macrhs_3d

end module macrhs_module
