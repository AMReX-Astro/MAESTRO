! Create the righthand side to the elliptic equation that is solved in 
! the MAC project step, \beta0 * (S - \bar{S}).  For the MAC projection, 
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

  subroutine make_macrhs(macrhs,rho0,S_cc,delta_gamma1_term,Sbar,beta0,dx, &
                         gamma1bar,p0,delta_p_term,dt,delta_chi,is_predictor)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical

    type(multifab) , intent(inout) :: macrhs(:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    type(multifab) , intent(in   ) :: S_cc(:)
    type(multifab) , intent(in   ) :: delta_gamma1_term(:)
    real(kind=dp_t), intent(in   ) :: Sbar(:,0:)
    real(kind=dp_t), intent(in   ) :: beta0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    type(multifab) , intent(in   ) :: delta_p_term(:)
    type(multifab) , intent(inout) :: delta_chi(:)
    logical        , intent(in   ) :: is_predictor

    real(kind=dp_t), pointer:: mp(:,:,:,:),sp(:,:,:,:),gp(:,:,:,:)
    real(kind=dp_t), pointer:: pop(:,:,:,:),dp(:,:,:,:)

    integer :: lo(get_dim(macrhs(1))),hi(get_dim(macrhs(1))),dm,nlevs
    integer :: i,n,ng_rh,ng_sr,ng_dg,ng_dp,ng_d

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_macrhs")

    dm = get_dim(macrhs(1))
    nlevs = size(macrhs)

    ng_rh = nghost(macrhs(1))
    ng_sr = nghost(S_cc(1))
    ng_dg = nghost(delta_gamma1_term(1))
    ng_dp = nghost(delta_p_term(1))
    ng_d  = nghost(delta_chi(1))

    do n = 1, nlevs

       do i = 1, nfabs(S_cc(n))
          mp => dataptr(macrhs(n), i)
          sp => dataptr(S_cc(n), i)
          gp => dataptr(delta_gamma1_term(n), i)
          dp => dataptr(delta_chi(n), i)
          pop => dataptr(delta_p_term(n), i)
          lo =  lwb(get_box(S_cc(n), i))
          hi =  upb(get_box(S_cc(n), i))
          select case (dm)
          case (1)
             call make_macrhs_1d(n,lo,hi,mp(:,1,1,1),ng_rh,sp(:,1,1,1),ng_sr, &
                                 gp(:,1,1,1),ng_dg,Sbar(n,:),beta0(n,:), &
                                 gamma1bar(n,:),p0(n,:), &
                                 pop(:,1,1,1),ng_dp, &
                                 dp(:,1,1,1),ng_d,dt,is_predictor)
          case (2)
             call make_macrhs_2d(n,lo,hi,mp(:,:,1,1),ng_rh,sp(:,:,1,1),ng_sr, &
                                 gp(:,:,1,1),ng_dg,Sbar(n,:),beta0(n,:), &
                                 gamma1bar(n,:),p0(n,:), &
                                 pop(:,:,1,1),ng_dp,dp(:,:,1,1),ng_d, &
                                 dt,is_predictor)
          case (3)
             if (spherical .eq. 1) then
                call make_macrhs_3d_sphr(lo,hi,rho0(1,:),mp(:,:,:,1),ng_rh,sp(:,:,:,1), &
                                         ng_sr,gp(:,:,:,1),ng_dg,Sbar(1,:),beta0(1,:), &
                                         dx(n,:),gamma1bar(1,:), &
                                         p0(1,:),pop(:,:,:,1),ng_dp, &
                                         dp(:,:,:,1),ng_d,dt,is_predictor)
             else
                call make_macrhs_3d(n,lo,hi,mp(:,:,:,1),ng_rh,sp(:,:,:,1),ng_sr, &
                                    gp(:,:,:,1),ng_dg,Sbar(n,:),beta0(n,:), &
                                    gamma1bar(n,:), &
                                    p0(n,:),pop(:,:,:,1),ng_dp, &
                                    dp(:,:,:,1),ng_d,dt,is_predictor)
             end if
          end select
       end do

    enddo

    call destroy(bpt)

  end subroutine make_macrhs

  subroutine make_macrhs_1d(n,lo,hi,rhs,ng_rh,S_cc,ng_sr,delta_gamma1_term,ng_dg, &
                            Sbar,beta0,gamma1bar,p0, &
                            delta_p_term,ng_dp,delta_chi,ng_d,dt,is_predictor)

    use probin_module, only: dpdt_factor
    use geometry, only: base_cutoff_density_coord

    integer         , intent(in   ) :: n, lo(:), hi(:), ng_rh, ng_sr, ng_dg, ng_dp, ng_d
    real (kind=dp_t), intent(  out) ::               rhs(lo(1)-ng_rh:)
    real (kind=dp_t), intent(in   ) ::              S_cc(lo(1)-ng_sr:)
    real (kind=dp_t), intent(in   ) :: delta_gamma1_term(lo(1)-ng_dg:)
    real (kind=dp_t), intent(in   ) :: Sbar(0:)  
    real (kind=dp_t), intent(in   ) :: beta0(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) ::      delta_p_term(lo(1)-ng_dp:)
    real (kind=dp_t), intent(inout) ::         delta_chi(lo(1)-ng_d :)
    real (kind=dp_t), intent(in   ) :: dt
    logical         , intent(in   ) :: is_predictor

    !     Local variables
    integer :: i

    do i = lo(1),hi(1)
       rhs(i) = beta0(i) * (S_cc(i) - Sbar(i) + delta_gamma1_term(i))
    end do

    if (dpdt_factor .gt. 0.0d0) then

       if (is_predictor) &
          delta_chi = 0.d0

       do i = lo(1),hi(1)
          if (i .lt. base_cutoff_density_coord(n)) then
             delta_chi(i) = delta_chi(i) + dpdt_factor * delta_p_term(i) / (dt*gamma1bar(i)*p0(i))
             rhs(i) = rhs(i) + beta0(i) * delta_chi(i)
          end if
       end do

    end if

  end subroutine make_macrhs_1d

  subroutine make_macrhs_2d(n,lo,hi,rhs,ng_rh,S_cc,ng_sr,delta_gamma1_term,ng_dg, &
                            Sbar,beta0,gamma1bar,p0, &
                            delta_p_term,ng_dp,delta_chi,ng_d,dt,is_predictor)

    use probin_module, only: dpdt_factor
    use geometry, only: base_cutoff_density_coord

    integer         , intent(in   ) :: n, lo(:), hi(:), ng_rh, ng_sr, ng_dg, ng_dp, ng_d
    real (kind=dp_t), intent(  out) ::               rhs(lo(1)-ng_rh:,lo(2)-ng_rh:)  
    real (kind=dp_t), intent(in   ) ::              S_cc(lo(1)-ng_sr:,lo(2)-ng_sr:)  
    real (kind=dp_t), intent(in   ) :: delta_gamma1_term(lo(1)-ng_dg:,lo(2)-ng_dg:)  
    real (kind=dp_t), intent(in   ) :: Sbar(0:)  
    real (kind=dp_t), intent(in   ) :: beta0(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) ::      delta_p_term(lo(1)-ng_dp:,lo(2)-ng_dp:)
    real (kind=dp_t), intent(inout) ::         delta_chi(lo(1)-ng_d :,lo(2)-ng_d :)
    real (kind=dp_t), intent(in   ) :: dt
    logical         , intent(in   ) :: is_predictor

    !     Local variables
    integer :: i, j

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          rhs(i,j) = beta0(j) * (S_cc(i,j) - Sbar(j) + delta_gamma1_term(i,j))
       end do
    end do

    if (dpdt_factor .gt. 0.0d0) then

       if (is_predictor) &
          delta_chi = 0.d0

       do j = lo(2),hi(2)
          if (j .lt. base_cutoff_density_coord(n)) then
             do i = lo(1),hi(1)
                delta_chi(i,j) = delta_chi(i,j) + dpdt_factor * delta_p_term(i,j) / (dt*gamma1bar(j)*p0(j))
                rhs(i,j) = rhs(i,j) + beta0(j) * delta_chi(i,j)
             end do
          end if
       end do

    end if

  end subroutine make_macrhs_2d

  subroutine make_macrhs_3d(n,lo,hi,rhs,ng_rh,S_cc,ng_sr,delta_gamma1_term,ng_dg, &
                            Sbar,beta0,gamma1bar,p0, &
                            delta_p_term,ng_dp,delta_chi,ng_d,dt,is_predictor)

    use geometry, only: base_cutoff_density_coord
    use probin_module, only: dpdt_factor
    use fill_3d_module

    integer         , intent(in   ) :: n,lo(:),hi(:),ng_rh,ng_sr,ng_dg,ng_dp,ng_d
    real (kind=dp_t), intent(  out) ::            rhs(lo(1)-ng_rh:,lo(2)-ng_rh:,lo(3)-ng_rh:)
    real (kind=dp_t), intent(in   ) ::           S_cc(lo(1)-ng_sr:,lo(2)-ng_sr:,lo(3)-ng_sr:)
    real (kind=dp_t), intent(in) :: delta_gamma1_term(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:)
    real (kind=dp_t), intent(in   ) ::      Sbar(0:)  
    real (kind=dp_t), intent(in   ) :: beta0(0:)  
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) ::   delta_p_term(lo(1)-ng_dp:,lo(2)-ng_dp:,lo(3)-ng_dp:)
    real (kind=dp_t), intent(inout) ::      delta_chi(lo(1)-ng_d :,lo(2)-ng_d :,lo(3)-ng_d: )
    real (kind=dp_t), intent(in   ) :: dt
    logical         , intent(in   ) :: is_predictor

    !     Local variables
    integer :: i, j, k
    
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhs(i,j,k) = beta0(k) * &
                  (S_cc(i,j,k) - Sbar(k) + delta_gamma1_term(i,j,k))
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    if (dpdt_factor .gt. 0.0d0) then

       if (is_predictor) &
          delta_chi = 0.d0

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          if (k .lt. base_cutoff_density_coord(n)) then
             do j = lo(2),hi(2)                
                do i = lo(1),hi(1)
                   delta_chi(i,j,k) = delta_chi(i,j,k) + dpdt_factor * delta_p_term(i,j,k) / &
                           (dt*gamma1bar(k)*p0(k))
                   rhs(i,j,k) = rhs(i,j,k) + beta0(k) * delta_chi(i,j,k)
                end do
             end do
          end if
       end do
       !$OMP END PARALLEL DO

    end if
       
  end subroutine make_macrhs_3d

  subroutine make_macrhs_3d_sphr(lo,hi,rho0,rhs,ng_rh,S_cc,ng_sr,delta_gamma1_term,ng_dg, &
                                 Sbar,beta0,dx,gamma1bar,p0, &
                                 delta_p_term,ng_dp,delta_chi,ng_d,dt,is_predictor)

    use probin_module, only: dpdt_factor, base_cutoff_density
    use fill_3d_module

    integer         , intent(in   ) :: lo(:),hi(:),ng_rh,ng_sr,ng_dg,ng_dp,ng_d
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(  out) ::            rhs(lo(1)-ng_rh:,lo(2)-ng_rh:,lo(3)-ng_rh:)
    real (kind=dp_t), intent(in   ) ::           S_cc(lo(1)-ng_sr:,lo(2)-ng_sr:,lo(3)-ng_sr:)
    real (kind=dp_t), intent(in) :: delta_gamma1_term(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:)
    real (kind=dp_t), intent(in   ) ::      Sbar(0:)  
    real (kind=dp_t), intent(in   ) :: beta0(0:)  
    real (kind=dp_t), intent(in   ) :: dx(:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) ::   delta_p_term(lo(1)-ng_dp:,lo(2)-ng_dp:,lo(3)-ng_dp:)
    real (kind=dp_t), intent(inout) ::      delta_chi(lo(1)-ng_d :,lo(2)-ng_d :,lo(3)-ng_d: )
    real (kind=dp_t), intent(in   ) :: dt
    logical         , intent(in   ) :: is_predictor

    !     Local variables
    integer :: i, j, k
    real(kind=dp_t), allocatable ::       div_cart(:,:,:,:)
    real(kind=dp_t), allocatable ::      Sbar_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: gamma1bar_cart(:,:,:,:)
    real(kind=dp_t), allocatable ::        p0_cart(:,:,:,:)
    real(kind=dp_t), allocatable ::      rho0_cart(:,:,:,:)

    allocate(div_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,beta0,div_cart,lo,hi,dx,0)
    
    allocate(Sbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,Sbar,Sbar_cart,lo,hi,dx,0)
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhs(i,j,k) = div_cart(i,j,k,1) * (S_cc(i,j,k) - Sbar_cart(i,j,k,1) + &
                          delta_gamma1_term(i,j,k))
          end do
       end do
    end do
    
    deallocate(Sbar_cart)
    
    if (dpdt_factor .ge. 0.0d0) then
       
       allocate(gamma1bar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(.false.,.false.,gamma1bar, &
                                         gamma1bar_cart,lo,hi,dx,0)
       
       allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(.false.,.false.,p0,p0_cart,lo,hi,dx,0)
       
       allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(.false.,.false.,rho0,rho0_cart,lo,hi,dx,0)
       
       if (is_predictor) &
          delta_chi = 0.d0

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)                
             do i = lo(1),hi(1)
                if (rho0_cart(i,j,k,1) .gt. base_cutoff_density) then
                   delta_chi(i,j,k) = delta_chi(i,j,k) + dpdt_factor * delta_p_term(i,j,k) / &
                           (dt*gamma1bar_cart(i,j,k,1)*p0_cart(i,j,k,1))
                   rhs(i,j,k) = rhs(i,j,k) + div_cart(i,j,k,1) * delta_chi(i,j,k)                           
                end if
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       
       deallocate(gamma1bar_cart,p0_cart,rho0_cart)
       
    end if
    
    deallocate(div_cart)
    
  end subroutine make_macrhs_3d_sphr
  
end module macrhs_module
