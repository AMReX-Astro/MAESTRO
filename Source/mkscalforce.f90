module mkscalforce_module
  !
  ! this module contains the 2d and 3d routines that make the forcing
  ! terms for the scalar equations.
  ! 
  ! mkrhohforce computes the  wtilde dp/dr + psi source for rho*h evolution
  !
  ! mktempforce computes the source terms that appear in the
  ! temperature evolution equation.  Note, this formulation is only
  ! used in the prediction of the temperature on the edges, when
  ! enthalpy_pred_type = predict_T_then_rhohprime or predict_T_then_h.  
  ! The edge temperatures are then converted to enthalpy for the full update.
  !

  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use ml_restrict_fill_module

  implicit none

  private

  public :: mkrhohforce, mkhprimeforce, mktempforce

contains

  subroutine mkrhohforce(mla, scal_force, is_prediction, thermal, umac, &
                         p0_old, p0_new, rho0_old, rho0_new, &
                         psi, dx, add_thermal, the_bc_level)

    use bl_prof_module
    use variables, only: foextrap_comp, rhoh_comp
    use geometry, only: spherical, nr_fine, nlevs_radial, polar
    use fill_3d_module
    use make_grav_module
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    type(multifab) , intent(inout) :: scal_force(:)
    logical        , intent(in   ) :: is_prediction
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) ::   p0_old(:,0:)
    real(kind=dp_t), intent(in   ) ::   p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) ::      psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    logical        , intent(in   ) :: add_thermal
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer                  :: i,n,comp,ng_f,ng_um,ng_th,ng_p0,ng_pm
    integer                  :: lo(mla%dim),hi(mla%dim),dm,nlevs

    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: tp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    real(kind=dp_t), pointer :: pmx(:,:,:,:)
    real(kind=dp_t), pointer :: pmy(:,:,:,:)
    real(kind=dp_t), pointer :: pmz(:,:,:,:)

    type(multifab) :: p0_cart(mla%nlevel)
    type(multifab) :: p0mac(mla%nlevel,mla%dim)

    real(kind=dp_t) :: p0_nph(nlevs_radial,0:nr_fine-1)
    real(kind=dp_t) ::   rho0(nlevs_radial,0:nr_fine-1)
    real(kind=dp_t) ::   grav(nlevs_radial,0:nr_fine-1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mkrhohforce")

    dm = mla%dim
    nlevs = mla%nlevel

    ! if we are doing the prediction, then it only makes sense to be in
    ! this routine if the quantity we are predicting is (rho h)' or h
    if (is_prediction .AND. &
         (enthalpy_pred_type == predict_T_then_rhohprime .OR. &
          enthalpy_pred_type == predict_T_then_h)) then
       call bl_error("ERROR: should not call mkrhohforce when predicting T")
    endif

    if (spherical .eq. 1 .or. polar .eq. 1) then
       p0_nph = HALF * (p0_old + p0_new)
       do n = 1,nlevs
          call multifab_build(p0_cart(n),mla%la(n),1,1)
          do comp=1,dm
             call multifab_build_edge(p0mac(n,comp),mla%la(n),1,1,comp)
          end do
       end do

       call put_1d_array_on_cart(p0_nph,p0_cart,foextrap_comp,.false.,.false.,&
                                 dx,the_bc_level,mla)

       call make_s0mac(mla,p0_nph,p0mac,dx,foextrap_comp,the_bc_level)
    end if

    ng_f  = nghost(scal_force(1))
    ng_um = nghost(umac(1,1))
    ng_th = nghost(thermal(1))
    ng_pm = nghost(p0mac(1,1))
    ng_p0 = nghost(p0_cart(1))

    rho0 = HALF*(rho0_old + rho0_new)

    call make_grav_cell(grav,rho0)
    
    do n=1,nlevs

       do i=1, nfabs(scal_force(n))
          fp => dataptr(scal_force(n), i)
          tp  => dataptr(thermal(n),i)
          lo = lwb(get_box(scal_force(n),i))
          hi = upb(get_box(scal_force(n),i))
          select case (dm)
          case (1)
             ump => dataptr(umac(n,1),i)
             call mkrhohforce_1d(n,fp(:,1,1,rhoh_comp), ng_f, is_prediction, &
                                 ump(:,1,1,1), ng_um, &
                                 tp(:,1,1,1), ng_th, lo, hi, p0_old(n,:), p0_new(n,:), &
                                 rho0(n,:), grav(n,:), psi(n,:), add_thermal)
          case (2)
             if (polar .eq. 0) then
                vmp => dataptr(umac(n,2),i)
                call mkrhohforce_2d(n,fp(:,:,1,rhoh_comp), ng_f, is_prediction, &
                                    vmp(:,:,1,1), ng_um, &
                                    tp(:,:,1,1), ng_th, lo, hi, p0_old(n,:), p0_new(n,:), &
                                    rho0(n,:), grav(n,:), psi(n,:), add_thermal)
            else 
                ump => dataptr(umac(n,1),i)
                vmp => dataptr(umac(n,2),i)
                pp  => dataptr(p0_cart(n),i)
                pmx => dataptr(p0mac(n,1),i)
                pmy => dataptr(p0mac(n,2),i)
                call mkrhohforce_2d_polar(fp(:,:,1,rhoh_comp), ng_f, is_prediction, &
                                         ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                         tp(:,:,1,1), ng_th, pp(:,:,1,1), ng_p0, &
                                         pmx(:,:,1,1), pmy(:,:,1,1), ng_pm, &
                                         lo, hi, dx(n,:), psi(1,:), add_thermal)
            end if
          case(3)
             wmp => dataptr(umac(n,3), i)
             if (spherical .eq. 0) then
                call mkrhohforce_3d(n,fp(:,:,:,rhoh_comp), ng_f, is_prediction, &
                                    wmp(:,:,:,1), ng_um, &
                                    tp(:,:,:,1), ng_th, lo, hi, p0_old(n,:), p0_new(n,:), &
                                    rho0(n,:), grav(n,:), psi(n,:), add_thermal)
             else
                ump => dataptr(umac(n,1),i)
                vmp => dataptr(umac(n,2),i)
                wmp => dataptr(umac(n,3),i)
                pp  => dataptr(p0_cart(n),i)
                pmx => dataptr(p0mac(n,1),i)
                pmy => dataptr(p0mac(n,2),i)
                pmz => dataptr(p0mac(n,3),i)
                call mkrhohforce_3d_sphr(fp(:,:,:,rhoh_comp), ng_f, is_prediction, &
                                         ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                         tp(:,:,:,1), ng_th, pp(:,:,:,1), ng_p0, &
                                         pmx(:,:,:,1), pmy(:,:,:,1), pmz(:,:,:,1), ng_pm, &
                                         lo, hi, dx(n,:), psi(1,:), add_thermal)
             end if
          end select
       end do

    end do

    if (spherical .eq. 1 .or. polar .eq. 1) then
       do n=1,nlevs
          call destroy(p0_cart(n))
          do comp=1,dm
             call destroy(p0mac(n,comp))
          end do
       end do
    end if

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,scal_force,mla%mba%rr,the_bc_level, &
                              icomp=rhoh_comp, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=scal_force(1)%ng)

    call destroy(bpt)
    
  end subroutine mkrhohforce

  subroutine mkrhohforce_1d(n,rhoh_force,ng_f,is_prediction, &
                            umac,ng_um,thermal,ng_th,lo,hi, &
                            p0_old,p0_new,rho0,grav,psi,add_thermal)

    use geometry, only: dr, nr, base_cutoff_density_coord
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}
    
    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the rhoh_force for the update, they will be used to time-center.

    integer,         intent(in   ) :: n,lo(:),hi(:),ng_f,ng_um,ng_th
    logical,         intent(in   ) :: is_prediction
    real(kind=dp_t), intent(  out) :: rhoh_force(lo(1)-ng_f :)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    logical        , intent(in   ) :: add_thermal

    real(kind=dp_t) :: gradp0, uadv
    integer :: i

    ! Add wtilde d(p0)/dr
    do i = lo(1),hi(1)

       if (i .lt. base_cutoff_density_coord(n)) then
          gradp0 = rho0(i) * grav(i)
       else if (i.eq.nr(n)-1) then
          ! NOTE: this should be zero since p0 is constant up here
          gradp0 = HALF * ( p0_old(i  ) + p0_new(i  ) &
                           -p0_old(i-1) - p0_new(i-1) ) / dr(n)
       else
          ! NOTE: this should be zero since p0 is constant up here
          gradp0 = HALF * ( p0_old(i+1) + p0_new(i+1) &
                           -p0_old(i  ) - p0_new(i  ) ) / dr(n)
       end if

       uadv = HALF*(umac(i)+umac(i+1))
       rhoh_force(i) = uadv * gradp0           
    enddo

    ! psi should always be in the force if we are doing the final update
    ! For prediction, it should not be in the force if we are predicting
    ! (rho h)', but should be there if we are predicting h or rhoh
    if ((is_prediction .AND. enthalpy_pred_type == predict_h) .OR. &
        (is_prediction .AND. enthalpy_pred_type == predict_rhoh) .OR. & 
        (.NOT. is_prediction)) then
       do i = lo(1),hi(1)
          rhoh_force(i) =  rhoh_force(i)  + psi(i)
       end do
    endif

    if (add_thermal) then
       do i=lo(1),hi(1)
          rhoh_force(i) = rhoh_force(i) + thermal(i)
       end do
    end if

  end subroutine mkrhohforce_1d

  subroutine mkrhohforce_2d(n,rhoh_force,ng_f,is_prediction, &
                            vmac,ng_um,thermal,ng_th,lo,hi, &
                            p0_old,p0_new,rho0,grav,psi,add_thermal)

    use geometry, only: dr, nr, base_cutoff_density_coord
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}
    
    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the rhoh_force for the update, they will be used to time-center.

    integer,         intent(in   ) :: n,lo(:),hi(:),ng_f,ng_um,ng_th
    logical,         intent(in   ) :: is_prediction
    real(kind=dp_t), intent(  out) :: rhoh_force(lo(1)-ng_f :,lo(2)-ng_f :)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:,lo(2)-ng_th:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    logical        , intent(in   ) :: add_thermal

    real(kind=dp_t) :: gradp0, vadv
    integer :: i,j


    ! Add wtilde d(p0)/dr
    do j = lo(2),hi(2)

       if (j .lt. base_cutoff_density_coord(n)) then
          gradp0 = rho0(j) * grav(j)
       else if (j.eq.nr(n)-1) then
          ! NOTE: this should be zero since p0 is constant up here
          gradp0 = HALF * ( p0_old(j  ) + p0_new(j  ) &
                           -p0_old(j-1) - p0_new(j-1) ) / dr(n)
       else
          ! NOTE: this should be zero since p0 is constant up here
          gradp0 = HALF * ( p0_old(j+1) + p0_new(j+1) &
                           -p0_old(j  ) - p0_new(j  ) ) / dr(n)
       end if


       do i = lo(1),hi(1)
          vadv = HALF*(vmac(i,j)+vmac(i,j+1))
          rhoh_force(i,j) =  vadv * gradp0           
       end do
    enddo

    ! psi should always be in the force if we are doing the final update
    ! For prediction, it should not be in the force if we are predicting
    ! (rho h)', but should be there if we are predicting h or rhoh
    if ((is_prediction .AND. enthalpy_pred_type == predict_h) .OR. &
        (is_prediction .AND. enthalpy_pred_type == predict_rhoh) .OR. &
        (.NOT. is_prediction)) then
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhoh_force(i,j) =  rhoh_force(i,j)  + psi(j)
          end do
       end do
    endif

    if (add_thermal) then
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhoh_force(i,j) = rhoh_force(i,j) + thermal(i,j)
          end do
       end do
    end if

  end subroutine mkrhohforce_2d

subroutine mkrhohforce_2d_polar(rhoh_force,ng_f,is_prediction, &
                                 umac,vmac,ng_um,thermal,ng_th, &
                                 p0_cart,ng_p0, &
                                 p0macx,p0macy, &
                                 ng_pm,lo,hi,dx,psi,add_thermal)

    use fill_3d_module
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer,         intent(in   ) :: lo(:),hi(:),ng_f,ng_um,ng_th,ng_p0,ng_pm
    logical,         intent(in   ) :: is_prediction
    real(kind=dp_t), intent(  out) :: rhoh_force(lo(1)-ng_f :,lo(2)-ng_f :)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:,lo(2)-ng_th:)
    real(kind=dp_t), intent(in   ) ::    p0_cart(lo(1)-ng_p0:,lo(2)-ng_p0:)
    real(kind=dp_t), intent(in   ) ::     p0macx(lo(1)-ng_pm:,lo(2)-ng_pm:)
    real(kind=dp_t), intent(in   ) ::     p0macy(lo(1)-ng_pm:,lo(2)-ng_pm:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    logical        , intent(in   ) :: add_thermal

    real(kind=dp_t), allocatable :: psi_cart(:,:,:)

    real(kind=dp_t) :: divup, p0divu
    integer         :: i,j
    !
    ! Here we make u grad p = div (u p) - p div (u)
    !
    !$OMP PARALLEL DO PRIVATE(i,j,divup,p0divu)
    do j = lo(2),hi(2)
        do i = lo(1),hi(1)

            divup = (umac(i+1,j) * p0macx(i+1,j) - umac(i,j) * p0macx(i,j)) / dx(1) + &
                    (vmac(i,j+1) * p0macy(i,j+1) - vmac(i,j) * p0macy(i,j)) / dx(2)
                    
            p0divu = ( (umac(i+1,j) - umac(i,j)) / dx(1) + &
                    (vmac(i,j+1) - vmac(i,j)) / dx(2)) * p0_cart(i,j)

            rhoh_force(i,j) = divup - p0divu

        end do
    end do
    !$OMP END PARALLEL DO
    
    
    ! psi should always be in the force if we are doing the final update
    ! For prediction, it should not be in the force if we are predicting
    ! (rho h)', but should be there if we are predicting h
    !
    if ((is_prediction .AND. enthalpy_pred_type == predict_h) .OR. &
        (is_prediction .AND. enthalpy_pred_type == predict_rhoh) .OR. &
        (.NOT. is_prediction)) then

       allocate(psi_cart(lo(1):hi(1),lo(2):hi(2),1))

       call put_1d_array_on_cart_2d_polar(.false.,.false.,psi,psi_cart,lo,hi,dx,0)

       !$OMP PARALLEL DO PRIVATE(i,j)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                rhoh_force(i,j) = rhoh_force(i,j) + psi_cart(i,j,1)
            enddo
        enddo
       !$OMP END PARALLEL DO
       deallocate(psi_cart)
    endif

    if (add_thermal) then
       !$OMP PARALLEL DO PRIVATE(i,j)
        do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                rhoh_force(i,j) = rhoh_force(i,j) + thermal(i,j)
            end do
        end do
       !$OMP END PARALLEL DO
    end if

  end subroutine mkrhohforce_2d_polar  
  
  subroutine mkrhohforce_3d(n,rhoh_force,ng_f,is_prediction, &
                            wmac,ng_um,thermal,ng_th,lo,hi,&
                            p0_old,p0_new,rho0,grav,psi,add_thermal)

    use geometry, only: dr, nr, base_cutoff_density_coord
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer,         intent(in   ) :: n,lo(:),hi(:),ng_f,ng_um,ng_th
    logical,         intent(in   ) :: is_prediction
    real(kind=dp_t), intent(  out) :: rhoh_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
    real(kind=dp_t), intent(in   ) ::       wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:,lo(2)-ng_th:,lo(3)-ng_th:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    logical        , intent(in   ) :: add_thermal

    real(kind=dp_t) :: gradp0,wadv
    integer :: i,j,k
    !
    ! Add wtilde d(p0)/dr
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k,gradp0,wadv)
    do k = lo(3),hi(3)

       if (k .lt. base_cutoff_density_coord(n)) then
          gradp0 = rho0(k) * grav(k)
       else if (k.eq.nr(n)-1) then
          gradp0 = HALF * ( p0_old(k  ) + p0_new(k  ) &
                           -p0_old(k-1) - p0_new(k-1) ) / dr(n)
       else
          gradp0 = HALF * ( p0_old(k+1) + p0_new(k+1) &
                           -p0_old(k  ) - p0_new(k  ) ) / dr(n)
       end if

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             wadv = HALF*(wmac(i,j,k)+wmac(i,j,k+1))
             rhoh_force(i,j,k) = wadv * gradp0 
          end do
       end do
    enddo
    !$OMP END PARALLEL DO

    ! psi should always be in the force if we are doing the final update
    ! For prediction, it should not be in the force if we are predicting
    ! (rho h)', but should be there if we are predicting h
    if ((is_prediction .AND. enthalpy_pred_type == predict_h) .OR. &
        (is_prediction .AND. enthalpy_pred_type == predict_rhoh) .OR. &
        (.NOT. is_prediction)) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + psi(k)
             end do
          end do
       enddo
       !$OMP END PARALLEL DO
    endif
       

    if (add_thermal) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + thermal(i,j,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

  end subroutine mkrhohforce_3d

  subroutine mkrhohforce_3d_sphr(rhoh_force,ng_f,is_prediction, &
                                 umac,vmac,wmac,ng_um,thermal,ng_th, &
                                 p0_cart,ng_p0, &
                                 p0macx,p0macy,p0macz, &
                                 ng_pm,lo,hi,dx,psi,add_thermal)

    use fill_3d_module
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer,         intent(in   ) :: lo(:),hi(:),ng_f,ng_um,ng_th,ng_p0,ng_pm
    logical,         intent(in   ) :: is_prediction
    real(kind=dp_t), intent(  out) :: rhoh_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:,lo(2)-ng_th:,lo(3)-ng_th:)
    real(kind=dp_t), intent(in   ) ::    p0_cart(lo(1)-ng_p0:,lo(2)-ng_p0:,lo(3)-ng_p0:)
    real(kind=dp_t), intent(in   ) ::     p0macx(lo(1)-ng_pm:,lo(2)-ng_pm:,lo(3)-ng_pm:)
    real(kind=dp_t), intent(in   ) ::     p0macy(lo(1)-ng_pm:,lo(2)-ng_pm:,lo(3)-ng_pm:)
    real(kind=dp_t), intent(in   ) ::     p0macz(lo(1)-ng_pm:,lo(2)-ng_pm:,lo(3)-ng_pm:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    logical        , intent(in   ) :: add_thermal

    real(kind=dp_t), allocatable :: psi_cart(:,:,:,:)

    real(kind=dp_t) :: divup, p0divu
    integer         :: i,j,k
    !
    ! Here we make u grad p = div (u p) - p div (u)
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k,divup,p0divu)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             divup = (umac(i+1,j,k) * p0macx(i+1,j,k) - umac(i,j,k) * p0macx(i,j,k)) / dx(1) + &
                     (vmac(i,j+1,k) * p0macy(i,j+1,k) - vmac(i,j,k) * p0macy(i,j,k)) / dx(2) + &
                     (wmac(i,j,k+1) * p0macz(i,j,k+1) - wmac(i,j,k) * p0macz(i,j,k)) / dx(3)

             p0divu = ( (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                        (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                        (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3) ) * p0_cart(i,j,k)

             rhoh_force(i,j,k) = divup - p0divu

          end do
       end do
    end do
    !$OMP END PARALLEL DO
    !
    ! psi should always be in the force if we are doing the final update
    ! For prediction, it should not be in the force if we are predicting
    ! (rho h)', but should be there if we are predicting h
    !
    if ((is_prediction .AND. enthalpy_pred_type == predict_h) .OR. &
        (is_prediction .AND. enthalpy_pred_type == predict_rhoh) .OR. &
        (.NOT. is_prediction)) then

       allocate(psi_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

       call put_1d_array_on_cart_3d_sphr(.false.,.false.,psi,psi_cart,lo,hi,dx,0)

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + psi_cart(i,j,k,1)
             enddo
          enddo
       enddo       
       !$OMP END PARALLEL DO
       deallocate(psi_cart)
    endif

    if (add_thermal) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + thermal(i,j,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

  end subroutine mkrhohforce_3d_sphr

  subroutine mkhprimeforce(mla, sold, snew, scal_force, is_prediction, thermal, umac, &
                           p0_old, p0_new, h0_old, h0_new, psi, dx, add_thermal, the_bc_level)

    use bl_prof_module
    use variables, only: foextrap_comp, rhoh_comp, rho_comp
    use geometry, only: spherical, polar, nr_fine
    use fill_3d_module, only: put_1d_array_on_cart
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    type(multifab) , intent(inout) :: scal_force(:)
    type(multifab) , intent(in   ) :: sold(:), snew(:)
    logical        , intent(in   ) :: is_prediction
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: h0_old(:,0:), h0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    logical        , intent(in   ) :: add_thermal
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer                  :: i,n
    integer                  :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer                  :: ng_f,ng_um,ng_th,ng_s
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: tp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    real(kind=dp_t), pointer :: hp(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)

    type(multifab)  :: p0_cart(mla%nlevel)
    type(multifab)  :: h0_cart(mla%nlevel)

    real(kind=dp_t) :: p0_nph(mla%nlevel,0:nr_fine-1)
    real(kind=dp_t) :: h0_nph(mla%nlevel,0:nr_fine-1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mkhprimeforce")

    ! if we are doing the prediction, then it only makes sense to be in
    ! this routine if the quantity we are predicting is h'
    if (is_prediction .AND. enthalpy_pred_type .ne. predict_hprime) then
       call bl_error("ERROR: should not call mkhprimeforce if enthalpy_pred_type .ne. predict_hprime")
    endif

    dm = mla%dim
    nlevs = mla%nlevel

    ng_f  = nghost(scal_force(1))
    ng_um = nghost(umac(1,1))
    ng_th = nghost(thermal(1))
    ng_s  = nghost(sold(1))

    if (spherical .eq. 1 .or. polar .eq. 1) then
       do n = 1,nlevs
          p0_nph(n,:) = HALF * (p0_old(n,:) + p0_new(n,:))
          call multifab_build(p0_cart(n),mla%la(n),1,1)
          h0_nph(n,:) = HALF * (h0_old(n,:) + h0_new(n,:))
          call multifab_build(h0_cart(n),mla%la(n),1,1)
       end do
    end if

    do n=1,nlevs

       if (spherical .eq. 1 .or. polar .eq. 1) then
          call put_1d_array_on_cart(p0_nph,p0_cart,foextrap_comp,.false.,.false.,&
                                    dx,the_bc_level,mla)
          call put_1d_array_on_cart(h0_nph,h0_cart,foextrap_comp,.false.,.false.,&
                                    dx,the_bc_level,mla)
       end if

       do i=1, nfabs(scal_force(n))
          fp => dataptr(scal_force(n), i)
          ump => dataptr(umac(n,1),i)
          vmp => dataptr(umac(n,2),i)
          tp  => dataptr(thermal(n),i)
          sop => dataptr(sold(n), i)
          snp => dataptr(snew(n), i)
          lo = lwb(get_box(scal_force(n),i))
          hi = upb(get_box(scal_force(n),i))
          select case (dm)
          case (2)
             if (polar .eq. 0) then
                call bl_error("mkscalforce: mkhprime force not written for 2d")
             else
                pp  => dataptr(p0_cart(n),i)
                hp  => dataptr(h0_cart(n),i)
                call mkhprimeforce_2d_polar(fp(:,:,1,rhoh_comp), ng_f, is_prediction, &
                                           ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                           tp(:,:,1,1), ng_th, &
                                           sop(:,:,1,rho_comp), snp(:,:,1,rho_comp), ng_s, &
                                           pp(:,:,1,1), hp(:,:,1,1), &
                                           lo, hi, dx(n,:), &
                                           psi(1,:), add_thermal)
             end if
          case(3)
             wmp  => dataptr(umac(n,3), i)
             if (spherical .eq. 0) then
                call bl_error("mkscalforce: mkhprime force not written for 3d plane parallel")
             else
                pp  => dataptr(p0_cart(n),i)
                hp  => dataptr(h0_cart(n),i)
                call mkhprimeforce_3d_sphr(fp(:,:,:,rhoh_comp), ng_f, is_prediction, &
                                           ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                           tp(:,:,:,1), ng_th, &
                                           sop(:,:,:,rho_comp), snp(:,:,:,rho_comp), ng_s, &
                                           pp(:,:,:,1), hp(:,:,:,1), &
                                           lo, hi, dx(n,:), &
                                           psi(1,:), add_thermal)
             end if
          end select
       end do

    end do

    if (spherical .eq. 1 .or. polar .eq. 1) then
       do n = 1,nlevs
          call destroy(p0_cart(n))
          call destroy(h0_cart(n))
       end do
    end if

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,scal_force,mla%mba%rr,the_bc_level, &
                              icomp=rhoh_comp, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=scal_force(1)%ng)

    call destroy(bpt)
    
  end subroutine mkhprimeforce

  subroutine mkhprimeforce_2d_polar(rhoh_force,ng_f,is_prediction,umac,vmac,ng_um, &
                                   thermal,ng_th,rhoold,rhonew,ng_s,p0_cart,h0_cart,lo,hi, &
                                   dx,psi,add_thermal)

    use fill_3d_module
    use pred_parameters

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer,         intent(in   ) :: lo(:),hi(:),ng_f,ng_um,ng_th,ng_s
    logical,         intent(in   ) :: is_prediction
    real(kind=dp_t), intent(  out) :: rhoh_force(lo(1)-ng_f :,lo(2)-ng_f :)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:,lo(2)-ng_th:)
    real(kind=dp_t), intent(in   ) ::     rhoold(lo(1)-ng_s :,lo(2)-ng_s :)
    real(kind=dp_t), intent(in   ) ::     rhonew(lo(1)-ng_s :,lo(2)-ng_s :)
    real(kind=dp_t), intent(in   ) ::    p0_cart(lo(1)-   1 :,lo(2)-   1 :)
    real(kind=dp_t), intent(in   ) ::    h0_cart(lo(1)-   1 :,lo(2)-   1 :)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    logical        , intent(in   ) :: add_thermal

    real(kind=dp_t), allocatable :: psi_cart(:,:,:)

    real(kind=dp_t) :: p0_lox,p0_hix,p0_loy,p0_hiy
    real(kind=dp_t) :: h0_lox,h0_hix,h0_loy,h0_hiy
    real(kind=dp_t) :: rhoavg
    real(kind=dp_t) :: divup, p0divu
    real(kind=dp_t) :: divuh, h0divu
    integer         :: i,j
    !
    ! Here we make u grad p = div (u p) - p div (u)
    !
    !$OMP PARALLEL DO PRIVATE(i,j,p0_lox,p0_hix,p0_loy,p0_hiy,divup,p0divu,rhoavg)
    do j = lo(2),hi(2)
        do i = lo(1),hi(1)

            p0_lox = HALF * (p0_cart(i,j) + p0_cart(i-1,j))
            p0_hix = HALF * (p0_cart(i,j) + p0_cart(i+1,j))
            p0_loy = HALF * (p0_cart(i,j) + p0_cart(i,j-1))
            p0_hiy = HALF * (p0_cart(i,j) + p0_cart(i,j+1))

            divup = (umac(i+1,j) * p0_hix - umac(i,j) * p0_lox) / dx(1) + &
                    (vmac(i,j+1) * p0_hiy - vmac(i,j) * p0_loy) / dx(2) 
                    
            p0divu = ( (umac(i+1,j) - umac(i,j)) / dx(1) + &
                    (vmac(i,j+1) - vmac(i,j)) / dx(2) ) * p0_cart(i,j)

            rhoavg = 0.5d0* (rhoold(i,j) + rhonew(i,j))

            rhoh_force(i,j) = (divup - p0divu) / rhoavg

        end do
    end do
    !$OMP END PARALLEL DO
    
    !
    ! Here we make u grad h_0 = div (u h_0) - h_0 div (u)
    !
    !$OMP PARALLEL DO PRIVATE(i,j,h0_lox,h0_hix,h0_loy,h0_hiy,divuh,h0divu)
    do j = lo(2),hi(2)
        do i = lo(1),hi(1)

            h0_lox = HALF * (h0_cart(i,j) + h0_cart(i-1,j))
            h0_hix = HALF * (h0_cart(i,j) + h0_cart(i+1,j))
            h0_loy = HALF * (h0_cart(i,j) + h0_cart(i,j-1))
            h0_hiy = HALF * (h0_cart(i,j) + h0_cart(i,j+1))
    
            divuh = (umac(i+1,j) * h0_hix - umac(i,j) * h0_lox) / dx(1) + &
                    (vmac(i,j+1) * h0_hiy - vmac(i,j) * h0_loy) / dx(2) 
                    
            h0divu = ( (umac(i+1,j) - umac(i,j)) / dx(1) + &
                    (vmac(i,j+1) - vmac(i,j)) / dx(2) ) * h0_cart(i,j)

            rhoh_force(i,j) = rhoh_force(i,j) - divuh + h0divu

        end do
    end do

    ! psi should always be in the force if we are doing the final update
    ! For prediction, it should not be in the force if we are predicting
    ! (rho h)', but should be there if we are predicting h
    if (.NOT. is_prediction) then

       allocate(psi_cart(lo(1):hi(1),lo(2):hi(2),1))
       call put_1d_array_on_cart_2d_polar(.false.,.false.,psi,psi_cart,lo,hi,dx,0)

       !$OMP PARALLEL DO PRIVATE(i,j)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                rhoh_force(i,j) = rhoh_force(i,j) + psi_cart(i,j,1)
            enddo
        enddo
       !$OMP END PARALLEL DO
       deallocate(psi_cart)
    endif

    if (add_thermal) then
       !$OMP PARALLEL DO PRIVATE(i,j)
        do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                rhoh_force(i,j) = rhoh_force(i,j) + thermal(i,j)
            end do
        end do
       !$OMP END PARALLEL DO
    end if

  end subroutine mkhprimeforce_2d_polar
  
  
  subroutine mkhprimeforce_3d_sphr(rhoh_force,ng_f,is_prediction,umac,vmac,wmac,ng_um, &
                                   thermal,ng_th,rhoold,rhonew,ng_s,p0_cart,h0_cart,lo,hi, &
                                   dx,psi,add_thermal)

    use fill_3d_module
    use pred_parameters

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer,         intent(in   ) :: lo(:),hi(:),ng_f,ng_um,ng_th,ng_s
    logical,         intent(in   ) :: is_prediction
    real(kind=dp_t), intent(  out) :: rhoh_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:,lo(2)-ng_th:,lo(3)-ng_th:)
    real(kind=dp_t), intent(in   ) ::     rhoold(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::     rhonew(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::    p0_cart(lo(1)-   1 :,lo(2)-   1 :,lo(3)-   1 :)
    real(kind=dp_t), intent(in   ) ::    h0_cart(lo(1)-   1 :,lo(2)-   1 :,lo(3)-   1 :)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    logical        , intent(in   ) :: add_thermal

    real(kind=dp_t), allocatable :: psi_cart(:,:,:,:)

    real(kind=dp_t) :: p0_lox,p0_hix,p0_loy,p0_hiy,p0_loz,p0_hiz
    real(kind=dp_t) :: h0_lox,h0_hix,h0_loy,h0_hiy,h0_loz,h0_hiz
    real(kind=dp_t) :: rhoavg
    real(kind=dp_t) :: divup, p0divu
    real(kind=dp_t) :: divuh, h0divu
    integer         :: i,j,k
    !
    ! Here we make u grad p = div (u p) - p div (u)
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k,p0_lox,p0_hix,p0_loy,p0_hiy,p0_loz,p0_hiz,divup,p0divu,rhoavg)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             p0_lox = HALF * (p0_cart(i,j,k) + p0_cart(i-1,j,k))
             p0_hix = HALF * (p0_cart(i,j,k) + p0_cart(i+1,j,k))
             p0_loy = HALF * (p0_cart(i,j,k) + p0_cart(i,j-1,k))
             p0_hiy = HALF * (p0_cart(i,j,k) + p0_cart(i,j+1,k))
             p0_loz = HALF * (p0_cart(i,j,k) + p0_cart(i,j,k-1))
             p0_hiz = HALF * (p0_cart(i,j,k) + p0_cart(i,j,k+1))

             divup = (umac(i+1,j,k) * p0_hix - umac(i,j,k) * p0_lox) / dx(1) + &
                     (vmac(i,j+1,k) * p0_hiy - vmac(i,j,k) * p0_loy) / dx(2) + &
                     (wmac(i,j,k+1) * p0_hiz - wmac(i,j,k) * p0_loz) / dx(3)

             p0divu = ( (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                        (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                        (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3) ) * p0_cart(i,j,k)

             rhoavg = 0.5d0* (rhoold(i,j,k) + rhonew(i,j,k))

             rhoh_force(i,j,k) = (divup - p0divu) / rhoavg

          end do
       end do
    end do
    !$OMP END PARALLEL DO
    !
    ! Here we make u grad h_0 = div (u h_0) - h_0 div (u)
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k,h0_lox,h0_hix,h0_loy,h0_hiy,h0_loz,h0_hiz,divuh,h0divu)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             h0_lox = HALF * (h0_cart(i,j,k) + h0_cart(i-1,j,k))
             h0_hix = HALF * (h0_cart(i,j,k) + h0_cart(i+1,j,k))
             h0_loy = HALF * (h0_cart(i,j,k) + h0_cart(i,j-1,k))
             h0_hiy = HALF * (h0_cart(i,j,k) + h0_cart(i,j+1,k))
             h0_loz = HALF * (h0_cart(i,j,k) + h0_cart(i,j,k-1))
             h0_hiz = HALF * (h0_cart(i,j,k) + h0_cart(i,j,k+1))

             divuh = (umac(i+1,j,k) * h0_hix - umac(i,j,k) * h0_lox) / dx(1) + &
                     (vmac(i,j+1,k) * h0_hiy - vmac(i,j,k) * h0_loy) / dx(2) + &
                     (wmac(i,j,k+1) * h0_hiz - wmac(i,j,k) * h0_loz) / dx(3)

             h0divu = ( (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                        (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                        (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3) ) * h0_cart(i,j,k)

             rhoh_force(i,j,k) = rhoh_force(i,j,k) - divuh + h0divu

          end do
       end do
    end do

    ! psi should always be in the force if we are doing the final update
    ! For prediction, it should not be in the force if we are predicting
    ! (rho h)', but should be there if we are predicting h
    if (.NOT. is_prediction) then

       allocate(psi_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(.false.,.false.,psi,psi_cart,lo,hi,dx,0)

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + psi_cart(i,j,k,1)
             enddo
          enddo
       enddo       
       !$OMP END PARALLEL DO
       deallocate(psi_cart)
    endif

    if (add_thermal) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + thermal(i,j,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

  end subroutine mkhprimeforce_3d_sphr

  subroutine mktempforce(mla,temp_force,umac,s,thermal,p0_old,p0_new,psi,dx,the_bc_level)

    use bl_prof_module
    use variables, only: foextrap_comp, temp_comp
    use geometry, only: spherical, polar, nr_fine
    use fill_3d_module, only: put_1d_array_on_cart

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: temp_force(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: thermal(:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer         :: i,n,ng_f,ng_um,ng_s,ng_th
    type(multifab)  :: p0_cart(mla%nlevel)
    real(kind=dp_t) :: p0_nph(mla%nlevel,0:nr_fine-1)
    integer         :: lo(mla%dim),hi(mla%dim),dm,nlevs

    real(kind=dp_t), pointer :: tp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mktempforce")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_f  = nghost(temp_force(1))
    ng_um = nghost(umac(1,1))
    ng_s  = nghost(s(1))
    ng_th = nghost(thermal(1))

    if (spherical .eq. 1 .or. polar .eq. 1) then
       p0_nph(1,:) = HALF * (p0_old(1,:) + p0_new(1,:))
       do n = 1,nlevs
          call multifab_build(p0_cart(n),mla%la(n),1,1)
       end do
       call put_1d_array_on_cart(p0_nph,p0_cart,foextrap_comp,.false.,.false.,&
                                 dx,the_bc_level,mla)
    end if

    do n=1,nlevs

       do i=1, nfabs(temp_force(n))
          fp  => dataptr(temp_force(n),i)
          ump => dataptr(umac(n,1),i)
          sp  => dataptr(s(n),i)
          tp  => dataptr(thermal(n),i)

          lo  =  lwb(get_box(s(n),i))
          hi  =  upb(get_box(s(n),i))

          select case (dm)
          case (1)
             call mktempforce_1d(n, fp(:,1,1,temp_comp), ng_f, sp(:,1,1,:), ng_s, &
                                 ump(:,1,1,1), ng_um, tp(:,1,1,1), ng_th, lo, hi, &
                                 p0_old(n,:), p0_new(n,:), psi(n,:))

          case (2)
             vmp => dataptr(umac(n,2),i)
             if (polar .eq. 1) then
                pp  => dataptr(p0_cart(n),i)
                call mktempforce_2d_polar(fp(:,:,1,temp_comp), ng_f, sp(:,:,1,:), ng_s, &
                                         ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                         tp(:,:,1,1), ng_th, lo, hi, &
                                         pp(:,:,1,1), psi(1,:), dx(n,:))
             
             else
                call mktempforce_2d(n, fp(:,:,1,temp_comp), ng_f, sp(:,:,1,:), ng_s, &
                                    vmp(:,:,1,1), ng_um, tp(:,:,1,1), ng_th, lo, hi, &
                                    p0_old(n,:), p0_new(n,:), psi(n,:))
             end if
          case(3)
             vmp => dataptr(umac(n,2),i)
             wmp => dataptr(umac(n,3),i)
             if (spherical .eq. 1) then
                pp  => dataptr(p0_cart(n),i)
                call mktempforce_3d_sphr(fp(:,:,:,temp_comp), ng_f, sp(:,:,:,:), ng_s, &
                                         ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                         tp(:,:,:,1), ng_th, lo, hi, &
                                         pp(:,:,:,1), psi(1,:), dx(n,:))
             else
                call mktempforce_3d(n, fp(:,:,:,temp_comp), ng_f, sp(:,:,:,:), ng_s, &
                                    wmp(:,:,:,1), ng_um, tp(:,:,:,1), ng_th, lo, hi, &
                                    p0_old(n,:), p0_new(n,:), psi(n,:))
             end if
          end select
       end do
   
    end do

    if (spherical .eq. 1 .or. polar .eq. 1) then
       do n = 1,nlevs
          call destroy(p0_cart(n))
       end do
    end if

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,temp_force,mla%mba%rr,the_bc_level, &
                              icomp=temp_comp, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=temp_force(1)%ng)

    call destroy(bpt)

  end subroutine mktempforce

  subroutine mktempforce_1d(n, temp_force, ng_f, s, ng_s, wmac, ng_um, thermal, ng_th, &
                            lo, hi, p0_old, p0_new, psi)

    use geometry, only: dr, nr
    use variables, only: temp_comp, rho_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec

    ! compute the source terms for temperature

    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the temp_force for the update, they will be used to time-center.

    integer,         intent(in   ) :: lo(:),hi(:),ng_f,ng_s,ng_um,ng_th
    integer,         intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: temp_force(lo(1)-ng_f :)
    real(kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s :,:)
    real(kind=dp_t), intent(in   ) ::       wmac(lo(1)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), p0_new(0:), psi(0:)

    integer :: i

    real(kind=dp_t) :: gradp0, wadv, dhdp

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    do i = lo(1),hi(1)

       if (i.eq.0) then
          gradp0 = HALF * ( p0_old(i+1) + p0_new(i+1) &
                           -p0_old(i  ) - p0_new(i  ) ) / dr(n)
       else if (i.eq.nr(n)-1) then
          gradp0 = HALF * ( p0_old(i  ) + p0_new(i  ) &
                           -p0_old(i-1) - p0_new(i-1) ) / dr(n)
       else
          gradp0 = FOURTH * ( p0_old(i+1) + p0_new(i+1) &
                             -p0_old(i-1) - p0_new(i-1) ) / dr(n)
       end if

       eos_state%T     = s(i,temp_comp)
       eos_state%rho   = s(i,rho_comp)
       eos_state%xn(:) = s(i,spec_comp:spec_comp+nspec-1) / s(i,rho_comp)

       pt_index(:) = (/i, -1, -1/)

       ! dens, temp, xmass inputs
       call eos(eos_input_rt, eos_state, pt_index)

       dhdp = ONE / s(i,rho_comp) + ( s(i,rho_comp) * eos_state%dedr - &
                                      eos_state%p / s(i,rho_comp) ) &
                                     / ( s(i,rho_comp) * eos_state%dpdr )

       wadv = HALF*(wmac(i)+wmac(i+1))

       temp_force(i) =  thermal(i) + (ONE - s(i,rho_comp) * dhdp) * &
                          (wadv * gradp0 + psi(i))
       temp_force(i) = temp_force(i) / (eos_state%cp * s(i,rho_comp))

    end do

  end subroutine mktempforce_1d

  subroutine mktempforce_2d(n, temp_force, ng_f, s, ng_s, wmac, ng_um, thermal, ng_th, &
                            lo, hi, p0_old, p0_new, psi)

    use geometry, only: dr, nr
    use variables, only: temp_comp, rho_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec

    ! compute the source terms for temperature

    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the temp_force for the update, they will be used to time-center.

    integer,         intent(in   ) :: lo(:),hi(:),ng_f,ng_s,ng_um,ng_th
    integer,         intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: temp_force(lo(1)-ng_f :,lo(2)-ng_f :)
    real(kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real(kind=dp_t), intent(in   ) ::       wmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:,lo(2)-ng_th:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), p0_new(0:), psi(0:)

    integer :: i,j

    real(kind=dp_t) :: gradp0, wadv, dhdp

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    do j = lo(2),hi(2)

       if (j.eq.0) then
          gradp0 = HALF * ( p0_old(j+1) + p0_new(j+1) &
                           -p0_old(j  ) - p0_new(j  ) ) / dr(n)
       else if (j.eq.nr(n)-1) then
          gradp0 = HALF * ( p0_old(j  ) + p0_new(j  ) &
                           -p0_old(j-1) - p0_new(j-1) ) / dr(n)
       else
          gradp0 = FOURTH * ( p0_old(j+1) + p0_new(j+1) &
                             -p0_old(j-1) - p0_new(j-1) ) / dr(n)
       end if

       do i = lo(1),hi(1)

          eos_state%T     = s(i,j,temp_comp)
          eos_state%rho   = s(i,j,rho_comp)
          eos_state%xn(:) = s(i,j,spec_comp:spec_comp+nspec-1) / s(i,j,rho_comp)

          pt_index(:) = (/i, j, -1/)

          ! dens, temp, xmass inputs
         call eos(eos_input_rt, eos_state, pt_index)

         dhdp = ONE / s(i,j,rho_comp) + ( s(i,j,rho_comp) * eos_state%dedr - &
                                          eos_state%p / s(i,j,rho_comp) ) &
                                        / ( s(i,j,rho_comp) * eos_state%dpdr )

         wadv = HALF*(wmac(i,j)+wmac(i,j+1))

         temp_force(i,j) =  thermal(i,j) + (ONE - s(i,j,rho_comp) * dhdp) * &
                            (wadv * gradp0 + psi(j))
         temp_force(i,j) = temp_force(i,j) / (eos_state%cp * s(i,j,rho_comp))

       end do
    end do

  end subroutine mktempforce_2d

  subroutine mktempforce_2d_polar(temp_force, ng_f, s, ng_s, umac, vmac, ng_um, &
                                 thermal, ng_th, lo, hi, p0_cart, psi, dx)

    use fill_3d_module
    use variables, only: temp_comp, rho_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    use pred_parameters

    ! compute the source terms for temperature

    integer,         intent(in   ) :: lo(:),hi(:),ng_f,ng_s,ng_um,ng_th
    real(kind=dp_t), intent(  out) :: temp_force(lo(1)-ng_f :,lo(2)-ng_f :)
    real(kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:,lo(2)-ng_th:)
    real(kind=dp_t), intent(in   ) ::    p0_cart(lo(1)-1    :,lo(2)-1    :)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer :: i,j
    real(kind=dp_t) :: p0_lox,p0_hix,p0_loy,p0_hiy
    real(kind=dp_t) :: divup,p0divu,ugradp,dhdp
    real(kind=dp_t), allocatable :: psi_cart(:,:,:)

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    allocate(psi_cart(lo(1):hi(1),lo(2):hi(2),1))

    call put_1d_array_on_cart_2d_polar(.false.,.false.,psi,psi_cart,lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,dhdp,p0_lox,p0_hix,p0_loy,p0_hiy,divup,p0divu,ugradp,eos_state,pt_index)
    do j = lo(2),hi(2)
        do i = lo(1),hi(1)
            
            eos_state%T     = s(i,j,temp_comp)
            eos_state%rho   = s(i,j,rho_comp)
            eos_state%xn(:) = s(i,j,spec_comp:spec_comp+nspec-1) / s(i,j,rho_comp)

            pt_index(:) = (/i, j, -1/)
            
            ! dens, temp, xmass inputs
            call eos(eos_input_rt, eos_state, pt_index)
            
            dhdp = ONE / s(i,j,rho_comp) + ( s(i,j,rho_comp) * eos_state%dedr - &
                                            eos_state%p / s(i,j,rho_comp) ) &
                                            / ( s(i,j,rho_comp) * eos_state%dpdr )

            p0_lox = HALF * (p0_cart(i,j) + p0_cart(i-1,j)) 
            p0_hix = HALF * (p0_cart(i,j) + p0_cart(i+1,j)) 
            p0_loy = HALF * (p0_cart(i,j) + p0_cart(i,j-1)) 
            p0_hiy = HALF * (p0_cart(i,j) + p0_cart(i,j+1))  
            
            divup = (umac(i+1,j) * p0_hix - umac(i,j) * p0_lox) / dx(1) + &
                    (vmac(i,j+1) * p0_hiy - vmac(i,j) * p0_loy) / dx(2)

            p0divu = ( (umac(i+1,j) - umac(i,j)) / dx(1) + &
                    (vmac(i,j+1) - vmac(i,j)) / dx(2) ) * p0_cart(i,j)
            
            ugradp = divup - p0divu
            
            temp_force(i,j) =  thermal(i,j) + &
                (ONE - s(i,j,rho_comp) * dhdp) * (ugradp + psi_cart(i,j,1))

            temp_force(i,j) = temp_force(i,j) / (eos_state%cp * s(i,j,rho_comp))
            
        end do
    end do
    !$OMP END PARALLEL DO

    deallocate(psi_cart)

  end subroutine mktempforce_2d_polar
  
  
  subroutine mktempforce_3d(n, temp_force, ng_f, s, ng_s, wmac, ng_um, thermal, ng_th, &
                            lo, hi, p0_old, p0_new, psi)

    use geometry,  only: dr, nr
    use variables, only: temp_comp, rho_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec

    ! compute the source terms for temperature

    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the temp_force for the update, they will be used to time-center.

    integer,         intent(in   ) :: lo(:),hi(:),n,ng_f,ng_s,ng_um,ng_th
    real(kind=dp_t), intent(  out) :: temp_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
    real(kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real(kind=dp_t), intent(in   ) ::       wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:,lo(2)-ng_th:,lo(3)-ng_th:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), p0_new(0:), psi(0:)

    integer         :: i,j,k
    real(kind=dp_t) :: dhdp, gradp0, wadv

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    !$OMP PARALLEL DO PRIVATE(i,j,k,gradp0,eos_state,pt_index,dhdp,wadv)
    do k = lo(3),hi(3)
       if (k.eq.0) then
          gradp0 = HALF * ( p0_old(k+1) + p0_new(k+1) &
                           -p0_old(k  ) - p0_new(k  ) ) / dr(n)
       else if (k.eq.nr(n)-1) then
          gradp0 = HALF * ( p0_old(k  ) + p0_new(k  ) &
                           -p0_old(k-1) - p0_new(k-1) ) / dr(n)
       else
          gradp0 = FOURTH * ( p0_old(k+1) + p0_new(k+1) &
                             -p0_old(k-1) - p0_new(k-1) ) / dr(n)
       end if

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             eos_state%T     = s(i,j,k,temp_comp)
             eos_state%rho   = s(i,j,k,rho_comp)
             eos_state%xn(:) = s(i,j,k,spec_comp:spec_comp+nspec-1) / s(i,j,k,rho_comp)

             pt_index(:) = (/i, j, k/)
             
             ! dens, temp, xmass inputs
             call eos(eos_input_rt, eos_state, pt_index)
             
             dhdp = ONE / s(i,j,k,rho_comp) + ( s(i,j,k,rho_comp) * eos_state%dedr - &
                  eos_state%p / s(i,j,k,rho_comp) ) / ( s(i,j,k,rho_comp) * eos_state%dpdr )
             
             wadv = HALF * (wmac(i,j,k+1) + wmac(i,j,k))
             
             temp_force(i,j,k) =  thermal(i,j,k) + &
                  (ONE - s(i,j,k,rho_comp) * dhdp) * (wadv * gradp0 + psi(k))
             
             temp_force(i,j,k) = temp_force(i,j,k) / (eos_state%cp * s(i,j,k,rho_comp))
             
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine mktempforce_3d

  subroutine mktempforce_3d_sphr(temp_force, ng_f, s, ng_s, umac, vmac, wmac, ng_um, &
                                 thermal, ng_th, lo, hi, p0_cart, psi, dx)

    use fill_3d_module
    use variables, only: temp_comp, rho_comp, spec_comp
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    use pred_parameters

    ! compute the source terms for temperature

    integer,         intent(in   ) :: lo(:),hi(:),ng_f,ng_s,ng_um,ng_th
    real(kind=dp_t), intent(  out) :: temp_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
    real(kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real(kind=dp_t), intent(in   ) ::       umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::       wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::    thermal(lo(1)-ng_th:,lo(2)-ng_th:,lo(3)-ng_th:)
    real(kind=dp_t), intent(in   ) ::    p0_cart(lo(1)-1    :,lo(2)-1    :,lo(3)-1    :)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer :: i,j,k
    real(kind=dp_t) :: p0_lox,p0_hix,p0_loy,p0_hiy,p0_loz,p0_hiz
    real(kind=dp_t) :: divup,p0divu,ugradp,dhdp
    real(kind=dp_t), allocatable :: psi_cart(:,:,:,:)

    type (eos_t) :: eos_state
    integer :: pt_index(MAX_SPACEDIM)

    allocate(psi_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,psi,psi_cart,lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k,dhdp,p0_lox,p0_hix,p0_loy,p0_hiy,p0_loz,p0_hiz,divup,p0divu,ugradp,eos_state,pt_index)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             eos_state%T     = s(i,j,k,temp_comp)
             eos_state%rho   = s(i,j,k,rho_comp)
             eos_state%xn(:) = s(i,j,k,spec_comp:spec_comp+nspec-1) / s(i,j,k,rho_comp)

             pt_index(:) = (/i, j, k/)
             
             ! dens, temp, xmass inputs
             call eos(eos_input_rt, eos_state, pt_index)
             
             dhdp = ONE / s(i,j,k,rho_comp) + ( s(i,j,k,rho_comp) * eos_state%dedr - &
                                                eos_state%p / s(i,j,k,rho_comp) ) &
                                                / ( s(i,j,k,rho_comp) * eos_state%dpdr )

             p0_lox = HALF * (p0_cart(i,j,k) + p0_cart(i-1,j,k)) 
             p0_hix = HALF * (p0_cart(i,j,k) + p0_cart(i+1,j,k)) 
             p0_loy = HALF * (p0_cart(i,j,k) + p0_cart(i,j-1,k)) 
             p0_hiy = HALF * (p0_cart(i,j,k) + p0_cart(i,j+1,k)) 
             p0_loz = HALF * (p0_cart(i,j,k) + p0_cart(i,j,k-1)) 
             p0_hiz = HALF * (p0_cart(i,j,k) + p0_cart(i,j,k+1)) 
             
             divup = (umac(i+1,j,k) * p0_hix - umac(i,j,k) * p0_lox) / dx(1) + &
                     (vmac(i,j+1,k) * p0_hiy - vmac(i,j,k) * p0_loy) / dx(2) + &
                     (wmac(i,j,k+1) * p0_hiz - wmac(i,j,k) * p0_loz) / dx(3)

             p0divu = ( (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                        (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                        (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3) ) * p0_cart(i,j,k)
             
             ugradp = divup - p0divu
             
             temp_force(i,j,k) =  thermal(i,j,k) + &
                  (ONE - s(i,j,k,rho_comp) * dhdp) * (ugradp + psi_cart(i,j,k,1))

             temp_force(i,j,k) = temp_force(i,j,k) / (eos_state%cp * s(i,j,k,rho_comp))
             
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(psi_cart)

  end subroutine mktempforce_3d_sphr

end module mkscalforce_module
