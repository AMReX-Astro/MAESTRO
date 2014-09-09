module update_vel_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private
  public :: update_velocity

contains

  subroutine update_velocity(uold,unew,umac,uedge,force,w0,w0mac,dx,dt, &
                             sponge,mla,the_bc_level)

    use bl_prof_module
    use bl_constants_module
    use ml_restrict_fill_module
    use geometry, only: spherical

    type(multifab)    , intent(in   ) :: uold(:)
    type(multifab)    , intent(inout) :: unew(:)
    type(multifab)    , intent(in   ) :: umac(:,:)
    type(multifab)    , intent(in   ) :: uedge(:,:)
    type(multifab)    , intent(in   ) :: force(:)
    real (kind = dp_t), intent(in   ) :: w0(:,0:)
    type(multifab)    , intent(in   ) :: w0mac(:,:)
    real (kind = dp_t), intent(in   ) :: dx(:,:)
    real (kind = dp_t), intent(in   ) :: dt
    type(multifab)    , intent(in   ) :: sponge(:)
    type(ml_layout)   , intent(inout) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: i,n,n_1d
    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: ng_uo,ng_un,ng_um,ng_ue,ng_sp,ng_f,ng_w0

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    real(kind=dp_t), pointer:: unp(:,:,:,:)
    real(kind=dp_t), pointer:: ump(:,:,:,:)
    real(kind=dp_t), pointer:: vmp(:,:,:,:)
    real(kind=dp_t), pointer:: wmp(:,:,:,:)
    real(kind=dp_t), pointer:: uepx(:,:,:,:)
    real(kind=dp_t), pointer:: uepy(:,:,:,:)
    real(kind=dp_t), pointer:: uepz(:,:,:,:)
    real(kind=dp_t), pointer:: spp(:,:,:,:)
    real(kind=dp_t), pointer:: fp(:,:,:,:)
    real(kind=dp_t), pointer:: w0xp(:,:,:,:)
    real(kind=dp_t), pointer:: w0yp(:,:,:,:)
    real(kind=dp_t), pointer:: w0zp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "update_velocity")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_uo = nghost(uold(1))
    ng_un = nghost(unew(1))
    ng_um = nghost(umac(1,1)) ! note we are assuming that ng is the same for all directions
    ng_ue = nghost(uedge(1,1))
    ng_sp = nghost(sponge(1))
    ng_f  = nghost(force(1))
    ng_w0 = nghost(w0mac(1,1))

    do n = 1, nlevs

       do i = 1, nfabs(uold(n))
          uop  => dataptr(uold(n),i)
          unp  => dataptr(unew(n),i)
          ump  => dataptr(umac(n,1),i)
          uepx => dataptr(uedge(n,1),i)
          spp  => dataptr(sponge(n),i)
          fp   =>  dataptr(force(n),i)
          lo = lwb(get_box(uold(n),i))
          hi = upb(get_box(uold(n),i))
          select case (dm)
          case (1)
             call update_velocity_1d(uop(:,1,1,:), ng_uo, unp(:,1,1,:), ng_un, &
                                     ump(:,1,1,1),  ng_um, &
                                     uepx(:,1,1,:), ng_ue, &
                                     fp(:,1,1,:), ng_f, w0(n,:), &
                                     lo, hi, dx(n,:), dt, spp(:,1,1,1), ng_sp)
          case (2)
             vmp  => dataptr(umac(n,2),i)
             uepy => dataptr(uedge(n,2),i)
             call update_velocity_2d(uop(:,:,1,:), ng_uo, unp(:,:,1,:), ng_un, &
                                     ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                     uepx(:,:,1,:), uepy(:,:,1,:), ng_ue, &
                                     fp(:,:,1,:), ng_f, w0(n,:), &
                                     lo, hi, dx(n,:), dt, spp(:,:,1,1), ng_sp)
          case (3)
             vmp  => dataptr(umac(n,2),i)
             wmp   => dataptr(umac(n,3),i)
             uepy => dataptr(uedge(n,2),i)
             uepz  => dataptr(uedge(n,3),i)
             w0xp   => dataptr(w0mac(n,1),i)
             w0yp   => dataptr(w0mac(n,2),i)
             w0zp   => dataptr(w0mac(n,3),i)
             if (spherical .eq. 1) then
                n_1d = 1
             else
                n_1d = n
             end if
             call update_velocity_3d(uop(:,:,:,:), ng_uo, unp(:,:,:,:), ng_un, &
                                     ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                     uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), ng_ue, &
                                     fp(:,:,:,:), ng_f, w0(n_1d,:), &
                                     w0xp(:,:,:,1), w0yp(:,:,:,1), w0zp(:,:,:,1), &
                                     ng_w0, lo, hi, dx(n,:), dt, spp(:,:,:,1), ng_sp)
          end select
       end do

    enddo

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,unew,mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=1, &
                              nc=dm, &
                              ng=unew(1)%ng)

    call destroy(bpt)

  end subroutine update_velocity

  subroutine update_velocity_1d(uold,ng_uo,unew,ng_un,umac,ng_um,uedgex,ng_ue, &
                                force,ng_f,w0,lo,hi,dx,dt,sponge,ng_sp)

    use bl_constants_module
    use probin_module, only: do_sponge

    integer, intent(in) :: lo(:), hi(:), ng_uo, ng_un, ng_um, ng_ue, ng_f, ng_sp
    real (kind = dp_t), intent(in   ) ::   uold(lo(1)-ng_uo:,:)  
    real (kind = dp_t), intent(  out) ::   unew(lo(1)-ng_un:,:)  
    real (kind = dp_t), intent(in   ) ::   umac(lo(1)-ng_um:)  
    real (kind = dp_t), intent(in   ) :: uedgex(lo(1)-ng_ue:,:)  
    real (kind = dp_t), intent(in   ) ::  force(lo(1)-ng_f :,:)  
    real (kind = dp_t), intent(in   ) :: sponge(lo(1)-ng_sp:)
    real (kind = dp_t), intent(in   ) ::     w0(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: dt

    integer :: i
    real (kind = dp_t) ubar
    real (kind = dp_t) ugradu

    do i = lo(1), hi(1)

       ! create cell-centered Utilde
       ubar = HALF*(umac(i) + umac(i+1))

       ! create (Utilde dot grad) Utilde
       ugradu = ubar*(uedgex(i+1,1) - uedgex(i,1))/dx(1) 

       ! update with (Utilde dot grad) Utilde and force
       unew(i,1) = uold(i,1) - dt * ugradu + dt * force(i,1)

       ! subtract (w0 dot grad) Utilde term
       ubar = HALF*(w0(i) + w0(i+1))
       unew(i,1) = unew(i,1) - dt * ubar*(uedgex(i+1,1) - uedgex(i,1))/dx(1)

       ! Add the sponge
       if (do_sponge) unew(i,:) = unew(i,:) * sponge(i)

    enddo

  end subroutine update_velocity_1d

  subroutine update_velocity_2d(uold,ng_uo,unew,ng_un,umac,vmac,ng_um,uedgex,uedgey,ng_ue, &
                                force,ng_f,w0,lo,hi,dx,dt,sponge,ng_sp)

    use bl_constants_module
    use probin_module, only: do_sponge

    integer, intent(in) :: lo(:), hi(:), ng_uo, ng_un, ng_um, ng_ue, ng_f, ng_sp
    real (kind = dp_t), intent(in   ) ::   uold(lo(1)-ng_uo:,lo(2)-ng_uo:,:)  
    real (kind = dp_t), intent(  out) ::   unew(lo(1)-ng_un:,lo(2)-ng_un:,:)  
    real (kind = dp_t), intent(in   ) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:)  
    real (kind = dp_t), intent(in   ) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:)  
    real (kind = dp_t), intent(in   ) :: uedgex(lo(1)-ng_ue:,lo(2)-ng_ue:,:)  
    real (kind = dp_t), intent(in   ) :: uedgey(lo(1)-ng_ue:,lo(2)-ng_ue:,:)  
    real (kind = dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,:)  
    real (kind = dp_t), intent(in   ) :: sponge(lo(1)-ng_sp:,lo(2)-ng_sp:)
    real (kind = dp_t), intent(in   ) ::     w0(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: dt

    integer :: i, j
    real (kind = dp_t) ubar,vbar
    real (kind = dp_t) ugradu,ugradv

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! create cell-centered Utilde
          ubar = HALF*(umac(i,j) + umac(i+1,j))
          vbar = HALF*(vmac(i,j) + vmac(i,j+1))

          ! create (Utilde dot grad) Utilde
          ugradu = ubar*(uedgex(i+1,j,1) - uedgex(i,j,1))/dx(1) + &
               vbar*(uedgey(i,j+1,1) - uedgey(i,j,1))/dx(2)

          ugradv = ubar*(uedgex(i+1,j,2) - uedgex(i,j,2))/dx(1) + &
               vbar*(uedgey(i,j+1,2) - uedgey(i,j,2))/dx(2)

          ! update with (Utilde dot grad) Utilde and force
          unew(i,j,1) = uold(i,j,1) - dt * ugradu + dt * force(i,j,1)
          unew(i,j,2) = uold(i,j,2) - dt * ugradv + dt * force(i,j,2)

          ! subtract (w0 dot grad) Utilde term
          vbar = HALF*(w0(j) + w0(j+1))
          unew(i,j,:) = unew(i,j,:) - dt * vbar*(uedgey(i,j+1,:) - uedgey(i,j,:))/dx(2)

          ! Add the sponge
          if (do_sponge) unew(i,j,:) = unew(i,j,:) * sponge(i,j)

       enddo
    enddo

  end subroutine update_velocity_2d

  subroutine update_velocity_3d(uold,ng_uo,unew,ng_un,umac,vmac,wmac,ng_um, &
                                uedgex,uedgey,uedgez,ng_ue,force,ng_f, &
                                w0,w0macx,w0macy,w0macz,ng_w0,lo,hi,dx,dt, &
                                sponge,ng_sp)

    use fill_3d_module
    use geometry, only: spherical
    use bl_constants_module
    use probin_module, only: do_sponge

    integer, intent(in) :: lo(:), hi(:)
    integer, intent(in) :: ng_uo, ng_un, ng_um, ng_ue, ng_f, ng_w0, ng_sp
    real (kind = dp_t), intent(in   ) ::    uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
    real (kind = dp_t), intent(  out) ::    unew(lo(1)-ng_un:,lo(2)-ng_un:,lo(3)-ng_un:,:)
    real (kind = dp_t), intent(in   ) ::    umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real (kind = dp_t), intent(in   ) ::    vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real (kind = dp_t), intent(in   ) ::    wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real (kind = dp_t), intent(in   ) ::  uedgex(lo(1)-ng_ue:,lo(2)-ng_ue:,lo(3)-ng_ue:,:)
    real (kind = dp_t), intent(in   ) ::  uedgey(lo(1)-ng_ue:,lo(2)-ng_ue:,lo(3)-ng_ue:,:)
    real (kind = dp_t), intent(in   ) ::  uedgez(lo(1)-ng_ue:,lo(2)-ng_ue:,lo(3)-ng_ue:,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real (kind = dp_t), intent(in   ) ::  sponge(lo(1)-ng_sp:,lo(2)-ng_sp:,lo(3)-ng_sp:) 
    real (kind = dp_t), intent(in   ) ::      w0(0:)
    real (kind = dp_t), intent(in   ) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real (kind = dp_t), intent(in   ) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real (kind = dp_t), intent(in   ) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real (kind = dp_t), intent(in   ) ::      dx(:)
    real (kind = dp_t), intent(in   ) :: dt

    integer :: i, j, k
    real (kind = dp_t) ubar,vbar,wbar
    real (kind = dp_t) ugradu,ugradv,ugradw
    real (kind = dp_t) :: gradux,graduy,graduz
    real (kind = dp_t) :: gradvx,gradvy,gradvz
    real (kind = dp_t) :: gradwx,gradwy,gradwz
    real (kind = dp_t) :: w0_gradur,w0_gradvr,w0_gradwr

    ! 1) Subtract (Utilde dot grad) Utilde term from old Utilde
    ! 2) Add forcing term to new Utilde
    !$OMP PARALLEL DO PRIVATE(i,j,k,ubar,vbar,wbar,ugradu,ugradv,ugradw)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! create cell-centered Utilde
             ubar = HALF*(umac(i,j,k) + umac(i+1,j,k))
             vbar = HALF*(vmac(i,j,k) + vmac(i,j+1,k))
             wbar = HALF*(wmac(i,j,k) + wmac(i,j,k+1))

             ! create (Utilde dot grad) Utilde
             ugradu = ubar*(uedgex(i+1,j,k,1) - uedgex(i,j,k,1))/dx(1) + &
                  vbar*(uedgey(i,j+1,k,1) - uedgey(i,j,k,1))/dx(2) + &
                  wbar*(uedgez(i,j,k+1,1) - uedgez(i,j,k,1))/dx(3)

             ugradv = ubar*(uedgex(i+1,j,k,2) - uedgex(i,j,k,2))/dx(1) + &
                  vbar*(uedgey(i,j+1,k,2) - uedgey(i,j,k,2))/dx(2) + &
                  wbar*(uedgez(i,j,k+1,2) - uedgez(i,j,k,2))/dx(3)

             ugradw = ubar*(uedgex(i+1,j,k,3) - uedgex(i,j,k,3))/dx(1) + &
                  vbar*(uedgey(i,j+1,k,3) - uedgey(i,j,k,3))/dx(2) + &
                  wbar*(uedgez(i,j,k+1,3) - uedgez(i,j,k,3))/dx(3)

             ! update with (Utilde dot grad) Utilde and force
             unew(i,j,k,1) = uold(i,j,k,1) - dt * ugradu + dt * force(i,j,k,1)
             unew(i,j,k,2) = uold(i,j,k,2) - dt * ugradv + dt * force(i,j,k,2)
             unew(i,j,k,3) = uold(i,j,k,3) - dt * ugradw + dt * force(i,j,k,3)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    if (spherical .eq. 0) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,wbar)
       do k = lo(3), hi(3)

          ! subtract (w0 dot grad) Utilde term
          wbar = HALF*(w0(k) + w0(k+1))
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                unew(i,j,k,:) = unew(i,j,k,:) - dt * wbar*(uedgez(i,j,k+1,:) &
                     - uedgez(i,j,k,:))/dx(3)

                ! Add the sponge
                if (do_sponge) unew(i,j,k,:) = unew(i,j,k,:) * sponge(i,j,k)

             enddo
          enddo

       enddo
       !$OMP END PARALLEL DO

    else

       !$OMP PARALLEL DO PRIVATE(i,j,k,gradux,gradvx,gradwx,graduy,gradvy,gradwy) &
       !$OMP PRIVATE(graduz,gradvz,gradwz,w0_gradur,w0_gradvr,w0_gradwr)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Subtract (w0 dot grad) Utilde term from new Utilde
                gradux = (uedgex(i+1,j,k,1) - uedgex(i,j,k,1))/dx(1)
                gradvx = (uedgex(i+1,j,k,2) - uedgex(i,j,k,2))/dx(1)
                gradwx = (uedgex(i+1,j,k,3) - uedgex(i,j,k,3))/dx(1)

                graduy = (uedgey(i,j+1,k,1) - uedgey(i,j,k,1))/dx(2)
                gradvy = (uedgey(i,j+1,k,2) - uedgey(i,j,k,2))/dx(2)
                gradwy = (uedgey(i,j+1,k,3) - uedgey(i,j,k,3))/dx(2)

                graduz = (uedgez(i,j,k+1,1) - uedgez(i,j,k,1))/dx(3)
                gradvz = (uedgez(i,j,k+1,2) - uedgez(i,j,k,2))/dx(3)
                gradwz = (uedgez(i,j,k+1,3) - uedgez(i,j,k,3))/dx(3)

                w0_gradur = gradux * HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)) &
                          + graduy * HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)) &
                          + graduz * HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))

                w0_gradvr = gradvx * HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)) &
                          + gradvy * HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)) &
                          + gradvz * HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))

                w0_gradwr = gradwx * HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)) &
                          + gradwy * HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)) &
                          + gradwz * HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))

                unew(i,j,k,1) = unew(i,j,k,1) - dt * w0_gradur
                unew(i,j,k,2) = unew(i,j,k,2) - dt * w0_gradvr
                unew(i,j,k,3) = unew(i,j,k,3) - dt * w0_gradwr

                ! Add the sponge
                if (do_sponge) unew(i,j,k,:) = unew(i,j,k,:) * sponge(i,j,k)

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    end if

  end subroutine update_velocity_3d

end module update_vel_module
