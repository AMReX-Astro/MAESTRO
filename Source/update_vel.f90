module update_vel_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private
  public :: update_velocity

contains

  subroutine update_velocity(nlevs,uold,unew,umac,uedge,force,w0,w0_cart,w0_force, &
                             w0_force_cart,dx,dt,sponge,do_sponge,mla,the_bc_level)

    use bl_prof_module
    use bl_constants_module
    use ml_restriction_module
    use multifab_fill_ghost_module
    use multifab_physbc_module

    implicit none

    integer           , intent(in   ) :: nlevs
    type(multifab)    , intent(in   ) :: uold(:)
    type(multifab)    , intent(inout) :: unew(:)
    type(multifab)    , intent(in   ) :: umac(:,:)
    type(multifab)    , intent(in   ) :: uedge(:,:)
    type(multifab)    , intent(in   ) :: force(:)
    real (kind = dp_t), intent(in   ) :: w0(:,0:)
    type(multifab)    , intent(in   ) :: w0_cart(:)
    real (kind = dp_t), intent(in   ) :: w0_force(:,0:)
    type(multifab)    , intent(in   ) :: w0_force_cart(:)
    real (kind = dp_t), intent(in   ) :: dx(:,:)
    real (kind = dp_t), intent(in   ) :: dt
    type(multifab)    , intent(in   ) :: sponge(:)
    logical           , intent(in   ) :: do_sponge
    type(ml_layout)   , intent(inout) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: i,dm,ng,n
    integer :: lo(uold(1)%dim),hi(uold(1)%dim)

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
    real(kind=dp_t), pointer:: w0p(:,:,:,:)
    real(kind=dp_t), pointer:: w0fp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "update_velocity")

    dm = uold(1)%dim
    ng = uold(1)%ng

    do n = 1, nlevs

       do i = 1, uold(n)%nboxes
          if ( multifab_remote(uold(n),i) ) cycle
          uop  => dataptr(uold(n),i)
          unp  => dataptr(unew(n),i)
          ump  => dataptr(umac(n,1),i)
          vmp  => dataptr(umac(n,2),i)
          uepx => dataptr(uedge(n,1),i)
          uepy => dataptr(uedge(n,2),i)
          spp  => dataptr(sponge(n),i)
          fp   =>  dataptr(force(n),i)
          lo = lwb(get_box(uold(n),i))
          hi = upb(get_box(uold(n),i))
          select case (dm)
          case (2)
             call update_velocity_2d(uop(:,:,1,:), unp(:,:,1,:), &
                                     ump(:,:,1,1), vmp(:,:,1,1), &
                                     uepx(:,:,1,:), uepy(:,:,1,:), &
                                     fp(:,:,1,:), w0(n,:), w0_force(n,:), &
                                     lo, hi, ng, dx(n,:), dt, spp(:,:,1,1), &
                                     do_sponge)
          case (3)
             wmp   => dataptr(umac(n,3),i)
             uepz  => dataptr(uedge(n,3),i)
             w0p   => dataptr(w0_cart(n),i)
             w0fp  => dataptr(w0_force_cart(n),i)
             call update_velocity_3d(uop(:,:,:,:), unp(:,:,:,:), &
                                     ump(:,:,:,1), vmp(:,:,:,1), &
                                     wmp(:,:,:,1), uepx(:,:,:,:), &
                                     uepy(:,:,:,:), uepz(:,:,:,:), &
                                     fp(:,:,:,:), w0(n,:), w0p(:,:,:,:), &
                                     w0_force(n,:), w0fp(:,:,:,:), &
                                     lo, hi, ng, dx(n,:), dt, spp(:,:,:,1), &
                                     do_sponge)
          end select
       end do

       call multifab_fill_boundary(unew(n))
       call multifab_physbc(unew(n),1,1,dm,dx(n,:),the_bc_level(n))

    enddo ! end loop over levels

    do n=nlevs,2,-1
       call ml_cc_restriction(unew(n-1),unew(n),mla%mba%rr(n-1,:))
       call multifab_fill_ghost_cells(unew(n),unew(n-1), &
                                      ng,mla%mba%rr(n-1,:), &
                                      the_bc_level(n-1),the_bc_level(n), &
                                      1,1,dm)
    enddo

    call destroy(bpt)

  end subroutine update_velocity

  subroutine update_velocity_2d(uold,unew,umac,vmac,uedgex,uedgey,force,w0,w0_force, &
                                lo,hi,ng,dx,dt,sponge,do_sponge)

    use bl_constants_module

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) ::     uold(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(  out) ::     unew(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) ::     umac(lo(1)- 1:,lo(2)- 1:)  
    real (kind = dp_t), intent(in   ) ::     vmac(lo(1)- 1:,lo(2)- 1:)  
    real (kind = dp_t), intent(in   ) ::   uedgex(lo(1)   :,lo(2)   :,:)  
    real (kind = dp_t), intent(in   ) ::   uedgey(lo(1)   :,lo(2)   :,:)  
    real (kind = dp_t), intent(in   ) ::    force(lo(1)- 1:,lo(2)- 1:,:)  
    real (kind = dp_t), intent(in   ) ::   sponge(lo(1)   :,lo(2)   :  )
    real (kind = dp_t), intent(in   ) ::       w0(0:)
    real (kind = dp_t), intent(in   ) :: w0_force(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: dt
    logical           , intent(in   ) :: do_sponge

    integer :: i, j
    real (kind = dp_t) ubar,vbar
    real (kind = dp_t) ugradu,ugradv

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ubar = HALF*(umac(i,j) + umac(i+1,j))
          vbar = HALF*(vmac(i,j) + vmac(i,j+1))

          ugradu = ubar*(uedgex(i+1,j,1) - uedgex(i,j,1))/dx(1) + &
               vbar*(uedgey(i,j+1,1) - uedgey(i,j,1))/dx(2)

          ugradv = ubar*(uedgex(i+1,j,2) - uedgex(i,j,2))/dx(1) + &
               vbar*(uedgey(i,j+1,2) - uedgey(i,j,2))/dx(2)

          unew(i,j,1) = uold(i,j,1) - dt * ugradu + dt * force(i,j,1)
          unew(i,j,2) = uold(i,j,2) - dt * ugradv + dt * force(i,j,2)

          ! Add w dot grad w0 term to w.
          unew(i,j,2) = unew(i,j,2) - dt * vbar*(w0(j+1) - w0(j))/dx(2)

          ! Add w0 dot grad u term to u and w.
          vbar = HALF*(w0(j) + w0(j+1))
          unew(i,j,:) = unew(i,j,:) - dt * vbar*(uedgey(i,j+1,:) - uedgey(i,j,:))/dx(2)

          ! Add in the pi0 term.
          unew(i,j,2) = unew(i,j,2) - dt * w0_force(j)

          ! Add the sponge
          if (do_sponge) unew(i,j,:) = unew(i,j,:) * sponge(i,j)

       enddo
    enddo

  end subroutine update_velocity_2d

  subroutine update_velocity_3d(uold,unew,umac,vmac,wmac,uedgex,uedgey,uedgez, &
                                force,w0,w0_cart,w0_force,w0_force_cart,lo,hi,ng,dx,dt, &
                                sponge,do_sponge)

    use geometry, only: spherical
    use bl_constants_module

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) ::     uold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) ::     unew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::     umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
    real (kind = dp_t), intent(in   ) ::     vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
    real (kind = dp_t), intent(in   ) ::     wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
    real (kind = dp_t), intent(in   ) ::   uedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::   uedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::   uedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::    force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   sponge(lo(1)   :,lo(2)   :,lo(3)   :  ) 
    real (kind = dp_t), intent(in   ) ::       w0(0:)
    real (kind = dp_t), intent(in   ) ::  w0_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) :: w0_force(0:)
    real (kind = dp_t), intent(in   ) ::  w0_force_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: dt
    logical           , intent(in   ) :: do_sponge

    integer :: i, j, k
    real (kind = dp_t) ubar,vbar,wbar
    real (kind = dp_t) ugradu,ugradv,ugradw
    real (kind = dp_t) :: gradux,graduy,graduz
    real (kind = dp_t) :: gradvx,gradvy,gradvz
    real (kind = dp_t) :: gradwx,gradwy,gradwz
    real (kind = dp_t) :: w0_gradur,w0_gradvr,w0_gradwr
    real (kind = dp_t) :: gradw0

    ! 1) Subtract (Utilde dot grad) Utilde term from old Utilde
    ! 2) Add forcing term to new Utilde
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ubar = HALF*(umac(i,j,k) + umac(i+1,j,k))
             vbar = HALF*(vmac(i,j,k) + vmac(i,j+1,k))
             wbar = HALF*(wmac(i,j,k) + wmac(i,j,k+1))

             ugradu = ubar*(uedgex(i+1,j,k,1) - uedgex(i,j,k,1))/dx(1) + &
                  vbar*(uedgey(i,j+1,k,1) - uedgey(i,j,k,1))/dx(2) + &
                  wbar*(uedgez(i,j,k+1,1) - uedgez(i,j,k,1))/dx(3)

             ugradv = ubar*(uedgex(i+1,j,k,2) - uedgex(i,j,k,2))/dx(1) + &
                  vbar*(uedgey(i,j+1,k,2) - uedgey(i,j,k,2))/dx(2) + &
                  wbar*(uedgez(i,j,k+1,2) - uedgez(i,j,k,2))/dx(3)

             ugradw = ubar*(uedgex(i+1,j,k,3) - uedgex(i,j,k,3))/dx(1) + &
                  vbar*(uedgey(i,j+1,k,3) - uedgey(i,j,k,3))/dx(2) + &
                  wbar*(uedgez(i,j,k+1,3) - uedgez(i,j,k,3))/dx(3)

             unew(i,j,k,1) = uold(i,j,k,1) - dt * ugradu + dt * force(i,j,k,1)
             unew(i,j,k,2) = uold(i,j,k,2) - dt * ugradv + dt * force(i,j,k,2)
             unew(i,j,k,3) = uold(i,j,k,3) - dt * ugradw + dt * force(i,j,k,3)

          enddo
       enddo
    enddo

    ! A) Subtract (Utilde dot er) dot grad w0 term from new Utilde.
    ! B) Subtract w0 dot grad U term from new Utilde
    if (spherical .eq. 0) then

       do k = lo(3), hi(3)

          gradw0 = (w0(k+1) - w0(k)) /dx(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                wbar = HALF*(wmac(i,j,k) + wmac(i,j,k+1))
                unew(i,j,k,3) = unew(i,j,k,3) - dt * wbar * gradw0
             enddo
          enddo

          wbar = HALF*(w0(k) + w0(k+1))
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                unew(i,j,k,:) = unew(i,j,k,:) - dt * wbar*(uedgez(i,j,k+1,:) &
                     - uedgez(i,j,k,:))/dx(3)

                ! Add in the pi0 term.
                unew(i,j,k,3) = unew(i,j,k,3) - dt * w0_force(k)

                ! Add the sponge
                if (do_sponge) unew(i,j,k,:) = unew(i,j,k,:) * sponge(i,j,k)

             enddo
          enddo

       enddo

    else

       ! A) Subtract (Utilde dot er) dot grad w0 term from new Utilde.

       ! B) Subtract (w0 dot grad) U term from new Utilde

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                gradux = (uedgex(i+1,j,k,1) - uedgex(i,j,k,1))/dx(1)
                gradvx = (uedgex(i+1,j,k,2) - uedgex(i,j,k,2))/dx(1)
                gradwx = (uedgex(i+1,j,k,3) - uedgex(i,j,k,3))/dx(1)
                graduy = (uedgey(i,j+1,k,1) - uedgey(i,j,k,1))/dx(2)
                gradvy = (uedgey(i,j+1,k,2) - uedgey(i,j,k,2))/dx(2)
                gradwy = (uedgey(i,j+1,k,3) - uedgey(i,j,k,3))/dx(2)
                graduz = (uedgez(i,j,k+1,1) - uedgez(i,j,k,1))/dx(3)
                gradvz = (uedgez(i,j,k+1,2) - uedgez(i,j,k,2))/dx(3)
                gradwz = (uedgez(i,j,k+1,3) - uedgez(i,j,k,3))/dx(3)

                w0_gradur = gradux * w0_cart(i,j,k,1) &
                     + graduy * w0_cart(i,j,k,2) + graduz * w0_cart(i,j,k,3)
                w0_gradvr = gradvx * w0_cart(i,j,k,1) &
                     + gradvy * w0_cart(i,j,k,2) + gradvz * w0_cart(i,j,k,3)
                w0_gradwr = gradwx * w0_cart(i,j,k,1) &
                     + gradwy * w0_cart(i,j,k,2) + gradwz * w0_cart(i,j,k,3)

                unew(i,j,k,1) = unew(i,j,k,1) - dt * w0_gradur
                unew(i,j,k,2) = unew(i,j,k,2) - dt * w0_gradvr
                unew(i,j,k,3) = unew(i,j,k,3) - dt * w0_gradwr

                ! Add in the pi0 term.
                unew(i,j,k,:) = unew(i,j,k,:) - dt * w0_force_cart(i,j,k,:)

                ! Add the sponge
                if (do_sponge) unew(i,j,k,:) = unew(i,j,k,:) * sponge(i,j,k)

             enddo
          enddo
       enddo

    end if

  end subroutine update_velocity_3d

end module update_vel_module
