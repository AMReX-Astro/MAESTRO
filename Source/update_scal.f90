module update_scal_module

  use bl_types
  use multifab_module
  use bl_constants_module
  use fill_3d_module
  use addw0_module
  use geometry
  use variables
  use network
  use probin_module, ONLY: evolve_base_state

  implicit none

  private
  public :: update_scal
  
contains

  subroutine update_scal(which_step,nstart,nstop,sold,snew,umac,w0,w0_cart_vec,eta,sedge, &
                         scal_force,s0_old,s0_edge_old,s0_new,s0_edge_new, &
                         s0_old_cart,s0_new_cart,domlo,domhi,dx,dt,evolve_base_state)

    implicit none

    integer           , intent(in   ) :: which_step, nstart, nstop
    type(multifab)    , intent(in   ) :: sold
    type(multifab)    , intent(inout) :: snew
    type(multifab)    , intent(in   ) :: umac(:)
    real(kind=dp_t)   , intent(in   ) :: w0(0:)
    type(multifab)    , intent(in   ) :: w0_cart_vec
    real(kind=dp_t)   , intent(inout) :: eta(0:,:)
    type(multifab)    , intent(in   ) :: sedge(:)
    type(multifab)    , intent(in   ) :: scal_force
    real(kind = dp_t) , intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind = dp_t) , intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)
    type(multifab)    , intent(in   ) :: s0_old_cart
    type(multifab)    , intent(in   ) :: s0_new_cart
    integer           , intent(in   ) :: domlo(:),domhi(:)
    real(kind = dp_t) , intent(in   ) :: dx(:),dt
    logical           , intent(in   ) :: evolve_base_state
    
    ! local
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)
    real(kind=dp_t), pointer :: s0op(:,:,:,:)
    real(kind=dp_t), pointer :: s0np(:,:,:,:)

    integer :: lo(sold%dim),hi(sold%dim)
    integer :: i,ng,dm

    dm = sold%dim
    ng = sold%ng

    do i = 1, sold%nboxes
       if ( multifab_remote(sold,i) ) cycle
       sop => dataptr(sold,i)
       snp => dataptr(snew,i)
       ump => dataptr(umac(1),i)
       vmp => dataptr(umac(2),i)
       sepx => dataptr(sedge(1),i)
       sepy => dataptr(sedge(2),i)
       fp => dataptr(scal_force,i)
       lo =  lwb(get_box(sold,i))
       hi =  upb(get_box(sold,i))
       select case (dm)
       case (2)
          call update_scal_2d(which_step, nstart, nstop, &
                              sop(:,:,1,:), snp(:,:,1,:), &
                              ump(:,:,1,1), vmp(:,:,1,1), w0, eta, &
                              sepx(:,:,1,:), sepy(:,:,1,:), fp(:,:,1,:), &
                              s0_old(:,:), s0_edge_old(:,:), &
                              s0_new(:,:), s0_edge_new(:,:), &
                              lo, hi, ng, dx, dt, evolve_base_state)
       case (3)
          wmp => dataptr(umac(3),i)
          sepz => dataptr(sedge(3),i)
          w0p => dataptr(w0_cart_vec,i)
          if (spherical .eq. 0) then
             call update_scal_3d_cart(which_step, nstart, nstop, &
                                      sop(:,:,:,:), snp(:,:,:,:), &
                                      ump(:,:,:,1), vmp(:,:,:,1), &
                                      wmp(:,:,:,1), w0, eta, &
                                      sepx(:,:,:,:), sepy(:,:,:,:), &
                                      sepz(:,:,:,:), fp(:,:,:,:), &
                                      s0_old(:,:), s0_edge_old(:,:), &
                                      s0_new(:,:), s0_edge_new(:,:), &
                                      lo, hi, ng, dx, dt, evolve_base_state)
          else
             s0op => dataptr(s0_old_cart, i)
             s0np => dataptr(s0_new_cart, i)
             call update_scal_3d_sphr(which_step, nstart, nstop, &
                                      sop(:,:,:,:), snp(:,:,:,:), &
                                      ump(:,:,:,1), vmp(:,:,:,1), &
                                      wmp(:,:,:,1), w0p(:,:,:,:), &
                                      sepx(:,:,:,:), sepy(:,:,:,:), &
                                      sepz(:,:,:,:), fp(:,:,:,:), &
                                      s0_old(:,:), s0_new(:,:), &
                                      s0op(:,:,:,:), s0np(:,:,:,:), &
                                      lo, hi, domlo, domhi, ng, dx, dt, &
                                      evolve_base_state)
          end if
       end select
    end do

  end subroutine update_scal
  
  subroutine update_scal_2d(which_step,nstart,nstop,sold,snew,umac,vmac,w0,eta, &
                            sedgex,sedgey,force,base_old,base_old_edge, &
                            base_new,base_new_edge,lo,hi,ng,dx,dt,evolve_base_state)

    ! update each scalar in time.  Here, it is assumed that the edge
    ! states (sedgex and sedgey) are for the perturbational quantities.
    
    implicit none
    
    integer           , intent(in   ) :: which_step, nstart, nstop, lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)
    real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)
    real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)
    real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   base_old(0:,:), base_old_edge(0:,:)
    real (kind = dp_t), intent(in   ) ::   base_new(0:,:), base_new_edge(0:,:)
    real (kind = dp_t), intent(in   ) :: w0(0:)
    real (kind = dp_t), intent(inout) :: eta(0:,:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: evolve_base_state
    
    integer :: i, j, comp, comp2
    real (kind = dp_t) :: divsu,divbaseu,delta_base
    real (kind = dp_t) :: delta,frac,sum,fac
    real (kind = dp_t), allocatable :: smin(:),smax(:)
    
    fac = ONE / dble(hi(1)-lo(1)+1)
    
    if (evolve_base_state) then
      if (which_step .eq. 1) then
       do comp = nstart, nstop
          eta(:,comp) = ZERO
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)
                eta(j,comp) = eta(j,comp) + vmac(i,j)*sedgey(i,j,comp)
!               eta(j,comp) = 0.d0
             end do
             eta(j,comp) = eta(j,comp) * fac
          end do
       end do
      end if
    end if
    
    do comp = nstart, nstop
       do j = lo(2), hi(2)
          
          delta_base = base_new(j,comp) - base_old(j,comp)
          
          do i = lo(1), hi(1)
             
             divsu = (umac(i+1,j) * sedgex(i+1,j,comp) &
                  -umac(i  ,j) * sedgex(i  ,j,comp) ) / dx(1) + &
                  ((vmac(i,j+1)+w0(j+1)) * sedgey(i,j+1,comp) &
                  -(vmac(i,j  )+w0(j  )) * sedgey(i,j  ,comp) ) / dx(2)
             
             if (which_step .eq. 1) then
                divbaseu = (umac(i+1,j) - umac(i,j) ) * base_old(j,comp) / dx(1) &
                     +( vmac(i,j+1) * base_old_edge(j+1,comp) &
                     -vmac(i,j  ) * base_old_edge(j  ,comp) ) / dx(2)
             else
                divbaseu = HALF * (umac(i+1,j) - umac(i,j) ) * &
                     (base_old(j,comp)+base_new(j,comp)) / dx(1) &
                     +HALF * (vmac(i,j+1) * (base_old_edge(j+1,comp) &
                     + base_new_edge(j+1,comp)) &
                     -vmac(i,j  ) * (base_old_edge(j  ,comp) &
                     + base_new_edge(j  ,comp)) ) / dx(2)
             end if
             
             snew(i,j,comp) = sold(i,j,comp) + delta_base &
                  - dt * (divsu + divbaseu) + dt * force(i,j,comp)
             
             if (evolve_base_state) then
               if (which_step .eq. 2) then
                snew(i,j,comp) = snew(i,j,comp) + dt / dx(2) * (eta(j+1,comp)-eta(j,comp))
               end if
             end if
             
          enddo
       enddo
    enddo
    
    ! Define the update to rho as the sum of the updates to (rho X)_i
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       
       allocate(smin(nstart:nstop),smax(nstart:nstop))
       smin(:) =  1.e20
       smax(:) = -1.e20
       
       snew(:,:,rho_comp) = sold(:,:,rho_comp)
       
       do comp = nstart, nstop
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                snew(i,j,rho_comp) = snew(i,j,rho_comp) + (snew(i,j,comp)-sold(i,j,comp))
                smin(comp) = min(smin(comp),snew(i,j,comp))
                smax(comp) = max(smax(comp),snew(i,j,comp))
             enddo
          enddo
       enddo
       
       if (which_step .eq. 1) then
          eta(:,rho_comp) = ZERO
          do comp = nstart, nstop
             do j = lo(2), hi(2)+1
                eta(j,rho_comp) = eta(j,rho_comp) + eta(j,comp)
             end do
          end do
       end if
       
    end if
    
    ! Do not allow the species to leave here negative.
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       do comp = nstart, nstop
          if (smin(comp) .lt. ZERO) then 
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (snew(i,j,comp) .lt. ZERO) then
                      delta = -snew(i,j,comp)
                      sum = ZERO 
                      do comp2 = nstart, nstop
                         if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) sum = sum + snew(i,j,comp2)
                      enddo
                      do comp2 = nstart, nstop
                         if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                            frac = snew(i,j,comp2) / sum
                            snew(i,j,comp2) = snew(i,j,comp2) - frac * delta
                         end if
                      enddo
                      snew(i,j,comp) = ZERO
                   end if
                enddo
             enddo
          end if
       enddo
       deallocate(smin,smax)
    end if
    
  end subroutine update_scal_2d
  
  subroutine update_scal_3d_cart(which_step,nstart,nstop,sold,snew,umac,vmac,wmac,w0, &
                                 eta,sedgex,sedgey,sedgez,force,base_old,base_old_edge, &
                                 base_new,base_new_edge,lo,hi,ng,dx,dt,evolve_base_state)
    
    implicit none
    
    integer           , intent(in   ) :: which_step, nstart, nstop, lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(inout) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(inout) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(inout) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   base_old(0:,:), base_old_edge(0:,:)
    real (kind = dp_t), intent(in   ) ::   base_new(0:,:), base_new_edge(0:,:)
    real (kind = dp_t), intent(in   ) :: w0(0:)
    real (kind = dp_t), intent(inout) :: eta(0:,:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: evolve_base_state
    
    integer :: i, j, k, comp, comp2
    real (kind = dp_t) :: divsu,divbaseu
    real (kind = dp_t) :: delta,frac,sum,delta_base,fac
    real (kind = dp_t), allocatable :: smin(:),smax(:)
    
    fac = ONE / dble( (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1) )
    
    if (which_step .eq. 1) then
       do comp = nstart, nstop
          do k = lo(3), hi(3)+1
             eta(k,comp) = ZERO
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
!               eta(k,comp) = eta(k,comp) + wmac(i,j,k)*sedgez(i,j,k,comp)
                   eta(k,comp) = 0.d0
                end do
             end do
             eta(k,comp) = eta(k,comp) * fac
          end do
       end do
    end if
    
    do comp = nstart, nstop
       
       do k = lo(3), hi(3)
          
          delta_base = base_new(k,comp) - base_old(k,comp)
          
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                
                divsu = (umac(i+1,j,k) * sedgex(i+1,j,k,comp) &
                     -umac(i  ,j,k) * sedgex(i  ,j,k,comp) ) / dx(1) + &
                     (vmac(i,j+1,k) * sedgey(i,j+1,k,comp) &
                     -vmac(i,j  ,k) * sedgey(i,j  ,k,comp) ) / dx(2) + &
                     ((wmac(i,j,k+1)+w0(k+1)) * sedgez(i,j,k+1,comp) &
                     -(wmac(i,j,k  )+w0(k  )) * sedgez(i,j,k  ,comp) ) / dx(3)
                
                if (which_step .eq. 1) then
                   divbaseu = base_old(k,comp) * &
                        ( (umac(i+1,j,k) - umac(i,j,k) ) / dx(1) &
                        +(vmac(i,j+1,k) - vmac(i,j,k) ) / dx(2) ) &
                        + ( base_old_edge(k+1,comp) * wmac(i,j,k+1) &
                        -base_old_edge(k  ,comp) * wmac(i,j,k  ) ) / dx(3)
                else
                   divbaseu = HALF * (base_old(k,comp) + base_new(k,comp)) * &
                        ( (umac(i+1,j,k) - umac(i,j,k) ) / dx(1) &
                        +(vmac(i,j+1,k) - vmac(i,j,k) ) / dx(2) ) &
                        + ( HALF * (base_old_edge(k+1,comp) + base_new_edge(k+1,comp)) &
                        * wmac(i,j,k+1) &
                        - HALF * (base_old_edge(k  ,comp) + base_new_edge(k  ,comp)) &
                        * wmac(i,j,k  ) ) / dx(3)
                end if
                
                snew(i,j,k,comp) = sold(i,j,k,comp) + delta_base &
                     - dt * (divsu + divbaseu) + dt * force(i,j,k,comp)
                
                if (which_step .eq. 2) then
                   snew(i,j,k,comp) = snew(i,j,k,comp) + dt / dx(3) * (eta(k+1,comp)-eta(k,comp))
                end if
                
             enddo
          enddo
       enddo
    end do
    
    ! Define the update to rho as the sum of the updates to (rho X)_i
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       
       allocate(smin(nstart:nstop),smax(nstart:nstop))
       smin(:) =  1.e20
       smax(:) = -1.e20
       
       snew(:,:,:,rho_comp) = sold(:,:,:,rho_comp)
       
       do comp = nstart, nstop
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) + (snew(i,j,k,comp)-sold(i,j,k,comp))
                   smin(comp) = min(smin(comp),snew(i,j,k,comp))
                   smax(comp) = max(smax(comp),snew(i,j,k,comp))
                enddo
             enddo
          enddo
       enddo
    end if
    
    ! Do not allow the species to leave here negative.
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       do comp = nstart, nstop
          if (smin(comp) .lt. ZERO) then
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      if (snew(i,j,k,comp) .lt. ZERO) then
                         delta = -snew(i,j,k,comp)
                         sum = ZERO
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) &
                                 sum = sum + snew(i,j,k,comp2)
                         enddo
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                               frac = snew(i,j,k,comp2) / sum
                               snew(i,j,k,comp2) = snew(i,j,k,comp2) - frac * delta
                            end if
                         enddo
                         snew(i,j,k,comp) = ZERO
                      end if
                   enddo
                enddo
             enddo
          end if
       enddo
       deallocate(smin,smax)
    end if
    
  end subroutine update_scal_3d_cart

  subroutine update_scal_3d_sphr(which_step,nstart,nstop,sold,snew,umac,vmac,wmac, &
                                 w0_cart,sedgex,sedgey,sedgez,force,base_old,base_new, &
                                 base_old_cart,base_new_cart,lo,hi,domlo,domhi,ng,dx,dt, &
                                 evolve_base_state)
    
    implicit none
    
    integer           , intent(in   ) :: which_step, nstart, nstop
    integer           , intent(in   ) :: lo(:), hi(:), domlo(:), domhi(:), ng
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(inout) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(inout) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(inout) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   base_old(0:,:)
    real (kind = dp_t), intent(in   ) ::   base_new(0:,:)
    real (kind = dp_t), intent(in   ) ::   base_old_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) ::   base_new_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) :: w0_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: evolve_base_state
    
    integer :: i, j, k, comp, comp2
    real (kind = dp_t) :: divsu,divbaseu,mult
    real (kind = dp_t) :: delta,frac,sum
    real (kind = dp_t) :: bc_lox,bc_loy,bc_loz
    real (kind = dp_t) :: bc_hix,bc_hiy,bc_hiz
    real (kind = dp_t), allocatable :: smin(:),smax(:)
    
    ! is spherical
    
    do comp = nstart, nstop
       
       ! Note the umac here does NOT have w0 in it
       if (which_step .eq. 1) then
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   
                   bc_lox = (base_old_cart(i,j,k,comp)+base_old_cart(i-1,j,k,comp)) * HALF
                   bc_loy = (base_old_cart(i,j,k,comp)+base_old_cart(i,j-1,k,comp)) * HALF
                   bc_loz = (base_old_cart(i,j,k,comp)+base_old_cart(i,j,k-1,comp)) * HALF
                   bc_hix = (base_old_cart(i,j,k,comp)+base_old_cart(i+1,j,k,comp)) * HALF
                   bc_hiy = (base_old_cart(i,j,k,comp)+base_old_cart(i,j+1,k,comp)) * HALF
                   bc_hiz = (base_old_cart(i,j,k,comp)+base_old_cart(i,j,k+1,comp)) * HALF
                   
                   if (i.eq.domlo(1)) bc_lox = base_old_cart(i,j,k,comp)
                   if (j.eq.domlo(2)) bc_loy = base_old_cart(i,j,k,comp)
                   if (k.eq.domlo(3)) bc_loz = base_old_cart(i,j,k,comp)
                   if (i.eq.domhi(1)) bc_hix = base_old_cart(i,j,k,comp)
                   if (j.eq.domhi(2)) bc_hiy = base_old_cart(i,j,k,comp)
                   if (k.eq.domhi(3)) bc_hiz = base_old_cart(i,j,k,comp)
                   
                   divbaseu =  &
                        ( bc_hix * umac(i+1,j,k) - bc_lox * umac(i,j,k) ) / dx(1) &
                        +( bc_hiy * vmac(i,j+1,k) - bc_loy * vmac(i,j,k) ) / dx(2) &
                        +( bc_hiz * wmac(i,j,k+1) - bc_loz * wmac(i,j,k) ) / dx(3)
                   
                   snew(i,j,k,comp) = sold(i,j,k,comp) + &
                        (base_new_cart(i,j,k,comp)-base_old_cart(i,j,k,comp)) - &
                        dt * divbaseu + dt * force(i,j,k,comp)
                   
                enddo
             enddo
          enddo
          
       else if (which_step .eq. 2) then
          
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   
                   bc_lox = (base_old_cart(i,j,k,comp)+base_old_cart(i-1,j,k,comp) &
                        +base_new_cart(i,j,k,comp)+base_new_cart(i-1,j,k,comp) ) * FOURTH
                   bc_loy = (base_old_cart(i,j,k,comp)+base_old_cart(i,j-1,k,comp) &
                        +base_new_cart(i,j,k,comp)+base_new_cart(i,j-1,k,comp) ) * FOURTH
                   bc_loz = (base_old_cart(i,j,k,comp)+base_old_cart(i,j,k-1,comp) &
                        +base_new_cart(i,j,k,comp)+base_new_cart(i,j,k-1,comp) ) * FOURTH
                   bc_hix = (base_old_cart(i,j,k,comp)+base_old_cart(i+1,j,k,comp) &
                        +base_new_cart(i,j,k,comp)+base_new_cart(i+1,j,k,comp) ) * FOURTH
                   bc_hiy = (base_old_cart(i,j,k,comp)+base_old_cart(i,j+1,k,comp) &
                        +base_new_cart(i,j,k,comp)+base_new_cart(i,j+1,k,comp) ) * FOURTH
                   bc_hiz = (base_old_cart(i,j,k,comp)+base_old_cart(i,j,k+1,comp) &
                        +base_new_cart(i,j,k,comp)+base_new_cart(i,j,k+1,comp) ) * FOURTH
                   
                   if (i.eq.domlo(1)) bc_lox = &
                        HALF * (base_old_cart(i,j,k,comp)+base_new_cart(i,j,k,comp))
                   if (j.eq.domlo(2)) bc_loy = &
                        HALF * (base_old_cart(i,j,k,comp)+base_new_cart(i,j,k,comp))
                   if (k.eq.domlo(3)) bc_loz = &
                        HALF * (base_old_cart(i,j,k,comp)+base_new_cart(i,j,k,comp))
                   if (i.eq.domhi(1)) bc_hix = &
                        HALF * (base_old_cart(i,j,k,comp)+base_new_cart(i,j,k,comp))
                   if (j.eq.domhi(2)) bc_hiy = &
                        HALF * (base_old_cart(i,j,k,comp)+base_new_cart(i,j,k,comp))
                   if (k.eq.domhi(3)) bc_hiz = &
                        HALF * (base_old_cart(i,j,k,comp)+base_new_cart(i,j,k,comp))
                   
                   divbaseu =  &
                        ( bc_hix * umac(i+1,j,k) - bc_lox * umac(i,j,k) ) / dx(1) &
                        +( bc_hiy * vmac(i,j+1,k) - bc_loy * vmac(i,j,k) ) / dx(2) &
                        +( bc_hiz * wmac(i,j,k+1) - bc_loz * wmac(i,j,k) ) / dx(3)
                   
                   snew(i,j,k,comp) = sold(i,j,k,comp) &
                        + (base_new_cart(i,j,k,comp)-base_old_cart(i,j,k,comp)) &
                        - dt * divbaseu + dt * force(i,j,k,comp)
                   
                enddo
             enddo
          enddo
          
       end if
       
    end do
    
    mult = ONE
    call addw0_3d_sphr(umac,vmac,wmac,w0_cart,lo,hi,mult)
    
    do comp = nstart, nstop
       
       ! Note the umac here DOES have w0 in it
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                
                divsu = (umac(i+1,j,k) * sedgex(i+1,j,k,comp) &
                     -umac(i  ,j,k) * sedgex(i  ,j,k,comp) ) / dx(1) + &
                     (vmac(i,j+1,k) * sedgey(i,j+1,k,comp) &
                     -vmac(i,j  ,k) * sedgey(i,j  ,k,comp) ) / dx(2) + &
                     (wmac(i,j,k+1) * sedgez(i,j,k+1,comp) &
                     -wmac(i,j,k  ) * sedgez(i,j,k  ,comp) ) / dx(3)
                
                snew(i,j,k,comp) = snew(i,j,k,comp) - dt * divsu
                
             enddo
          enddo
       enddo
    enddo
    
    mult = -ONE
    call addw0_3d_sphr(umac,vmac,wmac,w0_cart,lo,hi,mult)
    
    ! Define the update to rho as the sum of the updates to (rho X)_i
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       
       allocate(smin(nstart:nstop),smax(nstart:nstop))
       smin(:) =  1.e20
       smax(:) = -1.e20
       
       snew(:,:,:,rho_comp) = sold(:,:,:,rho_comp)
       
       do comp = nstart, nstop
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) + (snew(i,j,k,comp)-sold(i,j,k,comp))
                   smin(comp) = min(smin(comp),snew(i,j,k,comp))
                   smax(comp) = max(smax(comp),snew(i,j,k,comp))
                enddo
             enddo
          enddo
       enddo
    end if
    
    ! Do not allow the species to leave here negative.
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       do comp = nstart, nstop
          if (smin(comp) .lt. ZERO) then
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      if (snew(i,j,k,comp) .lt. ZERO) then
                         delta = -snew(i,j,k,comp)
                         sum = ZERO
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) &
                                 sum = sum + snew(i,j,k,comp2)
                         enddo
                         do comp2 = nstart, nstop
                            if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                               frac = snew(i,j,k,comp2) / sum
                               snew(i,j,k,comp2) = snew(i,j,k,comp2) - frac * delta
                            end if
                         enddo
                         snew(i,j,k,comp) = ZERO
                      end if
                   enddo
                enddo
             enddo
          end if
       enddo
       deallocate(smin,smax)
    end if
    
  end subroutine update_scal_3d_sphr
  
end module update_scal_module
