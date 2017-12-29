module rhoh_vs_t_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use ml_restrict_fill_module

  implicit none

  private

  public :: makeHfromRhoT_edge, makeTfromRhoH, makeTfromRhoP, makePfromRhoH
  
contains

  !============================================================================
  ! makeHfromRhoT_edge
  !============================================================================
  subroutine makeHfromRhoT_edge(u,sedge, &
                                rho0_old,rhoh0_old,t0_old, &
                                rho0_edge_old,rhoh0_edge_old,t0_edge_old, &
                                rho0_new,rhoh0_new,t0_new, &
                                rho0_edge_new,rhoh0_edge_new,t0_edge_new, &
                                the_bc_level,dx,mla)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical, nr_fine
    use variables
    use network
    use fill_3d_module
    use pred_parameters
    
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: t0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_edge_old(:,0:)
    real(kind=dp_t), intent(in   ) :: t0_edge_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: t0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_edge_new(:,0:)
    real(kind=dp_t), intent(in   ) :: t0_edge_new(:,0:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(in   ) :: mla
    
    ! local
    integer :: i,r,n,ng_u,ng_se,ng_r0,ng_rh0,ng_t0,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer ::   rp(:,:,:,:)
    real(kind=dp_t), pointer ::  rhp(:,:,:,:)
    real(kind=dp_t), pointer ::   tp(:,:,:,:)

    real(kind=dp_t) ::  rho0_halftime(1,0:nr_fine-1)
    real(kind=dp_t) :: rhoh0_halftime(1,0:nr_fine-1)
    real(kind=dp_t) ::    t0_halftime(1,0:nr_fine-1)
    type(multifab)  ::  rho0_cart(mla%nlevel)
    type(multifab)  :: rhoh0_cart(mla%nlevel)
    type(multifab)  ::    t0_cart(mla%nlevel)
    type(layout)    :: la

    type(bl_prof_timer), save :: bpt

    call build(bpt, "makeHfromRhoT_edge")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_u  = nghost(u(1))
    ng_se = nghost(sedge(1,1))

    if (spherical .eq. 1) then

       do n=1,nlevs
          la=get_layout(u(n))
          call multifab_build( rho0_cart(n),la,1,2)
          call multifab_build(rhoh0_cart(n),la,1,2)
          call multifab_build(   t0_cart(n),la,1,2)
       end do

       ng_r0  = nghost(rho0_cart(1))
       ng_rh0 = nghost(rhoh0_cart(1))
       ng_t0  = nghost(t0_cart(1))

       !$OMP PARALLEL DO PRIVATE(r)
       do r=0,nr_fine-1
          rho0_halftime(1,r)  = HALF * (rho0_old(1,r)  + rho0_new(1,r)  )
          rhoh0_halftime(1,r) = HALF * (rhoh0_old(1,r) + rhoh0_new(1,r) )
          t0_halftime(1,r)    = HALF * (t0_old(1,r)    + t0_new(1,r) )
       end do
       !$OMP END PARALLEL DO

       call put_1d_array_on_cart(rho0_halftime,rho0_cart,dm+rho_comp,.false.,.false.,dx, &
                                 the_bc_level,mla)
       call put_1d_array_on_cart(rhoh0_halftime,rhoh0_cart,dm+rhoh_comp,.false.,.false.,dx, &
                                 the_bc_level,mla)
       call put_1d_array_on_cart(t0_halftime,t0_cart,dm+temp_comp,.false.,.false.,dx, &
                                 the_bc_level,mla)
   endif

   do n=1,nlevs

       do i=1,nfabs(u(n))
          sepx => dataptr(sedge(n,1), i)
          lo = lwb(get_box(u(n),i))
          hi = upb(get_box(u(n),i))
          select case (dm)
          case (1)
             call makeHfromRhoT_edge_1d(sepx(:,1,1,:), ng_se, &
                                        rho0_edge_old(n,:), rhoh0_edge_old(n,:), &
                                        t0_edge_old(n,:), &
                                        rho0_edge_new(n,:), &
                                        rhoh0_edge_new(n,:), t0_edge_new(n,:), lo, hi)

          case (2)
             sepy => dataptr(sedge(n,2), i)
             call makeHfromRhoT_edge_2d(sepx(:,:,1,:), sepy(:,:,1,:), ng_se, &
                                        rho0_old(n,:), rhoh0_old(n,:), t0_old(n,:), &
                                        rho0_edge_old(n,:), rhoh0_edge_old(n,:), &
                                        t0_edge_old(n,:), rho0_new(n,:), rhoh0_new(n,:), &
                                        t0_new(n,:), rho0_edge_new(n,:), &
                                        rhoh0_edge_new(n,:), t0_edge_new(n,:), lo, hi)
          case (3)
             sepy => dataptr(sedge(n,2), i)
             sepz => dataptr(sedge(n,3), i)
             if (spherical .eq. 1) then
               rp   => dataptr( rho0_cart(n), i)
               rhp  => dataptr(rhoh0_cart(n), i)
               tp   => dataptr(   t0_cart(n), i)
               call makeHfromRhoT_edge_3d_sphr(sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                               ng_se, rp(:,:,:,1), ng_r0, rhp(:,:,:,1), &
                                               ng_rh0, tp(:,:,:,1), ng_t0, lo, hi)
             else
               call makeHfromRhoT_edge_3d_cart(sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                               ng_se, rho0_old(n,:), rhoh0_old(n,:), &
                                               t0_old(n,:), rho0_edge_old(n,:), &
                                               rhoh0_edge_old(n,:), t0_edge_old(n,:), &
                                               rho0_new(n,:), rhoh0_new(n,:), &
                                               t0_new(n,:), rho0_edge_new(n,:), &
                                               rhoh0_edge_new(n,:), t0_edge_new(n,:), &
                                               lo, hi)
             end if
          end select
       end do

    end do

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy( rho0_cart(n))
          call destroy(rhoh0_cart(n))
          call destroy(   t0_cart(n))
       end do
    end if

    call destroy(bpt)
    
  end subroutine makeHfromRhoT_edge

  !----------------------------------------------------------------------------
  ! makeHfromRhoT_edge_1d
  !----------------------------------------------------------------------------
  subroutine makeHfromRhoT_edge_1d(sx,ng_se, &
                                   rho0_edge_old,rhoh0_edge_old,t0_edge_old, &
                                   rho0_edge_new,rhoh0_edge_new,t0_edge_new, &
                                   lo,hi)

    use bl_constants_module
    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use eos_module,    only: eos, eos_input_rt
    use eos_type_module
    use network,       only: nspec
    use probin_module, only: enthalpy_pred_type, species_pred_type, small_temp
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:),ng_se
    real(kind=dp_t), intent(inout) :: sx(lo(1)-ng_se:,:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(0:),rhoh0_edge_old(0:),t0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(0:),rhoh0_edge_new(0:),t0_edge_new(0:)
 
    integer :: i
    real(kind=dp_t) :: t0_edge
    
    integer :: pt_index(MAX_SPACEDIM)
    type(eos_t) :: eos_state


    do i = lo(1), hi(1)+1

       ! get edge-centered temperature
       if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
          t0_edge = HALF*(t0_edge_old(i)+t0_edge_new(i))
          eos_state%T = max(sx(i,temp_comp)+t0_edge,small_temp)
       else
          eos_state%T = max(sx(i,temp_comp),small_temp)
       end if

       ! get edge-centered density and species
       if (species_pred_type .eq. predict_rhoprime_and_X) then
          
          ! interface states are rho' and X
          eos_state%rho = sx(i,rho_comp) + &
               HALF * (rho0_edge_old(i) + rho0_edge_new(i))
          
          eos_state%xn(:) = sx(i,spec_comp:spec_comp+nspec-1)
          
       else if (species_pred_type .eq. predict_rhoX) then
          
          ! interface states are rho and (rho X)
          eos_state%rho = sx(i,rho_comp)
          
          eos_state%xn(:) = sx(i,spec_comp:spec_comp+nspec-1)/eos_state%rho
          
       else if (species_pred_type .eq. predict_rho_and_X) then
          
          ! interface states are rho and X
          eos_state%rho = sx(i,rho_comp)

          eos_state%xn(:) = sx(i,spec_comp:spec_comp+nspec-1)

       endif

       pt_index(:) = (/i, -1, -1/)

       call eos(eos_input_rt, eos_state, pt_index)

       if (enthalpy_pred_type .eq. predict_T_then_h .or. &
           enthalpy_pred_type .eq. predict_Tprime_then_h) then
          sx(i,rhoh_comp) = eos_state%h

       else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
          sx(i,rhoh_comp) = eos_state%rho*eos_state%h &
               - HALF*(rhoh0_edge_old(i) + rhoh0_edge_new(i))
       end if

    enddo

  end subroutine makeHfromRhoT_edge_1d

  !----------------------------------------------------------------------------
  ! makeHfromRhoT_edge_2d
  !----------------------------------------------------------------------------
  subroutine makeHfromRhoT_edge_2d(sx,sy,ng_se, &
                                   rho0_old,rhoh0_old,t0_old, &
                                   rho0_edge_old,rhoh0_edge_old,t0_edge_old, &
                                   rho0_new,rhoh0_new,t0_new, &
                                   rho0_edge_new,rhoh0_edge_new,t0_edge_new, &
                                   lo,hi)
    
    use bl_constants_module
    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use eos_module,    only: eos, eos_input_rt
    use eos_type_module
    use network,       only: nspec
    use probin_module, only: enthalpy_pred_type, species_pred_type, small_temp
    use pred_parameters

    integer        , intent(in   ) :: lo(:),hi(:),ng_se
    real(kind=dp_t), intent(inout) :: sx(lo(1)-ng_se:,lo(2)-ng_se:,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1)-ng_se:,lo(2)-ng_se:,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:),rhoh0_old(0:),t0_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(0:),rhoh0_edge_old(0:),t0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:),rhoh0_new(0:),t0_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(0:),rhoh0_edge_new(0:),t0_edge_new(0:)

    integer :: i,j
    real(kind=dp_t) :: t0_edge

    integer :: pt_index(MAX_SPACEDIM)
    type(eos_t) :: eos_state

    ! x-edge
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1

          ! get edge-centered temperature
          if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
             t0_edge = HALF*(t0_old(j)+t0_new(j))
             eos_state%T = max(sx(i,j,temp_comp)+t0_edge,small_temp)
          else
             eos_state%T = max(sx(i,j,temp_comp),small_temp)
          end if
          
          ! get edge-centered density and species
          if (species_pred_type .eq. predict_rhoprime_and_X) then
             
             ! interface states are rho' and X
             eos_state%rho  = sx(i,j,rho_comp) + &
                  HALF * (rho0_old(j) + rho0_new(j))
             
             eos_state%xn(:) = sx(i,j,spec_comp:spec_comp+nspec-1)
             
          else if (species_pred_type .eq. predict_rhoX) then
             
             ! interface states are rho and (rho X)
             eos_state%rho = sx(i,j,rho_comp)
              
             eos_state%xn(:) = sx(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho
              
          else if (species_pred_type .eq. predict_rho_and_X) then
              
             ! interface states are rho and X
             eos_state%rho = sx(i,j,rho_comp)
             
             eos_state%xn(:) = sx(i,j,spec_comp:spec_comp+nspec-1)
              
          endif
           
          pt_index(:) = (/i, j, -1/)
          
          call eos(eos_input_rt, eos_state, pt_index)
           
          if (enthalpy_pred_type .eq. predict_T_then_h .or. &
               enthalpy_pred_type .eq. predict_Tprime_then_h) then
             sx(i,j,rhoh_comp) = eos_state%h
              
          else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
             sx(i,j,rhoh_comp) = eos_state%rho*eos_state%h - HALF*(rhoh0_old(j)+rhoh0_new(j))
          end if
          
       enddo
    enddo
     
    ! y-edge
    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
           
          ! get edge-centered temperature
          if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
             t0_edge = HALF*(t0_edge_old(j)+t0_edge_new(j))
             eos_state%T = max(sy(i,j,temp_comp)+t0_edge,small_temp)
          else
             eos_state%T = max(sy(i,j,temp_comp),small_temp)
          end if
          
          ! get edge-centered density and species
          if (species_pred_type .eq. predict_rhoprime_and_X) then
             
             ! interface states are rho' and X
             eos_state%rho  = sy(i,j,rho_comp) + &
                  HALF * (rho0_edge_old(j) + rho0_edge_new(j))
             
             eos_state%xn(:) = sy(i,j,spec_comp:spec_comp+nspec-1)

          else if (species_pred_type .eq. predict_rhoX) then
              
             ! interface states are rho and (rho X)
             eos_state%rho = sy(i,j,rho_comp)
              
             eos_state%xn(:) = sy(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho
              
          else if (species_pred_type .eq. predict_rho_and_X) then
              
             ! interface states are rho and X
             eos_state%rho = sy(i,j,rho_comp)
              
             eos_state%xn(:) = sy(i,j,spec_comp:spec_comp+nspec-1)
              
          endif


          pt_index(:) = (/i, j, -1/)
          
          call eos(eos_input_rt, eos_state, pt_index)
           
          if (enthalpy_pred_type .eq. predict_T_then_h .or. &
               enthalpy_pred_type .eq. predict_Tprime_then_h) then
             sy(i,j,rhoh_comp) = eos_state%h 
              
          else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
             sy(i,j,rhoh_comp) = eos_state%rho*eos_state%h &
                  - HALF*(rhoh0_edge_old(j) + rhoh0_edge_new(j))
          end if
           
       enddo
    enddo
    
  end subroutine makeHfromRhoT_edge_2d
  
  !----------------------------------------------------------------------------
  ! makeHfromRhoT_edge_3d_cart
  !----------------------------------------------------------------------------
  subroutine makeHfromRhoT_edge_3d_cart(sx,sy,sz,ng_se, &
                                        rho0_old,rhoh0_old,t0_old, &
                                        rho0_edge_old,rhoh0_edge_old,t0_edge_old, &
                                        rho0_new,rhoh0_new,t0_new, &
                                        rho0_edge_new,rhoh0_edge_new,t0_edge_new, &
                                        lo,hi)

    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use eos_module,    only: eos, eos_input_rt
    use eos_type_module
    use network,       only: nspec
    use probin_module, only: species_pred_type, enthalpy_pred_type, small_temp
    use pred_parameters
    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_se
    real(kind=dp_t), intent(inout) :: sx(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) :: sz(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(in   ) :: rho0_old(0:),rhoh0_old(0:),t0_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_old(0:),rhoh0_edge_old(0:),t0_edge_old(0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(0:),rhoh0_new(0:),t0_new(0:)
    real(kind=dp_t), intent(in   ) :: rho0_edge_new(0:),rhoh0_edge_new(0:),t0_edge_new(0:)
    
    integer         :: i,j,k
    real(kind=dp_t) :: t0_edge

    integer :: pt_index(MAX_SPACEDIM)
    type(eos_t) :: eos_state

    ! x-edge
    !$OMP PARALLEL DO PRIVATE(i,j,k,t0_edge,eos_state,pt_index)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1

             ! get edge-centered temperature
             if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                t0_edge = HALF*(t0_old(k)+t0_new(k))
                eos_state%T = max(sx(i,j,k,temp_comp)+t0_edge,small_temp)
             else
                eos_state%T = max(sx(i,j,k,temp_comp),small_temp)
             end if

             ! get edge-centered density and species
             if (species_pred_type .eq. predict_rhoprime_and_X) then

                ! interface states are rho' and X
                eos_state%rho = sx(i,j,k,rho_comp) + &
                     HALF * (rho0_old(k) + rho0_new(k))

                eos_state%xn(:) = sx(i,j,k,spec_comp:spec_comp+nspec-1)

             else if (species_pred_type .eq. predict_rhoX) then
                
                ! interface states are rho and (rho X)
                eos_state%rho = sx(i,j,k,rho_comp)

                eos_state%xn(:) = sx(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             else if (species_pred_type .eq. predict_rho_and_X) then

                ! interface states are rho and X
                eos_state%rho = sx(i,j,k,rho_comp)
                
                eos_state%xn(:) = sx(i,j,k,spec_comp:spec_comp+nspec-1)

             endif

             pt_index(:) = (/i, j, k/)
             
             call eos(eos_input_rt, eos_state, pt_index)
             
             if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                 enthalpy_pred_type .eq. predict_Tprime_then_h) then
                sx(i,j,k,rhoh_comp) = eos_state%h

             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                sx(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - HALF*(rhoh0_old(k)+rhoh0_new(k))
             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! y-edge
    !$OMP PARALLEL DO PRIVATE(i,j,k,t0_edge,eos_state,pt_index)    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             
             ! get edge-centered temperature
             if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                t0_edge = HALF*(t0_old(k)+t0_new(k))
                eos_state%T = max(sx(i,j,k,temp_comp)+t0_edge,small_temp)
             else
                eos_state%T = max(sy(i,j,k,temp_comp),small_temp)
             end if
       
             ! get edge-centered density and species
             if (species_pred_type .eq. predict_rhoprime_and_X) then

                ! interface states are rho' and X
                eos_state%rho  = sy(i,j,k,rho_comp) + &
                     HALF * (rho0_old(k) + rho0_new(k))

                eos_state%xn(:) = sy(i,j,k,spec_comp:spec_comp+nspec-1)

             else if (species_pred_type .eq. predict_rhoX) then

                ! interface states are rho and (rho X)
                eos_state%rho = sy(i,j,k,rho_comp) 

                eos_state%xn(:) = sy(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             else if (species_pred_type .eq. predict_rho_and_X) then

                ! interface states are rho and X
                eos_state%rho = sy(i,j,k,rho_comp) 

                eos_state%xn(:) = sy(i,j,k,spec_comp:spec_comp+nspec-1)

             endif                

             pt_index(:) = (/i, j, k/)
             
             call eos(eos_input_rt, eos_state, pt_index)
             
             if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                 enthalpy_pred_type .eq. predict_Tprime_then_h) then
                sy(i,j,k,rhoh_comp) = eos_state%h

             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                sy(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - HALF*(rhoh0_old(k)+rhoh0_new(k))
             end if

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! z-edge
    !$OMP PARALLEL DO PRIVATE(i,j,k,t0_edge,eos_state,pt_index) 
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! get edge-centered temperature
             if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                t0_edge = HALF*(t0_edge_old(k)+t0_edge_new(k))
                eos_state%T = max(sx(i,j,k,temp_comp)+t0_edge,small_temp)
             else
                eos_state%T = max(sz(i,j,k,temp_comp),small_temp)
             end if

             ! get edge-centered density and species
             if (species_pred_type .eq. predict_rhoprime_and_X) then

                ! interface states are rho' and X
                eos_state%rho = sz(i,j,k,rho_comp) + &
                     HALF * (rho0_edge_old(k) + rho0_edge_new(k))

                eos_state%xn(:) = sz(i,j,k,spec_comp:spec_comp+nspec-1)

             else if (species_pred_type .eq. predict_rhoX) then

                ! interface states are rho and (rho X)
                eos_state%rho = sz(i,j,k,rho_comp)  

                eos_state%xn(:) = sz(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             else if (species_pred_type .eq. predict_rho_and_X) then

                ! interface states are rho and X
                eos_state%rho = sz(i,j,k,rho_comp)                 

                eos_state%xn(:) = sz(i,j,k,spec_comp:spec_comp+nspec-1)

             endif

             pt_index(:) = (/i, j, k/)

             call eos(eos_input_rt, eos_state, pt_index)
             
             if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                 enthalpy_pred_type .eq. predict_Tprime_then_h) then
                sz(i,j,k,rhoh_comp) = eos_state%h

             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                sz(i,j,k,rhoh_comp) =  eos_state%rho*eos_state%h &
                     - HALF*(rhoh0_edge_old(k)+rhoh0_edge_new(k))
             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine makeHfromRhoT_edge_3d_cart

  !----------------------------------------------------------------------------
  ! makeHfromRhoT_edge_3d_sphr
  !----------------------------------------------------------------------------
  subroutine makeHfromRhoT_edge_3d_sphr(sx,sy,sz,ng_se,rho0_cart,ng_r0,rhoh0_cart,ng_rh0, &
                                        t0_cart,ng_t0,lo,hi)

    use variables,     only: rho_comp, temp_comp, spec_comp, rhoh_comp
    use eos_module,    only: eos, eos_input_rt
    use eos_type_module
    use network,       only: nspec
    use probin_module, only: specieS_pred_type, enthalpy_pred_type, small_temp
    use pred_parameters
    use bl_constants_module

    integer        , intent(in   ) :: ng_se,ng_r0,ng_rh0,ng_t0
    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) :: sz(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::  rho0_cart(lo(1)-ng_r0:,lo(2)-ng_r0:,lo(3)-ng_r0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_cart(lo(1)-ng_rh0:,lo(2)-ng_rh0:,lo(3)-ng_rh0:)
    real(kind=dp_t), intent(in   ) ::    t0_cart(lo(1)-ng_t0 :,lo(2)-ng_t0 :,lo(3)-ng_t0 :)
    
    ! Local variables
    integer :: i, j, k
    real(kind=dp_t) rho0_edge, rhoh0_edge, t0_edge

    integer :: pt_index(MAX_SPACEDIM)
    type(eos_t) :: eos_state
    
    ! x-edge
    !$OMP PARALLEL DO PRIVATE(i,j,k,t0_edge,rho0_edge,rhoh0_edge,eos_state,pt_index)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             
             ! get edge-centered temperature
             if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                t0_edge = HALF* (t0_cart(i-1,j,k) + t0_cart(i,j,k))
                eos_state%T = max(sx(i,j,k,temp_comp)+t0_edge,small_temp)
             else
                eos_state%T = max(sx(i,j,k,temp_comp),small_temp)
             end if

             ! get edge-centered density and species
             if (species_pred_type .eq. predict_rhoprime_and_X) then
                
                ! interface states are rho' and X
                rho0_edge = HALF * (rho0_cart(i-1,j,k) + rho0_cart(i,j,k))
                eos_state%rho = sx(i,j,k,rho_comp) + rho0_edge

                eos_state%xn(:) = sx(i,j,k,spec_comp:spec_comp+nspec-1)

             else if (species_pred_type .eq. predict_rhoX) then

                ! interface states are rho and (rho X)
                eos_state%rho = sx(i,j,k,rho_comp) 

                eos_state%xn(:) = sx(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             else if (species_pred_type .eq. predict_rho_and_X) then

                ! interface states are rho and X
                eos_state%rho = sx(i,j,k,rho_comp) 

                eos_state%xn(:) = sx(i,j,k,spec_comp:spec_comp+nspec-1)

             endif

             pt_index(:) = (/i, j, k/)

             call eos(eos_input_rt, eos_state, pt_index)
             
             if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                 enthalpy_pred_type .eq. predict_Tprime_then_h) then
                sx(i,j,k,rhoh_comp) = eos_state%h

             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                rhoh0_edge = HALF * (rhoh0_cart(i-1,j,k) + rhoh0_cart(i,j,k))
                sx(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - rhoh0_edge
             end if

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO


    ! y-edge
    !$OMP PARALLEL DO PRIVATE(i,j,k,t0_edge,rho0_edge,rhoh0_edge,eos_state,pt_index)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             
             ! get edge-centered temperature
             if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                t0_edge = HALF * (t0_cart(i,j-1,k) + t0_cart(i,j,k))
                eos_state%T = max(sy(i,j,k,temp_comp)+t0_edge,small_temp)
             else
                eos_state%T = max(sy(i,j,k,temp_comp),small_temp)
             end if

             ! get edge-centered density and species
             if (species_pred_type .eq. predict_rhoprime_and_X) then
                
                ! interface states are rho' and X
                rho0_edge = HALF * (rho0_cart(i,j-1,k) + rho0_cart(i,j,k))
                eos_state%rho = sy(i,j,k,rho_comp) + rho0_edge

                eos_state%xn(:) = sy(i,j,k,spec_comp:spec_comp+nspec-1) 

             else if (species_pred_type .eq. predict_rhoX) then

                ! interface states are rho and (rho X)
                eos_state%rho = sy(i,j,k,rho_comp)

                eos_state%xn(:) = sy(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             else if (species_pred_type .eq. predict_rho_and_X) then

                ! interface states are rho and X
                eos_state%rho = sy(i,j,k,rho_comp)

                eos_state%xn(:) = sy(i,j,k,spec_comp:spec_comp+nspec-1)

             endif

             pt_index(:) = (/i, j, k/)

             call eos(eos_input_rt, eos_state, pt_index)
             
             if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                 enthalpy_pred_type .eq. predict_Tprime_then_h) then
                sy(i,j,k,rhoh_comp) = eos_state%h

             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                rhoh0_edge = HALF * (rhoh0_cart(i,j-1,k) + rhoh0_cart(i,j,k))
                sy(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - rhoh0_edge
             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! z-edge
    !$OMP PARALLEL DO PRIVATE(i,j,k,t0_edge,rho0_edge,rhoh0_edge,eos_state,pt_index)
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
                 
             ! get edge-centered temperature
             if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                t0_edge = HALF * (t0_cart(i,j,k-1) + t0_cart(i,j,k))
                eos_state%T = max(sz(i,j,k,temp_comp)+t0_edge,small_temp)
             else
                eos_state%T = max(sz(i,j,k,temp_comp),small_temp)
             end if

             ! get edge-centered density and species
             if (species_pred_type .eq. predict_rhoprime_and_X) then

                ! interface states are rho' and X
                rho0_edge = HALF * (rho0_cart(i,j,k-1) + rho0_cart(i,j,k))
                eos_state%rho = sz(i,j,k,rho_comp) + rho0_edge
             
                eos_state%xn(:) = sz(i,j,k,spec_comp:spec_comp+nspec-1) 

             else if (species_pred_type .eq. predict_rhoX) then

                ! interface states are rho and (rho X)
                eos_state%rho = sz(i,j,k,rho_comp)

                eos_state%xn(:) = sz(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             else if (species_pred_type .eq. predict_rho_and_X) then

                ! interface states are rho and X
                eos_state%rho = sz(i,j,k,rho_comp)

                eos_state%xn(:) = sz(i,j,k,spec_comp:spec_comp+nspec-1)

             endif

             pt_index(:) = (/i, j, k/)

             call eos(eos_input_rt, eos_state, pt_index)
             
             if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                 enthalpy_pred_type .eq. predict_Tprime_then_h) then
                sz(i,j,k,rhoh_comp) = eos_state%h
             else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                rhoh0_edge = HALF * (rhoh0_cart(i,j,k-1) + rhoh0_cart(i,j,k))
                sz(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - rhoh0_edge
             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine makeHfromRhoT_edge_3d_sphr
  

  !============================================================================
  ! makeTfromRhoH
  !============================================================================
  subroutine makeTfromRhoH(state,p0,mla,the_bc_level,dx)

    use variables,             only: temp_comp
    use bl_prof_module
    use geometry, only: spherical

    type(multifab)    , intent(inout) :: state(:)
    real (kind = dp_t), intent(in   ) :: p0(:,0:)
    type(ml_layout)   , intent(in   ) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t)   , intent(in   ) :: dx(:,:)

    ! local
    integer                  :: i,ng,n
    integer                  :: lo(mla%dim),hi(mla%dim),dm,nlevs
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "makeTfromRhoH")

    dm = mla%dim
    nlevs = mla%nlevel

    ng = nghost(state(1))

    do n=1,nlevs

       do i=1,nfabs(state(n))
          sp => dataptr(state(n),i)
          lo = lwb(get_box(state(n),i))
          hi = upb(get_box(state(n),i))
          select case (dm)
          case (1)
             call makeTfromRhoH_1d(sp(:,1,1,:), lo, hi, ng, p0(n,:))
          case (2)
             call makeTfromRhoH_2d(sp(:,:,1,:), lo, hi, ng, p0(n,:))
          case (3)
             if (spherical .eq. 1) then
                call makeTfromRhoH_3d_sphr(sp(:,:,:,:), lo, hi, ng, p0(1,:), dx(n,:))
             else
                call makeTfromRhoH_3d(sp(:,:,:,:), lo, hi, ng, p0(n,:))
             endif
          end select
       end do

    end do

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,state,mla%mba%rr,the_bc_level, &
                              icomp=temp_comp, &
                              bcomp=dm+temp_comp, &
                              nc=1, &
                              ng=state(1)%ng)

    call destroy(bpt)

  end subroutine makeTfromRhoH

  !----------------------------------------------------------------------------
  ! makeTfromRhoH_1d
  !----------------------------------------------------------------------------
  subroutine makeTfromRhoH_1d(state,lo,hi,ng,p0)

    use variables, only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use network, only: nspec
    use eos_module, only: eos_input_re, eos_input_rh, eos
    use eos_type_module
    use probin_module, only: use_eos_e_instead_of_h
    
    integer, intent(in)               :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: state(lo(1)-ng:,:)
    real (kind = dp_t), intent(in   ) :: p0(0:)
    
    
    ! Local variables
    integer :: i
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state

    if (use_eos_e_instead_of_h) then

       do i = lo(1), hi(1)

          ! (rho, (h->e)) --> T, p

          eos_state%rho   = state(i,rho_comp)
          eos_state%T     = state(i,temp_comp)
          eos_state%xn(:) = state(i,spec_comp:spec_comp+nspec-1)/eos_state%rho

          ! e = h - p/rho
          eos_state%e = state(i,rhoh_comp) / state(i,rho_comp) - &
               p0(i) / state(i,rho_comp)

          pt_index(:) = (/i, -1, -1/)

          call eos(eos_input_re, eos_state, pt_index)
          
          state(i,temp_comp) = eos_state%T

       enddo

    else

       do i = lo(1), hi(1)

          ! (rho, h) --> T, p

          eos_state%rho   = state(i,rho_comp)
          eos_state%T     = state(i,temp_comp)
          eos_state%xn(:) = state(i,spec_comp:spec_comp+nspec-1)/eos_state%rho

          eos_state%h = state(i,rhoh_comp) / state(i,rho_comp)

          pt_index(:) = (/i, -1, -1/)

          call eos(eos_input_rh, eos_state, pt_index)
          
          state(i,temp_comp) = eos_state%T

       enddo

    endif

  end subroutine makeTfromRhoH_1d

  !----------------------------------------------------------------------------
  ! makeTfromRhoH_2d
  !----------------------------------------------------------------------------
  subroutine makeTfromRhoH_2d(state,lo,hi,ng,p0)

    use variables, only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use network, only: nspec
    use eos_module, only: eos_input_re, eos_input_rh, eos
    use eos_type_module
    use probin_module, only: use_eos_e_instead_of_h

    integer, intent(in)               :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:)
    
    ! Local variables
    integer :: i, j
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state

    if (use_eos_e_instead_of_h) then
    
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, (h->e)) --> T, p
          
             eos_state%rho   = state(i,j,rho_comp)
             eos_state%T     = state(i,j,temp_comp)
             eos_state%xn(:) = state(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho

             ! e = h - p/rho             
             eos_state%e = state(i,j,rhoh_comp) / state(i,j,rho_comp) - &
                  p0(j) / state(i,j,rho_comp)

             pt_index(:) = (/i, j, -1/)
          
             call eos(eos_input_re, eos_state, pt_index)
          
             state(i,j,temp_comp) = eos_state%T
          
          enddo
       enddo
    
    else

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, h) --> T, p
          
             eos_state%rho   = state(i,j,rho_comp)
             eos_state%T     = state(i,j,temp_comp)
             eos_state%xn(:) = state(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho
             
             eos_state%h = state(i,j,rhoh_comp) / state(i,j,rho_comp)

             pt_index(:) = (/i, j, -1/)
          
             call eos(eos_input_rh, eos_state, pt_index)
          
             state(i,j,temp_comp) = eos_state%T
          
          enddo
       enddo

    endif

  end subroutine makeTfromRhoH_2d

  !----------------------------------------------------------------------------
  ! makeTfromRhoH_3d
  !----------------------------------------------------------------------------
  subroutine makeTfromRhoH_3d(state,lo,hi,ng,p0)

    use variables,     only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use eos_module,    only: eos_input_re, eos_input_rh, eos
    use eos_type_module
    use network,       only: nspec
    use probin_module, only: use_eos_e_instead_of_h

    integer          , intent(in)    :: lo(:), hi(:), ng
    real(kind=dp_t)  , intent(inout) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind = dp_t), intent(in   ) ::  p0(0:)

    ! Local variables
    integer :: i, j, k
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state

    if (use_eos_e_instead_of_h) then

       !$OMP PARALLEL DO PRIVATE(i,j,k, eos_state, pt_index)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
             
                ! (rho, (h->e)) --> T, p
             
                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                ! e = h - p/rho
                eos_state%e = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - &
                     p0(k) / state(i,j,k,rho_comp)
                
                pt_index(:) = (/i, j, k/)
             
                call eos(eos_input_re, eos_state, pt_index)
             
                state(i,j,k,temp_comp) = eos_state%T
             
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    else

       !$OMP PARALLEL DO PRIVATE(i,j,k, eos_state, pt_index)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
             
                ! (rho, h) --> T, p
             
                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                eos_state%h = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
                
                pt_index(:) = (/i, j, k/)
             
                call eos(eos_input_rh, eos_state, pt_index)
             
                state(i,j,k,temp_comp) = eos_state%T
             
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    endif

  end subroutine makeTfromRhoH_3d

  !----------------------------------------------------------------------------
  ! makeTfromRhoH_3d_sphr
  !----------------------------------------------------------------------------
  subroutine makeTfromRhoH_3d_sphr(state,lo,hi,ng,p0,dx)

    use variables,     only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use eos_module,    only: eos_input_re, eos_input_rh, eos
    use eos_type_module
    use network,       only: nspec
    use probin_module, only: use_eos_e_instead_of_h
    use fill_3d_module

    integer          , intent(in)    :: lo(:), hi(:), ng
    real(kind=dp_t)  , intent(inout) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t)  , intent(in   ) ::  p0(0:)
    real(kind=dp_t)  , intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j, k
    real(kind=dp_t), allocatable :: p0_cart(:,:,:,:)
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state

    if (use_eos_e_instead_of_h) then

       allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_3d_sphr(.false.,.false.,p0,p0_cart,lo,hi,dx,0)

       !$OMP PARALLEL DO PRIVATE(i,j,k, eos_state, pt_index)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
             
                ! (rho, (h->e)) --> T, p
             
                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                ! e = h - p/rho
                eos_state%e = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - &
                     p0_cart(i,j,k,1) / state(i,j,k,rho_comp)
                
                pt_index(:) = (/i, j, k/)
             
                call eos(eos_input_re, eos_state, pt_index)
             
                state(i,j,k,temp_comp) = eos_state%T
             
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

       deallocate(p0_cart)

    else

       !$OMP PARALLEL DO PRIVATE(i,j,k, eos_state, pt_index)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
             
                ! (rho, h) --> T, p
             
                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                eos_state%h = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
                
                pt_index(:) = (/i, j, k/)
             
                call eos(eos_input_rh, eos_state, pt_index)
             
                state(i,j,k,temp_comp) = eos_state%T
             
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    endif

  end subroutine makeTfromRhoH_3d_sphr



  !============================================================================
  ! makeTfromRhoP
  !============================================================================
  subroutine makeTfromRhoP(state,p0,mla,the_bc_level,dx,updateRhoH_in)

    use variables,             only: temp_comp, rhoh_comp
    use bl_prof_module
    use geometry, only: spherical

    type(multifab)    , intent(inout) :: state(:)
    real (kind = dp_t), intent(in   ) :: p0(:,0:)
    type(ml_layout)   , intent(in   ) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t)   , intent(in   ) :: dx(:,:)
    integer, intent(in), optional     :: updateRhoH_in

    ! local
    integer                  :: i,ng,n,updateRhoH
    integer                  :: lo(mla%dim),hi(mla%dim),dm,nlevs
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "makeTfromRhoP")

    updateRhoH = 0
    if (present(updateRhoH_in)) updateRhoH = updateRhoH_in

    dm = mla%dim
    nlevs = mla%nlevel

    ng = nghost(state(1))

    do n=1,nlevs

       do i=1,nfabs(state(n))
          sp => dataptr(state(n),i)
          lo = lwb(get_box(state(n),i))
          hi = upb(get_box(state(n),i))
          select case (dm)
          case (1)
             call makeTfromRhoP_1d(sp(:,1,1,:),lo,hi,ng,p0(n,:),updateRhoH)
          case (2)
             call makeTfromRhoP_2d(sp(:,:,1,:),lo,hi,ng,p0(n,:),updateRhoH)
          case (3)
             if (spherical .eq. 1) then
                call makeTfromRhoP_3d_sphr(sp(:,:,:,:),lo,hi,ng,p0(1,:),dx(n,:),updateRhoH)
             else
                call makeTfromRhoP_3d(sp(:,:,:,:),lo,hi,ng,p0(n,:),updateRhoH)
             end if
          end select
       end do

    end do

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,state,mla%mba%rr,the_bc_level, &
                              icomp=temp_comp, &
                              bcomp=dm+temp_comp, &
                              nc=1, &
                              ng=state(1)%ng)

    if (updateRhoH .eq. 1) then
       call ml_restrict_and_fill(nlevs,state,mla%mba%rr,the_bc_level, &
                                 icomp=rhoh_comp, &
                                 bcomp=dm+rhoh_comp, &
                                 nc=1, &
                                 ng=state(1)%ng)

    end if

    call destroy(bpt)

  end subroutine makeTfromRhoP

  !----------------------------------------------------------------------------
  ! makeTfromRhoP_1d
  !----------------------------------------------------------------------------
  subroutine makeTfromRhoP_1d(state,lo,hi,ng,p0,updateRhoH)

    use variables,     only: rho_comp, spec_comp, temp_comp, pi_comp, rhoh_comp
    use eos_module,    only: eos_input_rp, eos
    use eos_type_module
    use network,       only: nspec
    use probin_module, only: use_pprime_in_tfromp

    integer, intent(in) :: lo(:), hi(:), ng, updateRhoH
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:)
    
    ! Local variables
    integer :: i
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state

    do i = lo(1), hi(1)

       ! (rho, p) --> T

       eos_state%rho   = state(i,rho_comp)
       eos_state%T     = state(i,temp_comp)
       eos_state%xn(:) = state(i,spec_comp:spec_comp+nspec-1)/eos_state%rho

       if (use_pprime_in_tfromp) then
          eos_state%p = p0(i) + state(i,pi_comp)
       else
          eos_state%p = p0(i)
       endif

       pt_index(:) = (/i, -1, -1/)

       call eos(eos_input_rp, eos_state, pt_index)

       state(i,temp_comp) = eos_state%T
       if (updateRhoH .eq. 1) then
          state(i,rhoh_comp) = eos_state%rho*eos_state%h
       end if

    enddo

  end subroutine makeTfromRhoP_1d

  !----------------------------------------------------------------------------
  ! makeTfromRhoP_2d
  !----------------------------------------------------------------------------
  subroutine makeTfromRhoP_2d(state,lo,hi,ng,p0,updateRhoH)

    use variables,     only: rho_comp, spec_comp, temp_comp, pi_comp, rhoh_comp
    use eos_module,    only: eos_input_rp, eos
    use eos_type_module
    use network,       only: nspec
    use probin_module, only: use_pprime_in_tfromp

    integer, intent(in) :: lo(:), hi(:), ng, updateRhoH
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:)
    
    ! Local variables
    integer :: i, j
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          
          ! (rho, p) --> T
          
          eos_state%rho  = state(i,j,rho_comp)
          eos_state%T = state(i,j,temp_comp)
          eos_state%xn(:) = state(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho
          
          if (use_pprime_in_tfromp) then
             eos_state%p = p0(j) + state(i,j,pi_comp)
          else
             eos_state%p = p0(j)
          endif

          pt_index(:) = (/i, j, -1/)
          
          call eos(eos_input_rp, eos_state, pt_index)
          
          state(i,j,temp_comp) = eos_state%T
          if (updateRhoH .eq. 1) then
             state(i,j,rhoh_comp) = eos_state%rho*eos_state%h
          end if
          
       enddo
    enddo
    
  end subroutine makeTfromRhoP_2d

  !----------------------------------------------------------------------------
  ! makeTfromRhoP_3d
  !----------------------------------------------------------------------------
  subroutine makeTfromRhoP_3d(state,lo,hi,ng,p0,updateRhoH)

    use variables,     only: rho_comp, spec_comp, temp_comp, pi_comp, rhoh_comp
    use eos_module,    only: eos_input_rp, eos
    use eos_type_module
    use network,       only: nspec
    use fill_3d_module
    use probin_module, only: use_pprime_in_tfromp

    integer, intent(in) :: lo(:), hi(:), ng, updateRhoH
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:)

    ! Local variables
    integer :: i, j, k
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state

    !$OMP PARALLEL DO PRIVATE(i,j,k, eos_state, pt_index)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             ! (rho, p) --> T
             
             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             if (use_pprime_in_tfromp) then
                eos_state%p     = p0(k) + state(i,j,k,pi_comp)
             else
                eos_state%p     = p0(k)
             endif

             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)
             
             call eos(eos_input_rp, eos_state, pt_index)
             
             state(i,j,k,temp_comp) = eos_state%T
             if (updateRhoH .eq. 1) then
                state(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h
             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine makeTfromRhoP_3d

  !----------------------------------------------------------------------------
  ! makeTfromRhoP_3d_sphr
  !----------------------------------------------------------------------------
  subroutine makeTfromRhoP_3d_sphr(state,lo,hi,ng,p0,dx,updateRhoH)

    use variables,     only: rho_comp, spec_comp, temp_comp, pi_comp, rhoh_comp
    use eos_module,    only: eos_input_rp, eos
    use eos_type_module
    use network,       only: nspec
    use fill_3d_module
    use probin_module, only: use_pprime_in_tfromp

    integer, intent(in) :: lo(:), hi(:), ng, updateRhoH
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j, k
    real(kind=dp_t), allocatable :: p0_cart(:,:,:,:)
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,p0,p0_cart,lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k, eos_state, pt_index)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             ! (rho, p) --> T
             
             eos_state%rho  = state(i,j,k,rho_comp)
             eos_state%T = state(i,j,k,temp_comp)
             if (use_pprime_in_tfromp) then
                eos_state%p = p0_cart(i,j,k,1) + state(i,j,k,pi_comp)
             else
                eos_state%p = p0_cart(i,j,k,1)
             endif
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)
             
             call eos(eos_input_rp, eos_state, pt_index)
             
             state(i,j,k,temp_comp) = eos_state%T
             if (updateRhoH .eq. 1) then
                state(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h
             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
    deallocate(p0_cart)

  end subroutine makeTfromRhoP_3d_sphr


  !============================================================================
  ! makePfromRhoH
  !============================================================================
  subroutine makePfromRhoH(state,sold,peos,mla,the_bc_level)

    use variables,             only: foextrap_comp, temp_comp
    use bl_prof_module

    type(multifab)    , intent(in   ) :: state(:)
    type(multifab)    , intent(in   ) :: sold(:)
    type(multifab)    , intent(inout) :: peos(:)
    type(ml_layout)   , intent(inout) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    integer                  :: i,ng_s,ng_so,ng_p,n
    integer                  :: lo(mla%dim),hi(mla%dim),dm,nlevs
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: pnp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "makePfromRhoH")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s  = nghost(state(1))
    ng_so = nghost(sold(1))
    ng_p  = nghost(peos(1))

    do n=1,nlevs

       do i=1,nfabs(state(n))
          snp => dataptr(state(n),i)
          sop => dataptr(sold(n),i)
          pnp => dataptr(peos(n),i)
          lo = lwb(get_box(state(n),i))
          hi = upb(get_box(state(n),i))
          select case (dm)
          case (1)
             call makePfromRhoH_1d(snp(:,1,1,:), sop(:,1,1,temp_comp), pnp(:,1,1,1), &
                                   lo, hi, ng_s, ng_so, ng_p)
          case (2)
             call makePfromRhoH_2d(snp(:,:,1,:), sop(:,:,1,temp_comp), pnp(:,:,1,1), &
                                   lo, hi, ng_s, ng_so, ng_p)
          case (3)
             call makePfromRhoH_3d(snp(:,:,:,:), sop(:,:,:,temp_comp), pnp(:,:,:,1), &
                                   lo, hi, ng_s, ng_so, ng_p)
          end select
       end do

    end do

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,peos,mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=peos(1)%ng)

    call destroy(bpt)

  end subroutine makePfromRhoH

  !----------------------------------------------------------------------------
  ! makePfromRhoH_1d
  !----------------------------------------------------------------------------
  subroutine makePfromRhoH_1d(state,temp_old,peos,lo,hi,ng_s,ng_so,ng_p)

    use variables,  only: rho_comp, spec_comp, rhoh_comp
    use eos_module, only: eos_input_rh, eos
    use eos_type_module
    use network,    only: nspec

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_so, ng_p
    real (kind = dp_t), intent(in   ) ::    state(lo(1)-ng_s :,:)
    real (kind = dp_t), intent(in   ) :: temp_old(lo(1)-ng_so:)    
    real (kind = dp_t), intent(inout) ::     peos(lo(1)-ng_p :)
    
    ! Local variables
    integer :: i
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state
    
    do i = lo(1), hi(1)

       ! (rho, H) --> T, p
       eos_state%rho   = state(i,rho_comp)
       eos_state%T     = temp_old(i)
       eos_state%xn(:) = state(i,spec_comp:spec_comp+nspec-1)/eos_state%rho
       
       eos_state%h = state(i,rhoh_comp) / state(i,rho_comp)

       pt_index(:) = (/i, -1, -1/)
       
       call eos(eos_input_rh, eos_state, pt_index)
       
       peos(i) = eos_state%p
       
    enddo

  end subroutine makePfromRhoH_1d

  !----------------------------------------------------------------------------
  ! makePfromRhoH_2d
  !----------------------------------------------------------------------------
  subroutine makePfromRhoH_2d(state,temp_old,peos,lo,hi,ng_s,ng_so,ng_p)

    use variables,  only: rho_comp, spec_comp, rhoh_comp
    use eos_module, only: eos_input_rh, eos
    use eos_type_module
    use network,    only: nspec

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_so, ng_p
    real (kind = dp_t), intent(in   ) ::    state(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real (kind = dp_t), intent(in   ) :: temp_old(lo(1)-ng_so:,lo(2)-ng_so:)    
    real (kind = dp_t), intent(inout) ::     peos(lo(1)-ng_p :,lo(2)-ng_p :)
    
    ! Local variables
    integer :: i, j
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! (rho, H) --> T, p
          eos_state%rho   = state(i,j,rho_comp)
          eos_state%T     = temp_old(i,j)
          eos_state%xn(:) = state(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho

          eos_state%h = state(i,j,rhoh_comp) / state(i,j,rho_comp)

          pt_index(:) = (/i, j, -1/)

          call eos(eos_input_rh, eos_state, pt_index)

          peos(i,j) = eos_state%p

       enddo
    enddo

  end subroutine makePfromRhoH_2d

  !----------------------------------------------------------------------------
  ! makePfromRhoH_3d
  !----------------------------------------------------------------------------
  subroutine makePfromRhoH_3d(state,temp_old,peos,lo,hi,ng_s,ng_so,ng_p)

    use variables,  only: rho_comp, spec_comp, rhoh_comp
    use eos_module, only: eos_input_rh, eos
    use eos_type_module
    use network,    only: nspec
    use fill_3d_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_s, ng_so, ng_p
    real (kind = dp_t), intent(in   ) ::    state(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real (kind = dp_t), intent(in   ) :: temp_old(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:)
    real (kind = dp_t), intent(inout) ::     peos(lo(1)-ng_p :,lo(2)-ng_p :,lo(3)-ng_p :)

    ! Local variables
    integer :: i, j, k
    integer :: pt_index(MAX_SPACEDIM)
    type (eos_t) :: eos_state

    !$OMP PARALLEL DO PRIVATE(i,j,k, eos_state, pt_index)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             ! (rho, H) --> T, p
             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = temp_old(i,j,k)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             eos_state%h = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)

             pt_index(:) = (/i, j, k/)
             
             call eos(eos_input_rh, eos_state, pt_index)
             
             peos(i,j,k) = eos_state%p
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine makePfromRhoH_3d

end module rhoh_vs_t_module
