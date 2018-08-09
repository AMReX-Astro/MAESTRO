module rk_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private
  public :: update_rk, update_thermo_rk

contains

  subroutine update_rk(snew,sold,stemp,szero,step,dt,mla,derivative_mode)

    use bl_prof_module
    use bl_constants_module
    
    integer, intent(in)  :: step
    type(multifab)    , intent(inout) :: snew(:)
    type(multifab)    , intent(inout) :: sold(:)
    type(multifab)    , intent(inout) :: stemp(:)
    type(multifab)    , intent(in   ) :: szero(:)    
    type(ml_layout)   , intent(inout) :: mla
    real(kind=dp_t)   , intent(in   ) :: dt
    logical           , intent(in   ), optional :: derivative_mode

    ! local
    integer :: i,n
    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: ng_sn,ng_so,ng_st,ng_sz
    logical :: deriv
    real(kind=dp_t) :: fac, dtfac
    real(kind=dp_t), pointer:: snp(:,:,:,:)
    real(kind=dp_t), pointer:: sop(:,:,:,:)
    real(kind=dp_t), pointer:: stp(:,:,:,:)
    real(kind=dp_t), pointer:: szp(:,:,:,:)
        
    dm = mla%dim
    nlevs = mla%nlevel
    
    if (present(derivative_mode)) then
        deriv = derivative_mode
    else
        deriv = .false.
    endif
    
    
    select case(step)
        case(1)
                fac = SIXTH 
                dtfac = HALF
        case(2)
                fac = THIRD
                dtfac = HALF
        case(3)
                fac = THIRD
                dtfac = ONE
        case(4)
                fac = THIRD
                dtfac = ZERO
    end select    
    
    
    ng_sn = nghost(snew(1))
    ng_so = nghost(sold(1))
    ng_st = nghost(stemp(1))
    ng_sz = nghost(szero(1))    

    do n = 1, nlevs

       do i = 1, nfabs(snew(n))
          snp  => dataptr(snew(n),i)
          sop  => dataptr(sold(n),i)
          stp  => dataptr(stemp(n),i)
          szp  => dataptr(szero(n),i)

          lo = lwb(get_box(snew(n),i))
          hi = upb(get_box(snew(n),i))
          select case (dm)
          case (1)
             call rk_update_1d(snp(:,1,1,:),ng_sn,sop(:,1,1,:),ng_so, &
                                stp(:,1,1,:),ng_st,szp(:,1,1,:),ng_sz, &
                                fac,dtfac,lo,hi,dt,deriv)
          case (2)
             call rk_update_2d(snp(:,:,1,:),ng_sn,sop(:,:,1,:),ng_so, &
                                stp(:,:,1,:),ng_st,szp(:,:,1,:),ng_sz, &
                                fac,dtfac,lo,hi,dt,deriv)
          case (3)
             call rk_update_3d(snp(:,:,:,:),ng_sn,sop(:,:,:,:),ng_so, &
                                stp(:,:,:,:),ng_st,szp(:,:,:,:),ng_sz, &
                                fac,dtfac,lo,hi,dt,deriv)
          end select
       end do

    enddo

  end subroutine update_rk

  subroutine rk_update_1d(snew,ng_sn,sold,ng_so,stemp,ng_st,szero,ng_sz, &
                         fac,dtfac,lo,hi,dt,derivative)
    use probin_module, only: rk_timestep_fac 
    
    integer, intent(in) :: lo(:), hi(:), ng_sn, ng_so, ng_st, ng_sz  
    real (kind = dp_t), intent(inout) ::   snew(lo(1)-ng_sn:,:)  
    real (kind = dp_t), intent(inout) ::   sold(lo(1)-ng_so:,:)  
    real (kind = dp_t), intent(inout) ::   stemp(lo(1)-ng_st:,:)  
    real (kind = dp_t), intent(in   ) ::   szero(lo(1)-ng_sz:,:)
    real (kind = dp_t), intent(in   ) :: fac, dtfac, dt
    logical           , intent(in   ) :: derivative
    integer :: i

    do i = lo(1)-ng_sn, hi(1)+ng_sn
            if (derivative) then
                stemp(i,:) = stemp(i,:)  + snew(i,:) * fac * rk_timestep_fac * dt
                sold(i,:) = szero(i,:) + snew(i,:) * dtfac * rk_timestep_fac * dt     
            else
                snew(i,:) = snew(i,:) - sold(i,:)
                stemp(i,:) = stemp(i,:)  + snew(i,:) * fac * rk_timestep_fac 
                sold(i,:) = szero(i,:) + snew(i,:) * dtfac * rk_timestep_fac 
            endif
    enddo

  end subroutine rk_update_1d


  subroutine rk_update_2d(snew,ng_sn,sold,ng_so,stemp,ng_st,szero,ng_sz, &
                         fac,dtfac,lo,hi,dt,derivative)
    use probin_module, only: rk_timestep_fac

    integer, intent(in) :: lo(:), hi(:), ng_sn, ng_so, ng_st, ng_sz  
    real (kind = dp_t), intent(inout) ::   snew(lo(1)-ng_sn:,lo(2)-ng_sn:,:)  
    real (kind = dp_t), intent(inout) ::   sold(lo(1)-ng_so:,lo(2)-ng_so:,:)  
    real (kind = dp_t), intent(inout) ::   stemp(lo(1)-ng_st:,lo(2)-ng_st:,:)  
    real (kind = dp_t), intent(in   ) ::   szero(lo(1)-ng_sz:,lo(2)-ng_sz:,:)
    real (kind = dp_t), intent(in   ) :: fac, dtfac, dt
    logical           , intent(in   ) :: derivative

    integer :: i,j
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j = lo(2)-ng_sn, hi(2)+ng_sn
            do i = lo(1)-ng_sn, hi(1)+ng_sn
                    if (derivative) then
                        stemp(i,j,:) = stemp(i,j,:)  + snew(i,j,:) * fac * rk_timestep_fac * dt
                        sold(i,j,:) = szero(i,j,:) + snew(i,j,:) * dtfac * rk_timestep_fac * dt
                    else
                        snew(i,j,:) = snew(i,j,:) - sold(i,j,:)
                        stemp(i,j,:) = stemp(i,j,:)  + snew(i,j,:) * fac * rk_timestep_fac
                        sold(i,j,:) = szero(i,j,:) + snew(i,j,:) * dtfac * rk_timestep_fac
                    endif
            enddo
    enddo
    !$OMP END PARALLEL DO    
  end subroutine rk_update_2d

  subroutine rk_update_3d(snew,ng_sn,sold,ng_so,stemp,ng_st,szero,ng_sz, &
                         fac,dtfac,lo,hi,dt,derivative)
    use probin_module, only: rk_timestep_fac

    integer, intent(in) :: lo(:), hi(:), ng_sn, ng_so, ng_st, ng_sz  
    real (kind = dp_t), intent(inout) ::   snew(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:,:)  
    real (kind = dp_t), intent(inout) ::   sold(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)  
    real (kind = dp_t), intent(inout) ::   stemp(lo(1)-ng_st:,lo(2)-ng_st:,lo(3)-ng_st:,:)  
    real (kind = dp_t), intent(in   ) ::   szero(lo(1)-ng_sz:,lo(2)-ng_sz:,lo(3)-ng_sz:,:)
    real (kind = dp_t), intent(in   ) :: fac, dtfac, dt
    logical           , intent(in   ) :: derivative

    integer :: i,j,k
    
    
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3)-ng_sn, hi(3)+ng_sn
            do j = lo(2)-ng_sn, hi(2)+ng_sn
                    do i = lo(1)-ng_sn, hi(1)+ng_sn
                            if (derivative) then
                                stemp(i,j,k,:) = stemp(i,j,k,:)  + snew(i,j,k,:) &
                                                    * fac   * rk_timestep_fac * dt
                                sold(i,j,k,:) = szero(i,j,k,:) + snew(i,j,k,:) &
                                                    * dtfac * rk_timestep_fac * dt
                            else
                                snew(i,j,k,:) = snew(i,j,k,:) - sold(i,j,k,:)
                                stemp(i,j,k,:) = stemp(i,j,k,:)  + snew(i,j,k,:) &
                                                    * fac   * rk_timestep_fac
                                sold(i,j,k,:) = szero(i,j,k,:) + snew(i,j,k,:) &
                                                    * dtfac * rk_timestep_fac
                                
                            endif
                    enddo
            enddo
    enddo
    !$OMP END PARALLEL DO    
  end subroutine rk_update_3d
  
  
    subroutine update_thermo_rk(s,p0,rho0,dx,mla,the_bc_level)

    use bl_prof_module
    use bl_constants_module
    use ml_restrict_fill_module
    use fill_3d_module, only : put_1d_array_on_cart
    use geometry,       only: polar, spherical
    use variables,      only: foextrap_comp,nscal,nspec,&
                              rho_comp,rhoh_comp,temp_comp,spec_comp
    
    type(multifab)    , intent(inout) :: s(:)
    real(kind = dp_t) , intent(in   ) :: p0(:,0:)
    real(kind = dp_t) , intent(in   ) :: rho0(:,0:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:)   
    type(ml_layout)   , intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    
    
    
    ! local
    integer :: i,n
    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: ng_s,ng_pc, ng_rc
    
    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: pcp(:,:,:,:)
    real(kind=dp_t), pointer:: rcp(:,:,:,:)
    
    
    type(multifab) :: p0_cart(mla%nlevel)
    type(multifab) :: rho0_cart(mla%nlevel)

    dm = mla%dim
    nlevs = mla%nlevel
    
    ng_s = nghost(s(1))
    
    if (polar .eq. 1 .or. spherical .eq. 1) then
            do n=1,nlevs
               call build(rho0_cart(n), get_layout(s(n)), 1, 1)  
               call build(p0_cart(n), get_layout(s(n)), 1, 1)  
            end do

            ng_rc = nghost(rho0_cart(1))            
            ng_pc = nghost(p0_cart(1))
            call put_1d_array_on_cart(p0,p0_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
            call put_1d_array_on_cart(rho0,rho0_cart,dm+rho_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
    endif
    
    do n = 1, nlevs

       do i = 1, nfabs(s(n))
          sp  => dataptr(s(n),i)
          
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (1)
             call rk_update_thermo_1d(sp(:,1,1,:),ng_s,p0(n,:),rho0(n,:),lo,hi)
          case (2)
             if (polar .eq. 1) then
                    pcp  => dataptr(p0_cart(n),i)
                    rcp  => dataptr(rho0_cart(n),i)
                    call rk_update_thermo_2d_polar(sp(:,:,1,:),ng_s,&
                                   pcp(:,:,1,1),ng_pc,rcp(:,:,1,1),ng_rc,lo,hi)
             else
                    call rk_update_thermo_2d(sp(:,:,1,:),ng_s,p0(n,:),rho0(n,:),lo,hi)
             endif
          case (3)
             if (spherical .eq. 1) then
                    pcp  => dataptr(p0_cart(n),i)
                    rcp  => dataptr(rho0_cart(n),i)
                    call rk_update_thermo_3d_spherical(sp(:,:,:,:),ng_s,&
                                   pcp(:,:,:,1),ng_pc,rcp(:,:,:,1),ng_rc,lo,hi)
             else
                    call rk_update_thermo_3d(sp(:,:,:,:),ng_s,p0(n,:),rho0(n,:),lo,hi)
             endif
          end select
       end do

    enddo
    
    ! restrict and fill all ghostcells. we potentially change rho, rhoX, rhoH and temp 
    ! rhoX
    call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                              icomp=spec_comp, &
                              bcomp=dm+spec_comp, &
                              nc=nspec, &
                              ng=s(1)%ng)
    ! rho
    call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                              icomp=rho_comp, &
                              bcomp=dm+rho_comp, &
                              nc=1, &
                              ng=s(1)%ng)

    ! rhoH
    call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                              icomp=rhoh_comp, &
                              bcomp=dm+rhoh_comp, &
                              nc=1, &
                              ng=s(1)%ng)
    ! temp
    call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                              icomp=temp_comp, &
                              bcomp=dm+temp_comp, &
                              nc=1, &
                              ng=s(1)%ng)
    
    
    if (polar .eq. 1 .or. spherical .eq. 1) then
            do n=1,nlevs
               call destroy(rho0_cart(n)) 
               call destroy(p0_cart(n))  
            end do
    endif
    
                              
  end subroutine update_thermo_rk


  subroutine rk_update_thermo_1d(s,ng_s,p0,rho0,lo,hi)
    
    use bl_constants_module
    use network,       only: nspec
    use eos_module,    only: eos, eos_input_rp
    use eos_type_module
    use probin_module, only: do_eos_h_above_cutoff, base_cutoff_density, use_pprime_in_tfromp
    use variables,     only: spec_comp, rho_comp, rhoh_comp, temp_comp, pi_comp
    
    integer, intent(in) :: lo(:), hi(:), ng_s  
    real (kind = dp_t), intent(inout) ::   s(lo(1)-ng_s:,:) 
    real (kind = dp_t), intent(in   ) ::   p0(0:)
    real (kind = dp_t), intent(in   ) ::   rho0(0:)
    
    integer :: i,comp,comp2
    integer :: pt_index(MAX_SPACEDIM)
    real(kind=dp_t) :: delta, sumX, frac
    logical :: has_negative_species
    
    type(eos_t) :: eos_state
    
    

    do i = lo(1), hi(1)

        has_negative_species = .false.

        ! define the update to rho as the sum of the updates to (rho X)_i
        do comp = spec_comp, spec_comp+nspec-1
           if (s(i,comp) .lt. ZERO) has_negative_species = .true.
        enddo

        ! enforce a density floor
        if (rho0(i) .lt. 0.5*base_cutoff_density) then
           do comp = spec_comp, spec_comp+nspec-1
              s(i,comp) = s(i,comp) * 0.5*base_cutoff_density/s(i,rho_comp)
           end do
           s(i,rho_comp) = 0.5*base_cutoff_density
        end if

        ! do not allow the species to leave here negative.
        if (has_negative_species) then
           do comp = spec_comp, spec_comp+nspec-1
             if (s(i,comp) .lt. ZERO) then
                 delta = -s(i,comp)
                 sumX = ZERO 
                 do comp2 = spec_comp, spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s(i,comp2) .ge. ZERO) then
                       sumX = sumX + s(i,comp2)
                    end if
                 enddo
                 do comp2 = spec_comp, spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s(i,comp2) .ge. ZERO) then
                       frac = s(i,comp2) / sumX
                       s(i,comp2) = s(i,comp2) - frac * delta
                    end if
                 enddo
                 s(i,comp) = ZERO
              end if
           end do
        end if

    enddo
               
    
    
    do i = lo(1), hi(1)
            
        eos_state%rho = s(i,rho_comp)
        eos_state%T   = s(i,temp_comp)
        if (use_pprime_in_tfromp) then
            eos_state%p = p0(i) + s(i,pi_comp)
        else
            eos_state%p = p0(i)
        endif
        eos_state%xn  = s(i,spec_comp:spec_comp+nspec-1)/eos_state%rho            
        
        pt_index(:) = (/i, -1, -1/)
                 
        ! (rho,P) --> T,h
        call eos(eos_input_rp, eos_state, pt_index)
                 
        if (do_eos_h_above_cutoff .and. s(i,rho_comp) .le. base_cutoff_density) then
            s(i,rhoh_comp) = s(i,rho_comp) * eos_state%h
        end if
        
        s(i,temp_comp) = eos_state%T
                                     
    enddo

  end subroutine rk_update_thermo_1d
   

  subroutine rk_update_thermo_2d(s,ng_s,p0,rho0,lo,hi)

    use bl_constants_module
    use network,       only: nspec
    use eos_module,    only: eos, eos_input_rp
    use eos_type_module
    use probin_module, only: do_eos_h_above_cutoff, base_cutoff_density, use_pprime_in_tfromp
    use variables,     only: spec_comp, rho_comp, rhoh_comp, temp_comp, pi_comp
    
    integer, intent(in) :: lo(:), hi(:), ng_s  
    real (kind = dp_t), intent(inout) ::   s(lo(1)-ng_s:,lo(2)-ng_s:,:) 
    real (kind = dp_t), intent(in   ) ::   p0(0:)
    real (kind = dp_t), intent(in   ) ::   rho0(0:)
    
    integer :: i,j,comp,comp2
    integer :: pt_index(MAX_SPACEDIM)
    real(kind=dp_t) :: delta, sumX, frac
    logical :: has_negative_species
    
    type(eos_t) :: eos_state
    
    
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)

            has_negative_species = .false.

            ! define the update to rho as the sum of the updates to (rho X)_i
            do comp = spec_comp, spec_comp+nspec-1
                if (s(i,j,comp) .lt. ZERO) has_negative_species = .true.
            enddo

            ! enforce a density floor
            if (rho0(j) .lt. 0.5*base_cutoff_density) then
               do comp = spec_comp, spec_comp+nspec-1
                  s(i,j,comp) = s(i,j,comp) * 0.5*base_cutoff_density/s(i,j,rho_comp)
               end do
               s(i,j,rho_comp) = 0.5*base_cutoff_density
            end if

            ! do not allow the species to leave here negative.
            if (has_negative_species) then
              do comp = spec_comp, spec_comp+nspec-1
                if (s(i,j,comp) .lt. ZERO) then
                  delta = -s(i,j,comp)
                  sumX = ZERO 
                  do comp2 = spec_comp, spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s(i,j,comp2) .ge. ZERO) then
                      sumX = sumX + s(i,j,comp2)
                    end if
                  enddo
                  do comp2 = spec_comp, spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s(i,j,comp2) .ge. ZERO) then
                      frac = s(i,j,comp2) / sumX
                      s(i,j,comp2) = s(i,j,comp2) - frac * delta
                     end if
                  enddo
                  s(i,j,comp) = ZERO
                end if
              end do
            end if
        enddo
    enddo
               
    
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)
                    
            eos_state%rho = s(i,j,rho_comp)
            eos_state%T   = s(i,j,temp_comp)
            if (use_pprime_in_tfromp) then
                eos_state%p = p0(j) + s(i,j,pi_comp)
            else
                eos_state%p = p0(j)
            endif
            eos_state%xn  = s(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho            
                    
            pt_index(:) = (/i, j, -1/)
                             
            ! (rho,P) --> T,h
            call eos(eos_input_rp, eos_state, pt_index)
                             
            if (do_eos_h_above_cutoff .and. s(i,j,rho_comp) .le. base_cutoff_density) then
               s(i,j,rhoh_comp) = s(i,j,rho_comp) * eos_state%h
            end if
                    
            s(i,j,temp_comp) = eos_state%T
        enddo                         
    enddo

  end subroutine rk_update_thermo_2d


 subroutine rk_update_thermo_2d_polar(s,ng_s,p0_cart,ng_pc,rho0_cart,ng_rc,lo,hi)

    use bl_constants_module
    use network,       only: nspec
    use eos_module,    only: eos, eos_input_rp
    use eos_type_module
    use probin_module, only: do_eos_h_above_cutoff, base_cutoff_density, use_pprime_in_tfromp
    use variables,     only: spec_comp, rho_comp, rhoh_comp, temp_comp, pi_comp
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_pc, ng_rc  
    real (kind = dp_t), intent(inout) ::   s(lo(1)-ng_s:,lo(2)-ng_s:,:) 
    real (kind = dp_t), intent(in   ) ::   p0_cart(lo(1)-ng_pc:,lo(2)-ng_pc:) 
    real (kind = dp_t), intent(in   ) ::   rho0_cart(lo(1)-ng_rc:,lo(2)-ng_rc:)
        
    integer :: i,j,comp,comp2
    integer :: pt_index(MAX_SPACEDIM)
    real(kind=dp_t) :: delta, sumX, frac
    logical :: has_negative_species
    
    type(eos_t) :: eos_state
    
    
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)

            has_negative_species = .false.

            ! define the update to rho as the sum of the updates to (rho X)_i
            do comp = spec_comp, spec_comp+nspec-1
                if (s(i,j,comp) .lt. ZERO) has_negative_species = .true.
            enddo

            ! enforce a density floor
            if (rho0_cart(i,j) .lt. 0.5*base_cutoff_density) then
               do comp = spec_comp, spec_comp+nspec-1
                  s(i,j,comp) = s(i,j,comp) * 0.5*base_cutoff_density/s(i,j,rho_comp)
               end do
               s(i,j,rho_comp) = 0.5*base_cutoff_density
            end if

            ! do not allow the species to leave here negative.
            if (has_negative_species) then
              do comp = spec_comp, spec_comp+nspec-1
                if (s(i,j,comp) .lt. ZERO) then
                  delta = -s(i,j,comp)
                  sumX = ZERO 
                  do comp2 = spec_comp, spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s(i,j,comp2) .ge. ZERO) then
                      sumX = sumX + s(i,j,comp2)
                    end if
                  enddo
                  do comp2 = spec_comp, spec_comp+nspec-1
                    if (comp2 .ne. comp .and. s(i,j,comp2) .ge. ZERO) then
                      frac = s(i,j,comp2) / sumX
                      s(i,j,comp2) = s(i,j,comp2) - frac * delta
                     end if
                  enddo
                  s(i,j,comp) = ZERO
                end if
              end do
            end if
        enddo
    enddo
               
    
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)
                    
            eos_state%rho = s(i,j,rho_comp)
            eos_state%T   = s(i,j,temp_comp)
            if (use_pprime_in_tfromp) then
                eos_state%p = p0_cart(i,j) + s(i,j,pi_comp)
            else
                eos_state%p = p0_cart(i,j)
            endif
            eos_state%xn  = s(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho            
                    
            pt_index(:) = (/i, j, -1/)
                             
            ! (rho,P) --> T,h
            call eos(eos_input_rp, eos_state, pt_index)
                             
            if (do_eos_h_above_cutoff .and. s(i,j,rho_comp) .le. base_cutoff_density) then
               s(i,j,rhoh_comp) = s(i,j,rho_comp) * eos_state%h
            end if
                    
            s(i,j,temp_comp) = eos_state%T
        enddo                         
    enddo

  end subroutine rk_update_thermo_2d_polar

  subroutine rk_update_thermo_3d(s,ng_s,p0,rho0,lo,hi)

    use bl_constants_module
    use network,       only: nspec
    use eos_module,    only: eos, eos_input_rp
    use eos_type_module
    use probin_module, only: do_eos_h_above_cutoff, base_cutoff_density, use_pprime_in_tfromp
    use variables,     only: spec_comp, rho_comp, rhoh_comp, temp_comp, pi_comp
    
    integer, intent(in) :: lo(:), hi(:), ng_s  
    real (kind = dp_t), intent(inout) ::   s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:) 
    real (kind = dp_t), intent(in   ) ::   p0(0:)
    real (kind = dp_t), intent(in   ) ::   rho0(0:)
    
    integer :: i,j,k,comp,comp2
    integer :: pt_index(MAX_SPACEDIM)
    real(kind=dp_t) :: delta, sumX, frac
    logical :: has_negative_species
    
    type(eos_t) :: eos_state
    
    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)

                has_negative_species = .false.

                ! define the update to rho as the sum of the updates to (rho X)_i
                do comp = spec_comp, spec_comp+nspec-1
                    if (s(i,j,k,comp) .lt. ZERO) has_negative_species = .true.
                enddo

                ! enforce a density floor
                if (rho0(k) .lt. 0.5*base_cutoff_density) then
                   do comp = spec_comp, spec_comp+nspec-1
                      s(i,j,k,comp) = s(i,j,k,comp) * 0.5*base_cutoff_density/s(i,j,k,rho_comp)
                   end do
                   s(i,j,k,rho_comp) = 0.5*base_cutoff_density
                end if

                ! do not allow the species to leave here negative.
                if (has_negative_species) then
                   do comp = spec_comp, spec_comp+nspec-1
                      if (s(i,j,k,comp) .lt. ZERO) then
                         delta = -s(i,j,k,comp)
                         sumX = ZERO 
                         do comp2 = spec_comp, spec_comp+nspec-1
                            if (comp2 .ne. comp .and. s(i,j,k,comp2) .ge. ZERO) then
                               sumX = sumX + s(i,j,k,comp2)
                            end if
                         enddo
                         do comp2 = spec_comp, spec_comp+nspec-1
                            if (comp2 .ne. comp .and. s(i,j,k,comp2) .ge. ZERO) then
                               frac = s(i,j,k,comp2) / sumX
                               s(i,j,k,comp2) = s(i,j,k,comp2) - frac * delta
                            end if
                         enddo
                         s(i,j,k,comp) = ZERO
                      end if
                   end do
                end if
            enddo
        enddo
    enddo
               
    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                        
                eos_state%rho = s(i,j,k,rho_comp)
                eos_state%T   = s(i,j,k,temp_comp)
                if (use_pprime_in_tfromp) then
                    eos_state%p = p0(k) + s(i,j,k,pi_comp)
                else
                    eos_state%p = p0(k)
                endif
                eos_state%xn  = s(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho            
                        
                pt_index(:) = (/i, j, k/)
                                 
                ! (rho,P) --> T,h
                call eos(eos_input_rp, eos_state, pt_index)
                                 
                if (do_eos_h_above_cutoff .and. s(i,j,k,rho_comp) .le. base_cutoff_density) then
                    s(i,j,k,rhoh_comp) = s(i,j,k,rho_comp) * eos_state%h
                end if
                        
                s(i,j,k,temp_comp) = eos_state%T
            enddo                         
        enddo
    enddo

  end subroutine rk_update_thermo_3d
  
  subroutine rk_update_thermo_3d_spherical(s,ng_s,p0_cart,ng_pc,rho0_cart,ng_rc,lo,hi)

    use bl_constants_module
    use network,       only: nspec
    use eos_module,    only: eos, eos_input_rp
    use eos_type_module
    use probin_module, only: do_eos_h_above_cutoff, base_cutoff_density, use_pprime_in_tfromp
    use variables,     only: spec_comp, rho_comp, rhoh_comp, temp_comp, pi_comp
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_pc, ng_rc
    real (kind = dp_t), intent(inout) ::   s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:) 
    real (kind = dp_t), intent(in   ) ::   p0_cart(lo(1)-ng_pc:,lo(2)-ng_pc:,lo(3)-ng_pc:) 
    real (kind = dp_t), intent(in   ) ::   rho0_cart(lo(1)-ng_rc:,lo(2)-ng_rc:,lo(3)-ng_rc:)
    
    integer :: i,j,k,comp,comp2
    integer :: pt_index(MAX_SPACEDIM)
    real(kind=dp_t) :: delta, sumX, frac
    logical :: has_negative_species
    
    type(eos_t) :: eos_state
    
    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)

                has_negative_species = .false.

                ! define the update to rho as the sum of the updates to (rho X)_i
                do comp = spec_comp, spec_comp+nspec-1
                    if (s(i,j,k,comp) .lt. ZERO) has_negative_species = .true.
                enddo

                ! enforce a density floor
                if (rho0_cart(i,j,k) .lt. 0.5*base_cutoff_density) then
                   do comp = spec_comp, spec_comp+nspec-1
                      s(i,j,k,comp) = s(i,j,k,comp) * 0.5*base_cutoff_density/s(i,j,k,rho_comp)
                   end do
                   s(i,j,k,rho_comp) = 0.5*base_cutoff_density
                end if

                ! do not allow the species to leave here negative.
                if (has_negative_species) then
                   do comp = spec_comp, spec_comp+nspec-1
                      if (s(i,j,k,comp) .lt. ZERO) then
                         delta = -s(i,j,k,comp)
                         sumX = ZERO 
                         do comp2 = spec_comp, spec_comp+nspec-1
                            if (comp2 .ne. comp .and. s(i,j,k,comp2) .ge. ZERO) then
                               sumX = sumX + s(i,j,k,comp2)
                            end if
                         enddo
                         do comp2 = spec_comp, spec_comp+nspec-1
                            if (comp2 .ne. comp .and. s(i,j,k,comp2) .ge. ZERO) then
                               frac = s(i,j,k,comp2) / sumX
                               s(i,j,k,comp2) = s(i,j,k,comp2) - frac * delta
                            end if
                         enddo
                         s(i,j,k,comp) = ZERO
                      end if
                   end do
                end if
            enddo
        enddo
    enddo
               
    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                        
                eos_state%rho = s(i,j,k,rho_comp)
                eos_state%T   = s(i,j,k,temp_comp)
                if (use_pprime_in_tfromp) then
                    eos_state%p = p0_cart(i,j,k) + s(i,j,k,pi_comp)
                else
                    eos_state%p = p0_cart(i,j,k)
                endif
                eos_state%xn  = s(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho            
                        
                pt_index(:) = (/i, j, k/)
                                 
                ! (rho,P) --> T,h
                call eos(eos_input_rp, eos_state, pt_index)
                                 
                if (do_eos_h_above_cutoff .and. s(i,j,k,rho_comp) .le. base_cutoff_density) then
                    s(i,j,k,rhoh_comp) = s(i,j,k,rho_comp) * eos_state%h
                end if
                        
                s(i,j,k,temp_comp) = eos_state%T
            enddo                         
        enddo
    enddo

  end subroutine rk_update_thermo_3d_spherical
     
end module rk_module
