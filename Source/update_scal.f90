module update_scal_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_scal

contains

  subroutine update_scal(mla,nstart,nstop,sold,snew,sflux,scal_force,p0_new,p0_new_cart, &
                         dx,dt,the_bc_level,derivative_mode)

    use bl_prof_module
    use bl_constants_module
    use geometry,  only: spherical, polar
    use variables, only: spec_comp, rho_comp
    use network,   only: nspec
    use ml_restrict_fill_module

    type(ml_layout)   , intent(inout) :: mla
    integer           , intent(in   ) :: nstart, nstop
    type(multifab)    , intent(in   ) :: sold(:)
    type(multifab)    , intent(inout) :: snew(:)
    type(multifab)    , intent(in   ) :: sflux(:,:)
    type(multifab)    , intent(in   ) :: scal_force(:)
    real(kind = dp_t) , intent(in   ) :: p0_new(:,0:)
    type(multifab)    , intent(in   ) :: p0_new_cart(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:),dt
    type(bc_level)    , intent(in   ) :: the_bc_level(:)
    logical           , intent(in), optional :: derivative_mode

    ! local
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: sfpx(:,:,:,:)
    real(kind=dp_t), pointer :: sfpy(:,:,:,:)
    real(kind=dp_t), pointer :: sfpz(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: p0p(:,:,:,:)
 
    logical :: derivative

    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: i,n
    integer :: ng_so,ng_sn,ng_sf,ng_f,ng_p

    type(bl_prof_timer), save :: bpt

    call build(bpt, "update_scal")
    
    if (present(derivative_mode)) then
        derivative = derivative_mode
    else
        derivative = .false.
    endif
    
    dm = mla%dim
    nlevs = mla%nlevel

    ng_so = nghost(sold(1))
    ng_sn = nghost(snew(1))
    ng_sf = nghost(sflux(1,1))
    ng_f  = nghost(scal_force(1))
    ng_p  = nghost(p0_new_cart(1))

    do n=1,nlevs

       do i = 1, nfabs(sold(n))
          sop => dataptr(sold(n),i)
          snp => dataptr(snew(n),i)
          sfpx => dataptr(sflux(n,1),i)
          fp => dataptr(scal_force(n),i)
          lo =  lwb(get_box(sold(n),i))
          hi =  upb(get_box(sold(n),i))
          select case (dm)
          case (1)
             call update_scal_1d(nstart, nstop, &
                                  sop(:,1,1,:), ng_so, &
                                  snp(:,1,1,:), ng_sn, &
                                 sfpx(:,1,1,:), ng_sf, &
                                   fp(:,1,1,:), ng_f, &
                                 p0_new(n,:), lo, hi, dx(n,:), dt, derivative)
          case (2)
             sfpy => dataptr(sflux(n,2),i)
             if (polar .eq. 0) then
                call update_scal_2d(nstart, nstop, &
                                    sop(:,:,1,:), ng_so, snp(:,:,1,:), ng_sn, &
                                    sfpx(:,:,1,:), sfpy(:,:,1,:), ng_sf, &
                                    fp(:,:,1,:), ng_f, &
                                    p0_new(n,:), lo, hi, dx(n,:), dt, derivative)
             else
                p0p => dataptr(p0_new_cart(n), i)
                call update_scal_2d_polar(nstart, nstop, &
                                         sop(:,:,1,:), ng_so, snp(:,:,1,:), ng_sn, &
                                         sfpx(:,:,1,:), sfpy(:,:,1,:), &
                                         ng_sf, fp(:,:,1,:), ng_f, &
                                         p0p(:,:,1,1), ng_p, lo, hi, dx(n,:), dt, derivative)
             end if
          case (3)
             sfpy => dataptr(sflux(n,2),i)
             sfpz => dataptr(sflux(n,3),i)
             if (spherical .eq. 0) then
                call update_scal_3d_cart(nstart, nstop, &
                                         sop(:,:,:,:), ng_so, snp(:,:,:,:), ng_sn, &
                                         sfpx(:,:,:,:), sfpy(:,:,:,:), sfpz(:,:,:,:), &
                                         ng_sf, fp(:,:,:,:), ng_f, &
                                         p0_new(n,:), lo, hi, dx(n,:), dt, derivative)
             else
                p0p => dataptr(p0_new_cart(n), i)
                call update_scal_3d_sphr(nstart, nstop, &
                                         sop(:,:,:,:), ng_so, snp(:,:,:,:), ng_sn, &
                                         sfpx(:,:,:,:), sfpy(:,:,:,:), sfpz(:,:,:,:), &
                                         ng_sf, fp(:,:,:,:), ng_f, &
                                         p0p(:,:,:,1), ng_p, lo, hi, dx(n,:), dt, derivative)
             end if
          end select
       end do

    end do

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,snew,mla%mba%rr,the_bc_level, &
                              icomp=nstart, &
                              bcomp=dm+nstart, &
                              nc=nstop-nstart+1, &
                              ng=snew(1)%ng)

    ! do the same for density if we updated the species
    if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
       call ml_restrict_and_fill(nlevs,snew,mla%mba%rr,the_bc_level, &
                                 icomp=rho_comp, &
                                 bcomp=dm+rho_comp, &
                                 nc=1, &
                                 ng=snew(1)%ng)

    end if

    call destroy(bpt)

  end subroutine update_scal

  subroutine update_scal_1d(nstart,nstop,sold,ng_so,snew,ng_sn,sfluxx,ng_sf, &
                            force,ng_f,p0_new,lo,hi,dx,dt,derivative)

    use network,       only: nspec
    use eos_module,    only: eos, eos_input_rp
    use eos_type_module
    use probin_module, only: do_eos_h_above_cutoff, base_cutoff_density
    use variables,     only: spec_comp, rho_comp, rhoh_comp, temp_comp
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:)
    integer           , intent(in   ) :: ng_so, ng_sn, ng_sf, ng_f
    real (kind = dp_t), intent(in   ) ::   sold(lo(1)-ng_so:,:)
    real (kind = dp_t), intent(  out) ::   snew(lo(1)-ng_sn:,:)
    real (kind = dp_t), intent(in   ) :: sfluxx(lo(1)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  force(lo(1)-ng_f :,:)
    real (kind = dp_t), intent(in   ) :: p0_new(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: derivative

    integer            :: i, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: delta,frac,sumX
    logical            :: has_negative_species

    integer :: pt_index(MAX_SPACEDIM)
    
    type(eos_t) :: eos_state
    
    if (derivative) then

            do comp = nstart, nstop

               do i=lo(1),hi(1)
                     
                  divterm = (sfluxx(i+1,comp) - sfluxx(i,comp))/dx(1)

                  snew(i,comp) = (-divterm + force(i,comp))
                  
               end do

            enddo

            ! update density
            if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
               
               snew(:,rho_comp) = ZERO
               
               do i = lo(1), hi(1)

                  ! define the update to rho as the sum of the updates to (rho X)_i
                  do comp = nstart, nstop
                     snew(i,rho_comp) = snew(i,rho_comp) + snew(i,comp)
                  enddo
                  
               enddo
               
            end if

    else
            do comp = nstart, nstop

               do i=lo(1),hi(1)
                     
                  divterm = (sfluxx(i+1,comp) - sfluxx(i,comp))/dx(1)

                  snew(i,comp) = sold(i,comp) + dt*(-divterm + force(i,comp))
                  
               end do

            enddo
            if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then
               
               do i = lo(1), hi(1)
                     
                  if (snew(i,rho_comp) .le. base_cutoff_density) then

                     eos_state%rho = snew(i,rho_comp)
                     eos_state%T   = sold(i,temp_comp)
                     eos_state%p   = p0_new(i)
                     eos_state%xn  =snew(i,spec_comp:spec_comp+nspec-1)/eos_state%rho
                     
                     pt_index(:) = (/i, -1, -1/)
                     
                     ! (rho,P) --> T,h
                     call eos(eos_input_rp, eos_state, pt_index)
                     
                     snew(i,rhoh_comp) = snew(i,rho_comp) * eos_state%h
                     
                  end if
                  
               enddo
               
            end if
            
            ! update density
            if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
               
               snew(:,rho_comp) = sold(:,rho_comp)
               
               do i = lo(1), hi(1)

                  has_negative_species = .false.

                  ! define the update to rho as the sum of the updates to (rho X)_i
                  do comp = nstart, nstop
                     snew(i,rho_comp) = snew(i,rho_comp) + (snew(i,comp)-sold(i,comp))
                     if (snew(i,comp) .lt. ZERO) has_negative_species = .true.
                  enddo

                  ! enforce a density floor
                  if (snew(i,rho_comp) .lt. 0.5d0*base_cutoff_density) then
                     do comp = nstart, nstop
                        snew(i,comp) = snew(i,comp) * 0.5d0*base_cutoff_density/snew(i,rho_comp)
                     end do
                     snew(i,rho_comp) = 0.5d0*base_cutoff_density
                  end if

                  ! do not allow the species to leave here negative.
                  if (has_negative_species) then
                     do comp = nstart, nstop
                        if (snew(i,comp) .lt. ZERO) then
                           delta = -snew(i,comp)
                           sumX = ZERO 
                           do comp2 = nstart, nstop
                              if (comp2 .ne. comp .and. snew(i,comp2) .ge. ZERO) then
                                 sumX = sumX + snew(i,comp2)
                              end if
                           enddo
                           do comp2 = nstart, nstop
                              if (comp2 .ne. comp .and. snew(i,comp2) .ge. ZERO) then
                                 frac = snew(i,comp2) / sumX
                                 snew(i,comp2) = snew(i,comp2) - frac * delta
                              end if
                           enddo
                           snew(i,comp) = ZERO
                        end if
                     end do
                  end if

               enddo
               
            end if
    endif


  end subroutine update_scal_1d

  subroutine update_scal_2d(nstart,nstop,sold,ng_so,snew,ng_sn,sfluxx,sfluxy,ng_sf, &
                            force,ng_f,p0_new,lo,hi,dx,dt, derivative)

    use network,       only: nspec
    use eos_module,    only: eos, eos_input_rp
    use eos_type_module
    use probin_module, only: do_eos_h_above_cutoff, base_cutoff_density
    use variables,     only: spec_comp, rho_comp, rhoh_comp, temp_comp
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:)
    integer           , intent(in   ) :: ng_so, ng_sn, ng_sf, ng_f
    real (kind = dp_t), intent(in   ) ::   sold(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real (kind = dp_t), intent(  out) ::   snew(lo(1)-ng_sn:,lo(2)-ng_sn:,:)
    real (kind = dp_t), intent(in   ) :: sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) :: sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real (kind = dp_t), intent(in   ) :: p0_new(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: derivative

    integer            :: i, j, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: delta,frac,sumX
    logical            :: has_negative_species

    integer :: pt_index(MAX_SPACEDIM)

    type(eos_t) :: eos_state


    if (derivative) then

            do comp = nstart, nstop

               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                       
                     divterm = (sfluxx(i+1,j,comp) - sfluxx(i,j,comp))/dx(1) &
                             + (sfluxy(i,j+1,comp) - sfluxy(i,j,comp))/dx(2)

                      snew(i,j,comp) = (-divterm + force(i,j,comp))
                   
                   end do              
               end do
                
            enddo

            ! update density
            if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
               
               snew(:,:,rho_comp) = ZERO
               
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                          ! define the update to rho as the sum of the updates to (rho X)_i
                          do comp = nstart, nstop
                             snew(i,j,rho_comp) = snew(i,j,rho_comp) + snew(i,j,comp)
                          enddo
                  enddo                  
               enddo
               
            end if

    else


            do comp = nstart, nstop

               do j=lo(2),hi(2)
                  do i=lo(1),hi(1)
                     
                     divterm = (sfluxx(i+1,j,comp) - sfluxx(i,j,comp))/dx(1) &
                             + (sfluxy(i,j+1,comp) - sfluxy(i,j,comp))/dx(2)

                     snew(i,j,comp) = sold(i,j,comp) + dt*(-divterm + force(i,j,comp))
                     
                  end do
               end do

            enddo

            if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then
               
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     
                     if (snew(i,j,rho_comp) .le. base_cutoff_density) then
                        eos_state%rho = snew(i,j,rho_comp)
                        eos_state%T   = sold(i,j,temp_comp)
                        eos_state%p   = p0_new(j)
                        eos_state%xn  = snew(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho

                        pt_index(:) = (/i, j, -1/)
                        
                        ! (rho,P) --> T,h
                        call eos(eos_input_rp, eos_state, pt_index)
                        
                        snew(i,j,rhoh_comp) = snew(i,j,rho_comp) * eos_state%h
                        
                     end if
                     
                  enddo
               enddo
               
            end if
            
            ! update density
            if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
               
               snew(:,:,rho_comp) = sold(:,:,rho_comp)
                   
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)

                     has_negative_species = .false.

                     ! define the update to rho as the sum of the updates to (rho X)_i  
                     do comp = nstart, nstop
                        snew(i,j,rho_comp) = snew(i,j,rho_comp) + (snew(i,j,comp)-sold(i,j,comp))
                        if (snew(i,j,comp) .lt. ZERO) has_negative_species = .true.
                     enddo

                     ! enforce a density floor
                     if (snew(i,j,rho_comp) .lt. 0.5d0*base_cutoff_density) then
                        do comp = nstart, nstop
                           snew(i,j,comp) = snew(i,j,comp) * &
                                0.5d0*base_cutoff_density/snew(i,j,rho_comp)
                        end do
                        snew(i,j,rho_comp) = 0.5d0*base_cutoff_density
                     end if

                     ! do not allow the species to leave here negative.
                     if (has_negative_species) then
                        do comp = nstart, nstop
                           if (snew(i,j,comp) .lt. ZERO) then
                              delta = -snew(i,j,comp)
                              sumX = ZERO 
                              do comp2 = nstart, nstop
                                 if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                                    sumX = sumX + snew(i,j,comp2)
                                 end if
                              enddo
                              do comp2 = nstart, nstop
                                 if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                                    frac = snew(i,j,comp2) / sumX
                                    snew(i,j,comp2) = snew(i,j,comp2) - frac * delta
                                 end if
                              enddo
                              snew(i,j,comp) = ZERO
                           end if
                        end do
                     end if

                  enddo
               enddo

            end if
    endif
     
  end subroutine update_scal_2d

  subroutine update_scal_2d_polar(nstart,nstop,sold,ng_so,snew,ng_sn,sfluxx,sfluxy, &
                                 ng_sf,force,ng_f,p0_new_cart,ng_p,lo,hi,dx,dt,derivative)

    use network,       only: nspec
    use eos_module,    only: eos, eos_input_rp
    use eos_type_module
    use probin_module, only: do_eos_h_above_cutoff, base_cutoff_density
    use variables,     only: spec_comp, rho_comp, rhoh_comp, temp_comp
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:)
    integer           , intent(in   ) :: ng_so, ng_sn, ng_sf, ng_f, ng_p
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_so:,lo(2)-ng_so:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_sn:,lo(2)-ng_sn:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real (kind = dp_t), intent(in   )::   p0_new_cart(lo(1)-ng_p :,lo(2)-ng_p :)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: derivative

    integer            :: i, j, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: delta,frac,sumX
    logical            :: has_negative_species

    integer :: pt_index(MAX_SPACEDIM)
    type(eos_t) :: eos_state



    if (derivative) then

            !$OMP PARALLEL PRIVATE(i,j,divterm,comp) 
            do comp = nstart, nstop
               !$OMP DO
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                       
                     divterm = (sfluxx(i+1,j,comp) - sfluxx(i,j,comp))/dx(1) &
                             + (sfluxy(i,j+1,comp) - sfluxy(i,j,comp))/dx(2)

                      snew(i,j,comp) = (-divterm + force(i,j,comp))
                   
                   end do              
               end do
               !$OMP END DO NOWAIT    
            enddo
            !$OMP END PARALLEL

            ! update density
            if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
               
               snew(:,:,rho_comp) = ZERO
               !$OMP PARALLEL DO PRIVATE(i,j,comp)             
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                          ! define the update to rho as the sum of the updates to (rho X)_i
                          do comp = nstart, nstop
                             snew(i,j,rho_comp) = snew(i,j,rho_comp) + snew(i,j,comp)
                          enddo
                  enddo                  
               enddo
               !$OMP END PARALLEL DO
            end if

    else


            !
            ! is spherical
            !
            !$OMP PARALLEL PRIVATE(i,j,divterm,comp) 
            do comp = nstart, nstop
               !$OMP DO
                do j = lo(2), hi(2)
                    do i = lo(1), hi(1)

                        divterm = (sfluxx(i+1,j,comp) - sfluxx(i,j,comp))/dx(1) &
                                + (sfluxy(i,j+1,comp) - sfluxy(i,j,comp))/dx(2)
                                
                        snew(i,j,comp) = sold(i,j,comp) + dt * (-divterm + force(i,j,comp))

                    enddo
                enddo
               !$OMP END DO NOWAIT
            end do
            !$OMP END PARALLEL

            if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then

               !$OMP PARALLEL DO PRIVATE(i,j,eos_state,pt_index)
                do j = lo(2), hi(2)
                    do i = lo(1), hi(1)

                        if (snew(i,j,rho_comp) .le. base_cutoff_density) then
                            eos_state%rho   = snew(i,j,rho_comp)
                            eos_state%T     = sold(i,j,temp_comp)
                            eos_state%p     = p0_new_cart(i,j)
                            eos_state%xn(:) = snew(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho

                            pt_index(:) = (/i, j, -1/)

                            ! (rho,P) --> T,h
                            call eos(eos_input_rp, eos_state, pt_index)

                            snew(i,j,rhoh_comp) = snew(i,j,rho_comp) * eos_state%h

                        end if

                    enddo
                enddo
               !$OMP END PARALLEL DO

            end if

            ! update density
            if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

               snew(:,:,rho_comp) = sold(:,:,rho_comp)

               !$OMP PARALLEL DO PRIVATE(i,j,has_negative_species,comp,delta,sumX,comp2,frac)
               do j = lo(2), hi(2)
                 do i = lo(1), hi(1)

                    has_negative_species = .false.

                    ! define the update to rho as the sum of the updates to (rho X)_i
                    do comp = nstart, nstop
                        snew(i,j,rho_comp) = snew(i,j,rho_comp) &
                            + (snew(i,j,comp)-sold(i,j,comp))
                        if (snew(i,j,comp) .lt. ZERO) has_negative_species = .true.
                    enddo

                    ! enforce a density floor
                    if (snew(i,j,rho_comp) .lt. 0.5d0*base_cutoff_density) then
                        do comp = nstart, nstop
                            snew(i,j,comp) = snew(i,j,comp) * &
                                0.5d0*base_cutoff_density/snew(i,j,rho_comp)
                        end do
                        snew(i,j,rho_comp) = 0.5d0*base_cutoff_density
                    end if

                    ! do not allow the species to leave here negative.
                    if (has_negative_species) then
                        do comp = nstart, nstop
                            if (snew(i,j,comp) .lt. ZERO) then
                                delta = -snew(i,j,comp)
                                sumX = ZERO 
                                do comp2 = nstart, nstop
                                    if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                                        sumX = sumX + snew(i,j,comp2)
                                    end if
                                enddo
                                do comp2 = nstart, nstop
                                    if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                                        frac = snew(i,j,comp2) / sumX
                                        snew(i,j,comp2) = snew(i,j,comp2) - frac * delta
                                    end if
                                enddo
                                snew(i,j,comp) = ZERO
                            end if
                        end do
                    end if

                 enddo
               enddo
               !$OMP END PARALLEL DO
            end if
      endif

  end subroutine update_scal_2d_polar  
  
  
  
  subroutine update_scal_3d_cart(nstart,nstop,sold,ng_so,snew,ng_sn,sfluxx,sfluxy,sfluxz, &
                                 ng_sf,force,ng_f,p0_new,lo,hi,dx,dt, derivative)

    use network,       only: nspec
    use eos_module,    only: eos, eos_input_rp
    use eos_type_module
    use probin_module, only: do_eos_h_above_cutoff, base_cutoff_density
    use variables,     only: spec_comp, rho_comp, rhoh_comp, temp_comp
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:)
    integer           , intent(in   ) :: ng_so, ng_sn, ng_sf, ng_f
    real (kind = dp_t), intent(in   ) ::   sold(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real (kind = dp_t), intent(  out) ::   snew(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:,:)
    real (kind = dp_t), intent(in   ) :: sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) :: sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) :: sfluxz(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real (kind = dp_t), intent(in   ) :: p0_new(0:)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: derivative

    integer            :: i, j, k, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: delta,frac,sumX
    logical            :: has_negative_species

    integer :: pt_index(MAX_SPACEDIM)

    type(eos_t) :: eos_state

    if (derivative) then

            !$OMP PARALLEL PRIVATE(i,j,divterm,comp) 
            do comp = nstart, nstop
               !$OMP DO
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        
                        divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                                + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                                + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)
           
                        snew(i,j,k,comp) = (-divterm + force(i,j,k,comp))
                     end do
                   end do              
               end do
               !$OMP END DO NOWAIT    
            enddo
            !$OMP END PARALLEL 

            ! update density
            if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
               
               snew(:,:,:,rho_comp) = ZERO
               !$OMP PARALLEL DO PRIVATE(i,j,comp)             
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                          ! define the update to rho as the sum of the updates to (rho X)_i
                          do comp = nstart, nstop
                             snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) + snew(i,j,k,comp)
                          enddo
                      enddo
                  enddo                  
               enddo
               !$OMP END PARALLEL DO
            end if

    else
            !$OMP PARALLEL PRIVATE(i,j,k,divterm,comp) 
            do comp = nstart, nstop
               !$OMP DO
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        
                        divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                                + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                                + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)
           
                        snew(i,j,k,comp) = sold(i,j,k,comp) + dt * (-divterm + force(i,j,k,comp))
                        
                     enddo
                  enddo
               enddo
               !$OMP END DO NOWAIT
            end do
            !$OMP END PARALLEL
            
            if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then

               !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state,pt_index)
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)

                        if (snew(i,j,k,rho_comp) .le. base_cutoff_density) then

                           eos_state%rho = snew(i,j,k,rho_comp)
                           eos_state%T   = sold(i,j,k,temp_comp)
                           eos_state%p   = p0_new(k)
                           eos_state%xn  = snew(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                           pt_index(:) = (/i, j, k/)

                           ! (rho,P) --> T,h
                           call eos(eos_input_rp, eos_state, pt_index)

                           snew(i,j,k,rhoh_comp) = snew(i,j,k,rho_comp) * eos_state%h

                        end if

                     enddo
                  enddo
               enddo
               !$OMP END PARALLEL DO

            end if


            ! update density
            if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

               snew(:,:,:,rho_comp) = sold(:,:,:,rho_comp)

               !$OMP PARALLEL DO PRIVATE(i,j,k,has_negative_species,comp,delta,sumX,comp2,frac)       
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)

                        has_negative_species = .false.

                        ! define the update to rho as the sum of the updates to (rho X)_i
                        do comp = nstart, nstop
                           snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) &
                                + (snew(i,j,k,comp)-sold(i,j,k,comp))
                           if (snew(i,j,k,comp) .lt. ZERO) has_negative_species = .true.
                        enddo

                        ! enforce a density floor
                        if (snew(i,j,k,rho_comp) .lt. 0.5d0*base_cutoff_density) then
                           do comp = nstart, nstop
                              snew(i,j,k,comp) = snew(i,j,k,comp) * &
                                   0.5d0*base_cutoff_density/snew(i,j,k,rho_comp)
                           end do
                           snew(i,j,k,rho_comp) = 0.5d0*base_cutoff_density
                        end if

                        ! do not allow the species to leave here negative.
                        if (has_negative_species) then
                           do comp = nstart, nstop
                              if (snew(i,j,k,comp) .lt. ZERO) then
                                 delta = -snew(i,j,k,comp)
                                 sumX = ZERO 
                                 do comp2 = nstart, nstop
                                    if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                                       sumX = sumX + snew(i,j,k,comp2)
                                    end if
                                 enddo
                                 do comp2 = nstart, nstop
                                    if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                                       frac = snew(i,j,k,comp2) / sumX
                                       snew(i,j,k,comp2) = snew(i,j,k,comp2) - frac * delta
                                    end if
                                 enddo
                                 snew(i,j,k,comp) = ZERO
                              end if
                           end do
                        end if

                     enddo
                  enddo
               enddo
               !$OMP END PARALLEL DO

            end if
      endif
  end subroutine update_scal_3d_cart

  subroutine update_scal_3d_sphr(nstart,nstop,sold,ng_so,snew,ng_sn,sfluxx,sfluxy,sfluxz, &
                                 ng_sf,force,ng_f,p0_new_cart,ng_p,lo,hi,dx,dt, derivative)

    use network,       only: nspec
    use eos_module,    only: eos, eos_input_rp
    use eos_type_module
    use probin_module, only: do_eos_h_above_cutoff, base_cutoff_density
    use variables,     only: spec_comp, rho_comp, rhoh_comp, temp_comp
    use pred_parameters
    use bl_constants_module

    integer           , intent(in   ) :: nstart, nstop, lo(:), hi(:)
    integer           , intent(in   ) :: ng_so, ng_sn, ng_sf, ng_f, ng_p
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_so:,lo(2)-ng_so:,lo(3)-ng_so:,:)
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_sn:,lo(2)-ng_sn:,lo(3)-ng_sn:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxx(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxy(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::  sfluxz(lo(1)-ng_sf:,lo(2)-ng_sf:,lo(3)-ng_sf:,:)
    real (kind = dp_t), intent(in   ) ::   force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real (kind = dp_t), intent(in   )::   p0_new_cart(lo(1)-ng_p :,lo(2)-ng_p :,lo(3)-ng_p :)
    real (kind = dp_t), intent(in   ) :: dt,dx(:)
    logical           , intent(in   ) :: derivative

    integer            :: i, j, k, comp, comp2
    real (kind = dp_t) :: divterm
    real (kind = dp_t) :: delta,frac,sumX
    logical            :: has_negative_species

    integer :: pt_index(MAX_SPACEDIM)
    type(eos_t) :: eos_state



    if (derivative) then

            !$OMP PARALLEL PRIVATE(i,j,divterm,comp) 
            do comp = nstart, nstop
               !$OMP DO
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        
                        divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                                + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                                + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)
           
                        snew(i,j,k,comp) = (-divterm + force(i,j,k,comp))
                      end do
                   end do              
               end do
               !$OMP END DO NOWAIT    
            enddo
            !$OMP END PARALLEL

            ! update density
            if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
               
               snew(:,:,:,rho_comp) = ZERO
               !$OMP PARALLEL DO PRIVATE(i,j,comp)             
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        
                          ! define the update to rho as the sum of the updates to (rho X)_i
                          do comp = nstart, nstop
                             snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) + snew(i,j,k,comp)
                          enddo
                      enddo
                  enddo                  
               enddo
               !$OMP END PARALLEL DO
            end if

    else

            !
            ! is spherical
            !
            !$OMP PARALLEL PRIVATE(i,j,k,divterm,comp) 
            do comp = nstart, nstop
               !$OMP DO
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)

                        divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                                + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                                + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)

                        snew(i,j,k,comp) = sold(i,j,k,comp) + dt * (-divterm + force(i,j,k,comp))

                     enddo
                  enddo
               enddo
               !$OMP END DO NOWAIT
            end do
            !$OMP END PARALLEL

            if ( do_eos_h_above_cutoff .and. (nstart .eq. rhoh_comp) ) then

               !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state,pt_index)
               do k = lo(3), hi(3) 
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)

                        if (snew(i,j,k,rho_comp) .le. base_cutoff_density) then
                           eos_state%rho   = snew(i,j,k,rho_comp)
                           eos_state%T     = sold(i,j,k,temp_comp)
                           eos_state%p     = p0_new_cart(i,j,k)
                           eos_state%xn(:) = snew(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                           pt_index(:) = (/i, j, k/)

                           ! (rho,P) --> T,h
                           call eos(eos_input_rp, eos_state, pt_index)

                           snew(i,j,k,rhoh_comp) = snew(i,j,k,rho_comp) * eos_state%h

                        end if

                     enddo
                  enddo
               enddo
               !$OMP END PARALLEL DO

            end if

            ! update density
            if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then

               snew(:,:,:,rho_comp) = sold(:,:,:,rho_comp)

               !$OMP PARALLEL DO PRIVATE(i,j,k,has_negative_species,comp,delta,sumX,comp2,frac)
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)

                        has_negative_species = .false.

                        ! define the update to rho as the sum of the updates to (rho X)_i
                        do comp = nstart, nstop
                           snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) &
                                + (snew(i,j,k,comp)-sold(i,j,k,comp))
                           if (snew(i,j,k,comp) .lt. ZERO) has_negative_species = .true.
                        enddo

                        ! enforce a density floor
                        if (snew(i,j,k,rho_comp) .lt. 0.5d0*base_cutoff_density) then
                           do comp = nstart, nstop
                              snew(i,j,k,comp) = snew(i,j,k,comp) * &
                                   0.5d0*base_cutoff_density/snew(i,j,k,rho_comp)
                           end do
                           snew(i,j,k,rho_comp) = 0.5d0*base_cutoff_density
                        end if

                        ! do not allow the species to leave here negative.
                        if (has_negative_species) then
                           do comp = nstart, nstop
                              if (snew(i,j,k,comp) .lt. ZERO) then
                                 delta = -snew(i,j,k,comp)
                                 sumX = ZERO 
                                 do comp2 = nstart, nstop
                                    if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                                       sumX = sumX + snew(i,j,k,comp2)
                                    end if
                                 enddo
                                 do comp2 = nstart, nstop
                                    if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                                       frac = snew(i,j,k,comp2) / sumX
                                       snew(i,j,k,comp2) = snew(i,j,k,comp2) - frac * delta
                                    end if
                                 enddo
                                 snew(i,j,k,comp) = ZERO
                              end if
                           end do
                        end if

                     enddo
                  enddo
               enddo
               !$OMP END PARALLEL DO
            end if
      endif
  end subroutine update_scal_3d_sphr

end module update_scal_module
