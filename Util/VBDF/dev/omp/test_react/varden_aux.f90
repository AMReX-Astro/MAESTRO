!This module provides auxiliary routines for the test_react varden().

module varden_aux
  !Global modules
  use BoxLib
  use bl_types
  use multifab_module
  use ml_layout_module
  use bl_constants_module
  implicit none

  !Global variables
  private
  type(multifab), allocatable, save :: react_s(:) !Multifab to hold react data

  integer, save :: rho_c, h_c, spec_c, t_c, p_c, omegadot_c, hnuc_c, &
                   lhnuc_c, hext_c, dxn_con_c, h_con_c, ncomps
  integer, save :: nlevs

  character(len=20), save, allocatable :: varnames(:)

  logical, save :: react_is_init = .false.

  !Subroutines
  public :: react_init, grid_init, therm_init, react_write, varden_close

contains

  !====== Public routines ======
  subroutine react_init(mla)
    !=== Data ===
    !Modules
    use network, only: nspec

    !Args
    type(ml_layout), intent(in) :: mla

    !Local variables
    integer :: n

    !=== Execution ===
    nlevs = mla%nlevel

    !Init react_s array keys and variable names array
    call keys_init()
    call names_init()

    !Init reaction multifab
    allocate(react_s(nlevs))
  
    do n = 1,nlevs
      call multifab_build(react_s(n), mla%la(n), ncomps, 1)
    end do
  
    react_is_init = .true.
  end subroutine react_init

  subroutine grid_init(mla, nlevs, bct, dx)
    !=== Data ===
    !Modules
    use network, only: nspec
    use probin_module, only: pmask, spherical_in, prob_lo, prob_hi, test_set
    use define_bc_module
    use geometry
    use box_util_module

    !Args
    type(ml_layout), intent(out)          :: mla
    integer        , intent(out)          :: nlevs
    type(bc_tower) , intent(out)          :: bct
    real(kind=dp_t), intent(out), pointer :: dx(:,:)

    !Local variables
    integer           :: n
    integer           :: dm
    type(ml_boxarray) :: mba

    !=== Execution ===
    !Grid/Geometry initializations (uses gr0_* file)
    call read_a_hgproj_grid(mba, test_set)
    print *, 'calling layout_build, ', pmask
    call ml_layout_build(mla,mba,pmask)
    
    !Check for proper nesting
    if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
       call bl_error('ERROR: fixed_grids not properly nested')
    end if
    
    !Initialize nlevs
    nlevs = mla%nlevel
    
    if (nlevs .ne. 1) then
      call bl_error('ERROR: only 1 level of refinement currently supported')
    end if

    !Initialize boundary conditions
    call initialize_bc(bct,nlevs,pmask)
    do n = 1,nlevs
       call bc_tower_level_build(bct,n,mla%la(n))
    end do


    !Initialize dm
    dm = mla%dim
    
    if (dm /= 3 .or. spherical_in /= 0) then
       call bl_error('ERROR: grid must be Cartesian and three-dimensional')
    endif

    !Initialize_dx
    call initialize_dx(dx,mba,nlevs)
 
    !Initialize base state grid properties
    !These need to be initialized before calling average()
    nr_fine = extent(mla%mba%pd(nlevs),dm)
    dr_fine = (prob_hi(dm)-prob_lo(dm)) / dble(nr_fine)

    !Clean up
    call destroy(mba)
  end subroutine grid_init

  subroutine therm_init(mla, s, tempbar, pbar)
    !=== Data ===
    !Modules
    use variables
    use probin_module, only: dens_max, temp_max, dens_min, temp_min, &
                             prob_hi, prob_lo 
    use network, only: nspec, spec_names
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use layout_module   , only: get_pd   

    !Args
    type(ml_layout),              intent(in   ) :: mla
    type(multifab) ,              intent(inout) :: s(:)
    real(kind=dp_t), allocatable, intent(  out) :: tempbar(:,:), pbar(:,:)

    !Local variables
    integer                     :: lo(mla%dim), hi(mla%dim)
    integer                     :: domlo(mla%dim), domhi(mla%dim)
    real(kind=dp_t), pointer    :: sp(:,:,:,:)
    integer                     :: n, i, ii, jj, kk
    integer                     :: res, nlevs
    real(kind=dp_t)             :: dlogrho, dlogT, dxn, summ
    real(kind=dp_t)             :: temp_zone, dens_zone
    real(kind=dp_t)             :: xn_zone(nspec, 0:extent(mla%mba%pd(1),3)-1)
    type (eos_t) :: eos_state

    !=== Execution ===

    ! get the lo and hi of the domain
    domlo = lwb(get_pd(get_layout(s(1))))
    domhi = upb(get_pd(get_layout(s(1))))

    !Initialize tempbar & pbar
    nlevs = mla%nlevel
    res   = domhi(3)-domlo(3)
    allocate(tempbar(nlevs,0:res))
    allocate(pbar(nlevs,0:res))
    tempbar = ZERO
    pbar    = ZERO

    !Density, temperature, and mass fraction increments
    dlogrho = (log10(dens_max) - log10(dens_min))/(extent(mla%mba%pd(1),1) - 1)
    dlogT   = (log10(temp_max) - log10(temp_min))/(extent(mla%mba%pd(1),2) - 1)
    call get_xn(xn_zone, domlo(3), domhi(3))

    !-- Initialization of the full thermodynamic grid --
    do n = 1, nlevs
       do i = 1, nfabs(s(n))
    
          sp  => dataptr(s(n), i)
 
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
        
          do kk = lo(3), hi(3)
             do jj = lo(2), hi(2)
                !Set the temperature
                temp_zone = 10.0_dp_t**(log10(temp_min) + dble(jj)*dlogT)

                do ii = lo(1), hi(1)
                   !Set the density
                   dens_zone = 10.0_dp_t**(log10(dens_min) + dble(ii)*dlogrho)
                   
                   !Call the EoS w/ rho, temp, & X as inputs
                   eos_state%T     = temp_zone
                   eos_state%rho   = dens_zone
                   eos_state%xn(:) = xn_zone(:,kk)
                   
                   call eos(eos_input_rt, eos_state, .false.)
    
                   !Initialize this element of the state
                   sp(ii,jj,kk,rho_comp)                    = dens_zone
                   sp(ii,jj,kk,rhoh_comp)                   = dens_zone * eos_state%h
                   sp(ii,jj,kk,spec_comp:spec_comp-1+nspec) = dens_zone * xn_zone(1:nspec,kk)
                   sp(ii,jj,kk,temp_comp)                   = temp_zone
                enddo
             enddo
          enddo

       enddo
    enddo
    
  end subroutine therm_init

  subroutine varden_close()
    !=== Data ===
    !Local variables
    integer :: n

    !=== Execution ===
    do n=1, nlevs
       call multifab_destroy(react_s(n))
    end do
    
    deallocate(react_s)
    deallocate(varnames)

    react_is_init = .false.
  end subroutine varden_close

  subroutine react_write(snew, sold, rho_omegadot, rho_Hnuc, rho_Hext, &
                        mla, dt, pref, the_bc_tower)
    !=== Data ===
    !Modules
    use variables
    use fabio_module
    use network, only: nspec
    use probin_module, only: run_prefix, &
                             prob_hi, prob_lo, &
                             do_heating, do_burning
    use define_bc_module, only: bc_tower

    !Args
    type(multifab)  , intent(in )   :: snew(:)
    type(multifab)  , intent(in )   :: sold(:)
    type(multifab)  , intent(in )   :: rho_omegadot(:)
    type(multifab)  , intent(in )   :: rho_Hnuc(:)
    type(multifab)  , intent(in )   :: rho_Hext(:)
    type(ml_layout) , intent(in )   :: mla
    real(kind=dp_t) , intent(in )   :: dt
    character(len=*), intent(in )   :: pref
    type(bc_tower)  , intent(in )   :: the_bc_tower

    !Local variables 
    real(kind=dp_t), pointer :: rsp(:,:,:,:), snp(:,:,:,:), sop(:,:,:,:), &
                                rwp(:,:,:,:), rnp(:,:,:,:), rep(:,:,:,:)
    real(kind=dp_t)          :: cur_rho, dxn
    integer                  :: n, i, j, ii, jj, kk
    integer, allocatable     :: lo(:), hi(:)

    !=== Execution ===
    if(.not. react_is_init) then
       call bl_error('ERROR: varden_aux reaction data not initialized')
    endif

    !Allocate bounding arrays.
    allocate(lo(mla%dim))
    allocate(hi(mla%dim))

    !Loop through cell by cell
    do n = 1, nlevs
       do i = 1, nfabs(react_s(n))
    
          rsp => dataptr(react_s(n),      i)
          snp => dataptr(snew(n),         i)
          sop => dataptr(sold(n),         i)
          rwp => dataptr(rho_omegadot(n), i)
          rnp => dataptr(rho_Hnuc(n),     i)
          rep => dataptr(rho_Hext(n),     i)
          
          lo = lwb(get_box(react_s(n), i))
          hi = upb(get_box(react_s(n), i))
          do kk = lo(3), hi(3)
             do jj = lo(2), hi(2)
                do ii = lo(1), hi(1)
                   cur_rho = snp(ii,jj,kk,rho_comp)

                   !Consistency checks
                   !1) Check that omegadot * dt = (change in mass fraction) 
                   !   for each species
                   do j=0, nspec-1
                      dxn = (snp(ii,jj,kk,spec_comp + j) / snp(ii,jj,kk,rho_comp)) - &
                           (sop(ii,jj,kk,spec_comp + j) / sop(ii,jj,kk,rho_comp))  
                      rsp(ii,jj,kk,dxn_con_c + j) = (rwp(ii,jj,kk,j + 1) / cur_rho) * dt - dxn
                   enddo
        
                   !2) Check that h_new = (h_old + Hnuc * dt + Hext * dt) 
                   rsp(ii,jj,kk,h_con_c) = snp(ii,jj,kk,rhoh_comp)   / cur_rho                   &
                                       - ((sop(ii,jj,kk,rhoh_comp) / sop(ii,jj,kk,rho_comp)) &
                                       + rnp(ii,jj,kk,1) / cur_rho * dt                           &
                                       + rep(ii,jj,kk,1) / cur_rho * dt)
           
                   !3) Check H_nuc and H_ext 
                   if(.not. do_heating .and. (rep(ii,jj,kk,1) /= ZERO) ) then
                      print *, 'ERROR: Non-zero H_ext with no heating'
                   endif
                   if(.not. do_burning .and. (rnp(ii,jj,kk,1) /= ZERO) ) then
                      print *, 'ERROR: Non-zero H_nuc with no burning'
                   endif
 
                   !Store reaction data
                   rsp(ii,jj,kk,rho_c) = cur_rho
                   rsp(ii,jj,kk,h_c)   = snp(ii,jj,kk,rhoh_comp) / cur_rho
                   do j=0, nspec-1
                      rsp(ii,jj,kk,spec_c + j)     = snp(ii,jj,kk,spec_comp + j) / cur_rho
                      rsp(ii,jj,kk,omegadot_c + j) = rwp(ii,jj,kk,j + 1)         / cur_rho
                   enddo
                   rsp(ii,jj,kk,t_c)    = snp(ii,jj,kk,temp_comp)
                   rsp(ii,jj,kk,hnuc_c) = rnp(ii,jj,kk,1) / cur_rho
                   if (rsp(ii,jj,kk,hnuc_c) >= ONE) then
                      rsp(ii,jj,kk,lhnuc_c) = log10(rsp(ii,jj,kk,hnuc_c))
                   else
                      rsp(ii,jj,kk,lhnuc_c) = ZERO
                   endif
                   rsp(ii,jj,kk,hext_c) = rep(ii,jj,kk,1) / cur_rho
                enddo
             enddo
          enddo

       enddo
    enddo

    !Write reaction data
    call fabio_ml_multifab_write_d(react_s,mla%mba%rr(:,1), &
                                   trim(run_prefix) // trim(pref), names=varnames, time=dt)

    call write_job_info(trim(run_prefix) // trim(pref), mla%mba, the_bc_tower, -ONE)
    
  end subroutine react_write

  !======= Private helper routines ======
  subroutine keys_init()
    use network, only: nspec

    rho_c      = 1 
    h_c        = 2
    spec_c     = 3
    t_c        = spec_c     + nspec 
    p_c        = t_c        + 1
    omegadot_c = p_c        + 1
    hnuc_c     = omegadot_c + nspec
    lhnuc_c    = hnuc_c     + 1
    hext_c     = lhnuc_c    + 1
    dxn_con_c  = hext_c     + 1
    h_con_c    = dxn_con_c  + nspec
    ncomps     = 8          + 3*nspec
    !Suggested components:
      !mode0_err_c (no  burn, no  heat)
      !mode1_err_c (no  burn, yes heat)
      !mode2_err_c (yes burn, no  heat)
      !mode3_err_c (yes burn, yes heat)
  end subroutine keys_init

  subroutine names_init()
    use network, only: nspec, spec_names
    integer :: i
    character(len=20) :: temp_buf
    
    allocate(varnames(ncomps))
    varnames(rho_c)        = 'density'
    varnames(h_c)          = 'enthalpy'
    do i = 0, nspec - 1
       write(temp_buf, *) i+1
       temp_buf = adjustl(temp_buf)
       varnames(spec_c     + i) = 'X_' // adjustl(trim(spec_names(i+1)))
       varnames(omegadot_c + i) = 'wdot(' // adjustl(trim(spec_names(i+1))) // ')'  
       varnames(dxn_con_c  + i) = adjustl(trim(spec_names(i+1))) // ' consistency'
    enddo
    varnames(t_c)          = 'temperature'
    varnames(p_c)          = 'pressure'
    varnames(hnuc_c)       = 'H_nuc'
    varnames(lhnuc_c)      = 'log10(H_nuc)'
    varnames(hext_c)       = 'H_ext'
    varnames(h_con_c)      = 'enthapy consistency'
  end subroutine names_init

  subroutine get_xn(xn_zone, lo, hi)
    !=== Data ===
    !Modules
    use network,       only: nspec, spec_names
    use probin_module, only: xin_file
    use bl_IO_module, only: unit_new

    !Args
    real(kind=dp_t), intent(  out) :: xn_zone(:,:)
    integer,         intent(in   ) :: lo, hi

    !Local data
    integer         :: un
    integer         :: i 
    real(kind=dp_t) :: summ, usr_in
    character (len=1024) :: line
    
    !=== Execution ===
    if (xin_file .eq. 'uniform') then
       summ = ZERO
       do i=1, nspec - 1
          print *, trim(adjustl(spec_names(i))) // ' mass fraction: '
          read(*,*) usr_in
          xn_zone(i,:) = usr_in
          summ = summ + usr_in
       enddo

       if(summ > ONE) then
          print *, 'ERROR: Mass fraction sum exceeds 1.0!'
          stop
       else
          xn_zone(nspec,:) = ONE - summ
       endif
    else
       ! read in an inputs file containing the mass fractions.
       ! each species is on its own line.
       ! Allow for comment lines with '#' in the first column
       un = unit_new()
       open(unit=un, file=xin_file, status='old')

       summ = ZERO

       !TODO: Add xn <= 1.0 error checking
       !TODO: Add proper cell count error checking
       i = 1
       do while (i <= nspec)
          ! read the line into a character buffer
          read (un,'(a)') line

          ! skip comments
          if (index(line, '#') == 1) cycle

          ! parse the line 
          read(line,*) xn_zone(i,:)

          i = i + 1
       enddo

       do i = 1, size(xn_zone,1)
          print *, i, sum(xn_zone(:,i))
       enddo

       close(un)
    endif
  end subroutine get_xn

end module varden_aux
