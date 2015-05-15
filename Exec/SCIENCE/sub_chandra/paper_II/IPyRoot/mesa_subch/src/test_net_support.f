! ***********************************************************************
!
!   Copyright (C) 2009  Bill Paxton, Frank Timmes
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful, 
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module test_net_support
      use chem_def
      use chem_lib
      use net_def
      use net_lib
      use const_def
      use rates_def
      
      implicit none
         
      logical, parameter :: extended_set = .true.
      logical, parameter :: sorted = .true.
      
      integer, parameter :: max_files = 20, max_cnt = 100000      
      
      logical :: qt
      
      character (len=256) :: eos_file_prefix, cache_suffix
      

      character (len=64) :: net_file
      
      integer :: handle
      type (Net_General_Info), pointer  :: g
      integer :: species, num_reactions
      
      integer, dimension(:), pointer :: net_iso, chem_id, isos_to_show

      integer :: which_rates_choice
      
      integer, pointer :: 
     >   reaction_id(:), reaction_table(:), rates_to_show(:), which_rates(:)

      real(dp) :: abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, eps_neu_total
      real(dp), dimension(:), pointer :: 
     >      xin, xin_copy, d_eps_nuc_dx, dxdt, d_dxdt_dRho, d_dxdt_dT
      real(dp), pointer :: d_dxdt_dx(:, :)  

      real(dp), dimension(:), pointer :: rho_vector, T_vector

      integer :: nrates_to_show, nisos_to_show

      real(dp) :: test_logT, test_logRho

      integer :: screening_mode
      real(dp), parameter :: theta_e_for_graboske_et_al = 1 ! for nondegenerate

      contains
      

      subroutine do_test_net(do_plots, symbolic)
         logical, intent(in) :: do_plots, symbolic
         call set_composition(g, species, xin)
         eta = 0
         if (do_plots) then
            call Create_Plot_Files(species, xin, g)
         else
            call Do_One_Net(species, g, symbolic)
         end if
      end subroutine do_test_net 
      
      
      subroutine test_net_setup(net_file_in)
         use rates_lib, only: rates_init
         character (len=*), intent(in) :: net_file_in

         integer :: info, i, ierr
         
         net_file = net_file_in

         allocate(net_iso(num_chem_isos), isos_to_show(num_chem_isos), chem_id(num_chem_isos))
        
         !Init net module 
         call net_init(ierr)
         if (ierr /= 0) stop 1
        
         !Handle is integer which is an index into the list of General_Net_Info
         !array.
         handle = alloc_net_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_net_handle failed'
            stop 2
         end if
         
         call net_start_def(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_start_def failed'
            stop 2
         end if
         
         call read_net_file(net_file, handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'read_net_file failed ', trim(net_file)
            stop 2
         end if
         
         call net_finish_def(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_finish_def failed'
            stop 2
         end if
     
        allocate(reaction_id(rates_reaction_id_max), reaction_table(rates_reaction_id_max))
        allocate(rates_to_show(rates_reaction_id_max), which_rates(rates_reaction_id_max))
      
         which_rates(:) = which_rates_choice
         
         !which_rates(irpp_to_he3) = 3

         call net_set_which_rates(handle, which_rates, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_set_which_rate_f17pg failed'
            stop 2
         end if
        
         cache_suffix = ''
         call net_setup_tables(handle, 'rate_tables', cache_suffix, info)
         if (ierr /= 0) then
            write(*,*) 'net_setup_tables failed'
            stop 2
         end if
         
         call net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_ptr failed'
            stop 2
         end if
         
         species = g % num_isos
         num_reactions = g% num_reactions
         
         call get_chem_id_table(handle, species, chem_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_chem_id_table failed'
            stop 2
         end if
         
         call get_net_iso_table(handle, net_iso, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_net_iso_table failed'
            stop 2
         end if
                  
         call get_reaction_id_table(handle, num_reactions, reaction_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_reaction_id_table failed'
            stop 2
         end if
         
         call get_net_reaction_table(handle, reaction_table, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_net_reaction_table failed'
            stop 2
         end if
         
         allocate(
     >      xin(species), xin_copy(species), d_eps_nuc_dx(species), 
     >      dxdt(species), d_dxdt_dRho(species), d_dxdt_dT(species), d_dxdt_dx(species, species))
     
      end subroutine test_net_setup
      

      subroutine Setup_eos(handle)
         ! allocate and load the eos tables
         use eos_def
         use eos_lib
         integer, intent(out) :: handle

         integer :: ierr
         logical, parameter :: use_cache = .true.
         
         eos_file_prefix = 'mesa'

         call eos_init(eos_file_prefix, '', '', use_cache, ierr)
         if (ierr /= 0) then
            write(*,*) 'eos_init failed in Setup_eos'
            stop 1
         end if
         
         handle = alloc_eos_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed trying to allocate eos handle'
            stop 1
         end if
      
      end subroutine Setup_eos
      
      
      subroutine test_net_cleanup
         integer :: ierr
         deallocate(xin)
         deallocate(d_eps_nuc_dx)
         deallocate(dxdt)
         deallocate(d_dxdt_dRho)
         deallocate(d_dxdt_dT)
         deallocate(d_dxdt_dx)
         call free_net_handle(handle)
      end subroutine test_net_cleanup
      
      
      subroutine change_net(new_net_file)
         character (len=*), intent(in) :: new_net_file
         call test_net_cleanup
         call test_net_setup(new_net_file)
      end subroutine change_net

      
      subroutine Do_One_Net(species, g, symbolic)
         use chem_lib, only:composition_info
         integer, intent(in) :: species
         type (Net_General_Info), pointer  :: g
         logical, intent(in) :: symbolic

         real(dp) :: logRho, logT, Rho, T, sum, mass_correction,
     >      eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, xh, xhe

         integer :: info, i, j, k, lwork, chem_id(species), num_reactions
         real(dp), dimension(species) :: dabar_dx, dzbar_dx, dmc_dx
         real(dp), dimension(:, :), pointer :: reaction_eps_nuc
         real(dp), dimension(:, :), pointer :: eps_nuc_categories
         real(dp), dimension(:, :), pointer :: rate_screened, rate_raw
         real(dp), pointer :: work(:), rate_factors(:), category_factors(:),
     >      actual_Qs(:), actual_neuQs(:)
         logical, pointer :: from_weaklib(:)
         
         num_reactions = g% num_reactions         

         logRho = test_logRho
         logT   = test_logT

         if (.not. qt) write(*,*)
         
         lwork = net_work_size(handle, info) 
         if (info /= 0) stop 1
         
         allocate(work(lwork), 
     >         rate_factors(num_reactions),  category_factors(num_categories),        
     >         rate_screened(num_rvs, num_reactions),         
     >         rate_raw(num_rvs, num_reactions), reaction_eps_nuc(num_rvs, num_reactions),
     >         eps_nuc_categories(num_rvs, num_categories),
     >         actual_Qs(num_reactions), actual_neuQs(num_reactions), from_weaklib(num_reactions),
     >         stat=info)
         if (info /= 0) stop 2
        
         call get_chem_id_table(handle, species, chem_id, info)
         if (info /= 0) stop 3

         call composition_info(
     >      species, chem_id, xin, xh, xhe, abar, zbar, z2bar, ye, mass_correction, sum, 
     >      dabar_dx, dzbar_dx, dmc_dx)

         Rho = 10**logRho
         T   = 10**logT
         
         rate_factors(:) = 1
         !Scale CAGO by 1.7
         i = reaction_table(ir_c12_ag_o16)
         rate_factors(i) = 1.7
         category_factors(:) = 1
         
         if (symbolic) then
            call net_get_symbolic_d_dxdt_dx(handle, species, num_reactions, 
     >            xin, T, logT, Rho, logRho, 
     >            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho,
     >            rate_factors, category_factors, 
     >            std_reaction_Qs, std_reaction_neuQs,
     >            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, 
     >            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, 
     >            screening_mode, theta_e_for_graboske_et_al,    
     >            rate_screened, rate_raw,
     >            reaction_eps_nuc, eps_nuc_categories, eps_neu_total,
     >            lwork, work, info)
         else
            call net_get_with_Qs(handle, species, num_reactions, 
     >            xin, T, logT, Rho, logRho, 
     >            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho,
     >            rate_factors, category_factors, 
     >            std_reaction_Qs, std_reaction_neuQs,
     >            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, 
     >            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, 
     >            screening_mode, theta_e_for_graboske_et_al,    
     >            rate_screened, rate_raw,
     >            reaction_eps_nuc, eps_nuc_categories, eps_neu_total,
     >            lwork, work, actual_Qs, actual_neuQs, from_weaklib, info)
         end if
         if (info /= 0) then
            write(*, *) 'bad return from net_get'
            stop 1
         end if
         
         if (symbolic .and..not. qt) then
            write(*,*) 'nonzero d_dxdt_dx entries'
            k = 0
            do j=1,species
               do i=1,species
                  if (d_dxdt_dx(i,j) /= 0) then
                     k = k + 1
                     write(*,'(a50,2i5)') 
     >                     trim(chem_isos% name(chem_id(i))) // 
     >                     ' ' // trim(chem_isos% name(chem_id(j))), i, j
                  end if
               end do
            end do
            write(*,*)
            write(*,'(a50,i5)') 'num non zeros', k
            write(*,*)
         else if (.not. qt) then
            call show_results(g, logT, logRho, species, num_reactions, xin, 
     >            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, 
     >            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, 
     >            rate_screened, rate_raw, reaction_eps_nuc,
     >            eps_nuc_categories, extended_set, sorted)         
         end if

         write(*,'(30x,4a20)') 'Q total', 'Q neutrino', 'Q total-neutrino'
         do i = 1, num_reactions
            if (from_weaklib(i)) then
               write(*,'(a30,99f20.10)') 'weaklib ' // trim(reaction_Name(reaction_id(i))), 
     >            actual_Qs(i), actual_neuQs(i), actual_Qs(i) - actual_neuQs(i)
            else
               write(*,'(a30,99f20.10)') trim(reaction_Name(reaction_id(i))), 
     >            actual_Qs(i), actual_neuQs(i), actual_Qs(i) - actual_neuQs(i)
            end if
         end do
         write(*,*)

         deallocate(work, rate_factors, category_factors, rate_screened, rate_raw, 
     >            eps_nuc_categories, reaction_eps_nuc, actual_Qs, actual_neuQs,
     >            from_weaklib)
 
         return

         write(*, *)
 1       format(a40, 6x, e25.10)
 2       format(a40, a6, e25.10)
         write(*, 1) 'abar', abar
         do i=1, species
            write(*, 2) 'dabar_dx', trim(chem_isos% name(chem_id(i))), dabar_dx(i)
         end do
         write(*, *)
         write(*, 1) 'zbar', zbar
         do i=1, species
            write(*, 2) 'dzbar_dx', trim(chem_isos% name(chem_id(i))), dzbar_dx(i)
         end do
         write(*, *)

      end subroutine Do_One_Net



      subroutine Create_Plot_Files(species, xin, g)
         integer, intent(in) :: species
         real(dp) :: xin(species)
         type (Net_General_Info), pointer  :: g

         integer:: i, j, k, lwork, ierr, num_reactions
         real(dp) :: T, logT, logT_min, logT_max
         real(dp) :: logRho_min, logRho_max, dlogT, dlogRho
         integer :: logT_points, logRho_points
         integer :: io, io_first, io_last, io_rho, io_tmp, io_params, num_out
         character (len=256) :: fname, dir

         real(dp), allocatable :: output_values(:, :, :)

         ! full range
         logT_max = 9.1d0
         logT_min = 6d0
         logRho_min = -3d0
         logRho_max = 10d0

         ! oxygen burning range
         logT_max = 9.5d0
         logT_min = 8d0
         logRho_min = -3d0
         logRho_max = 12d0

         ! test FL
         logT_max = 9d0
         logT_min = 7d0
         logRho_min = 5d0
         logRho_max = 10.2d0

         ! test FL
         logT_max = 8.4d0
         logT_min = 7.8d0
         logRho_min = 2d0
         logRho_max = 6d0

         ! test C+C
         logT_max = 10d0
         logT_min = 7.5d0
         logRho_min = 6d0
         logRho_max = 12d0

         logT_points = 251
         logRho_points = 251
         
         dir = 'plot_data'
         
         write(*, *) trim(dir)
         
 01   format(E30.22)

         dlogT = (logT_max-logT_min)/(logT_points-1)
         dlogRho = (logRho_max-logRho_min)/(logRho_points-1)         

         io_params = 40
         io_rho = 41
         io_tmp = 42
         io_first = 43

         fname = trim(dir) // '/' // 'params.data'
         open(unit=io_params, file=trim(fname))
         write(io_params, '(6f16.6)') 
     >         xin(net_iso(ih1)), xin(net_iso(ihe4)), xin(net_iso(ic12)),  
     >         xin(net_iso(in14)), xin(net_iso(io16))
         close(io_params)

         fname = trim(dir) // '/' // 'rho.data'
         open(unit=io_rho, file=trim(fname))

         fname = trim(dir) // '/' // 'tmp.data'
         open(unit=io_tmp, file=trim(fname))

         io = io_first-1
         io_last = Open_Files(io, dir)
         num_out = io_last - io_first + 1
         
         lwork = net_work_size(handle, ierr)
         if (ierr /= 0) stop 1
         num_reactions = g% num_reactions

         allocate(output_values(logRho_points, logT_points, num_out))
         
!xx$OMP PARALLEL DO PRIVATE(logT, T, j)
         do j=1, logT_points
            logT = logT_min + dlogT*(j-1)
            T = 10 ** logT

            call do_inner_loop(species, num_reactions, logT, T, j, output_values, g, xin, 
     >               logRho_points, logRho_min, dlogRho, lwork)
            
         end do
!xx$OMP END PARALLEL DO

         write(*, *) 'write the files'


         ! write out the results
         do j=1, logRho_points
            write(io_rho, 01) logRho_min + dlogRho*(j-1)
         end do
         close(io_rho)

         do i=1, logT_points
            write(io_tmp, 01) logT_min + dlogT*(i-1)
         enddo
         close(io_tmp)
         
!$OMP PARALLEL DO PRIVATE(k)
         do k = 1, num_out
            write(*, *) k
            write(io_first+k-1, '(e14.6)') output_values(:, :, k)
         end do
!$OMP END PARALLEL DO
         
         do io=io_first, io_last
            close(io)
         end do
         
      end subroutine Create_Plot_Files


      subroutine do_inner_loop(species, num_reactions, logT, T, j, output_values, g, xin, 
     >         logRho_points, logRho_min, dlogRho, lwork)
         integer, intent(in) :: species, num_reactions, lwork
         real(dp), intent(in) :: logT, T
         integer, intent(in) :: j, logRho_points
         type (Net_General_Info), pointer :: g
         real(dp), intent(OUT) :: output_values(:, :, :)
         real(dp), intent(in) :: xin(species), logRho_min, dlogRho

         integer :: i
         real(dp) :: logRho, Rho

         do i=1, logRho_points
            logRho = logRho_min + dlogRho*(i-1)
            Rho = 10 ** logRho
            call do_one_net_eval(species, num_reactions, logT, T, logRho, Rho, 
     >                  i, j, output_values, g, xin, lwork)
         enddo
         
      end subroutine do_inner_loop
      
      
      subroutine do_one_net_eval(species, num_reactions, logT, T, logRho, Rho, 
     >         i, j, output_values, g, xin, lwork)
         use chem_lib, only:composition_info
         integer, intent(in) :: species, num_reactions, lwork
         real(dp), intent(in) :: logT, T, logRho, Rho
         integer, intent(in) :: i, j
         type (Net_General_Info), pointer :: g
         real(dp), intent(OUT) :: output_values(:, :, :)
         real(dp), intent(in) :: xin(species)
      
         real(dp) :: abar, zbar, z2bar, ye, sum, mx
         real(dp) :: rate_factors(num_reactions), category_factors(num_categories)

         real(dp) :: eps_nuc
         real(dp) :: d_eps_nuc_dT
         real(dp) :: d_eps_nuc_dRho
         real(dp) :: d_eps_nuc_dx(species) 
            ! partial derivatives wrt mass fractions
         
         real(dp) :: dxdt(species)  
            ! rate of change of mass fractions caused by nuclear reactions
         real(dp) :: d_dxdt_dRho(species)
         real(dp) :: d_dxdt_dT(species)
         real(dp) :: d_dxdt_dx(species, species)  
            ! partial derivatives of rates wrt mass fractions
         real(dp) :: rate_screened(num_rvs, num_reactions)  
         real(dp) :: reaction_eps_nuc(num_rvs, num_reactions)  
         real(dp) :: rate_raw(num_rvs, num_reactions)  
         real(dp) :: eps_nuc_categories(num_rvs, num_categories)  

         integer :: info, k, h1, he4, chem_id(species)
         real(dp) :: xh, xhe, mass_correction
         real(dp), dimension(species) :: dabar_dx, dzbar_dx, dmc_dx
         real(dp), pointer :: work(:)
         real(dp), target :: work_ary(lwork)
         
         work => work_ary
         
         h1 = net_iso(ih1)
         he4 = net_iso(ihe4)
      
         call get_chem_id_table(handle, species, chem_id, info)
         if (info /= 0) stop 3

         call composition_info(
     >      species, chem_id, xin, xh, xhe, abar, zbar, z2bar, ye, 
     >      mass_correction, sum, dabar_dx, dzbar_dx, dmc_dx)
     
         rate_factors(:) = 1
         category_factors(:) = 1
         
         call net_get(handle, species, num_reactions, 
     >            xin, T, logT, Rho, logRho, 
     >            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho,
     >            rate_factors, category_factors, 
     >            std_reaction_Qs, std_reaction_neuQs,
     >            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, 
     >            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, 
     >            screening_mode, theta_e_for_graboske_et_al,    
     >            rate_screened, rate_raw,
     >            reaction_eps_nuc, eps_nuc_categories, eps_neu_total,
     >            lwork, work, info)
         if (info /= 0) then
            write(*, *) 'bad result from net_get'
            stop 1
         end if

         sum = 0d0
         mx = 0d0
         do k = species, 1, -1
            if (abs(dxdt(k)) > mx) mx = abs(dxdt(k))
            sum = sum + dxdt(k)
         end do
         if (mx <= 0) mx = 1d0
         if (sum-sum /= 0) then
            write(*, *) logRho, logT
            stop 1
         end if
         
         !eps_nuc = EPS3ALP(T,RHO,1d0,2d0) ! check F&L; use Y=1 and UE=2

         k =   1; output_values(i, j, k) = safe_log10(eps_nuc)
         
         if (.false.) then
            ! for testing, use Y=1 and UE=2
            !k = k+1; output_values(i, j, k) = safe_log10(eps_nuc/HK_3ALF(T,RHO,1d0))
            k = k+1; output_values(i, j, k) = eps_nuc/EPS3ALP(T,RHO,1d0,2d0)
         else
            k = k+1; output_values(i, j, k) =  sum / max(1d-20, mx)
         end if

         k = k+1; output_values(i, j, k) =  d_eps_nuc_dT * T / max(1d-20, eps_nuc)
         k = k+1; output_values(i, j, k) =  d_eps_nuc_dRho * Rho / max(1d-20, eps_nuc)
         k = k+1; output_values(i, j, k) =  safe_log10(abs(d_eps_nuc_dx(h1)))
         k = k+1; output_values(i, j, k) =  safe_log10(abs(d_eps_nuc_dx(he4)))
      
         if (k > max_files) then
            write(*, *) 'need to enlarge max_files'
            stop 1
         end if

      end subroutine do_one_net_eval
      
      
      FUNCTION EPS3ALP(T,RHO,Y,UE)
!
!     Gives the 3-alpha burning rate from Fushiki and Lamb (ap J,
!     317, 368). We found that our other routine, which
!     uses the screening of Salpeter and Van Horn is nearly
!     always in good agreement with this form. However, this
!     formula is also valid in the pycnonuclear regime. 
!
      IMPLICIT real(dp) (A-H,O-Z)
      include 'formats.dek'
      T6=T/1E6
      R6=RHO/1E6
      R6T=2.0*R6/UE
      R6T13=R6T**(1.0/3.0)
      R6T16=R6T**(1.0/6.0)
      T623=T6**(2.0/3.0)
      T632=T6**(3.0/2.0)
      T613=T6**(1.0/3.0)
      U=1.35*R6T13/T623
      IF (U.LT.1) THEN
         B1=(1-4.222E-2*T623)**2+2.643E-5*T6**(5.0/3.0)
         B2=(1-2.807E-2*T623)**2+2.704E-6*T6**(5.0/3.0)
         B1=16.16*EXP(-134.92/T6**(1.0/3.0))/(B1*T623)
         B2=244.6*(1.0+3.528E-3*T623)**5*EXP(-235.72/T613)/(B2*T623)         
         IF(5.458E3-R6T.GT.0) THEN
            B1=B1+EXP(-1065.1/T6)/T632
         ENDIF
         
         IF(1.836E4-R6T.GT.0) THEN
            B2=B2+EXP(-3336.4/T6)/T632
         ENDIF
         
         G1=B1*EXP(60.492*R6T13/T6)
         G2=B2*EXP(106.35*R6T13/T6)
      ELSE
         AF=(1.0/U**(1.5)+1.0)
         B1=(1.0-5.680E-2*R6T13)**2+8.815E-7*T6*T6
         B1=1.178*AF*EXP(-77.554/R6T16)/(B1*SQRT(T6))
         B2=(1.0-3.791E-2*R6T13)**2+5.162E-8*T6*T6
         B2=(1.0+5.070E-3*R6T13)**5*EXP(-135.08/R6T16)/B2
         B2=B2*13.48*AF/SQRT(T6)
         IF(5.458E3-R6T.GT.0) THEN
            G1=B1+EXP((60.492*R6T13-1065.1)/T6)/T632
         ELSE
            G1=B1
         ENDIF
         IF(1.836E4-R6T.GT.0) THEN
            G2=B2+EXP((106.35*R6T13-3336.4)/T6)/T632
         ELSE
            G2=B2
         ENDIF
       ENDIF
       EPS3ALP=5.120E29*G1*G2*Y**3*R6**2
       RETURN
       write(*,1) 'T', T
       write(*,1) 'RHO', RHO
       write(*,1) 'UE', UE
       write(*,1) 'Y', Y
       write(*,1) 'G1', G1
       write(*,1) 'G2', G2
       write(*,1) 'EPS3ALP', EPS3ALP
       write(*,*)
       stop
       END FUNCTION EPS3ALP
      
      
      FUNCTION HK_3ALF(T,RHO,Y) ! H&K eqn 6.80
      IMPLICIT real(dp) (A-H,O-Z)
      include 'formats.dek'
      T9 = T*1d-9
      HK_3ALF = 5.1d8*rho**2*Y**3*exp(-4.4027/T9)/T9**3 ! erg/g/s
      END FUNCTION HK_3ALF

         
      
      
      integer function Open_Files(io_start, dir)
         integer, intent(in) :: io_start
         character (len=256), intent(in) :: dir
         character (len=256) fname
         integer :: io
         io = io_start
         
         fname = trim(dir) // '/' // 'log_net_eps.data'
         io = io + 1; open(unit=io, file=trim(fname))

         fname = trim(dir) // '/' // 'sum_dxdt_div_max_dxdt.data'
         io = io+1; open(unit=io, file=trim(fname))

         fname = trim(dir) // '/' // 'd_lneps_dlnT.data'
         io = io + 1; open(unit=io, file=trim(fname))

         fname = trim(dir) // '/' // 'd_lneps_dlnRho.data'
         io = io + 1; open(unit=io, file=trim(fname))

         fname = trim(dir) // '/' // 'd_eps_nuc_dxh1.data'
         io = io + 1; open(unit=io, file=trim(fname))

         fname = trim(dir) // '/' // 'd_eps_nuc_dxhe4.data'
         io = io + 1; open(unit=io, file=trim(fname))
         
         Open_Files = io
         
      end function Open_Files
         

      real(dp) function safe_log10(x)
         real(dp), intent(in) :: x
         if (x <= 0) then
            safe_log10 = -99d0
         else
            safe_log10 = log10(x)
         end if
      end function safe_log10


      subroutine show_results(
     >         g, logT, logRho, species, num_reactions, xin, 
     >         eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, 
     >         dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, 
     >         rate_screened, rate_raw, reaction_eps_nuc,
     >         eps_nuc_categories, 
     >         extended_set, sorted)
         type (Net_General_Info), pointer  :: g
         real(dp), intent(in) :: logT, logRho
         integer, intent(in) :: species, num_reactions
         real(dp), intent(in) :: xin(species)
         real(dp), intent(in) :: eps_nuc
         real(dp), intent(in) :: d_eps_nuc_dT
         real(dp), intent(in) :: d_eps_nuc_dRho
         real(dp), intent(in) :: d_eps_nuc_dx(species) 
         real(dp), intent(in) :: dxdt(species) 
         real(dp), intent(in) :: d_dxdt_dRho(species) 
         real(dp), intent(in) :: d_dxdt_dT(species) 
         real(dp), intent(in) :: d_dxdt_dx(species, species) 
         real(dp), intent(in), dimension(num_rvs, num_categories) :: 
     >         eps_nuc_categories
         real(dp), intent(in), dimension(num_rvs, num_reactions) :: 
     >         rate_screened, rate_raw, reaction_eps_nuc 
         logical, intent(in) :: extended_set
         logical, intent(in) :: sorted

         write(*, *)
         call show_summary_results(logT, logRho,
     >         eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx) 
         
         if (extended_set) then
            write(*, *)
            call show_all_rates(
     >          g, num_reactions, rate_raw, rate_screened, reaction_eps_nuc, logT, logRho, sorted)
         end if
         
         write(*, *)
         call show_by_category(
     >            g, num_reactions, rate_screened, rate_raw, 
     >            eps_nuc_categories, 
     >            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, 
     >            sorted)
         
         if (.not. extended_set) return
         
         write(*, *)
         call show_dx_dt(g, species, xin, dxdt, sorted)
         
         write(*, *)
         call show_d_eps_nuc_dx(g, species, xin, d_eps_nuc_dx, sorted)

         write(*, *)

      end subroutine show_results

      
      subroutine show_summary_results(logT, logRho, 
     >         eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx) 
         real(dp), intent(in) :: logT, logRho
         real(dp), intent(in) :: eps_nuc
         real(dp), intent(in) :: d_eps_nuc_dT
         real(dp), intent(in) :: d_eps_nuc_dRho
         real(dp), intent(in) :: d_eps_nuc_dx(species) 

         real(dp) :: T, Rho, eps, d_eps_dt, d_eps_dd
         T = 10**logT; Rho = 10**logRho

         write(*, *)
         write(*, '(a40, f20.9)') 'log temp', logT
         write(*, '(a40, f20.9)') 'log rho', logRho
         eps = eps_nuc
         d_eps_dt = d_eps_nuc_dT
         d_eps_dd = d_eps_nuc_dRho
         write(*, *)
         write(*, '(a40, f20.9)') 'log(eps_nuc)', safe_log10(eps_nuc)
         write(*, *)
         write(*, '(a40, e20.9)') 'eps_nuc', eps_nuc
         write(*, *)
         write(*, '(a40, f20.9)') 'd_lneps_dlnT', d_eps_dt * T / eps
         write(*, '(a40, f20.9)') 'd_lneps_dlnRho', d_eps_dd * Rho / eps
      
      end subroutine show_summary_results

      
      subroutine show_all_rates(
     >      g, num_reactions, rate_raw, rate_screened, reaction_eps_nuc, logT, logRho, sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: num_reactions
         real(dp), intent(in) :: logT, logRho
         real(dp), dimension(num_rvs, num_reactions), intent(in) ::
     >      rate_raw, rate_screened, reaction_eps_nuc
         logical, intent(in) :: sorted
         
         real(dp), dimension(num_rvs, num_reactions) :: rfact
         integer :: i
         real(dp) :: T, Rho
         T = 10**logT; Rho = 10**logRho

         write(*, *)
         write(*, *) 'summary of log raw rates'
         write(*, *)
         call show_log_rates(g, rate_raw, T, Rho, sorted)
         write(*, *)
         write(*, *) 'summary of screening factors'
         write(*, *)
         do i=1,num_reactions
            if (rate_raw(i_rate, i) > 1d-50) then
               rfact(i_rate, i) = rate_screened(i_rate, i) / rate_raw(i_rate, i)
            else
               rfact(i_rate, i) = 1
            end if
         end do
         call show_rates(g, rfact, T, Rho, sorted)
         write(*, *)
         write(*, *) 'summary of log screened rates (reactions/gm/sec)'
         write(*, *)
         call show_log_rates(g, rate_screened, T, Rho, sorted)
         write(*, *)
         write(*, *) 'summary of log energy generation (reactions/gm/sec)'
         write(*, *)
         call show_log_rates(g, reaction_eps_nuc, T, Rho, sorted)
         write(*, *)


      end subroutine show_all_rates

      subroutine show_by_category(
     >         g, num_reactions, rate_screened, rate_raw, 
     >         eps_nuc_categories,
     >         eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, 
     >         sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: num_reactions
         real(dp), intent(in), dimension(num_rvs, num_reactions) :: 
     >         rate_screened, rate_raw
         real(dp), intent(in), dimension(num_rvs, num_categories) :: 
     >         eps_nuc_categories
         real(dp), intent(in) :: 
     >         eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx(species)
         logical, intent(in) :: sorted

         real(dp) :: mx
         integer :: k, j, jmx
         logical :: flgs(rates_reaction_id_max)
   
         integer :: info
         
         write(*, *)
         write(*, *) 'energy generation by category'
         write(*, *)
         write(*, '(a40, 3x, a20)') 'category', 'log rate (erg/g/sec)'
         write(*, *)
         flgs = .false.
         do k = 1, num_categories
            if (.not. sorted) then
               jmx = k
            else
               mx = -99d99; jmx = -1
               do j = 1, num_categories
                  if ((.not. flgs(j)) .and. eps_nuc_categories(i_rate, j) > mx) then
                     mx = eps_nuc_categories(i_rate, j); jmx = j
                  end if
               end do
               if (jmx <= 0) exit
               if (mx < 1) exit ! FOR TEST OUTPUT
               flgs(jmx) = .true.
            end if
            write(*, '(a40, 2x, f15.6, e15.6)') 
     >            trim(category_name(jmx)), safe_log10(eps_nuc_categories(i_rate, jmx)),
     >            eps_nuc_categories(i_rate, jmx)     
         end do
         write(*, *)
         write(*, '(a40, 2x, f15.6, e15.6)') 
     >            'log10(-photodisintegration)', safe_log10(-eps_nuc_categories(i_rate, iphoto)),
     >            -eps_nuc_categories(i_rate, iphoto)
         
         write(*, *)
         
      end subroutine show_by_category

      subroutine show_rates(g, rts, T, Rho, sorted)
         type (Net_General_Info), pointer  :: g
         real(dp), intent(in) :: rts(num_rvs, rates_reaction_id_max), T, Rho
         logical, intent(in) :: sorted
         
         logical :: flgs(rates_reaction_id_max)
         real(dp) :: mx
         integer :: k, j, jmx
         
         flgs = .false.
         
         do k = 1, g% num_reactions
            if (.not. sorted) then
               jmx = k; mx = rts(i_rate, jmx)
            else
               mx = -99d99; jmx = -1
               do j = 1, g% num_reactions
                  if ((.not. flgs(j)) .and. rts(i_rate, j) > mx) then
                     mx = rts(i_rate, j); jmx = j
                  end if
               end do
               if (jmx <= 0) exit
               if (mx < 1d-60) exit
               flgs(jmx) = .true.
            end if
            if (mx == 1) cycle
            write(*, '(a40, e20.9, 2e17.6)') trim(reaction_name(reaction_id(jmx))), mx
         end do
         
      end subroutine show_rates


      subroutine show_log_rates(g, rts, T, Rho, sorted)
         type (Net_General_Info), pointer  :: g
         real(dp), intent(in) :: rts(num_rvs, rates_reaction_id_max), T, Rho
         logical, intent(in) :: sorted
         
         logical :: flgs(rates_reaction_id_max)
         real(dp) :: mx
         integer :: k, j, jmx
         
         flgs = .false.
         
         do k = 1, g% num_reactions
            if (.not. sorted) then
               jmx = k; mx = rts(i_rate, jmx)
            else
               mx = -99d99; jmx = -1
               do j = 1, g% num_reactions
                  if ((.not. flgs(j)) .and. rts(i_rate, j) > mx) then
                     mx = rts(i_rate, j); jmx = j
                  end if
               end do
               if (jmx <= 0) exit
               if (mx < 1d-40) exit
               flgs(jmx) = .true.
            end if
            if (mx == 1) cycle
            write(*, '(a40, f20.9, 2e17.6)') trim(reaction_name(reaction_id(jmx))), safe_log10(mx)
         end do
         
      end subroutine show_log_rates
      
      
      subroutine show_dx_dt(g, species, xin, dxdt, sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: species
         real(dp), intent(in) :: xin(species)
         real(dp), intent(in) :: dxdt(species)
         logical, intent(in) :: sorted

         write(*, *)
         write(*, *) 'summary of isotope mass abundance changes'
         write(*, *)
         write(*, '(a40, 2(a17))') 'isotope', 'x initial', 'dx_dt   '
         call show_partials(g, species, xin, dxdt, .true., sorted)
         
      end subroutine show_dx_dt
      
      
      subroutine show_d_eps_nuc_dx(g, species, xin, d_eps_nuc_dx, sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: species
         real(dp), intent(in) :: xin(species)
         real(dp), intent(in) :: d_eps_nuc_dx(species)
         logical, intent(in) :: sorted

         write(*, *)
         write(*, *) 'summary of d_eps_nuc_dx'
         write(*, *)
         write(*, '(a40, a17)') 'isotope', 'd_eps_nuc_dx'
         call show_partials(g, species, xin, d_eps_nuc_dx, .false., sorted)
         
      end subroutine show_d_eps_nuc_dx
      
      
      subroutine show_partials(g, species, xin, derivs, initX_flag, sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: species
         real(dp), intent(in) :: xin(species)
         real(dp), intent(in) :: derivs(species)
         logical, intent(in) :: initX_flag, sorted

         real(dp) :: mx
         integer :: k, j, jmx
         integer, pointer :: chem_id(:)
         logical :: iflgs(species)
         chem_id => g% chem_id
         write(*, *)
         iflgs = .false.
         do k = 1, species
            if (.not. sorted) then
               jmx = k
            else
               mx = -99d99; jmx = -1
               do j = 1, species
                  if ((.not. iflgs(j)) .and. abs(derivs(j)) > mx) then
                     mx = abs(derivs(j)); jmx = j
                  end if
               end do
               if (jmx <= 0) exit
               if (mx < 1d-40) exit
            end if
            if (initX_flag) then
               write(*, '(a40, 2e17.6)') trim(chem_isos% name(chem_id(jmx))), xin(jmx), derivs(jmx)
            else
               write(*, '(a40, e25.14)') trim(chem_isos% name(chem_id(jmx))), derivs(jmx)
            end if
            iflgs(jmx) = .true.
         end do
         
      end subroutine show_partials


      subroutine set_composition(g, species, xin)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: species
         real(dp), intent(OUT) :: xin(species)
   
         real(dp) :: sum
         integer :: i, adjustment_iso
         
         eta = 0
         
         adjustment_iso = net_iso(ih1)

         if (net_file == 'data/sub_chandra.net') then   
         
            adjustment_iso = net_iso(ih1)
           
            xin = 0
            xin(net_iso(ih1))   =  0.05000D+00    ! h1  
            xin(net_iso(ihe3))  =  0.00000D+00    ! he3  
            xin(net_iso(ihe4))  =  0.84645D+00    ! he4  
            xin(net_iso(ic12))  =  0.04700D+00    ! c12
            xin(net_iso(in13))  =  0.00000D+00    ! n13
            xin(net_iso(in14))  =  0.00770D+00    ! n14
            xin(net_iso(ic14))  =  0.00085D+00    ! c14
            xin(net_iso(io16))  =  0.04700D+00    ! o16 
            xin(net_iso(io18))  =  0.00000D+00    ! o18 
            xin(net_iso(if18))  =  0.00000D+00    ! f18 
            xin(net_iso(ine20)) =  0.00000D+00    ! ne20 
            xin(net_iso(ine21)) =  0.00000D+00    ! ne21 
            xin(net_iso(ine22)) =  0.00100D+00    ! ne22 
            xin(net_iso(img24)) =  0.00000D+00    ! mg24

         else
            
            write(*,*) 'net_file ' // trim(net_file)
            stop 'set_composition: do not recognize net_file'
      
         end if

         sum = 0d0
         do i=1, species
            if (xin(i) < 1d-99) xin(i) = 1d-99
            sum = sum + xin(i)
         end do
         if (abs(1d0-sum) > 1d-4) write(*, *) 'change abundance sum by', 1d0-sum
         xin(adjustment_iso) = xin(adjustment_iso) + (1d0 - sum)
         if (xin(adjustment_iso) < 0d0) stop 'error in sum of abundances'

      end subroutine set_composition
      
      
      subroutine read_test_data(filename, n, rho_vec, T_vec, ierr)
         ! the data files have columns of mass, radius, density, temp
         use utils_lib
         character (len=*), intent(in) :: filename
         integer, intent(out) :: n
         real(dp), dimension(:), pointer :: rho_vec, T_vec ! to be allocated and filled
         integer, intent(out) :: ierr
         
         integer :: iounit, i
         real(dp) :: junk
         
         ierr = 0
         iounit = alloc_iounit(ierr); if (ierr /= 0) return
         open(unit=iounit, file=trim(filename), action='read', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'failed to open ', trim(filename)
            return
         end if
         
         i = 0
         do
            read(unit=iounit, fmt=*, iostat=ierr) junk, junk, junk, junk
            if (ierr == 0) then
               i = i+1; cycle
            end if
            ierr = 0
            n = i
            exit
         end do
         rewind(iounit)
         
         allocate(rho_vec(n), T_vec(n), stat=ierr); if (ierr /= 0) return
         
         do i=1, n
            read(iounit, *) junk, junk, rho_vec(i), T_vec(i)
         end do
      
         close(iounit)
      
         call free_iounit(iounit)
         
      end subroutine read_test_data

      
      subroutine Do_One_Test(net_file, do_timing)
         use chem_lib, only:composition_info
         use rates_lib
         character (len=*), intent(in) :: net_file
         logical, intent(in) :: do_timing
         call Do_One_Testcase(net_file, do_timing, .false.)
      end subroutine Do_One_Test

      
      subroutine Do_One_Test_and_show_Qs(net_file, do_timing)
         use chem_lib, only:composition_info
         use rates_lib
         character (len=*), intent(in) :: net_file
         logical, intent(in) :: do_timing
         call Do_One_Testcase(net_file, do_timing, .true.)
      end subroutine Do_One_Test_and_show_Qs
      
      
      subroutine Do_One_Testcase(net_file, do_timing, show_Qs)
         use chem_lib, only:composition_info
         use rates_lib
         character (len=*), intent(in) :: net_file
         logical, intent(in) :: do_timing, show_Qs
         
         real(dp) :: logRho, logT, Rho, T, xsum, Q1, Q2,
     >     eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, theta_e,
     >     eps_nuc_categories(num_rvs, num_categories), xh, xhe, mass_correction !approx_abar, approx_zbar
         real(dp), dimension(:, :), pointer :: reaction_eps_nuc
         real(dp), dimension(:, :), pointer :: rate_screened, rate_raw
         integer :: i, j, k, info, ierr, nreps, rep, times_total, elapsed_time, clock_rate, time0, time1
         integer :: lwork, adjustment_iso, ir_c12_c12_to_he4_ne20, ir_he4_ne20_to_c12_c12
         real(dp), dimension(:), pointer :: d_eps_nuc_dx, dabar_dx, dzbar_dx, dmc_dx
         real(dp), pointer :: work(:), rate_factors(:), category_factors(:),
     >      actual_Qs(:), actual_neuQs(:)       
         logical, pointer :: from_weaklib(:)
         
         
         include 'formats.dek'
         
         write(*,*) 'Do_One_Test ' // trim(net_file)
                  
         theta_e = theta_e_for_graboske_et_al
         call test_net_setup(net_file)
            
         ierr = 0
         
         write(*,*) 'species', species
         
         allocate(d_eps_nuc_dx(species), dabar_dx(species), dzbar_dx(species), dmc_dx(species))
         
         info = 0
         lwork = net_work_size(handle, info) 
         if (info /= 0) stop 1
         
         allocate(work(lwork), 
     >         rate_factors(num_reactions), category_factors(num_categories),
     >         rate_screened(num_rvs, num_reactions), rate_raw(num_rvs, num_reactions), 
     >         actual_Qs(num_reactions), actual_neuQs(num_reactions), from_weaklib(num_reactions),
     >         stat=info)
         if (info /= 0) stop 2
         
         rate_factors(:) = 1
         category_factors(:) = 1
         
         if (.false.) then ! get neutrino Q
         
            Q1 = eval_neutrino_Q(img22, is30)
            write(*,1) 'Qneu mg22->s30', Q1
         
            Q1 = eval_neutrino_Q(is30, ini56)
            write(*,1) 'Qneu s30->ni56', Q1
         
            Q1 = eval_neutrino_Q(ica38, ini56)
            write(*,1) 'Qneu ca38->ni56', Q1
         
            Q1 = eval_neutrino_Q(ini56, ige64)
            write(*,1) 'Qneu ni56->ge64', Q1

            Q1 = eval_neutrino_Q(ige64, ise68)
            write(*,1) 'Qneu ge64->se68', Q1

            Q1 = eval_neutrino_Q(ise68, ikr72)
            write(*,1) 'Qneu se68->kr72', Q1

            Q1 = eval_neutrino_Q(ikr72, isr76)
            write(*,1) 'Qneu kr72->sr76', Q1

            Q1 = eval_neutrino_Q(isr76, isn104)
            write(*,1) 'Qneu sr76->sn104', Q1

            stop
         end if
         
         if (.false.) then ! get reaction Q
            ! co55 -> fe55
            Q1 = isoB(ife55) - isoB(ico55)
            write(*,1) 'Q co55->fe55', Q1
            stop
         end if

         xin = 0
         eta = 0

         if (net_file == 'basic.net') then
            
            nrates_to_show = 10
         
            rates_to_show(1:nrates_to_show) = (/ 
     >         irpp_to_he3,
     >         irpep_to_he3,
     >         ir_he3_he3_to_h1_h1_he4,
     >         ir34_pp2,
     >         ir34_pp3,
     >         ir_h1_he3_wk_he4,
     >         irc12_to_n14,
     >         irn14_to_c12,
     >         irn14_to_o16,
     >         iro16_to_n14  /)

                     xin(net_iso(ih1))=     3.499985d-01
                    xin(net_iso(ihe3))=     1.3914298d-05
                    xin(net_iso(ic12))=     1.721186d-05  
                    xin(net_iso(in14))=     5.084861d-03
                    xin(net_iso(io16))=     9.502718d-03
                   xin(net_iso(ine20))=     0
                    xin(net_iso(ihe4))=     1d0 - 0.02d0 - xin(net_iso(ih1)) - xin(net_iso(ihe3))
                   xin(net_iso(img24))=     1d0 - sum(xin(:))
                              
                              write(*,*) 'sum xin', sum(xin(:))

                                  logT =    7.160481D+00
                                logRho =    2.178040D+00
                                   eta =   0
                               theta_e =    1
                                   
                                   
               screening_mode = extended_screening
               
               call net_set_logTcut(handle, 0d0, 0d0, info)
               if (info /= 0) then
                  write(*,*) 'failed in net_set_logTcut'
                  stop 1
               end if

         else if (net_file == 'pp_extras.net') then
            
            nrates_to_show = 6
         
            rates_to_show(1:nrates_to_show) = (/ 
     >      ir_h1_h1_wk_h2,
     >      ir_h2_pg_he3,
     >      ir_be7_ec_li7,
     >      ir_b8_wk_he4_he4,
     >      ir_he3_he3_to_h1_h1_he4,
     >      ir_he3_ag_be7 /)
     
                     xin(net_iso(ih1))=     7.0999999999999996D-01
                     xin(net_iso(ih2))=     2.0000000000000002D-05
                    xin(net_iso(ihe3))=     2.0000000000000002D-05
                    xin(net_iso(ihe4))=     2.7000000000000002D-01
                    xin(net_iso(ili7))=     1.0514111184895020D-08
                    xin(net_iso(ibe7))=     1.0000000000000000D-99
                     xin(net_iso(ib8))=     1.0000000000000000D-99
                    xin(net_iso(ic12))=     3.4351722028723463D-03
                    xin(net_iso(in14))=     1.0065620355000817D-03
                    xin(net_iso(io16))=     9.3438352969958897D-03
                   xin(net_iso(ine20))=     2.0956982142222580D-03
                   xin(net_iso(img24))=     4.0787217362982173D-03
                              
                              write(*,*) 'sum xin', sum(xin(:))

                                  logT =    5.6864273893515023D+00
                                logRho =    2.0591020210828619D+00
                                   eta =   -1.4317150417353590D+01
                               theta_e =    1
                                   
                                   
               screening_mode = classic_screening
               
               call net_set_logTcut(handle, 0d0, 0d0, info)
               if (info /= 0) then
                  write(*,*) 'failed in net_set_logTcut'
                  stop 1
               end if
               
               if (.false.) then
                  if (.false.) then
                     rate_factors(:) = 0
                     i = reaction_table(ir_h2_pg_he3)
                     if (i == 0) stop 1
                     rate_factors(i) = 1
                  else
                     i = reaction_table(ir_ne20_ag_mg24)
                     if (i == 0) stop 1
                     rate_factors(i) = 0
                     i = reaction_table(ir_o16_ag_ne20)
                     if (i == 0) stop 1
                     rate_factors(i) = 0
                  end if
               end if

         

         else if (net_file == 'agb.net') then
            
            nrates_to_show = 4
         
            rates_to_show(1:nrates_to_show) = (/ 
     >      ir_h1_h1_wk_h2,
     >      ir_c13_an_o16,
     >      ir_f19_ap_ne22,
     >      ir_he3_ag_be7 /)
     
                     xin(net_iso(ih1))= 1
                              
                              write(*,*) 'sum xin', sum(xin(:))

                                  logT =    8.6864273893515023D+00
                                logRho =    2.0591020210828619D+00
                                   eta =   -1.4317150417353590D+01
                               theta_e =    1
                                   
                                   
               screening_mode = classic_screening
               
               call net_set_logTcut(handle, 0d0, 0d0, info)
               if (info /= 0) then
                  write(*,*) 'failed in net_set_logTcut'
                  stop 1
               end if
               
               if (.false.) then
                  if (.false.) then
                     rate_factors(:) = 0
                     i = reaction_table(ir_h2_pg_he3)
                     if (i == 0) stop 1
                     rate_factors(i) = 1
                  else
                     i = reaction_table(ir_ne20_ag_mg24)
                     if (i == 0) stop 1
                     rate_factors(i) = 0
                     i = reaction_table(ir_o16_ag_ne20)
                     if (i == 0) stop 1
                     rate_factors(i) = 0
                  end if
               end if

         

         else if (net_file == 'pp_and_cno_extras.net') then
            
            nrates_to_show = 8
         
            rates_to_show(1:nrates_to_show) = (/ 
     >      rates_reaction_id('r_n13_wk_c13'),               
     >      rates_reaction_id('r_o15_wk_n15'),               
     >      rates_reaction_id('r_f17_wk_o17'),               
     >      rates_reaction_id('r_f18_wk_o18'),               
     >      rates_reaction_id('r_o14_wk_n14'),               
     >      rates_reaction_id('r_ne18_wk_f18'),               
     >      rates_reaction_id('r_ne19_wk_f19'),               
     >      ir_he4_he4_he4_to_c12 /)
     
         xin = 0
                     xin(net_iso(ih1))=     7.2265805432969643D-01
                    xin(net_iso(ihe3))=     6.7801726921522655D-04
                    xin(net_iso(ihe4))=     2.6667042876319019D-01
                    xin(net_iso(ic12))=     1.9056849943017622D-03
                    xin(net_iso(in13))=     1.4791107757148081D-04
                    xin(net_iso(in14))=     6.0253770632534619D-04
                    xin(net_iso(in15))=     1.9263919190065488D-07
                    xin(net_iso(io14))=     9.5977897394247582D-09
                    xin(net_iso(io15))=     7.6666060610240002D-08
                    xin(net_iso(io16))=     5.5886173952684358D-03
                    xin(net_iso(io17))=     2.0449560316023342D-06
                    xin(net_iso(io18))=     1.1313355448541673D-05
                    xin(net_iso(if19))=     2.1343308067891499D-07
                   xin(net_iso(ine20))=     1.2563102679570999D-03
                   xin(net_iso(img24))=     4.7858754879924638D-04
 
                                     T =    8.5648111120065376D+06
                                  logT =    6.9327177894944265D+00
                                   rho =    7.8571498592117219D+00
                                logRho =    8.9526503651107059D-01
                                  abar =    1.2655060647252907D+00
                                  zbar =    1.0901664301076275D+00
                                 z2bar =    1.3036906023574921D+00
                                    ye =    8.6144702146826535D-01
                                   eta =   -3.4387570967781595D+00                                   
                                   
               screening_mode = extended_screening

         else if (net_file == 'approx29.net') then

            nrates_to_show = 11
         
            rates_to_show(1:nrates_to_show) = (/ 
     >      ir1212,           
     >      irc12_to_n14,
     >      ir_c12_ag_o16,
     >      ir1216_to_mg24,           
     >      ir1216_to_si28,           
     >      ir1616a,
     >      ir1616ppa,
     >      ir1616ppg,
     >      ir_o16_ag_ne20,
     >      ir_ne20_ag_mg24,
     >      iro16_to_n14
     >      /)
     
            
                                   xin(net_iso(io16))=     8.4078806793430971D-01
                                   xin(net_iso(ic12))=     1.5667831148778857D-01
                                  xin(net_iso(ine20))=     2.5324828950044945D-03
                                  xin(net_iso(img24))=     1.1374886503238785D-06
                                  xin(net_iso(isi28))=     1.9320589618754134D-10
                                   xin(net_iso(ihe4))=     1.0410585725431680D-12
                                   xin(net_iso(is32))=     9.7508528557715739D-18
                                  xin(net_iso(iar36))=     2.3719635461470835D-26
                                  xin(net_iso(ica40))=     8.9365928208049725D-37
                                  xin(net_iso(iti44))=     5.9447148967850597D-49
                                  xin(net_iso(ini58))=     3.4324388103659015D-56
                                  xin(net_iso(ini60))=     1.3592016146385322D-56
                                  xin(net_iso(ife58))=     2.5667097198284398D-57
                                  xin(net_iso(ini62))=     1.9271671069951470D-57
                                  xin(net_iso(ini64))=     5.0462911310261966D-58
                                  xin(net_iso(iti50))=     1.1414849804824868D-58
                                  xin(net_iso(icr48))=     4.2643755024382846D-63
                                  xin(net_iso(icr54))=     3.4371501049170877D-72
                                  xin(net_iso(ife52))=     1.2270370811284163D-78
                                  xin(net_iso(ife56))=     1.0673812072272265D-97
                                  xin(net_iso(iprot))=     5.2120663852112486D-99
                                   xin(net_iso(in14))=     1.5987902769634659D-99
                                  xin(net_iso(ini56))=     1.1993659789578421D-99
                                  xin(net_iso(ife54))=     1.0831894305010897D-99
                                    xin(net_iso(ih1))=     1.0227678792105984D-99
                                  xin(net_iso(icr56))=     1.0000006088916377D-99
                                  xin(net_iso(ife60))=     1.0000000000000006D-99
                                   xin(net_iso(ihe3))=     1.0000000000000006D-99
                                  xin(net_iso(ineut))=     0.0000000000000000D+00
             
            xin(species) = 1d0 - sum(xin(1:species-1))
                              
                              write(*,*) 'sum xin', sum(xin(:))

                                                    T =    6.8918604410301852D+08
                                                 logT =    8.8383364744776500D+00
                                                  rho =    2.4288677461994484D+04
                                               logRho =    4.3854038677653975D+00
                                                 abar =    1.5213185776121293D+01
                                                 zbar =    7.6065928880606464D+00
                                                z2bar =    5.8507728595911630D+01
                                                   ye =    5.0000000000000000D-01
                                                  eta =   -2.6646884260237336D+00
                        screening_mode = extended_screening
                                              theta_e =    0.0000000000000000D+00
               
         
         else if (net_file == 'wd_o_ne_ignite.net') then


            nrates_to_show = 6
         
         
            rates_to_show(1:nrates_to_show) = (/ 
     >         rates_reaction_id('r_mg24_np_na24'),               
     >         rates_reaction_id('r_mg24_wk_na24'),               
     >         rates_reaction_id('rmg24ap_to_si28'),               
     >         rates_reaction_id('r_mg24_ag_si28'),               
     >         rates_reaction_id('r_mg24_ga_ne20'),               
     >         rates_reaction_id('r_mg24_gp_na23') /)
     
                   xin(net_iso(ineut))=     0.0000000000000000D+00
                     xin(net_iso(ih1))=     0.0000000000000000D+00
                    xin(net_iso(ihe3))=     0.0000000000000000D+00
                    xin(net_iso(ihe4))=     1.5957030460155789D-15
                    xin(net_iso(ic12))=     6.5187788875984993D-03
                    xin(net_iso(in14))=     1.2243249970731128D-08
                    xin(net_iso(io16))=     5.6259113770759883D-01
                    xin(net_iso(io20))=     2.0328589644116813D-10
                    xin(net_iso(if20))=     9.3235184256934400D-12
                   xin(net_iso(ine20))=     3.9890879552663960D-01
                   xin(net_iso(ine24))=     5.4471392360365520D-04
                   xin(net_iso(ina23))=     7.7578429229253277D-05
                   xin(net_iso(ina24))=     2.9424728593175371D-10
                   xin(net_iso(img24))=     2.8608638006352681D-02
                   xin(net_iso(isi28))=     8.9096274199580966D-04
                    xin(net_iso(is32))=     4.1509540873424419D-04
                   xin(net_iso(iar36))=     8.1172438275055043D-05
                   xin(net_iso(ica40))=     6.2815722448016586D-05
                   xin(net_iso(iti44))=     1.0479767011423344D-15
                   xin(net_iso(icr48))=     0.0000000000000000D+00
                   xin(net_iso(icr56))=     0.0000000000000000D+00
                   xin(net_iso(ife52))=     0.0000000000000000D+00
                   xin(net_iso(ife54))=     7.4774179803179899D-05
                   xin(net_iso(ife56))=     1.2255242776117584D-03
                   xin(net_iso(ini56))=     0.0000000000000000D+00
                   xin(net_iso(iprot))=     0.0000000000000000D+00

                                     T =    1.6602302786935723D+08
                                  logT =    8.2201683301062278D+00
                                   rho =    1.1362219397057323D+10
                                logRho =    1.0055463170965831D+01
                                  abar =    1.7562161143883930D+01
                                  zbar =    8.7794607650170917D+00
                                 z2bar =    7.8422787604638103D+01
                                    ye =    4.9990776722115221D-01
                                   eta =    6.0822366730380656D+02
                               theta_e =    0.0000000000000000D+00
                        screening_mode = extended_screening
         
         else if (net_file == 'approx21.net') then

            ! TESTING
            nrates_to_show = 2
         
            rates_to_show(1:nrates_to_show) = (/ 
     >         irpp_to_he3,         ! p(p,e+nu)h2(p,g)he3       
     >         irpep_to_he3        ! p(e-p,nu)h2(p,g)he3     
     >      /)
     
                                  xin(net_iso(isi28))=     5.7736594568628330D-01
                                   xin(net_iso(is32))=     3.3161413041024329D-01
                                  xin(net_iso(iar36))=     4.5440174228012949D-02
                                  xin(net_iso(ica40))=     4.4326435432300693D-02
                                  xin(net_iso(ife56))=     1.2343187771546734D-03
                                   xin(net_iso(io16))=     1.0269660544212729D-05
                                  xin(net_iso(icr48))=     5.5479714118961378D-06
                                  xin(net_iso(iti44))=     2.6150279013440417D-06
                                  xin(net_iso(img24))=     4.9126585387509215D-07
                                  xin(net_iso(ineut))=     6.7235812080405344D-08
                                   xin(net_iso(ic12))=     4.2308061053742621D-09
                                  xin(net_iso(ine20))=     5.3544403631109066D-11
                                  xin(net_iso(iprot))=     1.9172559653136226D-11
                                   xin(net_iso(ihe4))=     9.5856466256448351D-13
                                  xin(net_iso(ini56))=     1.3263075666515746D-16
                                  xin(net_iso(icr56))=     1.0420775439667622D-18
                                  xin(net_iso(ife52))=     1.0863095786231788D-20
                                  xin(net_iso(ife54))=     8.3327063420412484D-22
                                   xin(net_iso(ihe3))=     3.8463147997752578D-46
                                   xin(net_iso(in14))=     1.0000000000000000D-99
                                    xin(net_iso(ih1))=     1.0000000000000000D-99

                                                    T =    2.2254570420573139D+09
                                                 logT =    9.3474192155237201D+00
                                                  rho =    3.7761775506596074D+07
                                               logRho =    7.5770524060335811D+00
                                                 abar =    2.9961210816200811D+01
                                                 zbar =    1.4979283626686712D+01
                                                z2bar =    2.2655888497684035D+02
                                                   ye =    4.9995588357821041D-01
                                                  eta =    4.5094610816657257D+00
                                             screening_mode = extended_screening
                                              theta_e =    0.0000000000000000D+00
                               

         else if (net_file == 'rp_si26.net') then

            xin = 0
            xin(net_iso(ih1))   =   1
 
 
                                     T =    9.0d8
                                  logT =    log10(T)
                                   rho =    4.5d5
                                logRho =    log10(rho)
                                  abar =    4.6439541133960681D+01
                                  zbar =    2.0327386428253433D+01
                                 z2bar =    5.0071386690625775D+02
                                    ye =    0.5
                                   eta =    6.3910099340916329D+00
                                   screening_mode = extended_screening
                               theta_e =    0.0000000000000000D+00

         else
            
            write(*, *) 'need to define setup for net_file ', trim(net_file)
            stop 'Do_One_Test'
         
         end if
         
         Rho = 10**logRho
         T = 10**logT
         
         write(*, *)
         write(*, *)
         
         info = 0
         allocate(rate_screened(num_rvs, num_reactions), reaction_eps_nuc(num_rvs, num_reactions),     
     >         rate_raw(num_rvs, num_reactions), stat=info)
         if (info /= 0) stop 2
         
         ierr = 0

         call composition_info(
     >      species, chem_id, xin, xh, xhe, abar, zbar, z2bar, ye, 
     >      mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
         

        do i = 1, species
           write(*,'(a40,i6,e26.16)')  'x ' // trim(chem_isos% name(chem_id(i))), i, xin(i)
        end do
        write(*,*)

         if (do_timing) then
            nreps = 10000
            call zero_net_timing(g)
            g% doing_timing = .true.
            call system_clock(time0)
         else
            nreps = 1
         end if

         do rep=1,nreps
            call net_get_with_Qs(handle, species, num_reactions, 
     >            xin, T, logT, Rho, logRho, 
     >            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho,
     >            rate_factors, category_factors, 
     >            std_reaction_Qs, std_reaction_neuQs,
     >            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, 
     >            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, 
     >            screening_mode, theta_e_for_graboske_et_al,     
     >            rate_screened, rate_raw, 
     >            reaction_eps_nuc, eps_nuc_categories, eps_neu_total, 
     >            lwork, work, actual_Qs, actual_neuQs, from_weaklib, info)
            if (info /= 0) then
               write(*,1) 'logT', logT
               write(*,1) 'logRho', logRho
               write(*, *) 'bad return from net_get'
               stop 1
            end if
         end do

         if (do_timing) then
            call system_clock(time1)
            elapsed_time = time1 - time0
            times_total = get_net_timing_total(g)
            if (times_total <= 0) then
               write(*,*) 'must set g% doing_timing = .true.'
            end if
            write(*,*)
            write(*,*) 'percentages of total net eval time'
            write(*,*)
            write(*,'(a30,f14.3)') 'clock_net_eval%', dble(100*g% clock_net_eval)/dble(elapsed_time)
            write(*,'(a30,f14.3)') 'clock_net_weak_rates%', dble(100*g% clock_net_weak_rates)/dble(elapsed_time)
            write(*,'(a30,f14.3)') 'clock_net_rate_tables%', dble(100*g% clock_net_rate_tables)/dble(elapsed_time)
            write(*,'(a30,f14.3)') 'clock_net_screen%', dble(100*g% clock_net_screen)/dble(elapsed_time)
            write(*,'(a30,f14.3)') 'clock_net_derivs%', dble(100*g% clock_net_derivs)/dble(elapsed_time)
            write(*,'(a30,f14.3)') 'other%', dble(100*(elapsed_time - times_total))/dble(elapsed_time)
            call system_clock(time0, clock_rate)
            write(*,*)
            write(*,'(a30,f14.3)') 'times_total', dble(times_total)/dble(clock_rate)
            write(*,'(a30,f14.3)') 'elapsed time', dble(elapsed_time)/dble(clock_rate)
            write(*,*)
            write(*,*)
            if (g% doing_derivs_timing) then
               write(*,*)
               write(*,'(a30,2f14.3)') 'clock_derivs_setup', 
     >            dble(g% clock_derivs_setup)/dble(clock_rate),
     >            100*dble(g% clock_derivs_setup)/dble(g% clock_net_derivs)
               write(*,'(a30,2f14.3)') 'clock_derivs_select', 
     >            dble(g% clock_derivs_select)/dble(clock_rate),
     >            100*dble(g% clock_derivs_select)/dble(g% clock_net_derivs)
               write(*,'(a30,2f14.3)') 'clock_derivs_general', 
     >            dble(g% clock_derivs_general)/dble(clock_rate),
     >            100*dble(g% clock_derivs_general)/dble(g% clock_net_derivs)
               times_total = g% clock_derivs_setup + g% clock_derivs_select + g% clock_derivs_general
               write(*,'(a30,2f14.3)') 'other derivs time', 
     >            dble(g% clock_net_derivs - times_total)/dble(clock_rate),
     >            100*dble(g% clock_net_derivs - times_total)/dble(g% clock_net_derivs)
               write(*,*)
               write(*,*)
               write(*,*)
            end if
            return
         end if

         if (show_Qs) then
            write(*,*)
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            write(*,*)
            write(*,'(30x,4a20)') 'Q total', 'Q neutrino', 'Q total-neutrino'
            do i = 1, num_reactions
               if (from_weaklib(i)) then
                  write(*,'(a30,99f20.10)') 'weaklib ' // trim(reaction_Name(reaction_id(i))), 
     >               actual_Qs(i), actual_neuQs(i), actual_Qs(i) - actual_neuQs(i)
               else
                  write(*,'(a30,99f20.10)') trim(reaction_Name(reaction_id(i))), 
     >               actual_Qs(i), actual_neuQs(i), actual_Qs(i) - actual_neuQs(i)
               end if
            end do
            write(*,*)
            stop
         end if
         
         
         call dealloc
         
         
         write(*,2) 'screening_mode', screening_mode
            
        if (.true.) then
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            if (.true.) then
              write(*,*) 'reaction_eps_nuc'
               do i=1,nrates_to_show
                  j = rates_to_show(i)
                  if (j == 0) cycle
                  write(*,1) 'eps_nuc ' // trim(reaction_Name(j)), reaction_eps_nuc(i_rate,reaction_table(j))
               end do
               write(*,*)
               !stop
            end if
            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'rate_raw ' // trim(reaction_Name(j)), 
     >                           rate_raw(i_rate,reaction_table(j))
            end do
            write(*,*)
            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'd_rate_raw_dT ' // trim(reaction_Name(j)), 
     >                           rate_raw(i_rate_dT,reaction_table(j))
            end do
            write(*,*)
            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'rate_screened ' // trim(reaction_Name(j)),
     >                           rate_screened(i_rate,reaction_table(j))
            end do
            write(*,*)
            write(*,*)
           do i = 1, species
              write(*,1)  'x ' // trim(chem_isos% name(chem_id(i))), xin(i)
           end do
           write(*,*)
           do i = 1, species
              write(*,1)  'dxdt ' // trim(chem_isos% name(chem_id(i))), dxdt(i)
           end do
           write(*,*)
           do i = 1, species
              if (-dxdt(i) > 1d-90) 
     >           write(*,1)  'x/dxdt ' // trim(chem_isos% name(chem_id(i))), xin(i)/dxdt(i)
           end do
           write(*,*)
           do i = 1, species
              write(*,1)  'd_dxdt_dlnRho ' // trim(chem_isos% name(chem_id(i))), d_dxdt_dRho(i)*rho
           end do
           write(*,*)
           do i = 1, species
              write(*,1)  'd_dxdt_dlnT ' // trim(chem_isos% name(chem_id(i))), d_dxdt_dT(i)*T
           end do
           write(*,*)
           do i = 1, species
              write(*,1)  'd_dxdt_dx(1,:) ' // trim(chem_isos% name(chem_id(i))), d_dxdt_dx(1,i)
           end do
           write(*,*)
           do i = 1, num_categories
              if (abs(eps_nuc_categories(i_rate,i)) < 1d-20) cycle
              write(*,1)  'eps_nuc_cat ' // trim(category_name(i)), eps_nuc_categories(i_rate,i)
           end do
           write(*,*)
           do i = 1, species
              write(*,1)  'd_eps_nuc_dx ' // trim(chem_isos% name(chem_id(i))), d_eps_nuc_dx(i)
           end do
           write(*,*)
           write(*,1) 'eps_neu_total', eps_neu_total
           write(*,1) 'eps_nuc', eps_nuc
           write(*,1) 'd_epsnuc_dlnd', d_eps_nuc_dRho*Rho
           write(*,1) 'd_epsnuc_dlnT', d_eps_nuc_dT*T
           stop
        end if

         write(*,*)
         write(*,*)
         write(*,*) 'net_name ', trim(net_file)
         write(*,*) 'species', species
         write(*,1) 'theta_e =', theta_e_for_graboske_et_al
         write(*,1) 'abar =', abar
         write(*,1) 'zbar =', zbar
         write(*,1) 'z2bar =', z2bar
         write(*,1) 'ye =', ye
         write(*, *)
         do i=1,nrates_to_show
            j = rates_to_show(i)
            if (j == 0) cycle
            if (reaction_table(j) == 0) then
               write(*,*) 'missing reaction_table(j) for ' // trim(reaction_Name(j))
               stop
            end if
            write(*,1) trim(reaction_Name(j)), reaction_eps_nuc(i_rate,reaction_table(j))
         end do
         write(*,*)
         write(*,1) 'eps_nuc', eps_nuc
         write(*,1) 'd_eps_nuc_dRho', d_eps_nuc_dRho
         write(*,1) 'd_eps_nuc_dT', d_eps_nuc_dT
         do j=1,species
            write(*,1) 'd_eps_nuc_dx ' // trim(chem_isos% name(chem_id(j))), d_eps_nuc_dx(j)
         end do
         write(*,*)
         do j=1,species
            write(*,1) 'dxdt ' // trim(chem_isos% name(chem_id(j))), dxdt(j)
         end do
         write(*,*)
         do j=1,species
            write(*,1) 'd_dxdt_dRho ' // trim(chem_isos% name(chem_id(j))), d_dxdt_dRho(j)
         end do
         write(*,*)
         do j=1,species
            write(*,1) 'd_dxdt_dT ' // trim(chem_isos% name(chem_id(j))), d_dxdt_dT(j)
         end do
         write(*,*)
         do j=1,species
            write(*,1) 'd_dxdt_dx(1,:) ' // trim(chem_isos% name(chem_id(j))), d_dxdt_dx(1,j)
         end do
         write(*, *)
         do j=1,num_categories
            write(*,1) trim(category_name(j)), eps_nuc_categories(i_rate, j)
         end do
         write(*,*)
         write(*,1) 'eta =', eta
         write(*,1) 'logT =', logT
         write(*,1) 'logRho =', logRho
         write(*,*) 'screening_mode =', screening_mode
         write(*,*)
         do j=1,species
            write(*,1) 'xin(net_iso(i' // trim(chem_isos% name(chem_id(j))) // '))= ', xin(j)
         end do
         write(*,*)
         write(*,1) 'sum(xin(1:species))', sum(xin(1:species))
         write(*,1) '1 - sum(xin(1:species))', 1 - sum(xin(1:species))
         write(*,*)
         write(*,1) 'eps_nuc', eps_nuc
         write(*,1) 'eps_nuc_neu_total', eps_neu_total




         deallocate(actual_Qs, actual_neuQs, from_weaklib, rate_screened, rate_raw, reaction_eps_nuc)
         
         
         contains
         
         subroutine dealloc
            deallocate(work, rate_factors, category_factors)
         end subroutine dealloc


      end subroutine Do_One_Testcase

      
      subroutine test_neutrino_Q
         real(dp), parameter :: Qnu_n13 = 0.714440d0 !..13n(e+nu)13c
         real(dp), parameter :: Qnu_o15 = 1.005513d0 !..15o(e+nu)15n
         real(dp), parameter :: Qnu_f17 = 1.009145d0 !..17f(e+nu)17o
         real(dp), parameter :: Qnu_f18 = 0.393075d0 !..18f(e+nu)18o   
         real(dp), parameter :: Qnu_o14 = 2.22d0 !..14o(e+nu)14n
         real(dp), parameter :: Qnu_ne18 = 1.87d0 !..18ne(e+nu)18f
         real(dp), parameter :: Qnu_ne19 = 1.25d0 !..19ne(e+nu)19f
         !real(dp), parameter :: Qnu_mg21 = 6.2d0 !..mg21(e+nu)na21
         real(dp), parameter :: Qnu_mg22 = 2.1d0 !..mg22(e+nu)na22
         
 1       format(a40, 1pe26.16)
         
         write(*, 1) 'expected Q for 13n(e+nu)13c', Qnu_n13
         write(*, 1) 'calculated Q for 13n(e+nu)13c', eval_neutrino_Q(in13, ic13)
         write(*, *)
         write(*, 1) 'expected Q for 15o(e+nu)15n', Qnu_o15
         write(*, 1) 'calculated Q for 15o(e+nu)15n', eval_neutrino_Q(io15, in15)
         write(*, *)
         write(*, 1) 'expected Q for 17f(e+nu)17o', Qnu_f17
         write(*, 1) 'calculated Q for 17f(e+nu)17o', eval_neutrino_Q(if17, io17)
         write(*, *)
         write(*, 1) 'expected Q for 18f(e+nu)18o', Qnu_f18
         write(*, 1) 'calculated Q for 18f(e+nu)18o', eval_neutrino_Q(if18, io18)
         write(*, *)
         write(*, 1) 'expected Q for 14o(e+nu)14n', Qnu_o14
         write(*, 1) 'calculated Q for 14o(e+nu)14n', eval_neutrino_Q(io14, in14)
         write(*, *)
         write(*, 1) 'expected Q for 18ne(e+nu)18f', Qnu_ne18
         write(*, 1) 'calculated Q for 18ne(e+nu)18f', eval_neutrino_Q(ine18, if18)
         write(*, *)
         write(*, 1) 'expected Q for 19ne(e+nu)19f', Qnu_ne19
         write(*, 1) 'calculated Q for 19ne(e+nu)19f', eval_neutrino_Q(ine19, if19)
         write(*, *)
         write(*, 1) 'expected Q for mg22(e+nu)na22', Qnu_mg22
         write(*, 1) 'calculated Q for mg22(e+nu)na22', eval_neutrino_Q(img22, ina22)
         write(*, *)
         stop
      
      end subroutine test_neutrino_Q
      

      
      end module test_net_support




