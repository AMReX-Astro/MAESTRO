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
      use crlibm_lib
      use net_def
      use net_lib
      use const_def
      use rates_def
      use test_net_do_one
      
      implicit none
      
      integer, parameter :: max_files = 20, max_cnt = 100000      
      
      
      character (len=256) :: eos_file_prefix, cache_suffix
      
      type (Net_General_Info), pointer  :: g
      integer :: num_reactions
      
      integer, dimension(:), pointer :: net_iso, chem_id, isos_to_show

      integer :: which_rates_choice
      
      integer, pointer :: reaction_table(:)
      integer, pointer :: rates_to_show(:)
      integer, pointer :: which_rates(:)


      real(dp), dimension(:), pointer :: rho_vector, T_vector

      integer :: nrates_to_show, nisos_to_show
      


      contains
      

      subroutine do_test_net(do_plots, symbolic)
         logical, intent(in) :: do_plots, symbolic
         call set_composition(g, species, xin)
         eta = 0
         if (do_plots) then
            call Create_Plot_Files(species, xin, g)
         else
            call Do_One_Net(symbolic)
         end if
      end subroutine do_test_net 
      
      
      subroutine test_net_setup(net_file_in)
         character (len=*), intent(in) :: net_file_in
         integer, pointer :: r_id(:)

         integer :: info, i, ierr
         
         include 'formats'
         
         net_file = net_file_in

      	allocate(net_iso(num_chem_isos), isos_to_show(num_chem_isos), chem_id(num_chem_isos))
         
         call net_init(ierr)
         if (ierr /= 0) stop 1
         
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
         
         species = g% num_isos
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
         
         call do_test_net_alloc(species)
#ifdef offload
         r_id => reaction_id
         !dir$ offload target(mic) in(r_id)
         call do_copy_reaction_id_to_coprocessor(r_id)
         !dir$ offload target(mic) in(species)
         call do_test_net_alloc(species)
#endif

      end subroutine test_net_setup
      
#ifdef offload
      !dir$ attributes offload: mic :: do_copy_reaction_id_to_coprocessor
      subroutine do_copy_reaction_id_to_coprocessor(r_id)
         integer, pointer, intent(in) :: r_id(:)
         integer :: i, sz
         sz = size(r_id,dim=1)
         if (associated(reaction_id)) deallocate(reaction_id)
         allocate(reaction_id(sz))
         do i=1,sz
            reaction_id(i) = r_id(i)
         end do
      end subroutine do_copy_reaction_id_to_coprocessor
#endif         
      
#ifdef offload
      !dir$ attributes offload: mic :: do_test_net_alloc
#endif         
      subroutine do_test_net_alloc(species)
         integer, intent(in) :: species
         allocate( &
            xin(species), xin_copy(species), d_eps_nuc_dx(species), dxdt(species), &
            d_dxdt_dRho(species), d_dxdt_dT(species), d_dxdt_dx(species, species))
      end subroutine do_test_net_alloc


      subroutine Do_One_Net(symbolic)
         logical, intent(in) :: symbolic
         integer :: i, id
         if (.false.) then
            do i = 1, g% num_reactions
               id = g% reaction_id(i)
               if (id > 0) write(*,'(a,2i6)') trim(reaction_Name(id)), i, id
            end do
         end if
#ifdef offload
         call setup_coprocessor
         !dir$ offload target(mic) in(symbolic)
#endif         
         call do1_net(symbolic)
      end subroutine Do_One_Net
      
#ifdef offload
      subroutine setup_coprocessor ! runs on host
         logical :: qt_in
         character (len=64) :: net_file_in
         integer :: screening_mode_in
         real(dp) :: theta_e_in
         real(dp) :: test_logT_in, test_logRho_in
         integer :: handle_in
         integer :: species_in
         real(dp) :: eta_in, d_eta_dlnT_in, d_eta_dlnRho_in
         logical :: reuse_rate_raw_in, reuse_rate_screened_in
         real(dp), pointer :: xin_in(:)
         qt_in = qt
         net_file_in = net_file
         screening_mode_in = screening_mode
         theta_e_in = theta_e
         test_logT_in = test_logT
         test_logRho_in = test_logRho
         handle_in = handle
         species_in = species
         eta_in = eta
         d_eta_dlnT_in = d_eta_dlnT
         d_eta_dlnRho_in = d_eta_dlnRho
         reuse_rate_raw_in = reuse_rate_raw
         reuse_rate_screened_in = reuse_rate_screened
         xin_in => xin

         !dir$ offload target(mic) in( &
            qt_in, &
            net_file_in, &
            screening_mode_in, &
            theta_e_in, &
            test_logT_in, &
            test_logRho_in, &
            handle_in, &
            species_in, &
            eta_in, &
            d_eta_dlnT_in, &
            d_eta_dlnRho_in, &
            reuse_rate_raw_in, &
            reuse_rate_screened_in, &
            xin_in)
         call do_setup_coprocessor( &
            qt_in, &
            net_file_in, &
            screening_mode_in, &
            theta_e_in, &
            test_logT_in, &
            test_logRho_in, &
            handle_in, &
            species_in, &
            eta_in, &
            d_eta_dlnT_in, &
            d_eta_dlnRho_in, &
            reuse_rate_raw_in, &
            reuse_rate_screened_in, &
            xin_in)
         
      end subroutine setup_coprocessor
         
      !dir$ attributes offload: mic :: do_setup_coprocessor
      subroutine do_setup_coprocessor( & ! runs on mic
            qt_in, &
            net_file_in, &
            screening_mode_in, &
            theta_e_in, &
            test_logT_in, &
            test_logRho_in, &
            handle_in, &
            species_in, &
            eta_in, &
            d_eta_dlnT_in, &
            d_eta_dlnRho_in, &
            reuse_rate_raw_in, &
            reuse_rate_screened_in, &
            xin_in)
         logical, intent(in) :: qt_in
         character (len=64), intent(in) :: net_file_in
         integer, intent(in) :: screening_mode_in
         real(dp), intent(in) :: theta_e_in
         real(dp), intent(in) :: test_logT_in, test_logRho_in
         integer, intent(in) :: handle_in
         integer, intent(in) :: species_in
         real(dp), intent(in) :: eta_in, d_eta_dlnT_in, d_eta_dlnRho_in
         logical, intent(in) :: reuse_rate_raw_in, reuse_rate_screened_in         
         real(dp), pointer, intent(in) :: xin_in(:)
         integer :: i      
         n => net_info_target
         qt = qt_in
         net_file = net_file_in
         screening_mode = screening_mode_in
         theta_e = theta_e_in
         test_logT = test_logT_in
         test_logRho = test_logRho_in
         handle = handle_in
         species = species_in
         eta = eta_in
         d_eta_dlnT = d_eta_dlnT_in
         d_eta_dlnRho = d_eta_dlnRho_in
         reuse_rate_raw = reuse_rate_raw_in
         reuse_rate_screened = reuse_rate_screened_in
         do i=1,species
            xin(i) = xin_in(i)
         end do
      end subroutine do_setup_coprocessor
#endif         

      subroutine Setup_eos(handle)
 ! allocate and load the eos tables
         use eos_def
         use eos_lib
         integer, intent(out) :: handle
         integer :: ierr
         logical, parameter :: use_cache = .true.
         eos_file_prefix = 'mesa'
         call eos_init(eos_file_prefix, '', '', '', use_cache, ierr)
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
         call do_test_net_cleanup
#ifdef offload
         !dir$ offload target(mic)
         call do_test_net_cleanup
#endif
         call free_net_handle(handle)
      end subroutine test_net_cleanup
      
#ifdef offload
      !dir$ attributes offload: mic :: do_test_net_cleanup
#endif         
      subroutine do_test_net_cleanup
         deallocate(xin)
         deallocate(d_eps_nuc_dx)
         deallocate(dxdt)
         deallocate(d_dxdt_dRho)
         deallocate(d_dxdt_dT)
         deallocate(d_dxdt_dx)
      end subroutine do_test_net_cleanup
      
      
      subroutine change_net(new_net_file)
         character (len=*), intent(in) :: new_net_file
         call test_net_cleanup
         call test_net_setup(new_net_file)
      end subroutine change_net


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
         write(io_params, '(6f16.6)')  &
               xin(net_iso(ih1)), xin(net_iso(ihe4)), xin(net_iso(ic12)),   &
               xin(net_iso(in14)), xin(net_iso(io16))
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
            T = exp10_cr(logT)

            call do_inner_loop(species, num_reactions, logT, T, j, output_values, g, xin,  &
                     logRho_points, logRho_min, dlogRho, lwork)
            
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


      subroutine do_inner_loop(species, num_reactions, logT, T, j, output_values, g, xin,  &
               logRho_points, logRho_min, dlogRho, lwork)
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
            Rho = exp10_cr(logRho)
            call do_one_net_eval(species, num_reactions, logT, T, logRho, Rho,  &
                        i, j, output_values, g, xin, lwork)
         enddo
         
      end subroutine do_inner_loop
      
      
      subroutine do_one_net_eval(species, num_reactions, logT, T, logRho, Rho,  &
               i, j, output_values, g, xin, lwork)
         use chem_lib, only:composition_info
         integer, intent(in) :: species, num_reactions, lwork
         real(dp), intent(in) :: logT, T, logRho, Rho
         integer, intent(in) :: i, j
         type (Net_General_Info), pointer :: g
         real(dp), intent(OUT) :: output_values(:, :, :)
         real(dp), intent(in) :: xin(species)
      
         real(dp) :: abar, zbar, z2bar, ye, sum, mx, weak_rate_factor
         real(dp), target :: rate_factors_a(num_reactions)
         real(dp), pointer :: rate_factors(:)

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
         real(dp) :: eps_nuc_categories(num_categories)  

         integer :: info, k, h1, he4, chem_id(species)
         real(dp) :: xh, xhe, mass_correction
         real(dp), dimension(species) :: dabar_dx, dzbar_dx, dmc_dx
         real(dp), pointer :: work(:)
         real(dp), target :: work_ary(lwork)
         logical :: skip_jacobian
         
         work => work_ary
         rate_factors => rate_factors_a
         
         h1 = net_iso(ih1)
         he4 = net_iso(ihe4)
      
         call get_chem_id_table(handle, species, chem_id, info)
         if (info /= 0) stop 3

         call composition_info( &
            species, chem_id, xin, xh, xhe, abar, zbar, z2bar, ye,  &
            mass_correction, sum, dabar_dx, dzbar_dx, dmc_dx)
     
         rate_factors(:) = 1
         weak_rate_factor = 1
         skip_jacobian = .false.
         
         call net_get(handle, skip_jacobian, n, species, num_reactions,  &
                  xin, T, logT, Rho, logRho,  &
                  abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
                  rate_factors, weak_rate_factor, &
                  std_reaction_Qs, std_reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
                  screening_mode, theta_e, &
                  eps_nuc_categories, eps_neu_total, &
                  lwork, work, info)
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

         k =   1; output_values(i, j, k) = safe_log10_cr(eps_nuc)
         
         k = k+1; output_values(i, j, k) =  sum / max(1d-20, mx)

         k = k+1; output_values(i, j, k) =  d_eps_nuc_dT * T / max(1d-20, eps_nuc)
         k = k+1; output_values(i, j, k) =  d_eps_nuc_dRho * Rho / max(1d-20, eps_nuc)
         k = k+1; output_values(i, j, k) =  safe_log10_cr(abs(d_eps_nuc_dx(h1)))
         k = k+1; output_values(i, j, k) =  safe_log10_cr(abs(d_eps_nuc_dx(he4)))
      
         if (k > max_files) then
            write(*, *) 'need to enlarge max_files'
            stop 1
         end if

      end subroutine do_one_net_eval
      
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


      subroutine set_composition(g, species, xin)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: species
         real(dp), intent(OUT) :: xin(species)
   
         real(dp) :: sum
         integer :: i, adjustment_iso
         
         eta = 0
         
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
         
         real(dp) :: logRho, logT, Rho, T, xsum, Q1, Q2, &
           eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, weak_rate_factor, &
           eps_nuc_categories(num_categories), xh, xhe, mass_correction !approx_abar, approx_zbar
         integer :: i, j, k, info, ierr, nreps, rep, times_total, elapsed_time, clock_rate, time0, time1
         integer :: lwork, adjustment_iso, ir_c12_c12_to_he4_ne20, ir_he4_ne20_to_c12_c12
         real(dp), dimension(:), pointer :: d_eps_nuc_dx, dabar_dx, dzbar_dx, dmc_dx
         real(dp), pointer :: work(:), rate_factors(:), &
            actual_Qs(:), actual_neuQs(:)       
         logical, pointer :: from_weaklib(:)
         logical :: skip_jacobian
         real(dp), dimension(:), pointer :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho         
         
         include 'formats.dek'
         
         write(*,*) 'Do_One_Test ' // trim(net_file)
                  
         call test_net_setup(net_file)
            
         ierr = 0
         
         write(*,*) 'species', species
         
         allocate(d_eps_nuc_dx(species), dabar_dx(species), dzbar_dx(species), dmc_dx(species))
         
         info = 0
         lwork = net_work_size(handle, info) 
         if (info /= 0) stop 1
         
         allocate(work(lwork),  &
               rate_factors(num_reactions), &
               actual_Qs(num_reactions), &
               actual_neuQs(num_reactions), &
               from_weaklib(num_reactions), &
               stat=info)
         if (info /= 0) stop 2
         
         rate_factors(:) = 1
         weak_rate_factor = 1
         
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

         if (net_file == 'mesa_201.net') then
            
            nrates_to_show = 2
         
            rates_to_show(1:nrates_to_show) = (/  &
               ir_ar36_ag_ca40, &
               ir_ca40_ga_ar36  /)

                                  xin(net_iso(ife56))=     6.3551174304779179D-01
                                  xin(net_iso(icr52))=     1.0518849309100423D-01
                                  xin(net_iso(ini60))=     6.8579809304547934D-02
                                  xin(net_iso(ife55))=     3.0800958229563902D-02
                                  xin(net_iso(imn55))=     2.1373990699858084D-02
                                  xin(net_iso(ico57))=     2.0785551124804885D-02
                                  xin(net_iso(ife57))=     1.6083543077536507D-02
                                  xin(net_iso(imn53))=     1.5959192021165327D-02
                                  xin(net_iso(ico59))=     1.3346191916961030D-02
                                  xin(net_iso(ife58))=     1.2810477227656202D-02
                                  xin(net_iso(ife54))=     1.2428871023625330D-02
                                  xin(net_iso(ini62))=     1.1831182533246063D-02
                                  xin(net_iso(icr53))=     8.1041573656346795D-03
                                  xin(net_iso(ini61))=     5.1462544919734727D-03
                                  xin(net_iso(imn54))=     4.8562102089927551D-03
                                  xin(net_iso(icr54))=     3.4271335428296490D-03
                                  xin(net_iso(ini59))=     2.9811021846265560D-03
                                  xin(net_iso(ini58))=     2.8713349650366189D-03
                                   xin(net_iso(iv51))=     2.0988150935152424D-03
                                  xin(net_iso(ico58))=     2.0282210582857861D-03
                                  xin(net_iso(icr51))=     1.2727926750247761D-03
                                  xin(net_iso(icr50))=     3.5421790633561727D-04
                                  xin(net_iso(icu63))=     3.0335040211002022D-04
                                  xin(net_iso(iti48))=     2.8512639104984498D-04
                                  xin(net_iso(iti50))=     2.8394752519368671D-04
                                  xin(net_iso(ico56))=     2.4310701594699283D-04
                                  xin(net_iso(ico60))=     1.5319574812323961D-04
                                   xin(net_iso(ihe4))=     1.3780017854631601D-04
                                  xin(net_iso(imn56))=     1.1066944504563417D-04
                                   xin(net_iso(iv50))=     8.7723703257792468D-05
                                   xin(net_iso(iv49))=     8.0539507322740820D-05
                                  xin(net_iso(icu61))=     6.0898214714637941D-05
                                  xin(net_iso(iti49))=     5.8026224632063015D-05
                                  xin(net_iso(imn52))=     4.2516037186521133D-05
                                  xin(net_iso(ico61))=     4.1332712278653746D-05
                                  xin(net_iso(ini63))=     3.9285182071536010D-05
                                  xin(net_iso(ico55))=     3.3984664667118780D-05
                                  xin(net_iso(ife59))=     3.1874001320244153D-05
                                  xin(net_iso(izn64))=     2.7904240321969555D-05
                                  xin(net_iso(izn66))=     2.3416686728782276D-05
                                  xin(net_iso(icu62))=     2.1144678609352833D-05
                                  xin(net_iso(ini57))=     1.1679099363693715D-05
                                   xin(net_iso(iv52))=     1.0016381776869645D-05
                                  xin(net_iso(ini64))=     8.9564547480773243D-06
                                  xin(net_iso(iti47))=     8.7562205406651916D-06
                                  xin(net_iso(icu65))=     8.3830901771893844D-06
                                  xin(net_iso(iti46))=     6.5019528109988749D-06
                                  xin(net_iso(icu64))=     6.1556773376380492D-06
                                  xin(net_iso(iar38))=     5.0886468271422554D-06
                                  xin(net_iso(ife53))=     4.4893983930970075D-06
                                   xin(net_iso(is34))=     4.2369862519232918D-06
                                  xin(net_iso(izn65))=     4.1809380230467465D-06
                                  xin(net_iso(icr55))=     3.5647078269917385D-06
                                  xin(net_iso(imn51))=     1.3849206094166064D-06
                                  xin(net_iso(ife60))=     1.0919469947619407D-06
                                    xin(net_iso(ih1))=     9.8211241034612710D-07
                                  xin(net_iso(ica44))=     6.4721481111312556D-07
                                  xin(net_iso(isi30))=     6.0482126630410841D-07
                                  xin(net_iso(ica42))=     5.8926684360031471D-07
                                  xin(net_iso(isi28))=     5.5342250413192408D-07
                                   xin(net_iso(iv48))=     5.4911502105605942D-07
                                   xin(net_iso(iv53))=     5.4571869339657211D-07
                                  xin(net_iso(icu60))=     4.5761289432308086D-07
                                  xin(net_iso(izn63))=     4.4652045082610858D-07
                                  xin(net_iso(isc47))=     4.4099241256913694D-07
                                  xin(net_iso(ini56))=     3.9547421607331126D-07
                                  xin(net_iso(iti51))=     3.6358450091693021D-07
                                  xin(net_iso(izn62))=     3.0268870369662295D-07
                                   xin(net_iso(is32))=     3.0184105591459131D-07
                                  xin(net_iso(icr49))=     2.7706803242086906D-07
                                  xin(net_iso(isc45))=     2.5466905171154329D-07
                                   xin(net_iso(ik39))=     1.9906885458096793D-07
                                  xin(net_iso(icl35))=     1.5826438708200070D-07
                                   xin(net_iso(is33))=     1.2794750349676913D-07
                                   xin(net_iso(ip31))=     1.1805533268957108D-07
                                  xin(net_iso(icl37))=     1.0108890802543261D-07
                                  xin(net_iso(ica43))=     9.5766190955127930D-08
                                  xin(net_iso(iar36))=     8.3287350202069049D-08
                                  xin(net_iso(isi29))=     7.6157023642649312D-08
                                  xin(net_iso(isc49))=     7.2645544618960897D-08
                                  xin(net_iso(icu59))=     6.9036259006119691D-08
                                  xin(net_iso(ica40))=     6.4387988887989985D-08
                                  xin(net_iso(isc46))=     6.1294200120692868D-08
                                  xin(net_iso(iar37))=     5.0049708949178339D-08
                                  xin(net_iso(ineut))=     4.7236564548991118D-08
                                  xin(net_iso(isc48))=     3.7392395931746378D-08
                                  xin(net_iso(ife52))=     3.0982331086000098D-08
                                  xin(net_iso(icr56))=     2.9477366114308358D-08
                                  xin(net_iso(ica41))=     2.7100443580454983D-08
                                   xin(net_iso(is35))=     2.6527306613430685D-08
                                  xin(net_iso(ica45))=     2.5532508173462626D-08
                                  xin(net_iso(ico62))=     2.3608043613616947D-08
                                  xin(net_iso(iar39))=     2.3503192609361190D-08
                                  xin(net_iso(ica46))=     2.2455595751043146D-08
                                   xin(net_iso(iv47))=     1.9789662501752072D-08
                                  xin(net_iso(icu66))=     1.9285467280104737D-08
                                  xin(net_iso(icl36))=     1.8710358753573163D-08
                                   xin(net_iso(is36))=     1.6433794658525574D-08
                                   xin(net_iso(ik41))=     1.1021812817088955D-08
                                  xin(net_iso(ini65))=     1.0432196346423500D-08
                                   xin(net_iso(ip33))=     8.9625871046594568D-09
                                   xin(net_iso(ik40))=     7.2106087354006743D-09
                                  xin(net_iso(iar40))=     7.1174073521523365D-09
                                  xin(net_iso(iti45))=     5.6870435962784455D-09
                                   xin(net_iso(ip32))=     4.9537776441010461D-09
                                  xin(net_iso(icr48))=     3.0709436939798925D-09
                                  xin(net_iso(isc44))=     2.7862949914482315D-09
                                  xin(net_iso(isc43))=     1.4283233379580965D-09
                                  xin(net_iso(isi31))=     1.3747008133379580D-09
                                  xin(net_iso(iti52))=     1.2866447687101849D-09
                                  xin(net_iso(ico63))=     1.0456195723572430D-09
                                  xin(net_iso(iti44))=     7.1579364938842059D-10
                                   xin(net_iso(io16))=     5.6818053126812287D-10
                                  xin(net_iso(ial27))=     5.6117840295305725D-10
                                  xin(net_iso(ica47))=     5.4798226854137171D-10
                                  xin(net_iso(img24))=     3.8679049700264050D-10
                                  xin(net_iso(ini66))=     2.8640604416758339D-10
                                  xin(net_iso(izn61))=     2.4995195933116682D-10
                                  xin(net_iso(ica48))=     1.9626419275146166D-10
                                  xin(net_iso(ife61))=     1.8711288748163173D-10
                                   xin(net_iso(ik42))=     1.4923884785773920D-10
                                  xin(net_iso(isi32))=     1.4443919428884894D-10
                                   xin(net_iso(ip30))=     1.3908837795296053D-10
                                   xin(net_iso(ic12))=     1.3185610207511085D-10
                                   xin(net_iso(ik43))=     1.1518951971920992D-10
                                   xin(net_iso(iv54))=     8.8440974843643898D-11
                                  xin(net_iso(img26))=     8.5260453209638504D-11
                                  xin(net_iso(icl38))=     2.5853881572793256D-11
                                  xin(net_iso(isc50))=     1.7292118708137885D-11
                                  xin(net_iso(iar41))=     1.0665700916421081D-11
                                  xin(net_iso(img25))=     9.2051027558630546D-12
                                  xin(net_iso(ial28))=     8.9602126799277990D-12
                                  xin(net_iso(izn60))=     7.0658622923470296D-12
                                   xin(net_iso(ip34))=     4.2003981867966896D-12
                                  xin(net_iso(icr57))=     2.8195149736768269D-12
                                  xin(net_iso(ine20))=     1.1829006560529443D-12
                                  xin(net_iso(ife62))=     1.0149415562457171D-12
                                   xin(net_iso(ik44))=     5.5654687896999145D-13
                                   xin(net_iso(is31))=     4.0845373817884839D-13
                                   xin(net_iso(iv55))=     2.8610531204397973D-13
                                  xin(net_iso(ina23))=     2.5122330767236196D-13
                                   xin(net_iso(is37))=     2.1529809224524208D-13
                                  xin(net_iso(iti53))=     1.4129020276964862D-13
                                  xin(net_iso(iar35))=     1.2934962278023034D-13
                                  xin(net_iso(ial26))=     1.2407881461251650D-13
                                   xin(net_iso(in15))=     1.1945765698183132D-13
                                  xin(net_iso(ico64))=     8.1590671587313667D-14
                                  xin(net_iso(img27))=     6.7969591881380018D-14
                                    xin(net_iso(ih2))=     5.6613884110319004D-14
                                  xin(net_iso(ini67))=     5.1506647175988674D-14
                                  xin(net_iso(ica39))=     3.8945192023978035D-14
                                  xin(net_iso(ini55))=     3.0883773851523300D-14
                                  xin(net_iso(ine22))=     1.3808913342478939D-14
                                  xin(net_iso(ica49))=     1.0542673783683999D-14
                                  xin(net_iso(isc51))=     9.2083686560605166D-15
                                  xin(net_iso(isi27))=     8.3464191225256989D-15
                                  xin(net_iso(ine21))=     6.6697557229646203D-15
                                  xin(net_iso(ife51))=     6.5787879505834196D-15
                                   xin(net_iso(io17))=     3.8661065919908615D-15
                                   xin(net_iso(ic13))=     2.2604411883637647D-15
                                  xin(net_iso(icr58))=     2.1539938862813166D-15
                                  xin(net_iso(isi33))=     1.4888165334138879D-15
                                  xin(net_iso(ina24))=     7.2312125147882022D-16
                                  xin(net_iso(ial25))=     5.3576177397850227D-16
                                  xin(net_iso(ico65))=     4.9566812556376405D-16
                                  xin(net_iso(icr47))=     3.0194388465619186D-16
                                   xin(net_iso(in14))=     2.6891785216435022D-16
                                  xin(net_iso(ina22))=     2.5931749891034273D-16
                                   xin(net_iso(io15))=     2.5261003710407275D-16
                                  xin(net_iso(ini68))=     2.2419710841076782D-16
                                   xin(net_iso(io18))=     1.8206042351246347D-16
                                  xin(net_iso(iti43))=     9.9895985562055471D-17
                                   xin(net_iso(if19))=     7.0445166521693986D-17
                                   xin(net_iso(ihe3))=     6.8591249543794866D-17
                                  xin(net_iso(iti54))=     3.3955810406819116D-17
                                  xin(net_iso(ife63))=     3.0222399040043984D-17
                                   xin(net_iso(in13))=     2.5097133179904447D-17
                                  xin(net_iso(izn59))=     2.3782224732332168D-17
                                  xin(net_iso(img23))=     2.2784900370197838D-17
                                   xin(net_iso(if17))=     1.0052207505944401D-17
                                   xin(net_iso(if18))=     6.3013547340360942D-18
                                   xin(net_iso(iv56))=     2.1677062516769839D-18
                                  xin(net_iso(ina21))=     1.9109919103480182D-18
                                  xin(net_iso(ine23))=     1.2143041459039416D-18
                                   xin(net_iso(ili6))=     1.4203908924439747D-19
                                   xin(net_iso(ib11))=     1.1320699694157424D-19
                                   xin(net_iso(if20))=     4.4936577179455616D-20
                                  xin(net_iso(ine19))=     2.6387620496069204D-20
                                  xin(net_iso(ife64))=     1.1841283303587735D-20
                                   xin(net_iso(ib10))=     1.0494357182634838D-20
                                   xin(net_iso(ibe9))=     7.8292004126082962D-21
                                  xin(net_iso(ico66))=     6.4514234785580282D-21
                                   xin(net_iso(ili7))=     6.4231341717191173D-21
                                   xin(net_iso(in16))=     4.7443318701592232D-21
                                   xin(net_iso(ibe7))=     3.3532048357084064D-22
                                   xin(net_iso(io19))=     3.2062683754970905D-22
                                  xin(net_iso(ibe10))=     6.5331631260779713D-23
                                  xin(net_iso(ico67))=     4.9629777598135805D-24
                                  xin(net_iso(ife65))=     3.7335326045384740D-26
                                  xin(net_iso(ife66))=     8.6772104841406824D-30
                                    xin(net_iso(ib8))=     4.1439050548710376D-31
                              
                              write(*,*) 'sum xin', sum(xin(:))

                                                 logT =    9.6532818288064650D+00
                                               logRho =    7.9479966082179185D+00
                                                  eta =    2.7403163311838425D+00
                                              theta_e =    0.0000000000000000D+00
                                   
                                   
               screening_mode = extended_screening
               
               call net_set_logTcut(handle, 0d0, 0d0, info)
               if (info /= 0) then
                  write(*,*) 'failed in net_set_logTcut'
                  stop 1
               end if

         else if (net_file == 'approx21_cr60_plus_co56.net') then
            
            nrates_to_show = 2
         
            rates_to_show(1:nrates_to_show) = (/  &
               ir_ar36_ag_ca40, &
               ir_ca40_ga_ar36  /)

                                  xin(net_iso(ife56))=     8.3990195536270140D-01
                                  xin(net_iso(ife54))=     9.1219636106456448D-02
                                   xin(net_iso(ihe4))=     5.4913235623813200D-02
                                  xin(net_iso(icr60))=     5.7674995592308142D-03
                                  xin(net_iso(ico56))=     4.9881560076575713D-03
                                  xin(net_iso(iprot))=     1.3632657569069438D-03
                                  xin(net_iso(isi28))=     7.8626218715835341D-04
                                   xin(net_iso(is32))=     4.1969480135437183D-04
                                  xin(net_iso(iar36))=     1.8526179220718401D-04
                                  xin(net_iso(ica40))=     1.2242616041375457D-04
                                  xin(net_iso(ineut))=     1.1980892341740248D-04
                                  xin(net_iso(ini56))=     9.7550916738808739D-05
                                  xin(net_iso(ife52))=     5.8905163394112936D-05
                                  xin(net_iso(icr48))=     1.8019113042269716D-05
                                   xin(net_iso(io16))=     1.3408668333884428D-05
                                   xin(net_iso(ic12))=     1.1016907822550631D-05
                                  xin(net_iso(img24))=     8.2286117186337854D-06
                                  xin(net_iso(iti44))=     5.4079788265411310D-06
                                  xin(net_iso(ine20))=     2.6035880565410054D-07
                                    xin(net_iso(ih1))=     4.9844580237752611D-20
                                   xin(net_iso(ihe3))=     8.8329616569722745D-21
                                   xin(net_iso(in14))=     1.8035412315021348D-22
                              
                              write(*,*) 'sum xin', sum(xin(:))

                                  logT =    7.2162656791795046D+00
                                logRho =    2.2626037247171089D+00
                                   eta =   -1.4737364370795314D+00
                               theta_e =    0
                                                 logT =    9.8477686652860221D+00
                                               logRho =    8.4951102889416124D+00
                                                  eta =    2.9827856303789755D+00
                                              theta_e =    0.0000000000000000D+00
                                   
                                   
               screening_mode = extended_screening
               
               call net_set_logTcut(handle, 0d0, 0d0, info)
               if (info /= 0) then
                  write(*,*) 'failed in net_set_logTcut'
                  stop 1
               end if

         else if (net_file == 'pp_extras.net') then
            
            nrates_to_show = 6
         
            rates_to_show(1:nrates_to_show) = (/  &
            ir_h1_h1_wk_h2, &
            ir_h2_pg_he3, &
            ir_be7_wk_li7, &
            ir_b8_wk_he4_he4, &
            ir_he3_he3_to_h1_h1_he4, &
            ir_he3_ag_be7 /)
     
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

         else if (net_file == 'agb.net') then
            
            nrates_to_show = 4
         
            rates_to_show(1:nrates_to_show) = (/  &
            ir_h1_h1_wk_h2, &
            ir_c13_an_o16, &
            ir_f19_ap_ne22, &
            ir_he3_ag_be7 /)
     
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
         
            rates_to_show(1:nrates_to_show) = (/  &
            rates_reaction_id('r_n13_wk_c13'),                &
            rates_reaction_id('r_o15_wk_n15'),                &
            rates_reaction_id('r_f17_wk_o17'),                &
            rates_reaction_id('r_f18_wk_o18'),                &
            rates_reaction_id('r_o14_wk_n14'),                &
            rates_reaction_id('r_ne18_wk_f18'),                &
            rates_reaction_id('r_ne19_wk_f19'),                &
            ir_he4_he4_he4_to_c12 /)
     
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
 
                                  logT =    9.6d0
                                logRho =    6.0d0
                                   rho =    7.8571498592117219D+00
                                     T =    8.5648111120065376D+06
                                  abar =    1.2655060647252907D+00
                                  zbar =    1.0901664301076275D+00
                                 z2bar =    1.3036906023574921D+00
                                    ye =    8.6144702146826535D-01
                                   eta =   -3.4387570967781595D+00                                   
                                   
               screening_mode = extended_screening

         else if (net_file == 'approx21.net' .or. &
                  net_file == 'approx21_plus_co56.net' .or. &
                  net_file == 'approx21_test.net' .or. &
                  net_file == 'approx21_old.net' .or. &
                  net_file == 'approx21_new.net') then
            

            nrates_to_show = 5
            
            rates_to_show(1:nrates_to_show) = (/  &
            irprot_to_neut,            &
            irneut_to_prot, &
            irni56ec_to_fe56, &
            ir_fe52_ag_ni56,            &
            ir_ni56_ga_fe52             &
            /)
                     xin = 0
            
                                  xin(net_iso(ife56))=     8.0387021484318166D-01
                                  xin(net_iso(ife54))=     1.6096648736760832D-01
                                  xin(net_iso(icr56))=     2.9480945535920525D-02
                                   xin(net_iso(ihe4))=     4.8624637161320565D-03
                                  xin(net_iso(ini56))=     2.8376270731890360D-04
                                  xin(net_iso(isi28))=     1.5018628906135739D-04
                                   xin(net_iso(is32))=     1.1613271635573457D-04
                                  xin(net_iso(iprot))=     1.1139431633673653D-04
                                  xin(net_iso(ica40))=     5.3688377473494185D-05
                                  xin(net_iso(iar36))=     5.2702831822567062D-05
                                  xin(net_iso(ife52))=     3.7866504131935185D-05
                                  xin(net_iso(icr48))=     5.8401123037667974D-06
                                  xin(net_iso(ineut))=     4.9141227703118397D-06
                                  xin(net_iso(iti44))=     1.5085038561154746D-06
                                   xin(net_iso(io16))=     9.5384019609049255D-07
                                  xin(net_iso(img24))=     6.3808207717725580D-07
                                   xin(net_iso(ic12))=     2.9048656673991868D-07
                                  xin(net_iso(ine20))=     9.6468865023609427D-09
                                   xin(net_iso(ihe3))=     4.6203603862263096D-80
                                   xin(net_iso(in14))=     7.5867472225841235D-99
                                    xin(net_iso(ih1))=     0 ! 9.9987777520212238-100
                              write(*,*) 'test case sum xin', sum(xin(1:species))
                              

                                                 logT =    9.6d0
                                               logRho =    6d0
                                                    T =    10**logT
                                                  rho =    10**logRho
                                                 abar =    5.2051025574883582D+01
                                                 zbar =    2.4269152265763136D+01
                                                z2bar =    6.2641754206392147D+02
                                                   ye =    4.6625694686549729D-01
                                                  eta =    4.3030680106736412D+00
                        screening_mode = extended_screening
                                              theta_e =    0.0000000000000000D+00

         else
            
            write(*, *) 'need to define setup for net_file ', trim(net_file)
            stop 'Do_One_Test'
         
         end if
         
         Rho = exp10_cr(logRho)
         T = exp10_cr(logT)
         
         write(*, *)
         write(*, *)
         
         info = 0
         
         ierr = 0
         call composition_info( &
            species, chem_id, xin, xh, xhe, abar, zbar, z2bar, ye,  &
            mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
         
      	write(*,'(a40,e26.16)') 'xh', xh
      	write(*,'(a40,e26.16)') 'xhe', xhe
      	write(*,'(a40,e26.16)') 'abar', abar
      	write(*,'(a40,e26.16)') 'zbar', zbar
      	write(*,'(a40,e26.16)') 'z2bar', z2bar
      	write(*,'(a40,e26.16)') 'ye', ye
      	do i = 1, species
      	   write(*,'(a40,i6,e26.16)')  'init x ' // trim(chem_isos% name(chem_id(i))), i, xin(i)
      	end do
      	write(*,*)
      	write(*,'(a40,e26.16)') 'logT', logT
      	write(*,'(a40,e26.16)') 'logRho', logRho
      	write(*,'(a40,e26.16)') 'eta', eta

         if (do_timing) then
            nreps = 100000
            call zero_net_timing(g)
            g% doing_timing = .true.
            call system_clock(time0)
         else
            nreps = 1
         end if
         
         skip_jacobian = .false.
         
         if (.false.) then
            write(*,*) 'call net_get_rates_only'
            call net_get_rates_only( &
               handle, n, species, num_reactions,  &
               xin, T, logT, Rho, logRho,  &
               abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
               rate_factors, weak_rate_factor, &
               std_reaction_Qs, std_reaction_neuQs, &
               screening_mode, theta_e_for_graboske_et_al,  &
               lwork, work, ierr)
            stop 'net_get_rates_only'
         end if

         do rep=1,nreps
            call net_get_with_Qs(handle, skip_jacobian, n, species, num_reactions,  &
                  xin, T, logT, Rho, logRho,  &
                  abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
                   rate_factors, weak_rate_factor, &
                  std_reaction_Qs, std_reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
                  screening_mode, theta_e,      &
                  eps_nuc_categories, eps_neu_total,  &
                  lwork, work, actual_Qs, actual_neuQs, from_weaklib, info)
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
                  write(*,'(a30,99f20.10)') 'weaklib ' // trim(reaction_Name(reaction_id(i))),  &
                     actual_Qs(i), actual_neuQs(i), actual_Qs(i) - actual_neuQs(i)
               else
                  write(*,'(a30,99f20.10)') trim(reaction_Name(reaction_id(i))),  &
                     actual_Qs(i), actual_neuQs(i), actual_Qs(i) - actual_neuQs(i)
               end if
            end do
            write(*,*)
            stop
         end if
         
         
         call get_net_rate_ptrs(g% handle, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_net_rate_ptrs'
            stop 1
         end if
         
         write(*,2) 'screening_mode', screening_mode
            
      	if (.true.) then
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
         	write(*,*)
         	
      	   write(*,1) 'eps_nuc', eps_nuc
      	   write(*,1) 'd_epsnuc_dlnd', d_eps_nuc_dRho*Rho
      	   write(*,1) 'd_epsnuc_dlnT', d_eps_nuc_dT*T
         	write(*,*)
         	
      	   write(*,1) 'log eps_nuc', log10(eps_nuc)
      	   write(*,1) 'd_lnepsnuc_dlnd', d_eps_nuc_dRho*Rho/eps_nuc
      	   write(*,1) 'd_lnepsnuc_dlnT', d_eps_nuc_dT*T/eps_nuc
         	write(*,*)

         	
      	   stop





         	do i = 1, species
         	   write(*,1)  'd_eps_nuc_dx ' // trim(chem_isos% name(chem_id(i))), d_eps_nuc_dx(i)
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


         	
         	if (.true.) then
            	do i = 1, species
            	   write(*,1)  'd_dxdt_dx(:,neut) ' // &
            	      trim(chem_isos% name(chem_id(i))), d_dxdt_dx(i, net_iso(ineut))
            	end do
            	write(*,*)
         	end if


         	do i = 1, species
         	   write(*,1)  'dxdt ' // trim(chem_isos% name(chem_id(i))), dxdt(i)
         	end do
         	write(*,1) 'sum(dxdt)', sum(dxdt(1:species))
            
            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'rate_raw ' // trim(reaction_Name(j)),  &
                                 rate_raw(reaction_table(j))
            end do
            write(*,*)
            
            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'd_rate_raw_dT ' // trim(reaction_Name(j)),  &
                                 rate_raw_dT(reaction_table(j))
            end do
            write(*,*)
            
            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'rate_screened ' // trim(reaction_Name(j)), &
                                 rate_screened(reaction_table(j))
            end do
            write(*,*)
            
         	do i = 1, species
         	   write(*,1)  'x ' // trim(chem_isos% name(chem_id(i))), xin(i)
         	end do
         	write(*,*)
         	
         	do i = 1, species
         	   if (-dxdt(i) > 1d-90)  &
               	write(*,1)  'x/dxdt ' // trim(chem_isos% name(chem_id(i))), xin(i)/dxdt(i)
         	end do
         	write(*,*)
         	
         	if (.false.) then
            	do i = 1, num_categories
            	   if (abs(eps_nuc_categories(i)) < 1d-20) cycle
            	   write(*,1)  'eps_nuc_cat ' // trim(category_name(i)), eps_nuc_categories(i)
            	end do
            	write(*,*)
         	end if
         	
         	stop
      	end if

         write(*,*)
         write(*,*)
         write(*,*) 'net_name ', trim(net_file)
         write(*,*) 'species', species
         write(*,1) 'theta_e =', theta_e
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
            write(*,1) trim(category_name(j)), eps_nuc_categories( j)
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


         deallocate(work, rate_factors, actual_Qs, actual_neuQs, from_weaklib)


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




