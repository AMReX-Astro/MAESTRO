! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module net_derivs
      use net_def
      use const_def
      use chem_def
      use crlibm_lib
      use net_derivs_support
      use rates_def

      implicit none
      
      real(dp), parameter :: tiny_rate = 1d-50
      
      contains

#ifdef offload
      !dir$ options /offload_attribute_target=mic
#endif

      subroutine get_derivs( &
            n, dydt, eps_nuc_MeV, eta, ye, logtemp, temp, den, abar, zbar, &
            reuse_rate_screened, num_reactions, rate_factors, &
            symbolic, just_dydt, ierr)
         type (Net_Info), pointer :: n
         real(qp), pointer, intent(out) :: dydt(:,:)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         logical, intent(in) :: reuse_rate_screened
         integer, intent(in) :: num_reactions
         real(dp), intent(in) ::eta, ye, logtemp, temp, den, abar, zbar, &
            rate_factors(:)
         logical, intent(in) :: symbolic, just_dydt
         integer, intent(out) :: ierr

         logical :: all_okay, show_derivs_dydt
         integer, pointer, dimension(:) :: &
            reaction_kind, reverse_id, itab, rtab, reaction_id
         integer :: ir, r_ir, r_i, i, j, op_err, kind, jmax, icat_f, icat_r
         logical, target :: deriv_flgs_data(num_reactions)
         logical, pointer :: deriv_flgs(:)
         type (Net_General_Info), pointer  :: g
         real(dp), pointer :: y(:)
         integer :: time0, time1
         real(dp) :: T9, T932, eps_nuc_cancel_factor, eps_factor
         real(dp) :: eps_nuc_MeV_react, eps_nuc_MeV_old
         include 'formats'
         
         ierr = 0

         T9 = temp*1d-9
         T932 = T9*sqrt(T9)
         
         y => n% y
         g => n% g
         
         if (.true. .or. logtemp <= g% logT_lo_eps_nuc_cancel) then
            eps_nuc_cancel_factor = 1d0
         else if (logtemp >= g% logT_hi_eps_nuc_cancel) then
            eps_nuc_cancel_factor = 0d0
         else
            eps_nuc_cancel_factor = &
               (g% logT_hi_eps_nuc_cancel - logtemp)/&
               (g% logT_hi_eps_nuc_cancel - g% logT_lo_eps_nuc_cancel)
         end if
         
         !write(*,1) 'eps_nuc_cancel_factor', eps_nuc_cancel_factor, logtemp

         itab => g% net_iso
         rtab => g% net_reaction
         reaction_kind => g% reaction_reaclib_kind
         reverse_id => g% reverse_id_for_kind_ne_other
         reaction_id => g% reaction_id
         
         show_derivs_dydt = .false.
         if (show_dydt) then
            i = itab(io16)
            if (i > 0) &
               show_derivs_dydt = (abs(n% y(i) - show_dydt_y) < 1d-14)
         end if
         
         
         deriv_flgs => deriv_flgs_data
         if (checking_deriv_flags) deriv_flgs(:) = .false.

         dydt = 0
         if (.not. just_dydt) n% d_dydt_dy = 0
         
         ierr = 0
         eps_nuc_MeV = 0d0
         eps_nuc_MeV_react = 0d0
         eps_nuc_MeV_old = 0d0

         if (just_dydt) then
            jmax = 1 ! =  i_rate
         else
            jmax = num_rvs
         end if
      
         i = 1
         do while (i <= num_reactions)
      
            if (ierr /= 0) exit
            
            ir = reaction_id(i)
            icat_f = reaction_categories(ir)

            if (i < num_reactions) then
               r_i = i+1 ! reactions are ordered so reverse immediately follows forward
               r_ir = reaction_id(r_i)
               icat_r = reaction_categories(r_ir)
            else
               r_i = 0
               r_ir = 0
               icat_r = 0
            end if
            

            kind = reaction_kind(i)   
            if ( &
                !kind == ng_kind .or. &
                !kind == pn_kind .or. &
                !kind == pg_kind .or. &                
                !kind == ap_kind .or. & 
                !kind == an_kind .or. &                
                !kind == ag_kind .or. &                
                !kind == general_one_one_kind .or. &
                !kind == general_two_one_kind .or. &
                !kind == general_two_two_kind .or. &
                kind == -1 &
                ) kind = other_kind
            
            !kind = other_kind  ! TESTING
            
            !if (reaction_name(ir) == 'rfe52aprot_to_ni56') then
            !   write(*,2) 'rfe52aprot_to_ni56 kind', kind
            !   stop
            !end if
            
            eps_factor = 1d0
                        
            ir = reaction_id(i)         
            write(*,'(/,a,2i6)') ' reaction name <' // trim(reaction_Name(ir)) // '>', i, ir

            select case(kind)
               case (other_kind)
                  call get1_derivs( &
                     n, dydt, eps_nuc_MeV, i, eta, ye, temp, den, abar, zbar, &
                     reuse_rate_screened, num_reactions, rate_factors, rtab, itab, &
                     deriv_flgs, symbolic, just_dydt, &
                     ierr)
                  i = i+1
                  eps_nuc_MeV_react = eps_nuc_MeV(i_rate) - eps_nuc_MeV_old
                  write(*,*) 'eps_nuc: ', eps_nuc_MeV_react*Qconv
                  eps_nuc_MeV_old = eps_nuc_MeV(i_rate)
                  cycle
               case (ng_kind)
                  eps_factor = eps_nuc_cancel_factor
                  call get_basic_2_to_1_derivs(i,ierr)
               case (pn_kind)
                  eps_factor = eps_nuc_cancel_factor
                  call get_basic_2_to_2_derivs(i,ierr)
               case (pg_kind)
                  eps_factor = eps_nuc_cancel_factor
                  call get_basic_2_to_1_derivs(i,ierr)
               case (ap_kind)
                  call get_basic_2_to_2_derivs(i,ierr)
               case (an_kind)
                  call get_basic_2_to_2_derivs(i,ierr)
               case (ag_kind)
                  call get_basic_2_to_1_derivs(i,ierr)
               case (general_one_one_kind)
                  call get_general_1_to_1_derivs(i,ierr)
               case (general_two_one_kind)
                  call get_general_2_to_1_derivs(i,ierr)
               case (general_two_two_kind)
                  call get_general_2_to_2_derivs(i,ierr)
               case default
                  stop 'confusion in net wrt reaction kind'
            end select
            i = i+2
            ir = reaction_id(i-1)         
            write(*,'(a,2i6)') ' reaction name <' // trim(reaction_Name(ir)) // '>', i-1, ir
            eps_nuc_MeV_react = eps_nuc_MeV(i_rate) - eps_nuc_MeV_old
            write(*,*) 'net eps_nuc: ', eps_nuc_MeV_react*Qconv
            eps_nuc_MeV_old = eps_nuc_MeV(i_rate)
         end do

         
         contains

#ifdef offload
         !dir$ attributes offload: mic :: get_general_1_to_1_derivs
#endif
         subroutine get_general_1_to_1_derivs(i,ierr) ! e.g., 2 c12 -> mg24, 3 he4 -> c12
            integer, intent(in) :: i
            integer, intent(out) :: ierr
            
            real(qp) :: b, b_f, b_r, rate
            real(dp) :: d, d_f, d_r, d1, d2, Q, ys_f, ys_r, &
               d_ysf_dy1, d_ysr_dy2
            integer :: c1, c2, i1, i2, o1, o3
            
            include 'formats'
         
            ierr = 0            
            
            ! forward reaction is c1 i1 -> c2 i2
            c1 = reaction_inputs(1,ir)  
            i1 = itab(reaction_inputs(2,ir))    
            c2 = reaction_outputs(1,ir)
            i2 = itab(reaction_outputs(2,ir))

            if (symbolic) then
               n% d_dydt_dy(i1,i1) = 1
               n% d_dydt_dy(i1,i2) = 1
               n% d_dydt_dy(i2,i1) = 1
               n% d_dydt_dy(i2,i2) = 1
               return
            end if
         
            select case(c1)
               case (1)
                  ys_f = y(i1)
                  d_ysf_dy1 = 1d0
               case (2)
                  ys_f = y(i1)*y(i1)/2
                  d_ysf_dy1 = y(i1)
               case (3)
                  ys_f = y(i1)*y(i1)*y(i1)/6
                  d_ysf_dy1 = y(i1)*y(i1)/2
               case default
                  write(*,2) 'c1 bad for ' // trim(reaction_name(ir)), c1
                  stop 'get_general_1_to_1_derivs'
            end select
            d_f = ys_f
            
            select case(c2)
               case (1)
                  ys_r = y(i2)
                  d_ysr_dy2 = 1d0
               case (2)
                  ys_r = y(i2)*y(i2)/2
                  d_ysr_dy2 = y(i2)
               case (3)
                  ys_r = y(i2)*y(i2)*y(i2)/6
                  d_ysr_dy2 = y(i2)*y(i2)/2
               case default
                  write(*,2) 'c2 bad for ' // trim(reaction_name(ir)), c2
                  stop 'get_general_1_to_1_derivs'
            end select
            d_r = ys_r

            rate = n% rate_screened(i)
            b_f = d_f*rate
            rate = n% rate_screened(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate,i1) = dydt(i_rate,i1) - c1*b
            dydt(i_rate,i2) = dydt(i_rate,i2) + c2*b
            if (just_dydt) return

            Q = n% reaction_Qs(ir)*eps_factor
            b_f = Q*b_f
            b_r = -Q*b_r
            b = b_f + b_r
            eps_nuc_MeV(i_rate) = eps_nuc_MeV(i_rate) + b
            n% eps_nuc_categories(icat_f) = n% eps_nuc_categories(icat_f) + b_f
            n% eps_nuc_categories(icat_r) = n% eps_nuc_categories(icat_r) + b_r
            if (show_eps_nuc ) &
               write(*,1) trim(reaction_Name(ir)) // ' eps_nuc',  b, b_f, b_r
                        
            rate = n% rate_screened_dT(i)
            b_f = d_f*rate
            rate = n% rate_screened_dT(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate_dT,i1) = dydt(i_rate_dT,i1) - c1*b
            dydt(i_rate_dT,i2) = dydt(i_rate_dT,i2) + c2*b
            b_f = Q*b_f
            b_r = -Q*b_r
            b = b_f + b_r
            eps_nuc_MeV(i_rate_dT) = eps_nuc_MeV(i_rate_dT) + b
                        
            rate = n% rate_screened_dRho(i)
            b_f = d_f*rate
            rate = n% rate_screened_dRho(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate_dRho,i1) = dydt(i_rate_dRho,i1) - c1*b
            dydt(i_rate_dRho,i2) = dydt(i_rate_dRho,i2) + c2*b
            b_f = Q*b_f
            b_r = -Q*b_r
            b = b_f + b_r
            eps_nuc_MeV(i_rate_dRho) = eps_nuc_MeV(i_rate_dRho) + b
               
            if (checking_deriv_flags) then
               deriv_flgs(i) = .true.
               deriv_flgs(r_i) = .true.
            end if

            d1 = d_ysf_dy1*n% rate_screened(i) ! d(rate_f)/d(y1)
            d2 = d_ysr_dy2*n% rate_screened(r_i) ! d(rate_r)/d(y2)

            n% d_eps_nuc_dy(i1) = n% d_eps_nuc_dy(i1) + Q*d1
            n% d_eps_nuc_dy(i2) = n% d_eps_nuc_dy(i2) - Q*d2  
            
            n% d_dydt_dy(i1,i1) = n% d_dydt_dy(i1,i1) - c1*d1
            n% d_dydt_dy(i2,i1) = n% d_dydt_dy(i2,i1) + c2*d1

            n% d_dydt_dy(i1,i2) = n% d_dydt_dy(i1,i2) + c1*d2
            n% d_dydt_dy(i2,i2) = n% d_dydt_dy(i2,i2) - c2*d2
                        
         end subroutine get_general_1_to_1_derivs

#ifdef offload
         !dir$ attributes offload: mic :: get_general_2_to_1_derivs
#endif
         subroutine get_general_2_to_1_derivs(i,ierr) ! e.g., r_he4_si28_to_o16_o16
            integer, intent(in) :: i
            integer, intent(out) :: ierr
            
            real(qp) :: b, b_f, b_r, rate
            real(dp) :: d, d_f, d_r, d1, d2, d3, Q, ys_f, ys_r, &
               d_ysf_dy1, d_ysf_dy2, d_ysr_dy3, y1, y2, y3
            integer :: c1, c2, c3, i1, i2, i3, o1, o3
            
            include 'formats'
         
            ierr = 0            
            
            ! forward reaction is c1 i1 + c2 i2 -> c3 i3
            c1 = reaction_inputs(1,ir)
            i1 = itab(reaction_inputs(2,ir))    
            c2 = reaction_inputs(3,ir)   
            i2 = itab(reaction_inputs(4,ir))    
            c3 = reaction_outputs(1,ir)
            i3 = itab(reaction_outputs(2,ir))

            if (symbolic) then
               n% d_dydt_dy(i1,i1) = 1
               n% d_dydt_dy(i1,i2) = 1
               n% d_dydt_dy(i1,i3) = 1
               n% d_dydt_dy(i2,i1) = 1
               n% d_dydt_dy(i2,i2) = 1
               n% d_dydt_dy(i2,i3) = 1
               n% d_dydt_dy(i3,i1) = 1
               n% d_dydt_dy(i3,i2) = 1
               n% d_dydt_dy(i3,i3) = 1
               return
            end if
            
            y1 = y(i1)
            y2 = y(i2)
            y3 = y(i3)
         
            select case(c1)
               case (1)
                  ys_f = y1
                  d_ysf_dy1 = 1d0
               case (2)
                  ys_f = y1*y1/2d0
                  d_ysf_dy1 = y1
               case (3)
                  ys_f = y1*y1*y1/6d0
                  d_ysf_dy1 = y1*y1/2d0
               case default
                  write(*,2) 'c1 too big for ' // trim(reaction_name(ir))
                  stop 'get_general_2_to_1_derivs'
            end select
            ! at this point, ys_f and d_ysf_dy1 only have y1 terms
            ! now combine with y2 info
            select case(c2)
               case (1)
                  d_ysf_dy2 = ys_f
                  ys_f = ys_f*y2
                  d_ysf_dy1 = d_ysf_dy1*y2
               case (2)
                  d_ysf_dy2 = ys_f*y2
                  ys_f = ys_f*(y2*y2/2d0)
                  d_ysf_dy1 = d_ysf_dy1*(y2*y2/2d0)
               case (3)
                  d_ysf_dy2 = ys_f*(y2*y2/2d0)
                  ys_f = ys_f*(y2*y2*y2/6d0)
                  d_ysf_dy1 = d_ysf_dy1*(y2*y2*y2/6d0)
               case default
                  write(*,2) 'c1 too big for ' // trim(reaction_name(ir))
                  stop 'get_general_2_to_1_derivs'
            end select
                        
            select case(c3)
               case (1)
                  ys_r = y3
                  d_ysr_dy3 = 1d0
               case (2)
                  ys_r = y3*y3/2d0
                  d_ysr_dy3 = y3
               case (3)
                  ys_r = y3*y3*y3/6d0
                  d_ysr_dy3 = y3*y3/2d0
               case default
                  write(*,2) 'c3 too big for ' // trim(reaction_name(ir))
                  stop 'get_general_2_to_1_derivs'
            end select

            d_f = ys_f
            d_r = ys_r

            rate = n% rate_screened(i)
            b_f = d_f*rate
            rate = n% rate_screened(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate,i1) = dydt(i_rate,i1) - c1*b
            dydt(i_rate,i2) = dydt(i_rate,i2) - c2*b
            dydt(i_rate,i3) = dydt(i_rate,i3) + c3*b
            if (just_dydt) return

            Q = n% reaction_Qs(ir)*eps_factor
            b_f = Q*b_f
            b_r = -Q*b_r
            b = b_f + b_r
            eps_nuc_MeV(i_rate) = eps_nuc_MeV(i_rate) + b
            n% eps_nuc_categories(icat_f) = n% eps_nuc_categories(icat_f) + b_f
            n% eps_nuc_categories(icat_r) = n% eps_nuc_categories(icat_r) + b_r
            if (show_eps_nuc ) &
               write(*,1) trim(reaction_Name(ir)) // ' eps_nuc',  b, b_f, b_r
                        
            rate = n% rate_screened_dT(i)
            b_f = d_f*rate
            rate = n% rate_screened_dT(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate_dT,i1) = dydt(i_rate_dT,i1) - c1*b
            dydt(i_rate_dT,i2) = dydt(i_rate_dT,i2) - c2*b
            dydt(i_rate_dT,i3) = dydt(i_rate_dT,i3) + c3*b
            b_f = Q*b_f
            b_r = -Q*b_r
            b = b_f + b_r
            eps_nuc_MeV(i_rate_dT) = eps_nuc_MeV(i_rate_dT) + b
                        
            rate = n% rate_screened_dRho(i)
            b_f = d_f*rate
            rate = n% rate_screened_dRho(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate_dRho,i1) = dydt(i_rate_dRho,i1) - c1*b
            dydt(i_rate_dRho,i2) = dydt(i_rate_dRho,i2) - c2*b
            dydt(i_rate_dRho,i3) = dydt(i_rate_dRho,i3) + c3*b
            b_f = Q*b_f
            b_r = -Q*b_r
            b = b_f + b_r
            eps_nuc_MeV(i_rate_dRho) = eps_nuc_MeV(i_rate_dRho) + b
               
            if (checking_deriv_flags) then
               deriv_flgs(i) = .true.
               deriv_flgs(r_i) = .true.
            end if

            d1 = d_ysf_dy1*n% rate_screened(i) ! d(rate_f)/d(y1)
            d2 = d_ysf_dy2*n% rate_screened(i) ! d(rate_f)/d(y2)
            d3 = d_ysr_dy3*n% rate_screened(r_i) ! d(rate_r)/d(y3)

            n% d_eps_nuc_dy(i1) = n% d_eps_nuc_dy(i1) + Q*d1
            n% d_eps_nuc_dy(i2) = n% d_eps_nuc_dy(i2) + Q*d2  
            n% d_eps_nuc_dy(i3) = n% d_eps_nuc_dy(i3) - Q*d3  

!           dydt(1,i1) = dydt(1,i1) - c1*(ys_f*n% rate_screened(i) - ys_r*n% rate_screened(r_i))
            n% d_dydt_dy(i1,i1) = n% d_dydt_dy(i1,i1) - c1*d1
            n% d_dydt_dy(i1,i2) = n% d_dydt_dy(i1,i2) - c1*d2
            n% d_dydt_dy(i1,i3) = n% d_dydt_dy(i1,i3) + c1*d3

!           dydt(1,i2) = dydt(1,i2) - c2*(ys_f*n% rate_screened(i) - ys_r*n% rate_screened(r_i))
            n% d_dydt_dy(i2,i1) = n% d_dydt_dy(i2,i1) - c2*d1
            n% d_dydt_dy(i2,i2) = n% d_dydt_dy(i2,i2) - c2*d2
            n% d_dydt_dy(i2,i3) = n% d_dydt_dy(i2,i3) + c2*d3

!           dydt(1,i3) = dydt(1,i3) + c3*(ys_f*n% rate_screened(i) - ys_r*n% rate_screened(r_i))
            n% d_dydt_dy(i3,i1) = n% d_dydt_dy(i3,i1) + c3*d1
            n% d_dydt_dy(i3,i2) = n% d_dydt_dy(i3,i2) + c3*d2
            n% d_dydt_dy(i3,i3) = n% d_dydt_dy(i3,i3) - c3*d3
            
            if (.false. .and. reaction_name(ir) == 'r_he4_si28_to_o16_o16') then ! .and. &
                  !y1 > 1d-20 .and. y2 > 1d-20 .and. y3 > 1d-20) then
               write(*,*)
               write(*,2) trim(reaction_name(ir))
               write(*,2) trim(reaction_name(r_ir))
               write(*,*)
               write(*,2) 'c1', c1
               write(*,2) 'c2', c2
               write(*,2) 'c3', c3
               write(*,*)
               write(*,1) 'd1', d1
               write(*,1) 'd2', d2
               write(*,*)
               write(*,1) 'dr1', d_ysf_dy1
               write(*,1) 'dr2', d_ysf_dy2
               write(*,*)
               write(*,1) 'd_ysf_dy1', d_ysf_dy1
               write(*,1) 'd_ysf_dy2', d_ysf_dy2
               write(*,*)
               write(*,1) 'y1',  y1
               write(*,1) 'y2',  y2
               write(*,1) 'y3',  y3
               write(*,*)
               write(*,1) 'd_dydt_dy(i1,i1)',  - c1*d1
               write(*,1) 'd_dydt_dy(i1,i2)',  - c1*d2
               write(*,1) 'd_dydt_dy(i1,i3)',  c1*d3
               write(*,*)
               write(*,1) 'd_dydt_dy(i2,i1)',  - c2*d1
               write(*,1) 'd_dydt_dy(i2,i2)',  - c2*d2
               write(*,1) 'd_dydt_dy(i2,i3)',  c2*d3
               write(*,*)
               write(*,1) 'd_dydt_dy(i3,i1)',  c3*d1
               write(*,1) 'd_dydt_dy(i3,i2)',  c3*d2
               write(*,1) 'd_dydt_dy(i3,i3)',  -c3*d3
               write(*,*)
               stop 'get_general_2_to_1_derivs'
            end if
                        
         end subroutine get_general_2_to_1_derivs

#ifdef offload
         !dir$ attributes offload: mic :: get_general_2_to_2_derivs
#endif
         subroutine get_general_2_to_2_derivs(i,ierr)
            integer, intent(in) :: i
            integer, intent(out) :: ierr
            
            real(qp) :: b, b_f, b_r, rate
            real(dp) :: d, d_f, d_r, d1, d2, d3, d4, Q, ys_f, ys_r, &
               d_ysf_dy1, d_ysf_dy2, d_ysr_dy3, d_ysr_dy4, y1, y2, y3, y4
            integer :: c1, c2, c3, c4, i1, i2, i3, i4
            
            include 'formats'
         
            ierr = 0            
            
            ! forward reaction is c1 i1 + c2 i2 -> c3 i3 + c4 i4
            c1 = reaction_inputs(1,ir)
            i1 = itab(reaction_inputs(2,ir))    
            c2 = reaction_inputs(3,ir)   
            i2 = itab(reaction_inputs(4,ir))    
            c3 = reaction_outputs(1,ir)
            i3 = itab(reaction_outputs(2,ir))
            c4 = reaction_outputs(3,ir)
            i4 = itab(reaction_outputs(4,ir))

            if (symbolic) then
               n% d_dydt_dy(i1,i1) = 1
               n% d_dydt_dy(i1,i2) = 1
               n% d_dydt_dy(i1,i3) = 1
               n% d_dydt_dy(i1,i4) = 1
               
               n% d_dydt_dy(i2,i1) = 1
               n% d_dydt_dy(i2,i2) = 1
               n% d_dydt_dy(i2,i3) = 1
               n% d_dydt_dy(i2,i4) = 1
               
               n% d_dydt_dy(i3,i1) = 1
               n% d_dydt_dy(i3,i2) = 1
               n% d_dydt_dy(i3,i3) = 1
               n% d_dydt_dy(i3,i4) = 1
               
               n% d_dydt_dy(i4,i1) = 1
               n% d_dydt_dy(i4,i2) = 1
               n% d_dydt_dy(i4,i3) = 1
               n% d_dydt_dy(i4,i4) = 1
               return
            end if
            
            y1 = y(i1)
            y2 = y(i2)
            y3 = y(i3)
            y4 = y(i4)
         
            select case(c1)
               case (1)
                  ys_f = y1
                  d_ysf_dy1 = 1d0
               case (2)
                  ys_f = y1*y1/2d0
                  d_ysf_dy1 = y1
               case (3)
                  ys_f = y1*y1*y1/6d0
                  d_ysf_dy1 = y1*y1/2d0
               case default
                  write(*,2) 'c1 too big for ' // trim(reaction_name(ir)), c1
                  stop 'get_general_2_to_2_derivs'
            end select
            ! at this point, ys_f and d_ysf_dy1 only have y1 terms
            ! now combine with y2 info
            select case(c2)
               case (1)
                  d_ysf_dy2 = ys_f
                  ys_f = ys_f*y2
                  d_ysf_dy1 = d_ysf_dy1*y2
               case (2)
                  d_ysf_dy2 = ys_f*y2
                  ys_f = ys_f*(y2*y2/2d0)
                  d_ysf_dy1 = d_ysf_dy1*(y2*y2/2d0)
               case (3)
                  d_ysf_dy2 = ys_f*(y2*y2/2d0)
                  ys_f = ys_f*(y2*y2*y2/6d0)
                  d_ysf_dy1 = d_ysf_dy1*(y2*y2*y2/6d0)
               case default
                  write(*,2) 'c2 too big for ' // trim(reaction_name(ir)), c2
                  stop 'get_general_2_to_2_derivs'
            end select
                        
            select case(c3)
               case (1)
                  ys_r = y3
                  d_ysr_dy3 = 1d0
               case (2)
                  ys_r = y3*y3/2d0
                  d_ysr_dy3 = y3
               case (3)
                  ys_r = y3*y3*y3/6d0
                  d_ysr_dy3 = y3*y3/2d0
               case default
                  write(*,2) 'c3 too big for ' // trim(reaction_name(ir)), c3
                  stop 'get_general_2_to_2_derivs'
            end select
            ! at this point, ys_r and d_ysf_dy3 only have y3 terms
            ! now combine with y4 info
            select case(c4)
               case (1)
                  d_ysr_dy4 = ys_r
                  ys_r = ys_r*y4
                  d_ysr_dy3 = d_ysr_dy3*y4
               case (2)
                  d_ysr_dy4 = ys_r*y4
                  ys_r = ys_r*(y4*y4/2d0)
                  d_ysr_dy3 = d_ysr_dy3*(y4*y4/2d0)
               case (3)
                  d_ysr_dy4 = ys_r*(y4*y4/2d0)
                  ys_r = ys_r*(y4*y4*y4/6d0)
                  d_ysr_dy3 = d_ysr_dy3*(y4*y4*y4/6d0)
               case default
                  write(*,2) 'c4 too big for ' // trim(reaction_name(ir))
                  stop 'get_general_2_to_2_derivs'
            end select

            d_f = ys_f
            d_r = ys_r

            rate = n% rate_screened(i)
            b_f = d_f*rate
            rate = n% rate_screened(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate,i1) = dydt(i_rate,i1) - c1*b
            dydt(i_rate,i2) = dydt(i_rate,i2) - c2*b
            dydt(i_rate,i3) = dydt(i_rate,i3) + c3*b
            dydt(i_rate,i4) = dydt(i_rate,i4) + c4*b
            if (just_dydt) return

            Q = n% reaction_Qs(ir)*eps_factor
            b_f = Q*b_f
            b_r = -Q*b_r
            b = b_f + b_r
            eps_nuc_MeV(i_rate) = eps_nuc_MeV(i_rate) + b
            n% eps_nuc_categories(icat_f) = n% eps_nuc_categories(icat_f) + b_f
            n% eps_nuc_categories(icat_r) = n% eps_nuc_categories(icat_r) + b_r
            if (show_eps_nuc ) &
               write(*,1) trim(reaction_Name(ir)) // ' eps_nuc',  b, b_f, b_r
                        
            rate = n% rate_screened_dT(i)
            b_f = d_f*rate
            rate = n% rate_screened_dT(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate_dT,i1) = dydt(i_rate_dT,i1) - c1*b
            dydt(i_rate_dT,i2) = dydt(i_rate_dT,i2) - c2*b
            dydt(i_rate_dT,i3) = dydt(i_rate_dT,i3) + c3*b
            dydt(i_rate_dT,i4) = dydt(i_rate_dT,i4) + c4*b
            b_f = Q*b_f
            b_r = -Q*b_r
            b = b_f + b_r
            eps_nuc_MeV(i_rate_dT) = eps_nuc_MeV(i_rate_dT) + b
                        
            rate = n% rate_screened_dRho(i)
            b_f = d_f*rate
            rate = n% rate_screened_dRho(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate_dRho,i1) = dydt(i_rate_dRho,i1) - c1*b
            dydt(i_rate_dRho,i2) = dydt(i_rate_dRho,i2) - c2*b
            dydt(i_rate_dRho,i3) = dydt(i_rate_dRho,i3) + c3*b
            dydt(i_rate_dRho,i4) = dydt(i_rate_dRho,i4) + c4*b
            b_f = Q*b_f
            b_r = -Q*b_r
            b = b_f + b_r
            eps_nuc_MeV(i_rate_dRho) = eps_nuc_MeV(i_rate_dRho) + b
               
            if (checking_deriv_flags) then
               deriv_flgs(i) = .true.
               deriv_flgs(r_i) = .true.
            end if

            d1 = d_ysf_dy1*n% rate_screened(i) ! d(rate_f)/d(y1)
            d2 = d_ysf_dy2*n% rate_screened(i) ! d(rate_f)/d(y2)
            d3 = d_ysr_dy3*n% rate_screened(r_i) ! d(rate_r)/d(y3)
            d4 = d_ysr_dy4*n% rate_screened(r_i) ! d(rate_r)/d(y4)

            n% d_eps_nuc_dy(i1) = n% d_eps_nuc_dy(i1) + Q*d1
            n% d_eps_nuc_dy(i2) = n% d_eps_nuc_dy(i2) + Q*d2  
            n% d_eps_nuc_dy(i3) = n% d_eps_nuc_dy(i3) - Q*d3  
            n% d_eps_nuc_dy(i4) = n% d_eps_nuc_dy(i4) - Q*d4  

            n% d_dydt_dy(i1,i1) = n% d_dydt_dy(i1,i1) - c1*d1
            n% d_dydt_dy(i1,i2) = n% d_dydt_dy(i1,i2) - c1*d2
            n% d_dydt_dy(i1,i3) = n% d_dydt_dy(i1,i3) + c1*d3
            n% d_dydt_dy(i1,i4) = n% d_dydt_dy(i1,i4) + c1*d4

            n% d_dydt_dy(i2,i1) = n% d_dydt_dy(i2,i1) - c2*d1
            n% d_dydt_dy(i2,i2) = n% d_dydt_dy(i2,i2) - c2*d2
            n% d_dydt_dy(i2,i3) = n% d_dydt_dy(i2,i3) + c2*d3
            n% d_dydt_dy(i2,i4) = n% d_dydt_dy(i2,i4) + c2*d4

            n% d_dydt_dy(i3,i1) = n% d_dydt_dy(i3,i1) + c3*d1
            n% d_dydt_dy(i3,i2) = n% d_dydt_dy(i3,i2) + c3*d2
            n% d_dydt_dy(i3,i3) = n% d_dydt_dy(i3,i3) - c3*d3
            n% d_dydt_dy(i3,i4) = n% d_dydt_dy(i3,i4) - c3*d4

            n% d_dydt_dy(i4,i1) = n% d_dydt_dy(i4,i1) + c4*d1
            n% d_dydt_dy(i4,i2) = n% d_dydt_dy(i4,i2) + c4*d2
            n% d_dydt_dy(i4,i3) = n% d_dydt_dy(i4,i3) - c4*d3
            n% d_dydt_dy(i4,i4) = n% d_dydt_dy(i4,i4) - c4*d4

         end subroutine get_general_2_to_2_derivs

#ifdef offload
         !dir$ attributes offload: mic :: get_basic_2_to_2_derivs
#endif
         subroutine get_basic_2_to_2_derivs(i,ierr)
            integer, intent(in) :: i
            integer, intent(out) :: ierr
            
            real(qp) :: b, b_f, b_r, rate
            real(dp) :: d_f, d_r, e_f, e_r, d1, d2, d3, d4, &
               ys_f, ys_r, Q, y1, y2, y3, y4
            integer :: i1, i2, i3, i4, o1, o2, o3, o4
            
            include 'formats'
         
            ierr = 0            
            
            ! forward reaction is i1 + i2 -> i3 + i4
            i1 = itab(reaction_inputs(2,ir))    
            i2 = itab(reaction_inputs(4,ir))    
            i3 = itab(reaction_outputs(2,ir))
            i4 = itab(reaction_outputs(4,ir))

            if (symbolic) then
               n% d_dydt_dy(i1,i1) = 1
               n% d_dydt_dy(i1,i2) = 1
               n% d_dydt_dy(i1,i3) = 1
               n% d_dydt_dy(i1,i4) = 1
               n% d_dydt_dy(i2,i1) = 1
               n% d_dydt_dy(i2,i2) = 1
               n% d_dydt_dy(i2,i3) = 1
               n% d_dydt_dy(i2,i4) = 1
               n% d_dydt_dy(i3,i1) = 1
               n% d_dydt_dy(i3,i2) = 1
               n% d_dydt_dy(i3,i3) = 1
               n% d_dydt_dy(i3,i4) = 1
               n% d_dydt_dy(i4,i1) = 1
               n% d_dydt_dy(i4,i2) = 1
               n% d_dydt_dy(i4,i3) = 1
               n% d_dydt_dy(i4,i4) = 1
               return
            end if
            
            y1 = y(i1)
            y2 = y(i2)
            y3 = y(i3)
            y4 = y(i4)
         
            ys_f = y1*y2
            d_f = ys_f
            
            ys_r = y3*y4
            d_r = ys_r

            rate = n% rate_screened(i)
            b_f = d_f*rate
            rate = n% rate_screened(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate,i1) = dydt(i_rate,i1) - b
            dydt(i_rate,i2) = dydt(i_rate,i2) - b
            dydt(i_rate,i3) = dydt(i_rate,i3) + b
            dydt(i_rate,i4) = dydt(i_rate,i4) + b
            if (just_dydt) return

            Q = n% reaction_Qs(ir)*eps_factor
            eps_nuc_MeV(i_rate) = eps_nuc_MeV(i_rate) + b*Q
            n% eps_nuc_categories(icat_f) = n% eps_nuc_categories(icat_f) + b_f*Q
            n% eps_nuc_categories(icat_r) = n% eps_nuc_categories(icat_r) - b_r*Q
            if (show_eps_nuc ) &
               write(*,1) trim(reaction_Name(ir)) // ' eps_nuc',  b, b_f, b_r
                        
            rate = n% rate_screened_dT(i)
            b_f = d_f*rate
            rate = n% rate_screened_dT(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate_dT,i1) = dydt(i_rate_dT,i1) - b
            dydt(i_rate_dT,i2) = dydt(i_rate_dT,i2) - b
            dydt(i_rate_dT,i3) = dydt(i_rate_dT,i3) + b
            dydt(i_rate_dT,i4) = dydt(i_rate_dT,i4) + b
            eps_nuc_MeV(i_rate_dT) = eps_nuc_MeV(i_rate_dT) + b*Q
                        
            rate = n% rate_screened_dRho(i)
            b_f = d_f*rate
            rate = n% rate_screened_dRho(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate_dRho,i1) = dydt(i_rate_dRho,i1) - b
            dydt(i_rate_dRho,i2) = dydt(i_rate_dRho,i2) - b
            dydt(i_rate_dRho,i3) = dydt(i_rate_dRho,i3) + b
            dydt(i_rate_dRho,i4) = dydt(i_rate_dRho,i4) + b
            eps_nuc_MeV(i_rate_dRho) = eps_nuc_MeV(i_rate_dRho) + b*Q
            
            if (checking_deriv_flags) then
               deriv_flgs(i) = .true.
               deriv_flgs(r_i) = .true.
            end if
            
            e_f = n% rate_screened(i)
            e_r = n% rate_screened(r_i)

            d1 = y2*e_f ! d(rate_f)/d(y1)
            d2 = y1*e_f ! d(rate_f)/d(y2)
            d3 = y4*e_r ! d(rate_r)/d(y3)
            d4 = y3*e_r ! d(rate_r)/d(y4)

            n% d_eps_nuc_dy(i1) = n% d_eps_nuc_dy(i1) + Q*d1
            n% d_eps_nuc_dy(i2) = n% d_eps_nuc_dy(i2) + Q*d2
            n% d_eps_nuc_dy(i3) = n% d_eps_nuc_dy(i3) - Q*d3      
            n% d_eps_nuc_dy(i4) = n% d_eps_nuc_dy(i4) - Q*d4     
            
            n% d_dydt_dy(i1,i1) = n% d_dydt_dy(i1,i1) - d1
            n% d_dydt_dy(i2,i1) = n% d_dydt_dy(i2,i1) - d1
            n% d_dydt_dy(i3,i1) = n% d_dydt_dy(i3,i1) + d1
            n% d_dydt_dy(i4,i1) = n% d_dydt_dy(i4,i1) + d1
            
            n% d_dydt_dy(i1,i2) = n% d_dydt_dy(i1,i2) - d2
            n% d_dydt_dy(i2,i2) = n% d_dydt_dy(i2,i2) - d2
            n% d_dydt_dy(i3,i2) = n% d_dydt_dy(i3,i2) + d2
            n% d_dydt_dy(i4,i2) = n% d_dydt_dy(i4,i2) + d2

            n% d_dydt_dy(i1,i3) = n% d_dydt_dy(i1,i3) + d3
            n% d_dydt_dy(i2,i3) = n% d_dydt_dy(i2,i3) + d3
            n% d_dydt_dy(i3,i3) = n% d_dydt_dy(i3,i3) - d3
            n% d_dydt_dy(i4,i3) = n% d_dydt_dy(i4,i3) - d3

            n% d_dydt_dy(i1,i4) = n% d_dydt_dy(i1,i4) + d4
            n% d_dydt_dy(i2,i4) = n% d_dydt_dy(i2,i4) + d4
            n% d_dydt_dy(i3,i4) = n% d_dydt_dy(i3,i4) - d4
            n% d_dydt_dy(i4,i4) = n% d_dydt_dy(i4,i4) - d4
            
         end subroutine get_basic_2_to_2_derivs

#ifdef offload
         !dir$ attributes offload: mic :: get_basic_2_to_1_derivs
#endif
         subroutine get_basic_2_to_1_derivs(i,ierr)
            integer, intent(in) :: i
            integer, intent(out) :: ierr
            
            real(qp) :: b, b_f, b_r, rate
            real(dp) :: d_f, d_r, e_f, e_r, d1, d2, d3, &
               ys_f, ys_r, Q, y1, y2, y3
            integer :: i1, i2, i3, o1, o2, o3, k
            
            include 'formats'
         
            ierr = 0            
            
            ! forward reaction is i1 + i2 -> i3
            i1 = itab(reaction_inputs(2,ir))    
            i2 = itab(reaction_inputs(4,ir))    
            i3 = itab(reaction_outputs(2,ir))
            
!            if (reaction_inputs(1,ir) /= 1 .or. &
!                reaction_inputs(3,ir) /= 1 .or. &
!                reaction_outputs(1,ir) /= 1 .or. &
!                reaction_inputs(6,ir) /= 0 .or. &
!                reaction_outputs(4,ir) /= 0) then
!               write(*,*) 'bad reaction for _ag_ ' // trim(reaction_name(ir))
!               stop
!            end if

            if (symbolic) then
               n% d_dydt_dy(i1,i1) = 1
               n% d_dydt_dy(i1,i2) = 1
               n% d_dydt_dy(i1,i3) = 1
               n% d_dydt_dy(i2,i1) = 1
               n% d_dydt_dy(i2,i2) = 1
               n% d_dydt_dy(i2,i3) = 1
               n% d_dydt_dy(i3,i1) = 1
               n% d_dydt_dy(i3,i2) = 1
               n% d_dydt_dy(i3,i3) = 1
               return
            end if

            y1 = y(i1)
            y2 = y(i2)
            y3 = y(i3)
         
            ys_f = y1*y2
            d_f = ys_f
            
            ys_r = y3
            d_r = ys_r

            rate = n% rate_screened(i)
            b_f = d_f*rate
            rate = n% rate_screened(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate,i1) = dydt(i_rate,i1) - b
            dydt(i_rate,i2) = dydt(i_rate,i2) - b
            dydt(i_rate,i3) = dydt(i_rate,i3) + b
            if (just_dydt) return

            Q = n% reaction_Qs(ir)*eps_factor
            eps_nuc_MeV(i_rate) = eps_nuc_MeV(i_rate) + b*Q
            n% eps_nuc_categories(icat_f) = n% eps_nuc_categories(icat_f) + b_f*Q
            n% eps_nuc_categories(icat_r) = n% eps_nuc_categories(icat_r) - b_r*Q
            if (show_eps_nuc ) &
               write(*,1) trim(reaction_Name(ir)) // ' eps_nuc',  b, b_f, b_r
                        
            rate = n% rate_screened_dT(i)
            b_f = d_f*rate
            rate = n% rate_screened_dT(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate_dT,i1) = dydt(i_rate_dT,i1) - b
            dydt(i_rate_dT,i2) = dydt(i_rate_dT,i2) - b
            dydt(i_rate_dT,i3) = dydt(i_rate_dT,i3) + b
            eps_nuc_MeV(i_rate_dT) = eps_nuc_MeV(i_rate_dT) + b*Q
                        
            rate = n% rate_screened_dRho(i)
            b_f = d_f*rate
            rate = n% rate_screened_dRho(r_i)
            b_r = d_r*rate
            b = b_f - b_r
            dydt(i_rate_dRho,i1) = dydt(i_rate_dRho,i1) - b
            dydt(i_rate_dRho,i2) = dydt(i_rate_dRho,i2) - b
            dydt(i_rate_dRho,i3) = dydt(i_rate_dRho,i3) + b
            eps_nuc_MeV(i_rate_dRho) = eps_nuc_MeV(i_rate_dRho) + b*Q
            
            if (checking_deriv_flags) then
               deriv_flgs(i) = .true.
               deriv_flgs(r_i) = .true.
            end if
            
            e_f = n% rate_screened(i)
            e_r = n% rate_screened(r_i)

            d1 = y2*e_f ! d(rate_f)/d(y1)
            d2 = y1*e_f ! d(rate_f)/d(y2)
            d3 = e_r ! d(rate_r)/d(y3)

            n% d_eps_nuc_dy(i1) = n% d_eps_nuc_dy(i1) + Q*d1
            n% d_eps_nuc_dy(i2) = n% d_eps_nuc_dy(i2) + Q*d2
            n% d_eps_nuc_dy(i3) = n% d_eps_nuc_dy(i3) - Q*d3   
            
            n% d_dydt_dy(i1,i1) = n% d_dydt_dy(i1,i1) - d1
            n% d_dydt_dy(i2,i1) = n% d_dydt_dy(i2,i1) - d1
            n% d_dydt_dy(i3,i1) = n% d_dydt_dy(i3,i1) + d1
            
            n% d_dydt_dy(i1,i2) = n% d_dydt_dy(i1,i2) - d2
            n% d_dydt_dy(i2,i2) = n% d_dydt_dy(i2,i2) - d2
            n% d_dydt_dy(i3,i2) = n% d_dydt_dy(i3,i2) + d2

            n% d_dydt_dy(i1,i3) = n% d_dydt_dy(i1,i3) + d3
            n% d_dydt_dy(i2,i3) = n% d_dydt_dy(i2,i3) + d3
            n% d_dydt_dy(i3,i3) = n% d_dydt_dy(i3,i3) - d3
               
         end subroutine get_basic_2_to_1_derivs

#ifdef offload
         !dir$ attributes offload: mic :: Check
#endif         
         subroutine Check
            integer :: nrates
            nrates = n% g% num_reactions
            
            do ir = 1, nrates
               if (.not. deriv_flgs(ir)) then
                  all_okay = .false.
                  write(*,'(a,i4,2x,a)') 'missing derivs for ', ir, &
                        trim(reaction_Name(g% reaction_id(ir)))
               end if
            end do
         
         end subroutine Check
         
      
      end subroutine get_derivs


      subroutine get1_derivs( &
            n, dydt, eps_nuc_MeV, i, eta, ye, temp, den, abar, zbar, &
            reuse_rate_screened, num_reactions, rate_factors, rtab, itab, &
            deriv_flgs, symbolic, just_dydt, ierr)
         use rates_lib, only: eval_n14_electron_capture_rate
         type (Net_Info), pointer :: n
         integer, intent(in) :: i, num_reactions
         real(qp), pointer, intent(out) :: dydt(:,:)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         logical, intent(in) :: reuse_rate_screened
         real(dp), intent(in) :: eta, ye, temp, den, abar, zbar, rate_factors(:)
         integer, pointer, intent(in) :: rtab(:), itab(:)
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         integer, intent(out) :: ierr
         
         integer :: ir, j, prot, neut, h1, he4, c14, n14, ne20, ne22, &
               mg21, mg22, mg23, mg24, al23, al24, si24, si25, si26,  &
               s28, s29, s30, cl31, ar32, ar33, ar34, k35, ca36, ca37, ca38, &
               ti41, ti42, v43, cr44, cr45, cr46, cr56, &
               fe48, fe49, fe50, fe51, fe52, fe54, fe56, fe58, fe60, fe62, fe64, &
               ni52, ni53, ni54, ni55, ni56, ni58, ni60, ni62, &
               zn57, zn58, zn59, zn60, ge62, ge63, ge64, &
               se68, kr72, sr76, mo84, sn104
         integer :: weak_id, num_reaction_inputs, in1, in2, in3, in4, in5
         integer :: cin1, cin2, cin3, cin4, cin5
         real(dp) :: din1, din2, din3, din4, din5
         integer :: num_reaction_outputs, out1, out2, out3, out4, out5
         integer :: cout1, cout2, cout3, cout4, cout5
         real(dp) :: dout1, dout2, dout3, dout4, dout5
         type (Net_General_Info), pointer  :: g
         integer, pointer :: reaction_id(:)
         real(dp), pointer :: y(:)
         integer :: i1, i2, i3, idr1, idr2, idr3, o1, o2, o3
         real(dp) :: r, dr1, dr2, dr3, rate, d_rate_dlnT, d_rate_dlnRho, Q, Qneu, rn14ec
         
         logical :: done, has_prot, has_neut, has_h1, switch_to_prot
         integer :: time0, time1, max_Z, Z_plus_N_for_max_Z
         integer, parameter :: min_Z_for_switch_to_prot = 12

         logical, parameter :: dbg = .false.
         
         include 'formats.dek'
         
         y => n% y
         g => n% g
         reaction_id => g% reaction_id

         ierr = 0
         prot = itab(iprot)
         neut = itab(ineut)
         h1 = itab(ih1)
         he4 = itab(ihe4)

         ir = reaction_id(i)         

         if (reaction_outputs(1,ir) == 0) return ! skip aux reactions
         
         if (dbg) &
            write(*,'(/,a,2i6)') ' reaction name <' // trim(reaction_Name(ir)) // '>', i, ir
         
         max_Z = g% reaction_max_Z(i)
         Z_plus_N_for_max_Z = g% reaction_max_Z_plus_N_for_max_Z(i)

         select case(ir)

            case(ir_he4_he4_he4_to_c12) ! triple alpha
               if (g% which_rates(ir) == use_rate_3a_FL87) then 
                  call do_FL_3alf(i) ! Fushiki and Lamb, Apj, 317, 368-388, 1987
                  return
               end if

            case(ir1212) ! c12 + c12 -> ne20 + he4; rate includes -> na23 + h1
               if (g% which_rates(ir) == use_rate_1212_G05) then 
                  ! use Gasques, et al. Phys Review C, 72, 025806 (2005)
                  call do_G05_1212_to_ne20(i)
                  return
               end if
            
            case(irn14ag_lite) ! n14 + 1.5 alpha => ne20
               n14 = itab(in14)
               ne20 = itab(ine20)
               r = y(n14) * y(he4)
               dr1 = y(n14)
               dr2 = y(he4)               
               call do_two_one(n, dydt, eps_nuc_MeV, i, 1d0, n14, 1.5d0, he4, &
                  r, 1d0, ne20, he4, dr1, n14, dr2, &
                  deriv_flgs, symbolic, just_dydt)
               return
         
         end select

         num_reaction_inputs = get_num_reaction_inputs(ir)
         num_reaction_outputs = get_num_reaction_outputs(ir)
         
         if (dbg) write(*,*) 'num_reaction_inputs', num_reaction_inputs
         if (dbg) write(*,*) 'num_reaction_outputs', num_reaction_outputs
         if (dbg) write(*,*)
         
         switch_to_prot = .false.
         cout1 = 0; out1 = 0; o1 = 0
         cout2 = 0; out2 = 0; o2 = 0
         cout3 = 0; out3 = 0; o3 = 0
         cin1 = 0; in1 = 0; i1 = 0
         cin2 = 0; in2 = 0; i2 = 0
         cin3 = 0; in3 = 0; i3 = 0

         if (num_reaction_outputs >= 1) then
            cout1 = reaction_outputs(1,ir); dout1 = cout1
            out1 = reaction_outputs(2,ir)
            o1 = itab(out1)
            if (o1 == 0) then
               write(*,*) trim(reaction_Name(ir))
               stop 'get1_derivs: itab(out1) = 0'
            end if
         end if
            
         if (num_reaction_outputs >= 2) then
            cout2 = reaction_outputs(3,ir); dout2 = cout2
            out2 = reaction_outputs(4,ir)
            o2 = itab(out2)
            if (o2 == 0) then
               write(*,*) trim(reaction_Name(ir))
               stop 'get1_derivs: itab(out2) = 0'
            end if
         end if
            
         if (num_reaction_outputs >= 3) then
            cout3 = reaction_outputs(5,ir); dout3 = cout3
            out3 = reaction_outputs(6,ir)
            o3 = itab(out3)
            if (o3 == 0) then
               write(*,*) trim(reaction_Name(ir))
               stop 'get1_derivs: itab(out3) = 0'
            end if
         end if
            
         if (num_reaction_outputs >= 4) then
            write(*,*) trim(reaction_Name(ir))
            stop 'get1_derivs: num_reaction_outputs >= 4'
         end if
         
         if (num_reaction_inputs == 1) then            
            cin1 = reaction_inputs(1,ir); din1 = cin1
            in1 = reaction_inputs(2,ir)
            i1 = itab(in1)            
         else if (num_reaction_inputs == 2 .or. num_reaction_inputs == 3) then            
            cin1 = reaction_inputs(1,ir); din1 = cin1
            in1 = reaction_inputs(2,ir)
            i1 = itab(in1)                        
            cin2 = reaction_inputs(3,ir); din2 = cin2
            in2 = reaction_inputs(4,ir)
            i2 = itab(in2)            
            if (num_reaction_inputs == 3) then
               cin3 = reaction_inputs(5,ir); din3 = cin3
               in3 = reaction_inputs(6,ir)
               i3 = itab(in3)
            end if
         end if
            
         switch_to_prot = (prot /= 0) .and. (max_Z >= min_Z_for_switch_to_prot)
         if (switch_to_prot) then
            if (i1 == h1) then
               in1 = iprot
               i1 = prot
            else if (i2 == h1) then
               in2 = iprot
               i2 = prot
            else if (i3 == h1) then
               in3 = iprot
               i3 = prot
            end if
            if (o1 == h1) then
               out1 = iprot
               o1 = prot
            else if (o2 == h1) then
               out2 = iprot
               o2 = prot
            else if (o3 == h1) then
               out3 = iprot
               o3 = prot
            end if
         end if

         if (num_reaction_inputs == 1) then
            
            if (i1 == 0) then
               write(*,*) trim(reaction_Name(ir))
               write(*,2) 'num_reaction_inputs', num_reaction_inputs
               stop 'get1_derivs: itab(in1) = 0'
            end if
            
            if (cin1 == 1) then
               r = y(i1)
               idr1 = i1
               dr1 = 1
            else if (cin1 == 3 .and. in1 /= ih1) then ! 3 he4
               !write(*,'(/,a)') '1/6*r  reaction name <' // trim(reaction_Name(ir)) // '>'
               r = (1d0/6d0)*y(i1)*y(i1)*y(i1)
               idr1 = i1
               dr1 = 0.5d0*y(i1)*y(i1)
            else ! 2 body
               !write(*,'(/,a)') '1/2*r  reaction name <' // trim(reaction_Name(ir)) // '>'
               !write(*,'(i3,3x,99e20.10)') i, n% rate_raw(i), n% rate_screened(i)
               r = 0.5d0*y(i1)*y(i1)
               idr1 = i1
               dr1 = y(i1)
               !stop
            end if
            
         else if (num_reaction_inputs == 2 .or. num_reaction_inputs == 3) then
                        
            if (reaction_ye_rho_exponents(2,ir) == 0) then
               ! treat as 1 body reaction
               r = y(i1)
               idr1 = i1
               dr1 = 1
               idr2 = i2
               dr2 = 0
               !write(*,*) 'get1_derivs rho=0: ' // trim(reaction_Name(ir))
               !stop 'net_derivs'
            else if ((cin1 == 1 .and. cin2 == 1) .or. reaction_ye_rho_exponents(2,ir) == 1) then
               ! reaction_ye_rho_exponents(2,ir) == 1 for electron captures; treat as 2 body reaction
               r = y(i1)*y(i2)
               dr1 = y(i1)
               idr1 = i2
               dr2 = y(i2)
               idr2 = i1
            else if (cin1 == 2 .and. cin2 == 1) then 
               r = 0.5d0*y(i1)*y(i1)*y(i2)
               dr1 = 0.5d0*y(i1)*y(i1)
               idr1 = i2
               dr2 = y(i1)*y(i2)
               idr2 = i1
            else if (cin1 == 1 .and. cin2 == 2) then 
               ! e.g., rhe4p, r_neut_he4_he4_to_be9, r_neut_h1_h1_to_h1_h2
               r = y(i1)*0.5d0*y(i2)*y(i2)
               dr1 = y(i1)*y(i2)
               idr1 = i2
               dr2 = 0.5d0*y(i2)*y(i2)
               idr2 = i1
            else if (cin1 == 2 .and. cin2 == 2) then 
               ! e.g., r_neut_neut_he4_he4_to_h3_li7, r_h1_h1_he4_he4_to_he3_be7
               r = 0.5d0*y(i1)*y(i1)*0.5d0*y(i2)*y(i2)
               dr1 = 0.5d0*y(i1)*y(i1)*y(i2)
               idr1 = i2
               dr2 = y(i1)*0.5d0*y(i2)*y(i2)
               idr2 = i1
            else
               write(*,*) 'get1_derivs: ' // trim(reaction_Name(ir)) // ' invalid coefficient'
               stop 'get1_derivs'
            end if            
            
            if (num_reaction_inputs == 3) then
               ! we assume that the 3rd kind of input is just "along for the ride"
               ! e.g., some compound reactions such as r34_pp2 are in this category.
               dr3 = 0
               idr3 = i3
               if (i3 == 0) then
                  write(*,*) trim(reaction_Name(ir))
                  stop 'get1_derivs: itab(in3) = 0'
               end if
            end if

         else

            write(*,*) 'get1_derivs: ' // trim(reaction_Name(ir)) // ' invalid specification'
            stop 'get1_derivs'

         end if
         
         if (num_reaction_inputs == 1) then
         
            if (num_reaction_outputs == 1) then 
               ! reaction of form din1 in1 -> dout1 out1
               if (dbg) write(*,*) ' do_one_one din1', din1, trim(chem_isos% name(g% chem_id(i1)))
               if (dbg) write(*,*) 'do_one_one dout1', dout1, trim(chem_isos% name(g% chem_id(o1)))
               
               weak_id = g% weak_reaction_index(i)
               
               done = .false.
               if (weak_id > 0) then
                  if (weak_id > g% num_weaklib_rates) then
                     write(*,2) 'i', i
                     write(*,2) 'ir', ir
                     write(*,2) 'weak_id', weak_id
                     write(*,2) 'g% num_weaklib_rates', g% num_weaklib_rates
                     write(*,*) trim(trim(reaction_Name(ir)))
                     stop 'derivs'
                  end if
                  if (g% weaklib_ids(weak_id) > 0) then ! > 0 means included in weaklib
                     
                     if (.not. reuse_rate_screened) then
                        n% rate_screened(i) = n% lambda(weak_id)
                        n% rate_screened_dT(i) = 0
                        n% rate_screened_dRho(i) = 0
                     end if
                     
                     call do_one_one_neu( &
                        n, dydt, eps_nuc_MeV, i, din1, i1, r, dout1, o1, &
                        idr1, dr1, n% Q(weak_id), n% Qneu(weak_id), &
                        deriv_flgs, symbolic, just_dydt)
                     done = .true.

                  end if
               end if
               if (.not. done) then ! weak reaction not in weaklib

                     if (.false. .and. trim(reaction_Name(ir)) == 'r_o15_wk_n15') then
                     write(*,*)
                     write(*,'(2(i5,x),a,2x,99e20.10)') i, ir, &
                        'do_one_one_neu ' // trim(reaction_Name(ir)) // ' ' // &
                        trim(chem_isos% name(g% chem_id(i1))) // ' => ' // &
                        trim(chem_isos% name(g% chem_id(o1)))
                     stop 'weak reaction not in weaklib'
                     end if

                  call do_one_one( & ! this will use reaction_neuQs &
                     n, dydt, eps_nuc_MeV, i, din1, i1, r, dout1, o1, &
                     idr1, dr1, deriv_flgs, symbolic, just_dydt)
               end if
     
            else if (num_reaction_outputs == 2) then
               ! reaction of form cin1 in1 -> dout1 out1 + dout2 out2
               if (dbg) write(*,*) ' do_one_two din1', din1, trim(chem_isos% name(g% chem_id(i1)))
               if (dbg) write(*,*) 'do_one_two dout1', dout1, trim(chem_isos% name(g% chem_id(o1)))
               if (dbg) write(*,*) 'do_one_two dout2', dout2, trim(chem_isos% name(g% chem_id(o2)))

               if (.false. .and. reaction_Name(ir) == 'r_he4_he4_he4_to_h1_b11' .and. r > 0) then
                  write(*,'(3i6,3x,a,2x,99e20.10)') i, ir, &
                     reaction_ye_rho_exponents(2,ir), &
                     'do_one_two ' // trim(reaction_Name(ir)) // ' ' // &
                     trim(chem_isos% name(g% chem_id(i1))) // ' => ' // &
                     trim(chem_isos% name(g% chem_id(o1))) // ' + ' // &
                     trim(chem_isos% name(g% chem_id(o2))), &
                     n% rate_screened(i), r, dr1, y(i1)
                  stop
               end if

               call do_one_two( &
                  n, dydt, eps_nuc_MeV, i, din1, i1, r, dout1, o1, dout2, o2, &
                  idr1, dr1, deriv_flgs, &
                  symbolic, just_dydt)
                  
            else if (num_reaction_outputs == 3) then
               ! reaction of form cin1 in1 -> dout1 out1 + dout2 out2 + dout3 out3
               if (dbg) write(*,*) ' do_one_three din1', din1, trim(chem_isos% name(g% chem_id(i1)))
               if (dbg) write(*,*) 'do_one_three dout1', dout1, trim(chem_isos% name(g% chem_id(o1)))
               if (dbg) write(*,*) 'do_one_three dout2', dout2, trim(chem_isos% name(g% chem_id(o2)))
               if (dbg) write(*,*) 'do_one_three dout3', dout3, trim(chem_isos% name(g% chem_id(o3)))
               call do_one_three( &
                  n, dydt, eps_nuc_MeV, i, din1, i1, r, dout1, o1, dout2, o2, dout3, o3, &
                  idr1, dr1, deriv_flgs, &
                  symbolic, just_dydt)
            else
               write(*,*) trim(reaction_Name(ir))
               stop 'too many reaction_outputs for num_reaction_inputs == 1'
            end if
            
         else if (num_reaction_inputs == 2) then
         
            if (num_reaction_outputs == 1) then 
               ! reaction of form din1 in1 + din2 in2 -> dout1 out1
               if (dbg) write(*,*) ' do_two_one din1', din1, trim(chem_isos% name(g% chem_id(i1)))
               if (dbg) write(*,*) ' do_two_one din2', din2, trim(chem_isos% name(g% chem_id(i2)))
               if (dbg) write(*,*) 'do_two_one dout1', dout1, trim(chem_isos% name(g% chem_id(o1)))
               
               if (.false. .and. reaction_Name(ir) == 'r_neut_he4_he4_to_be9' .and. r > 0 .and. &
                     abs(y(i1) - 7.7763751756339478D-05) < 1d-20) then
                  write(*,'(i3,3x,a,2x,99e20.10)') i, &
                     'do_two_one ' // trim(reaction_Name(ir)) // ' ' // &
                     trim(chem_isos% name(g% chem_id(i1))) // ' + ' // &
                     trim(chem_isos% name(g% chem_id(i2))) // ' => ' // &
                     trim(chem_isos% name(g% chem_id(o1))), &
                     r, dr1, dr2, y(i1), y(i2)
                  !stop
               end if
               
               if (.false. .and. reaction_Name(ir) == 'r_he4_si28_to_o16_o16') then
                  write(*,2) 'y(i1)', i1, y(i1)
                  write(*,2) 'y(i2)', i2, y(i2)
                  write(*,1) 'r', r
                  write(*,1) 'rate screened', n% rate_screened(i)
                  write(*,1) 'r*y1*y2', y(i1)*y(i2)*n% rate_screened(i)
                  !stop
               end if
               
               call do_two_one( &
                  n, dydt, eps_nuc_MeV, i, din1, i1, din2, i2, r, dout1, o1, &
                  idr1, dr1, idr2, dr2, deriv_flgs, &
                  symbolic, just_dydt)
            else if (num_reaction_outputs == 2) then
               ! reaction of form din1 in1 + din2 in2 -> dout1 out1 + dout2 out2
               if (dbg) write(*,*) ' do_two_two din1', din1, trim(chem_isos% name(g% chem_id(i1)))
               if (dbg) write(*,*) ' do_two_two din2', din2, trim(chem_isos% name(g% chem_id(i2)))
               if (dbg) write(*,*) 'do_two_two dout1', dout1, trim(chem_isos% name(g% chem_id(o1)))
               if (dbg) write(*,*) 'do_two_two dout2', dout2, trim(chem_isos% name(g% chem_id(o2)))
               call do_two_two( &
                  n, dydt, eps_nuc_MeV, i, din1, i1, din2, i2, r, dout1, o1, dout2, o2, &
                  idr1, dr1, idr2, dr2, deriv_flgs, &
                  symbolic, just_dydt)
            else if (num_reaction_outputs == 3) then
               ! reaction of form din1 in1 + din2 in2 -> dout1 out1 + dout2 out2 + dout3 out3
               if (dbg) write(*,*) ' do_two_three din1', din1, trim(chem_isos% name(g% chem_id(i1)))
               if (dbg) write(*,*) ' do_two_three din2', din2, trim(chem_isos% name(g% chem_id(i2)))
               if (dbg) write(*,*) 'do_two_three dout1', dout1, trim(chem_isos% name(g% chem_id(o1)))
               if (dbg) write(*,*) 'do_two_three dout2', dout2, trim(chem_isos% name(g% chem_id(o2)))
               if (dbg) write(*,*) 'do_two_three dout3', dout3, trim(chem_isos% name(g% chem_id(o3)))
               
               call do_two_three( &
                  n, dydt, eps_nuc_MeV, i, din1, i1, din2, i2, r, dout1, o1, dout2, o2, dout3, o3, &
                  idr1, dr1, idr2, dr2, deriv_flgs, &
                  symbolic, just_dydt)
            else
               write(*,*) trim(reaction_Name(ir))
               stop 'too many reaction_outputs for num_reaction_inputs == 2'
            end if
            
         else if (num_reaction_inputs == 3) then

            if (num_reaction_outputs == 1) then 
               ! reaction of form din1 in1 + din2 in2 + din3 in3 -> dout1 out1
               if (dbg) write(*,*) ' do_three_one din1', din1, trim(chem_isos% name(g% chem_id(i1)))
               if (dbg) write(*,*) ' do_three_one din2', din2, trim(chem_isos% name(g% chem_id(i2)))
               if (dbg) write(*,*) ' do_three_one din3', din3, trim(chem_isos% name(g% chem_id(i3)))
               if (dbg) write(*,*) 'do_three_one dout1', dout1, trim(chem_isos% name(g% chem_id(o1)))
               call do_three_one( &
                  n, dydt, eps_nuc_MeV, i, din1, i1, din2, i2, din3, i3, r, dout1, o1, &
                  idr1, dr1, idr2, dr2, deriv_flgs, &
                  symbolic, just_dydt)
            else if (num_reaction_outputs == 2) then
               ! reaction of form din1 in1 + din2 in2 + din3 in3 -> dout1 out1 + dout2 out2
               if (dbg) write(*,*) ' do_three_two din1', din1, trim(chem_isos% name(g% chem_id(i1)))
               if (dbg) write(*,*) ' do_three_two din2', din2, trim(chem_isos% name(g% chem_id(i2)))
               if (dbg) write(*,*) ' do_three_two din3', din3, trim(chem_isos% name(g% chem_id(i3)))
               if (dbg) write(*,*) 'do_three_two dout1', dout1, trim(chem_isos% name(g% chem_id(o1)))
               if (dbg) write(*,*) 'do_three_two dout2', dout2, trim(chem_isos% name(g% chem_id(o2)))
               call do_three_two( &
                  n, dydt, eps_nuc_MeV, i, din1, i1, din2, i2, din3, i3, r, dout1, o1, dout2, o2, &
                  idr1, dr1, idr2, dr2, idr3, dr3, deriv_flgs, &
                  symbolic, just_dydt)
            else
               write(*,*) trim(reaction_Name(ir))
               stop 'too many reaction_outputs for num_reaction_inputs == 3'
            end if
            
         else
            stop 'too many reaction_inputs'
         end if
               
         
         contains


#ifdef offload
         !dir$ attributes offload: mic :: do_FL_3alf
#endif         
         subroutine do_FL_3alf(i) ! Fushiki and Lamb, Apj, 317, 368-388, 1987
            use rates_lib, only: eval_FL_epsnuc_3alf
            integer, intent(in) :: i
            integer :: he4, c12
            real(dp) :: UE, XHe4, YHe4, &
                  FLeps_nuc, dFLeps_nuc_dT, dFLeps_nuc_dRho, r, drdT, drdRho, conv
            include 'formats.dek'
            he4 = itab(ihe4)
            c12 = itab(ic12)
            UE = abar/zbar
            YHe4 = y(he4)
            XHe4 = 4*YHe4
            if (YHe4 < 1d-50) then
               n% rate_screened(i) = 0
               n% rate_screened_dT(i) = 0
               n% rate_screened_dRho(i) = 0
               r = 0
               dr1 = 0
            else
               call eval_FL_epsnuc_3alf( &
                        temp, den, XHe4, UE, FLeps_nuc, dFLeps_nuc_dT, dFLeps_nuc_dRho)
               conv = Qconv*n% reaction_Qs(ir)
               r = YHe4*YHe4*YHe4/6d0
               dr1 = 0.5d0*YHe4*YHe4
               n% rate_screened(i) = FLeps_nuc/r*rate_factors(i)/conv
               n% rate_screened_dT(i) = dFLeps_nuc_dT/r*rate_factors(i)/conv
               n% rate_screened_dRho(i) = dFLeps_nuc_dRho/r*rate_factors(i)/conv    
            end if
            call do_one_one(n, dydt, eps_nuc_MeV, i, 3d0, he4, r, 1d0, c12, he4, dr1, &
                     deriv_flgs, symbolic, just_dydt)
         end subroutine do_FL_3alf
         
#ifdef offload
         !dir$ attributes offload: mic :: do_G05_1212_to_ne20
#endif         
         subroutine do_G05_1212_to_ne20(i) ! use Gasques, et al. Phys Review C, 72, 025806 (2005)
            use rates_lib, only: eval_G05_epsnuc_CC
            integer, intent(in) :: i
            integer :: he4, c12, ne20
            real(dp) :: X12, Y12, eps, deps_dT, deps_dRho, r, drdT, drdRho, conv
            include 'formats.dek'
            he4 = itab(ihe4)
            c12 = itab(ic12)
            ne20 = itab(ine20)
            Y12 = y(c12)
            X12 = 12*Y12
            if (Y12 < 1d-50) then
               n% rate_screened(i) = 0
               n% rate_screened_dT(i) = 0
               n% rate_screened_dRho(i) = 0
               r = 0
               dr1 = 0
            else
               call eval_G05_epsnuc_CC(temp, den, X12, eps, deps_dT, deps_dRho)
               conv = Qconv*n% reaction_Qs(ir)
               r = Y12*Y12/2d0
               dr1 = Y12
               n% rate_screened(i) = eps/r*rate_factors(i)/conv
               if (.false. .and. abs(temp - 1.8975458896321258D+08) < 1d-5) then
                  write(*,1) 'temp', temp
                  write(*,1) 'den', den
                  write(*,1) 'x c12', X12
                  write(*,1) 'eps', eps
                  write(*,1) 'rate', n% rate_raw(i)
                  stop 'net_derivs use_1212_G05'
               end if
               n% rate_screened_dT(i) = deps_dT/r*rate_factors(i)/conv
               n% rate_screened_dRho(i) = deps_dRho/r*rate_factors(i)/conv         
            end if
            call do_one_two(n, dydt, eps_nuc_MeV, i, 2d0, c12, r, &
               1d0, ne20, 1d0, he4, c12, dr1, deriv_flgs, &
               symbolic, just_dydt)            
         end subroutine do_G05_1212_to_ne20

          
      end subroutine get1_derivs


      subroutine eval_ni56_ec_rate( &
            temp, den, ye, eta, zbar, weak_rate_factor, &
            rate, d_rate_dlnT, d_rate_dlnRho, Q, Qneu, reuse_rate_screened, ierr)       
         real(dp), intent(in) :: temp, den, ye, eta, zbar, weak_rate_factor
         real(dp), intent(out) :: rate, d_rate_dlnT, d_rate_dlnRho, Q, Qneu
         logical, intent(in) :: reuse_rate_screened
		   integer, intent(out) :: ierr
		   include 'formats'
         call eval1_weak_rate( &
            weak_rate_id_for_ni56_ec, temp, ye, den, eta, zbar, &
            weak_rate_factor, rate, d_rate_dlnT, d_rate_dlnRho, Q, Qneu, &
            reuse_rate_screened, ierr)
      end subroutine eval_ni56_ec_rate


      subroutine eval_co56_ec_rate( &
            temp, den, ye, eta, zbar, weak_rate_factor, &
            rate, d_rate_dlnT, d_rate_dlnRho, Q, Qneu, reuse_rate_screened, ierr)       
         real(dp), intent(in) :: temp, den, ye, eta, zbar, weak_rate_factor
         real(dp), intent(out) :: rate, d_rate_dlnT, d_rate_dlnRho, Q, Qneu
         logical, intent(in) :: reuse_rate_screened
		   integer, intent(out) :: ierr
		   include 'formats'
         call eval1_weak_rate( &
            weak_rate_id_for_co56_ec, temp, ye, den, eta, zbar, &
            weak_rate_factor, rate, d_rate_dlnT, d_rate_dlnRho, Q, Qneu, &
            reuse_rate_screened, ierr)
      end subroutine eval_co56_ec_rate


      subroutine eval1_weak_rate( &
            id, temp, ye, rho, eta, zbar, weak_rate_factor, &
            rate_out, d_rate_dlnT, d_rate_dlnRho, Q_out, Qneu_out, &
            reuse_rate_screened, ierr)       
         use rates_def, only: Coulomb_Info
         use rates_lib, only: eval_weak_reaction_info
		   integer, intent(in) :: id
         real(dp), intent(in) :: temp, ye, rho, eta, zbar, weak_rate_factor
         real(dp), intent(out) :: &
            rate_out, d_rate_dlnT, d_rate_dlnRho, Q_out, Qneu_out
         logical, intent(in) :: reuse_rate_screened
		   integer, intent(out) :: ierr
         
		   integer :: ids(1)
           type(Coulomb_Info), pointer :: cc
           type(Coulomb_Info), target :: cc_info
		   real(dp) :: T9, YeRho, d_eta_dlnT, d_eta_dlnRho
		   ! lambda = combined rate = 10**ldecay + 10**lcapture
		   ! Q and Qneu are for combined rate of beta decay and electron capture.
		   ! Q is total, so Q-Qneu is the actual thermal energy.
   	   ! note: lambdas include Ye Rho factors for electron captures.
   	   ! so treat the rates as if just beta decays
		   real(dp), dimension(1), target :: &
            ldecay_a, d_ldecay_dT9_a, d_ldecay_dlYeRho_a, &
            lcapture_a, d_lcapture_dT9_a, d_lcapture_dlYeRho_a, &
            lneutrino_a, d_lneutrino_dT9_a, d_lneutrino_dlYeRho_a, &
            lambda_a, dlambda_dlnT_a, dlambda_dlnRho_a, &
            Q_a, dQ_dlnT_a, dQ_dlnRho_a, &
            Qneu_a, dQneu_dlnT_a, dQneu_dlnRho_a
		   real(dp), dimension(:), pointer :: &
            ldecay, d_ldecay_dT9, d_ldecay_dlYeRho, &
            lcapture, d_lcapture_dT9, d_lcapture_dlYeRho, &
            lneutrino, d_lneutrino_dT9, d_lneutrino_dlYeRho, &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho
         
         include 'formats'
         
         ierr = 0
         if (id <= 0) then
            ierr = -1
            return
         end if
		   
	      ldecay => ldecay_a
	      d_ldecay_dT9 => d_ldecay_dT9_a
	      d_ldecay_dlYeRho => d_ldecay_dlYeRho_a
	      
	      lcapture => lcapture_a
	      d_lcapture_dT9 => d_lcapture_dT9_a
	      d_lcapture_dlYeRho => d_lcapture_dlYeRho_a
	      
	      lneutrino => lneutrino_a
	      d_lneutrino_dT9 => d_lneutrino_dT9_a
	      d_lneutrino_dlYeRho => d_lneutrino_dlYeRho_a
	      
	      lambda => lambda_a
	      dlambda_dlnT => dlambda_dlnT_a
	      dlambda_dlnRho => dlambda_dlnRho_a
	      
	      Q => Q_a
	      dQ_dlnT => dQ_dlnT_a
	      dQ_dlnRho => dQ_dlnRho_a
	      
	      Qneu => Qneu_a
	      dQneu_dlnT => dQneu_dlnT_a
	      dQneu_dlnRho => dQneu_dlnRho_a
	      
	      ids(1) = id
	      T9 = temp*1d-9
	      YeRho = ye*rho
	      d_eta_dlnT = 0
	      d_eta_dlnRho = 0
          cc => cc_info        
         call eval_weak_reaction_info( &
            1, ids, cc, T9, YeRho, zbar, &
            eta, d_eta_dlnT, d_eta_dlnRho, &
            ldecay, d_ldecay_dT9, d_ldecay_dlYeRho, &
            lcapture, d_lcapture_dT9, d_lcapture_dlYeRho, &
            lneutrino, d_lneutrino_dT9, d_lneutrino_dlYeRho, &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho, &
            ierr)
            
         if (ierr /= 0) then
            return
            
            write(*,*) 'failed in eval_weak_reaction_info'
            stop 'eval1_weak_rate'
         end if
         
         if (.not. reuse_rate_screened) then
            rate_out = lambda(1)*weak_rate_factor
            d_rate_dlnT = dlambda_dlnT(1)*weak_rate_factor
            d_rate_dlnRho = dlambda_dlnRho(1)*weak_rate_factor
         end if
         
         Q_out = Q(1) 
         Qneu_out = Qneu(1) 
         
      end subroutine eval1_weak_rate

#ifdef offload
      !dir$ end options
#endif

      end module net_derivs













