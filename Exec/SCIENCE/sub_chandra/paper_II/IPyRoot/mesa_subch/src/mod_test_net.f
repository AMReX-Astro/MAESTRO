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

!Annotated by Adam Jacobs
!
!This module is used by the driver program executed by the 'rn' script.
!Below I've annotated and explained what's going on.

      module mod_test_net
      
        implicit none
      
      contains
     
        !The main driving subroutine.  This is called by test_net.f.
        !quiet is a logical that determines if the test spits out output
        !or not. 
        subroutine test(quiet)
          !!! Modules
          !Includes MESA's chem, net, reaclib, weaklib, rates,
          !ecapture, screen, and const modules (see main mesa directory).
          !So even though this is the test for the net module, it requires
          !several other MESA modules.
          !
          !Locally, test_burn, test_burn_const_P, and test_net_support modules
          !are used.
          use chem_lib
          use net_def
          use net_lib
          use rates_lib, only: rates_init
          use const_lib
          use const_def, only: mesa_dir
          use test_net_support
     
          !!! Variable declaration, data, and initialization
          logical, intent(in) :: quiet
          character (len=64) :: my_mesa_dir
      
          integer :: ierr
      
          qt = quiet
      
          my_mesa_dir = '/home/ajacobs/Research/Codebase/mesa'
          call const_init(my_mesa_dir,ierr)     
          if (ierr /= 0) then
            write(*,*) 'const_init failed'
            stop 1
          end if        
      
          ierr = 0
          call chem_init('isotopes.data', ierr)
          if (ierr /= 0) then
            write(*,*) 'chem_init failed'
            stop 1
          end if
      
          call rates_init('reactions.list', '', ierr)
          if (ierr /= 0) stop 1
     
          !!! Execution 

          !Choose the rates to use, here setting a preference for the 
          !NACRE rates (Angulo, C., et al. 1999, Nucl. Phys. A, 656, 3) 
          !where available in the parameter space. We also establish the
          !temperature, density, and screening to be used.
          which_rates_choice = rates_NACRE_if_available
          
          !solar-like
          !test_logT = 7.5d0                   
          !test_logRho = 2d0                   
          
          !sub-Chandra 200 MK series
          !test_logT = 8.3d0  
          !test_logRho = 5d0  
          !test_logRho = 5.2d0
          !test_logRho = 5.4d0
          !test_logRho = 5.6d0
          !test_logRho = 5.8d0
          !test_logRho = 6.0d0
          !test_logRho = 6.2d0
          !test_logRho = 6.4d0
          !test_logRho = 6.6d0
          !test_logRho = 6.8d0
          !test_logRho = 7.0d0

          !sub-Chandra 300 MK series
          test_logT = 8.477d0                   
          !test_logT = 8.6d0                   
          !test_logRho = 5d0  
          !test_logRho = 5.2d0
          !test_logRho = 5.4d0
          test_logRho = 5.6d0
          !test_logRho = 5.8d0
          !test_logRho = 6.0d0
          !test_logRho = 6.2d0
          !test_logRho = 6.4d0
          !test_logRho = 6.6d0
          !test_logRho = 6.8d0
          !test_logRho = 7.0d0
          
          !test_logRho = 4.5d0

          screening_mode = extended_screening 
     
          !Below different networks are exercised and the results output. All of
          !the network files can be found in $MESA_DIR/data/net_data/nets
          
          !First we exercise the basic network.  It represents simple hydrogen
          !and helium burning. A test makes use of 
          !one_zone_burn.f's test_net_setup() and
          !test_net_support.f's do_test_net(), change_net() 
          !
          !The basic network includes h1, he3, he4 (a), c12, n14, o16, ne20, and mg24.
          !As for reactions, the network includes: 
          ! +PP chain: 
          !   p(p,e+nu)h2(p,g)he3 and p(e-p,nu)h2(p,g)he3
          !
          !   PPI, PPII, PPIII, and the weak reaction he3(p e+nu)he4
          !
          ! +CNO cycles 1 and 2
          !
          ! +Helium burning:
          !   Triple alpha: 
          !     he4(2he4 g)c12
          !   First bit of the alpha chain:
          !     c12(a g)o16 and c12(a p)n15(p g)o16
          !     o16(a g)ne20 and o16(a p)f19(p g)ne20
          !     ne20(a g)mg24 and ne20(a p)na23(p g)mg24
          !   Also some reaction I'm not familiar with:
          !     n14 + 1.5 alpha --> ne20
          !   
          !An excellent resource for (non-MESA-specific) reaction info and
          !details (including definition of some of the terms above) is 
          !Frank Timmes' website http://cococubed.asu.edu/code_pages/burn.shtml
          if (.not. qt) write(*,*) ' **************** Basic **************** '
      
          call test_net_setup('data/sub_chandra.net')
          call do_test_net(.false.,.false.)
      
          !!!Clean up 
          call test_net_cleanup
          call net_shutdown
          if (.not. qt) write(*,*)
      
        end subroutine test
      end module mod_test_net
