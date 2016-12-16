module make_hp_module
  
  use bl_types

  implicit none

  private

  public :: make_hp
  
  contains

  subroutine make_hp(hp,p0)

    use bl_constants_module
    use geometry, only: spherical, nr_fine,r_cc_loc, dr, dr_fine, &
	  nlevs_radial, nr, numdisjointchunks, & 
	  r_start_coord, r_end_coord
    use probin_module, only: do_2d_planar_octant, ref_ratio
    use restrict_base_module
    use parallel

    ! compute the pressure scale height
    ! The base state uses 0-based indexing, so hp_cell 
    ! does too.
    
    real(kind=dp_t), intent(  out) :: hp(:,0:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)

    ! Local variables
    integer                      :: r, n, i
    real(kind=dp_t)              :: dp, ddr

    if (spherical .eq. 0) then
      if (do_2d_planar_octant .eq. 1) then

          ! compute pressure scale height as in the planar case
	n = 1

	do r=0,nr_fine-1
	    !forward difference
	    if (r == 0) then
	      dp = p0(1,r+1)-p0(1,r)
	      ddr= dr(1)
	    !backward difference
	    else if (r == nr_fine-1) then
	      dp = p0(1,r)-p0(1,r-1)
	      ddr= dr(1)
	    else 
	    !centered difference
	      dp = p0(1,r+1)-p0(1,r-1)
	      ddr= 2 *dr(1)	     
	    endif
	    
	    if (dp == ZERO) then
	      hp(1,r) = 1e30
	    else
	      hp(n,r) = -(ddr * p0(1,r))/dp
	    endif
	enddo

          do n = 2, nlevs_radial
             do i=1,numdisjointchunks(n)
                
                ! forward difference
                if (r_start_coord(n,i) .eq. 0) then
                  r = 0
		  dp = p0(n,1)-p0(n,0)
		  ddr = dr(n)
		else 
		   r = r_start_coord(n,i)
		   dp = p0(n,r)-p0(n-1,r/ref_ratio)
		   ddr = r_cc_loc(n,r)-r_cc_loc(n-1,r/ref_ratio)
                end if
		
		if (dp == ZERO) then
		    hp(n,r) = -huge(ZERO)
		else
		    hp(n,r) = (ddr * p0(n,r))/dp
		endif


                do r=r_start_coord(n,i)+1,r_end_coord(n,i)
                  !backward difference
                  if (r .eq. r_end_coord(n,i)) then
		     dp = p0(n,r)-p0(n,r-1)
		     ddr= dr(n)                     
                  else
                  !centered difference
		    dp = p0(n,r+1)-p0(n,r-1)
		     ddr= 2 *dr(n)	     
                  endif
                   
                   
		  if (dp == ZERO) then
		      hp(n,r) = 1e30
		  else
		      hp(n,r) = -(ddr * p0(n,r))/dp
		  endif		    
                end do
             enddo
          end do

          call restrict_base(hp,.true.)
          call fill_ghost_base(hp,.true.)  
       endif

    else  ! spherical = 1

          
       do r=0,nr_fine-1
          !forward difference
          if (r == 0) then
            dp = p0(1,r+1)-p0(1,r)
            ddr = dr(1)
          !backward difference
          else if (r == nr_fine-1) then
            dp = p0(1,r)-p0(1,r-1)
            ddr = dr(1)
          else 
          !centered difference
            dp = p0(1,r+1)-p0(1,r-1)
            ddr = 2*dr(1)
          endif
          
          if (dp == ZERO) then
           hp(1,r) = -huge(ZERO)
          else
	    hp(1,r) = -(ddr * p0(1,r))/dp
	  endif
       enddo
    end if
    
  end subroutine make_hp

end module make_hp_module
  