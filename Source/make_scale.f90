module make_scale_module
  
  use bl_types

  implicit none

  private

  public :: make_scale
  
  contains

  subroutine make_scale(hq,q0)

    use bl_constants_module
    use geometry, only: spherical, nr_fine,r_cc_loc, dr, dr_fine, &
	  nlevs_radial, nr, numdisjointchunks, & 
	  r_start_coord, r_end_coord
    use probin_module, only: ref_ratio
    use restrict_base_module
    use parallel

    ! compute the q-scale height
    ! The base state uses 0-based indexing, so hq_cell 
    ! does too.
    
    real(kind=dp_t), intent(  out) :: hq(:,0:)
    real(kind=dp_t), intent(in   ) :: q0(:,0:)

    ! Local variables
    integer                      :: r, n, i
    real(kind=dp_t)              :: dq, ddr

    if (spherical .eq. 0) then

	! compute scale height as in the planar case
      n = 1

      do r=0,nr_fine-1
	  !forward difference
	  if (r == 0) then
	    dq = q0(1,r+1)-q0(1,r)
	    ddr= dr(1)
	  !backward difference
	  else if (r == nr_fine-1) then
	    dq = q0(1,r)-q0(1,r-1)
	    ddr= dr(1)
	  else 
	  !centered difference
	    dq = q0(1,r+1)-q0(1,r-1)
	    ddr= 2 *dr(1)	     
	  endif
	  
	  if (dq == ZERO) then
	    hq(1,r) = 1e30
	  else
	    hq(1,r) = -(ddr * q0(1,r))/dq
	  endif
      enddo

      do n = 2, nlevs_radial
	  do i=1,numdisjointchunks(n)
	    
	    ! forward difference
	    if (r_start_coord(n,i) .eq. 0) then
	      r = 0
	      dq = q0(n,1)-q0(n,0)
	      ddr = dr(n)
	    else 
		r = r_start_coord(n,i)
		dq = q0(n,r)-q0(n-1,r/ref_ratio)
		ddr = r_cc_loc(n,r)-r_cc_loc(n-1,r/ref_ratio)
	    end if
	    
	    if (dq == ZERO) then
		hq(n,r) = 1e30
	    else
		hq(n,r) = (ddr * q0(n,r))/dq
	    endif


	    do r=r_start_coord(n,i)+1,r_end_coord(n,i)
	      !backward difference
	      if (r .eq. r_end_coord(n,i)) then
		  dq = q0(n,r)-q0(n,r-1)
		  ddr= dr(n)                     
	      else
	      !centered difference
		dq = q0(n,r+1)-q0(n,r-1)
		  ddr= 2 *dr(n)	     
	      endif
		
		
	      if (dq == ZERO) then
		  hq(n,r) = 1e30
	      else
		  hq(n,r) = -(ddr * q0(n,r))/dq
	      endif		    
	    end do
	  enddo
      end do

      call restrict_base(hq,.true.)
      call fill_ghost_base(hq,.true.) 

    else  ! spherical = 1

       ! Computing only the finest level is enough for spherical symmetry,
       ! since we have to map to the 3D grid anyway, which is also possible with a finer hq          
       do r=0,nr_fine-1
          !forward difference
          if (r == 0) then
            dq = q0(1,r+1)-q0(1,r)
            ddr = dr(1)
          !backward difference
          else if (r == nr_fine-1) then
            dq = q0(1,r)-q0(1,r-1)
            ddr = dr(1)
          else 
          !centered difference
            dq = q0(1,r+1)-q0(1,r-1)
            ddr = 2*dr(1)
          endif
          
          if (dq == ZERO) then
           hq(1,r) = 1e30
          else
	   hq(1,r) = -(ddr * q0(1,r))/dq
	  endif
       enddo
    end if
    
  end subroutine make_scale

end module make_scale_module
  