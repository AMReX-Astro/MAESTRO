module make_hp_module
  
  use bl_types

  implicit none

  private

  public :: make_hp
  
  contains

  subroutine make_hp(hp,p0)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, dr_fine
    use restrict_base_module
    use parallel

    ! compute the pressure scale height
    ! The base state uses 0-based indexing, so hp_cell 
    ! does too.
    
    real(kind=dp_t), intent(  out) :: hp(:,0:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)

    ! Local variables
    integer                      :: r, n, i
    real(kind=dp_t)              :: dp, dr

    if (spherical .eq. 0) then

    else  ! spherical = 1

          
       do r=0,nr_fine-1
          !forward difference
          if (r == 0) then
            dp = p0(1,r+1)-p0(1,r)
            dr = dr_fine
          !backward difference
          else if (r == nr_fine-1) then
            dp = p0(1,r)-p0(1,r-1)
            dr = dr_fine
          else 
          !centered difference
            dp = p0(1,r+1)-p0(1,r-1)
            dr = 2*dr_fine
          endif
          
          if (dp == ZERO) then
           hp(1,r) = -huge(ZERO)
          else
	    hp(1,r) = (dr * p0(1,r))/dp
	  endif
       enddo
    end if
    
  end subroutine make_hp

end module make_hp_module
  