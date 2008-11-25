module react_base_module
  ! 
  ! This is React Base from paper III.
  ! 
  use bl_types

  implicit none
  
  private

  public :: react_base
  
contains
  
  subroutine react_base(rhoh0_in,rho_Hnucbar,rho_Hextbar,halfdt_in,rhoh0_out)

    use geometry, only: r_start_coord, r_end_coord, numdisjointchunks, nlevs_radial
    use bl_prof_module
    use restrict_base_module
     
    real(kind=dp_t), intent(in   ) :: rhoh0_in(:,0:)
    real(kind=dp_t), intent(in   ) :: rho_Hnucbar(:,0:)
    real(kind=dp_t), intent(in   ) :: rho_Hextbar(:,0:)
    real(kind=dp_t), intent(in   ) :: halfdt_in
    real(kind=dp_t), intent(  out) :: rhoh0_out(:,0:)
    
    integer :: n,r,i

    type(bl_prof_timer), save :: bpt

    call build(bpt, "react_base")
    
    do n=1,nlevs_radial
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             
             ! update enthalpy
             rhoh0_out(n,r) = rhoh0_in(n,r) &
                  + halfdt_in * rho_Hnucbar(n,r) + halfdt_in * rho_Hextbar(n,r)
             
          end do
       end do
    enddo ! end loop over levels

    call restrict_base(rhoh0_out,.true.)
    call fill_ghost_base(rhoh0_out,.true.)

    call destroy(bpt)

  end subroutine react_base

end module react_base_module
