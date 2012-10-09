module initialize_module

  use bl_constants_module
  use ml_boxarray_module

  implicit none

  private

  public :: initialize_dx

contains

  subroutine initialize_dx(dx,mba,num_levs)

    use probin_module, only: prob_lo, prob_hi

    real(dp_t)       , pointer     :: dx(:,:)
    type(ml_boxarray), intent(in ) :: mba
    integer          , intent(in ) :: num_levs
    
    integer :: n,d,dm

    dm = mba%dim
    
    allocate(dx(num_levs,dm))
    
    do d=1,dm
       dx(1,d) = (prob_hi(d)-prob_lo(d)) / real(extent(mba%pd(1),d),kind=dp_t)
    end do
    do n=2,num_levs
       dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
    end do

  end subroutine initialize_dx
  
end module initialize_module
