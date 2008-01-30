module cell_to_edge_module

  use bl_types

  implicit none
  
  private

  public :: cell_to_edge, cell_to_edge_allcomps
  
contains
  
  subroutine cell_to_edge(n,s0_cell,s0_edge)

    use bl_constants_module
    use geometry, only: nr

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(in   ) :: s0_cell(0:)
    real(kind=dp_t), intent(inout) :: s0_edge(0:)
    
    real(kind=dp_t)                ::  s0min,s0max,tmp
    integer                        ::  r
    
    s0_edge(    0) = s0_cell(      0)
    s0_edge(nr(n)) = s0_cell(nr(n)-1)
    
    s0_edge(      1) = HALF * (s0_cell(      0) + s0_cell(      1))
    s0_edge(nr(n)-1) = HALF * (s0_cell(nr(n)-1) + s0_cell(nr(n)-2))
    
    do r = 2,nr(n)-2
       tmp = 7.d0/12.d0 * (s0_cell(r  ) + s0_cell(r-1)) &
            -1.d0/12.d0 * (s0_cell(r+1) + s0_cell(r-2))
       s0min      = min(s0_cell(r),s0_cell(r-1))
       s0max      = max(s0_cell(r),s0_cell(r-1))
       s0_edge(r) = min(max(tmp,s0min),s0max)
    end do
    
  end subroutine cell_to_edge
  
  subroutine cell_to_edge_allcomps(n,s0_cell,s0_edge)
    
    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(in   ) :: s0_cell(0:,:)
    real(kind=dp_t), intent(inout) :: s0_edge(0:,:)
    
    ! local
    integer                        ::  comp
    
    do comp = 1,size(s0_edge,dim=2)
       call cell_to_edge(n,s0_cell(:,comp),s0_edge(:,comp))
    end do
    
  end subroutine cell_to_edge_allcomps
  
end module cell_to_edge_module
