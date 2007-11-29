module cell_to_edge_module

  use bl_types
  use bl_constants_module

  implicit none

  private
  public :: cell_to_edge, cell_to_edge_n

  contains

   subroutine cell_to_edge(s0_cell,s0_edge)

      real(kind=dp_t), intent(in   ) ::  s0_cell(0:)
      real(kind=dp_t), intent(inout) ::  s0_edge(0:)

      real(kind=dp_t)                ::  s0min,s0max
      integer                        ::  k,nr

      nr = size(s0_cell,dim=1)

      s0_edge( 0) = s0_cell(   0)
      s0_edge(nr) = s0_cell(nr-1)

      s0_edge(   1) = HALF * (s0_cell(   0) + s0_cell(   1))
      s0_edge(nr-1) = HALF * (s0_cell(nr-1) + s0_cell(nr-2))

      do k = 2, nr-2
         s0_edge(k) = 7.d0/12.d0 * (s0_cell(k  ) + s0_cell(k-1)) &
                     -1.d0/12.d0 * (s0_cell(k+1) + s0_cell(k-2))
         s0min = min(s0_cell(k),s0_cell(k-1))
         s0max = max(s0_cell(k),s0_cell(k-1))
         s0_edge(k) = max(s0_edge(k),s0min)
         s0_edge(k) = min(s0_edge(k),s0max)
      end do

   end subroutine cell_to_edge

   subroutine cell_to_edge_n(s0_cell,s0_edge)

      real(kind=dp_t), intent(in   ) ::  s0_cell(0:,:)
      real(kind=dp_t), intent(inout) ::  s0_edge(0:,:)
      integer                        ::  n

      do n = 1,size(s0_edge,dim=2)
        call cell_to_edge(s0_cell(:,n),s0_edge(:,n))
      end do

   end subroutine cell_to_edge_n

end module cell_to_edge_module
