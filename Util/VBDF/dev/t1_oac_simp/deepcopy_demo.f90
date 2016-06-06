!Make with
!$ pgf95 -acc -Minfo=acc -Mcuda=cuda7.0 -ta=nvidia:maxwell deepcopy_demo.f90

module my_mod
   type :: user_type
      real :: scalar
      real, allocatable :: mem_arr(0:)
   end type
end module my_mod

program deepcopy_demo
   use my_mod
   implicit none

   integer, parameter :: NCELLS = 10
   type(user_type)    :: ut(NCELLS)
   integer            :: i

   !$acc enter data create(ut)
   do i = 1,  NCELLS
      allocate(ut(i)%mem_arr(0:3))
      ut(i)%mem_arr = -1.0
      ut(i)%scalar  = 100.
       
      !$acc update device(ut(i)%scalar)
      !$acc enter data copyin(ut(i)%mem_arr(0:3))
   enddo
    
   !$acc parallel loop gang vector present(ut)
   do i = 1, NCELLS
      ut(i)%scalar     = 200
      ut(i)%mem_arr(0) = -200. !ut(i)%scalar
      ut(i)%mem_arr(1) = -300. !ut(i)%scalar
      ut(i)%mem_arr(2) = -400. !ut(i)%scalar
   end do
  
   print *, 'ut sc 1 b4 update: ', ut(1)%scalar
   print *, 'ut ma 1 b4 update: ', ut(1)%mem_arr
   !$acc update host(ut(1)%scalar)
   !$acc update host(ut(1)%mem_arr(0:3))
   print *, 'ut sc 1 af update: ', ut(1)%scalar
   print *, 'ut ma 1 af update: ', ut(1)%mem_arr
   print *, 'ut sc 2 af update: ', ut(2)%scalar
   print *, 'ut ma 2 af update: ', ut(2)%mem_arr
   do i=1, NCELLS
      !$acc exit data delete(ut(i)%mem_arr)
      deallocate(ut(i)%mem_arr)
   end do

   !$acc exit data delete(ut(:))
end program deepcopy_demo

