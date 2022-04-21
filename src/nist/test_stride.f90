program test_stride

   implicit none

   interface 
      subroutine get_stride(a,inca) bind(c,name='get_stride')
         real, intent(in) :: a(:)
         integer, intent(out) :: inca
      end subroutine
   end interface

   real, allocatable :: a(:)
   integer :: inca

   allocate(a(25))

   call get_stride(a,inca)
   print *, inca

   call get_stride(a(1:6:2),inca)
   print *, inca

end program