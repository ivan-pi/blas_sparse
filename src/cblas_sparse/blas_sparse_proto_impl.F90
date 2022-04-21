submodule (blas_sparse_proto) blas_sparse_proto_impl

   implicit none (type, external)

   integer, parameter :: blas_success = 0
   integer, parameter :: blas_failure = -1

contains

   module subroutine duscr_begin( m, n, a, istat )
      integer, intent(in) :: m, n
      integer, intent(out) :: a, istat

      interface
         integer(c_int) function blas_duscr_begin( m, n ) bind(c,name="BLAS_duscr_begin")
            use, intrinsic :: iso_c_binding, only: c_int
            integer(c_int), intent(in), value :: m, n
         end function blas_duscr_begin
      end interface   

      istat = blas_success
      
      a = blas_duscr_begin( m, n )
      
      if (a == blas_failure) then
         istat = blas_failure
      end if

   end subroutine duscr_begin

   module subroutine duscr_insert_entries( a, val, indx, jndx, istat )
      integer, intent(in) :: a, indx(:), jndx(:)
      real(dp), intent(in) :: val(:)
      integer, intent(out) :: istat

      interface
         integer(c_int) function blas_duscr_insert_entries( a, nnz, val, indx, jndx ) bind(c,name="BLAS_duscr_insert_entries")
            use, intrinsic :: iso_c_binding, only: c_double, c_int
            integer(c_int), intent(in), value :: a
            integer(c_int), intent(in), value :: nnz
            real(c_double), intent(in) :: val(*)
            integer(c_int), intent(in) :: indx(*), jndx(*)
         end function blas_duscr_insert_entries
      end interface

      integer :: nnz

      nnz = size(val)
      ! TODO: check sizes of val, indx, and jndx are matching
      !       if not, set istat and return   

      istat = blas_duscr_insert_entries( a, nnz, val, indx, jndx)

   end subroutine duscr_insert_entries

   module subroutine duscr_insert_row( a, i, val, indx, istat )
      integer, intent(in) :: a, i, indx(:)
      real(dp), intent(in) :: val(:)
      integer, intent(out) :: istat

      interface
         integer(c_int) function blas_duscr_insert_row( a, i, nnz, val, indx ) bind(c,name="BLAS_duscr_insert_row")
            use, intrinsic :: iso_c_binding, only: c_double, c_int
            integer(c_int), intent(in), value :: a
            integer(c_int), intent(in), value :: i, nnz
            real(c_double), intent(in) :: val(*)
            integer(c_int), intent(in) :: indx(*)
         end function blas_duscr_insert_row
      end interface

      integer :: nnz

      nnz = size(val)
      ! TODO: check sizes of val and indx are matching
      !       if not, set istat and return   

      istat = blas_duscr_insert_row( a, i, nnz, val, indx)

   end subroutine duscr_insert_row


   module subroutine uscr_end( a, istat )
      integer, intent(in) :: a
      integer, intent(out) :: istat

      interface
         integer(c_int) function blas_uscr_end( a ) bind(c,name="BLAS_uscr_end")
            use, intrinsic :: iso_c_binding, only: c_int
            integer(c_int), intent(in), value :: a
         end function
      end interface

      istat = blas_uscr_end( a )

   end subroutine uscr_end


   module subroutine usgp( a, pname, m )
      integer, intent(in) :: a
      integer, intent(in) :: pname
      integer, intent(out) :: m

      interface
         integer(c_int) function blas_usgp( a, pname ) bind(c,name="BLAS_usgp")
            use, intrinsic :: iso_c_binding, only: c_int
            integer(c_int), intent(in), value :: a
            integer(c_int), intent(in), value :: pname
         end function
      end interface      

      m = blas_usgp( a, pname )

   end subroutine usgp


   module subroutine ussp( a, pname, m )
      integer, intent(inout) :: a  ! TODO: wrong intent?
      integer, intent(in) :: pname
      integer, intent(out) :: m

      interface
         integer(c_int) function blas_ussp( a, pname ) bind(c,name="BLAS_ussp")
            use, intrinsic :: iso_c_binding, only: c_int
            integer(c_int), intent(in), value :: a
            integer(c_int), intent(in), value :: pname
         end function
      end interface      

      m = blas_ussp( a, pname )

   end subroutine ussp


   module subroutine usds( a, istat )
      integer, intent(in) :: a
      integer, intent(out) :: istat

      interface
         integer(c_int) function blas_usds( a ) bind(c,name="BLAS_usds")
            use, intrinsic :: iso_c_binding, only: c_int
            integer(c_int), intent(in), value :: a
         end function
      end interface

      istat = blas_usds( a )
   
   end subroutine usds


end submodule