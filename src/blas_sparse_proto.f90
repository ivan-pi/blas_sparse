module blas_sparse_proto

   use blas_sparse_namedconstants, only: blas_conj_type, blas_trans_type
   
   implicit none
   public

   private :: blas_conj_type, blas_trans_type

   integer, parameter, private :: sp = kind(1.0e0)
   integer, parameter, private :: dp = kind(1.0d0)
!
! 3.8.2 LEVEL 1 COMPUTATIONAL ROUTINES
!

   !
   ! usdot (sparse dot product)
   !
   interface usdot
      module function susdot( x, indx, y, conj ) result(r)
         integer, intent(in) :: indx(:)
         real(sp), intent(in) :: x(:), y(:)
         integer(blas_conj_type), intent(in), optional :: conj
         real(sp) :: r
      end function susdot
      module function dusdot( x, indx, y, conj ) result(r)
         integer, intent(in) :: indx(:)
         real(dp), intent(in) :: x(:), y(:)
         integer(blas_conj_type), intent(in), optional :: conj
         real(dp) :: r
      end function dusdot
      module function cusdot( x, indx, y, conj ) result(r)
         integer, intent(in) :: indx(:)
         complex(sp), intent(in) :: x(:), y(:)
         integer(blas_conj_type), intent(in), optional :: conj
         complex(sp) :: r
      end function cusdot
      module function zusdot( x, indx, y, conj ) result(r)
         integer, intent(in) :: indx(:)
         complex(dp), intent(in) :: x(:), y(:)
         integer(blas_conj_type), intent(in), optional :: conj
         complex(dp) :: r
      end function zusdot
   end interface 

   !
   ! usaxpy (sparse vector update)
   !
   interface usaxpy
      module subroutine susaxpy( x, indx, y, alpha )
         real(sp), intent(in) :: x(:)
         real(sp), intent(inout) :: y(:)
         integer, intent(in) :: indx(:)
         real(sp), intent(in), optional :: alpha
      end subroutine susaxpy
      module subroutine dusaxpy( x, indx, y, alpha )
         real(dp), intent(in) :: x(:)
         real(dp), intent(inout) :: y(:)
         integer, intent(in) :: indx(:)
         real(dp), intent(in), optional :: alpha
      end subroutine dusaxpy
      module subroutine cusaxpy( x, indx, y, alpha )
         complex(sp), intent(in) :: x(:)
         complex(sp), intent(inout) :: y(:)
         integer, intent(in) :: indx(:)
         complex(sp), intent(in), optional :: alpha
      end subroutine cusaxpy
      module subroutine zusaxpy( x, indx, y, alpha )
         complex(dp), intent(in) :: x(:)
         complex(dp), intent(inout) :: y(:)
         integer, intent(in) :: indx(:)
         complex(dp), intent(in), optional :: alpha
      end subroutine zusaxpy
   end interface 

   !
   ! usga (sparse gather into compressed form)
   !
   interface usga
      module subroutine susga( y, x, indx )
         real(sp), intent(in) :: y(:)
         real(sp), intent(out) :: x(:)
         integer, intent(in) :: indx(:)
      end subroutine susga
      module subroutine dusga( y, x, indx )
         real(dp), intent(in) :: y(:)
         real(dp), intent(out) :: x(:)
         integer, intent(in) :: indx(:)
      end subroutine dusga
      module subroutine cusga( y, x, indx )
         complex(sp), intent(in) :: y(:)
         complex(sp), intent(out) :: x(:)
         integer, intent(in) :: indx(:)
      end subroutine cusga
      module subroutine zusga( y, x, indx )
         complex(dp), intent(in) :: y(:)
         complex(dp), intent(out) :: x(:)
         integer, intent(in) :: indx(:)
      end subroutine zusga
   end interface 

   !
   ! usgz (sparse gather and zero)
   !
   interface usgz
      module subroutine susgz( y, x, indx )
         real(sp), intent(in) :: y(:)
         real(sp), intent(out) :: x(:)
         integer, intent(in) :: indx(:)
      end subroutine susgz
      module subroutine dusgz( y, x, indx )
         real(dp), intent(in) :: y(:)
         real(dp), intent(out) :: x(:)
         integer, intent(in) :: indx(:)
      end subroutine dusgz
      module subroutine cusgz( y, x, indx )
         complex(sp), intent(in) :: y(:)
         complex(sp), intent(out) :: x(:)
         integer, intent(in) :: indx(:)
      end subroutine cusgz
      module subroutine zusgz( y, x, indx )
         complex(dp), intent(in) :: y(:)
         complex(dp), intent(out) :: x(:)
         integer, intent(in) :: indx(:)
      end subroutine zusgz
   end interface 

   !
   ! ussc (sparse scatter)
   !
   interface ussc
      module subroutine sussc( x, y, indx )
         real(sp), intent(in) :: x(:)
         real(sp), intent(inout) :: y(:)
         integer, intent(in) :: indx(:)
      end subroutine sussc
      module subroutine dussc( x, y, indx )
         real(dp), intent(in) :: x(:)
         real(dp), intent(inout) :: y(:)
         integer, intent(in) :: indx(:)
      end subroutine dussc
      module subroutine cussc( x, y, indx )
         complex(sp), intent(in) :: x(:)
         complex(sp), intent(inout) :: y(:)
         integer, intent(in) :: indx(:)
      end subroutine cussc
      module subroutine zussc( x, y, indx )
         complex(dp), intent(in) :: x(:)
         complex(dp), intent(inout) :: y(:)
         integer, intent(in) :: indx(:)
      end subroutine zussc
   end interface 

!
! 3.8.3 LEVEL 2 COMPUTATIONAL ROUTINES
!

   !
   ! usmv (sparse matrix/vector multiply)
   !
   interface usmv
      module subroutine susmv( a, x, y, istat, transa, alpha )
         integer, intent(in) :: a
         real(sp), intent(in) :: x(:)
         real(sp), intent(inout) :: y(:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transa
         real(sp), intent(in), optional :: alpha
      end subroutine susmv
      module subroutine dusmv( a, x, y, istat, transa, alpha )
         integer, intent(in) :: a
         real(dp), intent(in) :: x(:)
         real(dp), intent(inout) :: y(:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transa
         real(dp), intent(in), optional :: alpha
      end subroutine dusmv
      module subroutine cusmv( a, x, y, istat, transa, alpha )
         integer, intent(in) :: a
         complex(sp), intent(in) :: x(:)
         complex(sp), intent(inout) :: y(:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transa
         complex(sp), intent(in), optional :: alpha
      end subroutine cusmv
      module subroutine zusmv( a, x, y, istat, transa, alpha )
         integer, intent(in) :: a
         complex(dp), intent(in) :: x(:)
         complex(dp), intent(inout) :: y(:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transa
         complex(dp), intent(in), optional :: alpha
      end subroutine zusmv
   end interface 

   !
   ! ussv (sparse triangular solve)
   !
   interface ussv
      module subroutine sussv( t, x, istat, transt, alpha )
         integer, intent(in) :: t
         real(sp), intent(inout) :: x(:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transt
         real(sp), intent(in), optional :: alpha
      end subroutine sussv
      module subroutine dussv( t, x, istat, transt, alpha )
         integer, intent(in) :: t
         real(dp), intent(inout) :: x(:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transt
         real(dp), intent(in), optional :: alpha
      end subroutine dussv
      module subroutine cussv( t, x, istat, transt, alpha )
         integer, intent(in) :: t
         complex(sp), intent(inout) :: x(:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transt
         complex(sp), intent(in), optional :: alpha
      end subroutine cussv
      module subroutine zussv( t, x, istat, transt, alpha )
         integer, intent(in) :: t
         complex(dp), intent(inout) :: x(:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transt
         complex(dp), intent(in), optional :: alpha
      end subroutine zussv
   end interface 

!
! 3.8.4 LEVEL 3 COMPUTATIONAL ROUTINES
!

   !
   ! usmm
   !
   interface usmm
      module subroutine susmm( a, b, c, istat, transa, alpha )
         integer, intent(in) :: a
         real(sp), intent(in) :: b(:,:)
         real(sp), intent(inout) :: c(:,:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transa
         real(sp), intent(in), optional :: alpha
      end subroutine susmm
      module subroutine dusmm( a, b, c, istat, transa, alpha )
         integer, intent(in) :: a
         real(dp), intent(in) :: b(:,:)
         real(dp), intent(inout) :: c(:,:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transa
         real(dp), intent(in), optional :: alpha
      end subroutine dusmm
      module subroutine cusmm( a, b, c, istat, transa, alpha )
         integer, intent(in) :: a
         complex(sp), intent(in) :: b(:,:)
         complex(sp), intent(inout) :: c(:,:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transa
         complex(sp), intent(in), optional :: alpha
      end subroutine cusmm
      module subroutine zusmm( a, b, c, istat, transa, alpha )
         integer, intent(in) :: a
         complex(dp), intent(in) :: b(:,:)
         complex(dp), intent(inout) :: c(:,:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transa
         complex(dp), intent(in), optional :: alpha
      end subroutine zusmm
   end interface 

   !
   ! ussm
   !
   interface ussm
      module subroutine sussm( t, b, istat, transt, alpha )
         integer, intent(in) :: t
         real(sp), intent(inout) :: b(:,:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transt
         real(sp), intent(in), optional :: alpha
      end subroutine sussm
      module subroutine dussm( t, b, istat, transt, alpha )
         integer, intent(in) :: t
         real(dp), intent(inout) :: b(:,:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transt
         real(dp), intent(in), optional :: alpha
      end subroutine dussm
      module subroutine cussm( t, b, istat, transt, alpha )
         integer, intent(in) :: t
         complex(sp), intent(inout) :: b(:,:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transt
         complex(sp), intent(in), optional :: alpha
      end subroutine cussm
      module subroutine zussm( t, b, istat, transt, alpha )
         integer, intent(in) :: t
         complex(dp), intent(inout) :: b(:,:)
         integer, intent(out) :: istat
         integer(blas_trans_type), intent(in), optional :: transt
         complex(dp), intent(in), optional :: alpha
      end subroutine zussm
   end interface 

!
! 3.8.6 CREATION ROUTINES
!

   !
   ! uscr_begin
   !
   interface
   module subroutine suscr_begin( m, n, a, istat )
      integer, intent(in) :: m, n
      integer, intent(out) :: a, istat
   end subroutine suscr_begin
   module subroutine duscr_begin( m, n, a, istat )
      integer, intent(in) :: m, n
      integer, intent(out) :: a, istat
   end subroutine duscr_begin
   module subroutine cuscr_begin( m, n, a, istat )
      integer, intent(in) :: m, n
      integer, intent(out) :: a, istat
   end subroutine cuscr_begin
   module subroutine zuscr_begin( m, n, a, istat )
      integer, intent(in) :: m, n
      integer, intent(out) :: a, istat
   end subroutine zuscr_begin
   end interface

   !
   ! uscr_block_begin
   !
   interface
   module subroutine suscr_block_begin( mb, nb, k, l, a, istat )
      integer, intent(in) :: Mb, Nb, k, l
      integer, intent(out) :: a, istat
   end subroutine suscr_block_begin
   module subroutine duscr_block_begin( mb, nb, k, l, a, istat )
      integer, intent(in) :: Mb, Nb, k, l
      integer, intent(out) :: a, istat
   end subroutine duscr_block_begin
   module subroutine cuscr_block_begin( mb, nb, k, l, a, istat )
      integer, intent(in) :: Mb, Nb, k, l
      integer, intent(out) :: a, istat
   end subroutine cuscr_block_begin
   module subroutine zuscr_block_begin( mb, nb, k, l, a, istat )
      integer, intent(in) :: Mb, Nb, k, l
      integer, intent(out) :: a, istat
   end subroutine zuscr_block_begin
   end interface

   !
   ! uscr_variable_block_begin
   !
   interface
   module subroutine suscr_variable_block_begin( mb, nb, k, l, a, istat )
      integer, intent(in) :: mb, nb, k, l
      integer, intent(out) :: a, istat
   end subroutine suscr_variable_block_begin
   module subroutine duscr_variable_block_begin( mb, nb, k, l, a, istat )
      integer, intent(in) :: mb, nb, k, l
      integer, intent(out) :: a, istat
   end subroutine duscr_variable_block_begin
   module subroutine cuscr_variable_block_begin( mb, nb, k, l, a, istat )
      integer, intent(in) :: mb, nb, k, l
      integer, intent(out) :: a, istat
   end subroutine cuscr_variable_block_begin
   module subroutine zuscr_variable_block_begin( mb, nb, k, l, a, istat )
      integer, intent(in) :: mb, nb, k, l
      integer, intent(out) :: a, istat
   end subroutine zuscr_variable_block_begin
   end interface

!
! 3.8.7 INSERTION ROUTINES
!

   !> uscr_insert_entry 
   !! 
   !! insert single value at coordinate (i,j)
   !!
   interface uscr_insert_entry
   module subroutine suscr_insert_entry( a, val, i, j, istat )
      integer, intent(in) :: a, i, j
      real(sp), intent(in) :: val
      integer, intent(out) :: istat
   end subroutine suscr_insert_entry
   module subroutine duscr_insert_entry( a, val, i, j, istat )
      integer, intent(in) :: a, i, j
      real(dp), intent(in) :: val
      integer, intent(out) :: istat
   end subroutine duscr_insert_entry
   module subroutine cuscr_insert_entry( a, val, i, j, istat )
      integer, intent(in) :: a, i, j
      complex(sp), intent(in) :: val
      integer, intent(out) :: istat
   end subroutine cuscr_insert_entry
   module subroutine zuscr_insert_entry( a, val, i, j, istat )
      integer, intent(in) :: a, i, j
      complex(dp), intent(in) :: val
      integer, intent(out) :: istat
   end subroutine zuscr_insert_entry
   end interface

   !> uscr_insert_entries 
   !!
   !! insert a list of values in coordinate form (val,i,j)
   !!
   interface uscr_insert_entries
   module subroutine suscr_insert_entries( a, val, indx, jndx, istat )
      integer, intent(in) :: a, indx(:), jndx(:)
      real(sp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine suscr_insert_entries
   module subroutine duscr_insert_entries( a, val, indx, jndx, istat )
      integer, intent(in) :: a, indx(:), jndx(:)
      real(dp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine duscr_insert_entries
   module subroutine cuscr_insert_entries( a, val, indx, jndx, istat )
      integer, intent(in) :: a, indx(:), jndx(:)
      complex(sp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine cuscr_insert_entries
   module subroutine zuscr_insert_entries( a, val, indx, jndx, istat )
      integer, intent(in) :: a, indx(:), jndx(:)
      complex(dp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine zuscr_insert_entries
   end interface

   !> uscr_insert_col 
   !!
   !! insert a compressed column
   !!
   interface uscr_insert_col
   module subroutine suscr_insert_col( a, j, val, indx, istat )
      integer, intent(in) :: a, j, indx(:)
      real(sp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine suscr_insert_col
   module subroutine duscr_insert_col( a, j, val, indx, istat )
      integer, intent(in) :: a, j, indx(:)
      real(dp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine duscr_insert_col
   module subroutine cuscr_insert_col( a, j, val, indx, istat )
      integer, intent(in) :: a, j, indx(:)
      complex(sp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine cuscr_insert_col
   module subroutine zuscr_insert_col( a, j, val, indx, istat )
      integer, intent(in) :: a, j, indx(:)
      complex(dp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine zuscr_insert_col
   end interface

   !> uscr_insert_row
   !!
   !! insert a compressed row
   !!
   interface uscr_insert_row
   module subroutine suscr_insert_row( a, i, val, indx, istat )
      integer, intent(in) :: a, i, indx(:)
      real(sp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine suscr_insert_row
   module subroutine duscr_insert_row( a, i, val, indx, istat )
      integer, intent(in) :: a, i, indx(:)
      real(dp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine duscr_insert_row
   module subroutine cuscr_insert_row( a, i, val, indx, istat )
      integer, intent(in) :: a, i, indx(:)
      complex(sp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine cuscr_insert_row
   module subroutine zuscr_insert_row( a, i, val, indx, istat )
      integer, intent(in) :: a, i, indx(:)
      complex(dp), intent(in) :: val(:)
      integer, intent(out) :: istat
   end subroutine zuscr_insert_row
   end interface

   !> uscr_insert_clique
   !!
   !! insert a dense matrix clique
   !!
   interface uscr_insert_clique
   module subroutine suscr_insert_clique( a, val, indx, jndx, istat )
      integer, intent(in) :: a, indx(:), jndx(:)
      real(sp), intent(in) :: val(:,:)
      integer, intent(out) :: istat
   end subroutine suscr_insert_clique
   module subroutine duscr_insert_clique( a, val, indx, jndx, istat )
      integer, intent(in) :: a, indx(:), jndx(:)
      real(dp), intent(in) :: val(:,:)
      integer, intent(out) :: istat
   end subroutine duscr_insert_clique
   module subroutine cuscr_insert_clique( a, val, indx, jndx, istat )
      integer, intent(in) :: a, indx(:), jndx(:)
      complex(sp), intent(in) :: val(:,:)
      integer, intent(out) :: istat
   end subroutine cuscr_insert_clique
   module subroutine zuscr_insert_clique( a, val, indx, jndx, istat )
      integer, intent(in) :: a, indx(:), jndx(:)
      complex(dp), intent(in) :: val(:,:)
      integer, intent(out) :: istat
   end subroutine zuscr_insert_clique
   end interface

   !> uscr_insert_block
   !!
   !! insert a block entry at block coordinate (bi, bj)
   !!
   interface uscr_insert_block
   module subroutine suscr_insert_block( a, val, bi, bj, istat )
      integer, intent(in) :: a, bi, bj
      real(sp), intent(in) :: val(:,:)
      integer, intent(out) :: istat
   end subroutine suscr_insert_block
   module subroutine duscr_insert_block( a, val, bi, bj, istat )
      integer, intent(in) :: a, bi, bj
      real(dp), intent(in) :: val(:,:)
      integer, intent(out) :: istat
   end subroutine duscr_insert_block
   module subroutine cuscr_insert_block( a, val, bi, bj, istat )
      integer, intent(in) :: a, bi, bj
      complex(sp), intent(in) :: val(:,:)
      integer, intent(out) :: istat
   end subroutine cuscr_insert_block
   module subroutine zuscr_insert_block( a, val, bi, bj, istat )
      integer, intent(in) :: a, bi, bj
      complex(dp), intent(in) :: val(:,:)
      integer, intent(out) :: istat
   end subroutine zuscr_insert_block
   end interface

!
! 3.8.8 COMPLETION OF CONSTRUCTION ROUTINE
!

   !
   ! uscr_end (entries completed; build valid matrix handle)
   !
   interface
   module subroutine uscr_end( a, istat )
      integer, intent(in) :: a
      integer, intent(out) :: istat
   end subroutine uscr_end
   end interface

!
! 3.8.9 MATRIX PROPERTY ROUTINES
!

   !
   ! usgp (get matrix property)
   !
   interface
   module subroutine usgp( a, pname, m )
      integer, intent(in) :: a
      integer, intent(in) :: pname
      integer, intent(out) :: m
   end subroutine usgp
   end interface

   !
   ! ussp (set matrix property)
   !
   interface
   module subroutine ussp( a, pname, m )
      integer, intent(inout) :: a
      integer, intent(in) :: pname
      integer, intent(out) :: m
   end subroutine ussp
   end interface

!
! 3.8.10 DESTRUCTION ROUTINE
!

   !
   ! usds (release matrix handle)
   !
   interface
   module subroutine usds( a, istat )
      integer, intent(in) :: a
      integer, intent(out) :: istat
   end subroutine usds
   end interface

end module