module mkl_sb

   use, intrinsic :: iso_fortran_env, only: int32
   use mkl_spblas

   implicit none
   private

   integer, parameter :: fsb_idx_ = int32

   ! TODO: ideally, a parameterized derived type of generic kind
   !       an alternative is a base matrix, and four children classes
   !       a further hierarchy based upon storage type is also possible
   !
   type :: fcsr_mat
      integer :: nrows, ncols, nnz
      integer(fsb_idx_), allocatable :: indx(:), jndx(:)
      real(dp), allocatable :: vals(:)
      integer :: params(20)
      logical :: ordered = .false.
   end type

   integer :: nmat = 0

   ! Currently, static array of matrix handles
   integer, parameter :: NMAX = 20

   ! The Fortran allocated matrices
   type(fcsr_mat), target :: cache_(NMAX)

   ! MKL matrix handles and descriptors
   type(sparse_matrix_t), target :: mkl_cache_(NMAX)
   type(matrix_descr), target :: mkl_descr_(NMAX)

   ! Initial number of nonzeros
   integer, parameter :: NNZ0 = 10000


contains

   !> Create a new handle
   integer function push_handle()
      nmat = nmat + 1
      push = nmat
   end function

   !> Remove a handle 
   subroutine pop_handle()
      nmat = nmat - 1
   end subroutine

   !> Check if a handle exists
   logical function is_handle_active(a)
      integer, intent(in) :: a
      is_active = 0 < a .and. a <= nmat
   end function

   module subroutine duscr_begin( m, n, a, istat )
      integer, intent(in) :: m, n
      integer, intent(out) :: a, istat

      a = push_handle()

      call fcsr_init(cache_(a),m,n,NNZ0,istat)

      ! TODO: if failed due to memory allocation
      !       pop handle

   end subroutine duscr_begin

   subroutine fcsr_init(mat,nrows,ncols,nnz,istat)
      type(fcsr_mat), intent(out) :: mat
      integer, intent(in) :: nrows, ncols, nnz
      integer, intent(out) :: istat

      integer :: alloc_err
      character(len=256) :: alloc_msg

      istat = 0

      mat%nrows = nrows
      mat%ncols = ncols
      mat%nnz = nnz

      allocate(mat%indx(nrows + 1), &
               mat%jndx(nnz), &
               mat%vals(nnz), &
               stat = alloc_err, &
               errmsg = alloc_msg)

      if (alloc_err /= 0) then
         write(*,*) "[fcsr_init] " // alloc_msg
         istat = -1
      end if
   end subroutine

   !> Release any memory held by the matrix
   subroutine fcsr_release(mat)
      ! Any allocatable members will be deallocated due to intent(out)
      type(fcsr_mat), intent(out) :: mat
      associate(n => mat%nrows)
         ! use associate to prevent false compiler warning
      end associate
   end subroutine

!
! 3.8.8 COMPLETION OF CONSTRUCTION ROUTINE
!

   !
   ! uscr_end (entries completed; build valid matrix handle)
   !
   module subroutine uscr_end( a, istat )
      integer, intent(in) :: a
      integer, intent(out) :: istat

      type(sparse_matrix_t) :: mkl_handle
      integer :: mkl_stat, prop
      integer :: indx_base

      ! TODO: check if `a` is a valid handle, if not return


      !
      ! Figure out indexing
      !
      indx_base = SPARSE_INDEX_BASE_ONE
      call usgp(a, pname, prop)
      if (prop == BLAS_ZERO_BASE) then
         indx_base = SPARSE_INDEX_BASE_ZERO
      end if

      ! TODO: pick correct routine, based on field type (real, complex)
      !       we should use a base class and dynamic polymorphism to 
      !       have a single cache_ of matrices, but of different type

      mkl_stat = mkl_sparse_d_create_csr( &
         mkl_handle, &
         indx_base, &
         cache_(a)%nrows, &
         cache_(a)%ncols, &
         cache_(a)%indx(1), &
         cache_(a)%indx(2), &
         cache_(a)%jndx(1), &
         cache_(a)%vals(2))

      if (mkl_stat /= SPARSE_STATUS_SUCCESS)
         ! handle error
      end if

      !
      ! Perform column ordering
      !
      mkl_stat = mkl_sparse_order(mkl_handle)

      if (mkl_stat == SPARSE_STATUS_SUCCESS)
         cache_(a)%ordered = .true.
      end if

      !
      ! Insert handle into cache_
      !
      mkl_cache_(a) = mkl_handle

   end subroutine uscr_end

!
! 3.8.9 MATRIX PROPERTY ROUTINES
!

   !
   ! usgp (get matrix property)
   !
   module subroutine usgp( a, pname, m )
      integer, intent(in) :: a
      integer, intent(in) :: pname
      integer, intent(out) :: m
   end subroutine usgp

   !
   ! ussp (set matrix property)
   !

   module subroutine ussp( a, pname, m )
      integer, intent(in) :: a
      integer, intent(in) :: pname
      integer, intent(out) :: m
   end subroutine ussp

!
! 3.8.10 DESTRUCTION ROUTINE
!

   !
   ! usds (release matrix handle)
   !
   subroutine usds(a,istat)
      integer, intent(in) :: a
      integer, intent(out) :: istat

      integer :: mkl_stat

      istat = 0

      mkl_stat = mkl_sparse_destroy(mkl_cache_(a))

      if (mkl_stat /= sparse_status_success) then
         ! handle error
         istat = -1
      end if

      call fcsr_release(cache_(a))

      call pop_handle()

   end subroutine

!
! 3.8.3 LEVEL 2 COMPUTATIONAL ROUTINES
!

   !
   ! usmv (sparse matrix/vector multiply)
   !
   module subroutine dusmv( a, x, y, istat, transa, alpha )

      use mkl_spblas, only: mkl_sparse_d_mv
      integer, intent(in) :: a
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: y(:)
      integer, intent(out) :: istat
      integer(blas_trans_type), intent(in), optional :: transa
      real(dp), intent(in), optional :: alpha

      real(dp) :: alpha_
      real(dp), parameter :: beta_ = 1.0_dp
      integer(c_int) :: op

      ! Assume success
      istat = 0

      ! Set optional alpha value
      alpha_ = 1.0_dp
      if (present(alpha)) then
         alpha_ = alpha
      end if

      ! Choose transpose operation
      op = SPARSE_OPERATION_NON_TRANSPOSE
      if (present(transa)) then
         select case(transa)
         case(blas_trans)
            op = SPARSE_OPERATION_TRANSPOSE
         case(blas_conj_trans)
            op = SPARSE_OPERATION_CONJUGATE_TRANSPOSE
         case default
            ! non-transpose
         end select
      end if
      ! TODO: recover from invalid input

      mkl_stat = mkl_sparse_d_mv( &
         op, &
         alpha_, &
         mkl_cache_(a), &
         mkl_descr_(a), &
         x, &
         beta_, &
         y)

      if (mkl_stat /= SPARSE_STATUS_SUCCESS) then
         istat = -1
      end if

   end subroutine dusmv

   !
   ! ussv (sparse triangular solve)
   !
   module subroutine dussv( t, x, istat, transt, alpha )
      integer, intent(in) :: t
      real(dp), intent(inout) :: x(:)
      integer, intent(out) :: istat
      integer(blas_trans_type), intent(in), optional :: transt
      real(dp), intent(in), optional :: alpha

      real(dp), allocatable :: y(:)

      real(dp) :: alpha_
      integer(c_int) :: op

      ! Assume success
      istat = 0

      ! Set optional alpha value
      alpha_ = 1.0_dp
      if (present(alpha)) then
         alpha_ = alpha
      end if

      ! Choose transpose operation
      op = SPARSE_OPERATION_NON_TRANSPOSE
      if (present(transa)) then
         select case(transa)
         case(blas_trans)
            op = SPARSE_OPERATION_TRANSPOSE
         case(blas_conj_trans)
            op = SPARSE_OPERATION_CONJUGATE_TRANSPOSE
         case default
            ! non-transpose
         end select
      end if

      ! TODO: recover from invalid input

      ! TODO: due to MKL API we need to create a copy,
      !       we should provide also an extension, where separate arrays
      !       are used for the RHS and solution vector. see also ussm
      y = x


      mkl_stat = mkl_sparse_d_trsv( &
         op, &
         alpha_, &
         mkl_cache_(t), &
         mkl_descr_(t), &
         y, &
         x)

      if (mkl_stat /= SPARSE_STATUS_SUCCESS) then
         istat = -1
      end if
      
   end subroutine dussv

!
! 3.8.4 LEVEL 3 COMPUTATIONAL ROUTINES
!

   !
   ! usmm
   !
   module subroutine dusmm( a, b, c, istat, transa, alpha )

      use mkl_spblas, only: mkl_sparse_d_mm
      integer, intent(in) :: a
      real(dp), intent(in) :: b(:,:)
      real(dp), intent(inout) :: c(:,:)
      integer, intent(out) :: istat
      integer(blas_trans_type), intent(in), optional :: transa
      real(dp), intent(in), optional :: alpha

      real(dp) :: alpha_
      real(dp), parameter :: beta_ = 1.0_dp
      integer(c_int) :: op, ncols, ldb, ldc

      ! Assume success
      istat = 0

      ! Set optional alpha value
      alpha_ = 1.0_dp
      if (present(alpha)) then
         alpha_ = alpha
      end if

      ! Choose transpose operation
      op = SPARSE_OPERATION_NON_TRANSPOSE
      if (present(transa)) then
         select case(transa)
         case(blas_trans)
            op = SPARSE_OPERATION_TRANSPOSE
         case(blas_conj_trans)
            op = SPARSE_OPERATION_CONJUGATE_TRANSPOSE
         case default
            ! non-transpose
         end select
      end if

      ! TODO: recover from invalid input

      ncols = size(c,2)
      ldb = size(b,1)
      ldc = size(c,1)

      ! TODO: check compatible dimensions

      mkl_stat = mkl_sparse_d_mm( &
         op, &
         alpha_, &
         mkl_cache_(a), &
         mkl_descr_(a), &
         SPARSE_LAYOUT_COLUMN_MAJOR,
         b, ncols, ldb, &
         beta_, &
         c, ldc)

      if (mkl_stat /= SPARSE_STATUS_SUCCESS) then
         istat = -1
      end if

   end subroutine dusmm

   !
   ! ussm
   !
   module subroutine dussm( t, b, istat, transt, alpha )
      integer, intent(in) :: t
      real(dp), intent(inout) :: b(:,:)
      integer, intent(out) :: istat
      integer(blas_trans_type), intent(in), optional :: transt
      real(dp), intent(in), optional :: alpha

      real(dp) :: alpha_
      integer(c_int) :: op, ncols, ldb, ldc

      real(dp), allocatable :: x(:,:)

      ! Assume success
      istat = 0

      ! Set optional alpha value
      alpha_ = 1.0_dp
      if (present(alpha)) then
         alpha_ = alpha
      end if

      ! Choose transpose operation
      op = SPARSE_OPERATION_NON_TRANSPOSE
      if (present(transt)) then
         select case(transt)
         case(blas_trans)
            op = SPARSE_OPERATION_TRANSPOSE
         case(blas_conj_trans)
            op = SPARSE_OPERATION_CONJUGATE_TRANSPOSE
         case default
            ! non-transpose
         end select
      end if

      ! TODO: recover from invalid input

      ncols = size(b,2)
      ldb = size(b,1)

      ! TODO: check compatible dimensions

      ! TODO: MKL doesn't overwrite the input matrix, so a copy needs to
      !       created internally, however if the user wants to preserve it,
      !       another copy will be made externally which is absurd. 
      !
      !       Ultimately, we should provide an extension for users who don't
      !       want to overwrite b

      ! Create copy of b
      ! via allocation on assignment
      x = b

      mkl_stat = mkl_sparse_d_trsm( &
         op, &
         alpha_, &
         mkl_cache_(t), &
         mkl_descr_(t), &
         SPARSE_LAYOUT_COLUMN_MAJOR,
         x, ncols, ldb, &
         beta_, &
         b, ldb)

      if (mkl_stat /= SPARSE_STATUS_SUCCESS) then
         istat = -1
      end if

   end subroutine dussm


!
!  ### Implementation-specific Extensions ###
!

   module subroutine dusmvdot( a, x, y, d, istat, transa, alpha )

      use mkl_spblas, only: mkl_sparse_d_dotmv
      integer, intent(in) :: a
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: y(:), d
      integer, intent(out) :: istat
      integer(blas_trans_type), intent(in), optional :: transa
      real(dp), intent(in), optional :: alpha

      real(dp) :: alpha_
      real(dp), parameter :: beta_ = 1.0_dp
      integer(c_int) :: op

      ! Assume success
      istat = 0

      ! Set optional alpha value
      alpha_ = 1.0_dp
      if (present(alpha)) then
         alpha_ = alpha
      end if

      ! Choose transpose operation
      op = SPARSE_OPERATION_NON_TRANSPOSE
      if (present(transa)) then
         select case(transa)
         case(blas_trans)
            op = SPARSE_OPERATION_TRANSPOSE
         case(blas_conj_trans)
            op = SPARSE_OPERATION_CONJUGATE_TRANSPOSE
         case default
            ! non-transpose
         end select
      end if

      ! TODO: recover from invalid input
      ! TODO: optionally, check matrix-vector sizes are compatible

      mkl_stat = mkl_sparse_d_dotmv( &
         op, &
         alpha_, &
         mkl_cache_(a), &
         mkl_descr_(a), &
         x, &
         beta_, &
         y, &
         d)

      if (mkl_stat /= SPARSE_STATUS_SUCCESS) then
         istat = -1
      end if

   end subroutine dusmv

end module