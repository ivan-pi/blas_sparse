module blas_sparse_namedconstants

   use, intrinsic :: iso_c_binding, only: c_int
   
   implicit none
   private :: c_int

   public

   integer, parameter :: blas_order_type = c_int
   enum, bind(c)
      enumerator :: &
         blas_rowmajor = 101, &
         blas_colmajor = 102
   end enum

   integer, parameter :: blas_trans_type = c_int
   enum, bind(c)
      enumerator :: &
         blas_no_trans   = 111, &
         blas_trans      = 112, &
         blas_conj_trans = 113
   end enum

   integer, parameter :: blas_uplo_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_upper = 121, &
         blas_lower = 122
   end enum

   integer, parameter :: blas_diag_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_non_unit_diag = 131, &
         blas_unit_diag     = 132
   end enum

   integer, parameter :: blas_side_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_left_side  = 141, &
         blas_right_side = 142
   end enum

   integer, parameter :: blas_cmach_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_base      = 151, &
         blas_t         = 152, &
         blas_rnd       = 153, &
         blas_ieee      = 154, &
         blas_emin      = 155, &
         blas_emax      = 156, &
         blas_eps       = 157, &
         blas_prec      = 158, &
         blas_underflow = 159, &
         blas_overflow  = 160, &
         blas_sfmin     = 161
   end enum

   integer, parameter :: blas_norm_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_one_norm       = 171, &
         blas_real_one_norm  = 172, &
         blas_two_norm       = 173, &
         blas_frobenius_norm = 174, &
         blas_inf_norm       = 175, &
         blas_real_inf_norm  = 176, &
         blas_max_norm       = 177, &
         blas_real_max_norm  = 178
   end enum

   integer, parameter :: blas_sort_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_increasing_order = 181, &
         blas_decreasing_order = 182
   end enum

   integer, parameter :: blas_conj_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_conj    = 191, &
         blas_no_conj = 192
   end enum

   integer, parameter :: blas_jrot_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_jrot_inner  = 201, &
         blas_jrot_outer  = 202, &
         blas_jrot_sorted = 203
   end enum

   integer, parameter :: blas_prec_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_prec_single     = 211, &
         blas_prec_double     = 212, &
         blas_prec_indigenous = 213, &
         blas_prec_extra      = 214
   end enum

   integer, parameter :: blas_base_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_zero_base = 221, &
         blas_one_base  = 222
   end enum

   integer, parameter :: blas_symmetry_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_general          = 231, &
         blas_symmetric        = 232, &
         blas_hermitian        = 233, &
         blas_triangular       = 234, &
         blas_lower_triangular = 235, &
         blas_upper_triangular = 236, &
         blas_lower_symmetric  = 237, &
         blas_upper_symmetric  = 238, &
         blas_lower_hermitian  = 239, &
         blas_upper_hermitian  = 240
   end enum

   integer, parameter :: blas_field_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_complex          = 241, &
         blas_real             = 242, &
         blas_double_precision = 243, &
         blas_single_precision = 244
   end enum

   integer, parameter :: blas_size_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_num_rows      = 251, &
         blas_num_cols      = 252, &
         blas_num_nonzeros  = 253
   end enum

   integer, parameter :: blas_handle_type = c_int
   enum, bind(c)
      enumerator :: &
         blas_invalid_handle = 261, &
         blas_new_handle     = 262, &
         blas_open_handle    = 263, &
         blas_valid_handle   = 264
   end enum

   integer, parameter :: blas_sparsity_optimization_type = c_int 
   enum, bind(c)
      enumerator :: &
         blas_regular       = 271, &
         blas_irregular     = 272, &
         blas_block         = 273, &
         blas_unassembled   = 274
   end enum

end module
