#include "ISO_Fortran_binding.h"

#include "blas_enum.h"
#include "blas_sparse_proto.h"

#define INC(X,d) X->dim[d].sm / X->elem_len
#define LD(X) X->dim[0].extent
 
void susmv( const int* a, 
            const CFI_cdesc_t *x,
            CFI_cdesc_t *y
            int *istat
            enum blas_trans_type transa,
            const float *alpha) {

   enum blas_trans_type transa_ = blas_no_trans;
   if (transa)
      transa_ = transa;

   float alpha_ = 1.0;
   if (alpha) 
      alpha_ = alpha;

   float *x_ = (float *) x->base_addr;
   int incx = INC(x,0);
   
   float *y_ = (float *) y->base_addr;
   int incy = INC(y,0);

   istat = BLAS_susmv( transa_, alpha_, a, x_, incx, y_, incy );

}


void sussv( const int* t, 
            CFI_cdesc_t *x,
            int *istat
            enum blas_trans_type transt,
            const float *alpha) {

   enum blas_trans_type transt_ = blas_no_trans;
   if (transt)
      transt_ = transt;

   float alpha_ = 1.0;
   if (alpha) 
      alpha_ = alpha;

   float *x_ = (float *) x->base_addr;
   int incx = INC(x,0);

   istat = BLAS_sussv( transt_, alpha_, t, x_, incx);

}



void susmm( const int* a, 
            const CFI_cdesc_t *b,
            CFI_cdesc_t *c
            int *istat
            enum blas_trans_type transa,
            const float *alpha) {

   enum blas_trans_type transa_ = blas_no_trans;
   if (transa)
      transa_ = transa;

   float alpha_ = 1.0;
   if (alpha) 
      alpha_ = alpha;

   float *b_ = (float *) b->base_addr;
   int ldb = LD(b);
   int nrhs = b->dim[1].extent;
   
   float *c_ = (float *) c->base_addr;
   int ldc = LD(c);

   istat = BLAS_susmm( blas_col_major, transa_, alpha_, a, b_, ldb, c_, ldc );

}


void sussm( const int* t, 
            CFI_cdesc_t *b,
            int *istat
            enum blas_trans_type transt,
            const float *alpha) {

   enum blas_trans_type transt_ = blas_no_trans;
   if (transt)
      transt_ = transt;

   float alpha_ = 1.0;
   if (alpha) 
      alpha_ = alpha;

   float *b_ = (float *) b->base_addr;
   int ldb = LD(b);
   int nrhs = b->dim[1].extent;

   istat = BLAS_susmm( blas_col_major, transt_, nrhs, alpha_, t, b_, ldb );

}

// build valid matrix handle
void uscr_end( const int *a, int *istat) {
   istat = BLAS_uscr_end( a );
}

// get matrix property
void usgp( const int *a, const int *pname, int *istat) {
   istat = BLAS_usgp( a, pname );
}

// set matrix property
void ussp( int *a, const int *pname, int *istat) {
   istat = BLAS_ussp( a, pname );
}

// usds (release matrix handle)
void usds( const int *a, int *istat) {
   istat = BLAS_usds( a );
}