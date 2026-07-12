#ifndef JXF_BLAS_HEADER
#define JXF_BLAS_HEADER

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * Prototypes
 *--------------------------------------------------------------------------*/

/* dasum.c */
JXF_Real
jxf_dasum( JXF_Int *n, JXF_Real *dx, JXF_Int *incx );

/* daxpy.c */
JXF_Int
jxf_daxpy( JXF_Int *n, JXF_Real *da, JXF_Real *dx, JXF_Int *incx, JXF_Real *dy, JXF_Int *incy );

/* dcopy.c */
JXF_Int
jxf_dcopy( JXF_Int *n, JXF_Real *dx, JXF_Int *incx, JXF_Real *dy, JXF_Int *incy );

/* ddot.c */
JXF_Real
jxf_ddot( JXF_Int *n, JXF_Real *dx, JXF_Int *incx, JXF_Real *dy, JXF_Int *incy );

/* dgemm.c */
JXF_Int
jxf_dgemm( const char *transa, const char *transb, JXF_Int *m, JXF_Int *n, JXF_Int *k,
JXF_Real *alpha, JXF_Real *a, JXF_Int *lda, JXF_Real *b, JXF_Int *ldb, JXF_Real *beta, JXF_Real *c, JXF_Int *ldc );

/* dgemv.c */
JXF_Int
jxf_dgemv( const char *trans, JXF_Int *m, JXF_Int *n, JXF_Real *alpha,
JXF_Real *a, JXF_Int *lda, JXF_Real *x, JXF_Int *incx, JXF_Real *beta, JXF_Real *y, JXF_Int *incy );

/* dger.c */
JXF_Int
jxf_dger( JXF_Int *m, JXF_Int *n, JXF_Real *alpha, JXF_Real *x, JXF_Int *incx, JXF_Real *y, JXF_Int *incy, JXF_Real *a, JXF_Int *lda );

/* dnrm2.c */
JXF_Real
jxf_dnrm2( JXF_Int *n, JXF_Real *dx, JXF_Int *incx );

/* drot.c */
JXF_Int
jxf_drot( JXF_Int *n, JXF_Real *dx, JXF_Int *incx, JXF_Real *dy, JXF_Int *incy, JXF_Real *c, JXF_Real *s );

/* dscal.c */
JXF_Int
jxf_dscal( JXF_Int *n, JXF_Real *da, JXF_Real *dx, JXF_Int *incx );

/* dswap.c */
JXF_Int
jxf_dswap( JXF_Int *n, JXF_Real *dx, JXF_Int *incx, JXF_Real *dy, JXF_Int *incy );

/* dsymm.c */
JXF_Int
jxf_dsymm( const char *side, const char *uplo, JXF_Int *m, JXF_Int *n, JXF_Real *alpha,
JXF_Real *a, JXF_Int *lda, JXF_Real *b, JXF_Int *ldb, JXF_Real *beta, JXF_Real *c__, JXF_Int *ldc );

/* dsymv.c */
JXF_Int
jxf_dsymv( const char *uplo, JXF_Int *n, JXF_Real *alpha, JXF_Real *a,
JXF_Int *lda, JXF_Real *x, JXF_Int *incx, JXF_Real *beta, JXF_Real *y, JXF_Int *incy );

/* dsyr2.c */
JXF_Int
jxf_dsyr2( const char *uplo, JXF_Int *n, JXF_Real *alpha, JXF_Real *x, JXF_Int *incx, JXF_Real *y, JXF_Int *incy, JXF_Real *a, JXF_Int *lda );

/* dsyr2k.c */
JXF_Int
jxf_dsyr2k( const char *uplo, const char *trans, JXF_Int *n, JXF_Int *k, JXF_Real *alpha,
JXF_Real *a, JXF_Int *lda, JXF_Real *b, JXF_Int *ldb, JXF_Real *beta, JXF_Real *c__, JXF_Int *ldc );

/* dsyrk.c */
JXF_Int
jxf_dsyrk( const char *uplo, const char *trans, JXF_Int *n, JXF_Int *k, JXF_Real *alpha,
JXF_Real *a, JXF_Int *lda, JXF_Real *beta, JXF_Real *c, JXF_Int *ldc );

/* dtrmm.c */
JXF_Int
jxf_dtrmm( const char *side, const char *uplo, const char *transa, const char *diag,
JXF_Int *m, JXF_Int *n, JXF_Real *alpha, JXF_Real *a, JXF_Int *lda, JXF_Real *b, JXF_Int *ldb );

/* dtrmv.c */
JXF_Int
jxf_dtrmv( const char *uplo, const char *trans, const char *diag, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *x, JXF_Int *incx );

/* dtrsm.c */
JXF_Int
jxf_dtrsm( const char *side, const char *uplo, const char *transa, const char *diag,
JXF_Int *m, JXF_Int *n, JXF_Real *alpha, JXF_Real *a, JXF_Int *lda, JXF_Real *b, JXF_Int *ldb );

/* dtrsv.c */
JXF_Int
jxf_dtrsv( const char *uplo, const char *trans, const char *diag, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *x, JXF_Int *incx );

/* idamax.c */
JXF_Int
jxf_idamax( JXF_Int *n, JXF_Real *dx, JXF_Int *incx );

#ifdef __cplusplus
}
#endif

#endif
