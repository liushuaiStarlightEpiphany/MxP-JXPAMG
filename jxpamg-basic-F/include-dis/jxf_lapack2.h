#ifndef JXF_LAPACK_HEADER
#define JXF_LAPACK_HEADER

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * Prototypes
 *--------------------------------------------------------------------------*/

/* dbdsqr.c */
JXF_Int
jxf_dbdsqr(const char *uplo, JXF_Int *n, JXF_Int *ncvt, JXF_Int *nru, JXF_Int *ncc,
JXF_Real *d__, JXF_Real *e, JXF_Real *vt, JXF_Int *ldvt, JXF_Real *u, JXF_Int *ldu, JXF_Real *c__, JXF_Int *ldc, JXF_Real *work, JXF_Int *info );

/* dgebd2.c */
JXF_Int
jxf_dgebd2( JXF_Int *m, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *d__,
JXF_Real *e, JXF_Real *tauq, JXF_Real *taup, JXF_Real *work, JXF_Int *info );

/* dgebrd.c */
JXF_Int
jxf_dgebrd( JXF_Int *m, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *d__,
JXF_Real *e, JXF_Real *tauq, JXF_Real *taup, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dgelq2.c */
JXF_Int
jxf_dgelq2( JXF_Int *m, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *info );

/* dgelqf.c */
JXF_Int
jxf_dgelqf( JXF_Int *m, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dgels.c */
JXF_Int
jxf_dgels( char *trans, JXF_Int *m, JXF_Int *n, JXF_Int *nrhs, JXF_Real *a,
JXF_Int *lda, JXF_Real *b, JXF_Int *ldb, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dgeqr2.c */
JXF_Int
jxf_dgeqr2( JXF_Int *m, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *info );

/* dgeqrf.c */
JXF_Int
jxf_dgeqrf( JXF_Int *m, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dgesvd.c */
JXF_Int
jxf_dgesvd( char *jobu, char *jobvt, JXF_Int *m, JXF_Int *n, JXF_Real *a, JXF_Int *lda,
JXF_Real *s, JXF_Real *u, JXF_Int *ldu, JXF_Real *vt, JXF_Int *ldvt, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dgetf2.c */
JXF_Int
jxf_dgetf2( JXF_Int *m, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Int *ipiv, JXF_Int *info );

/* dgetrf.c */
JXF_Int
jxf_dgetrf( JXF_Int *m, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Int *ipiv, JXF_Int *info );

/* dgetri.c */
JXF_Int
jxf_dgetri( JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Int *ipiv, JXF_Real *work, JXF_Int *lwork, JXF_Int *info);

/* dgetrs.c */
JXF_Int
jxf_dgetrs( const char *trans, JXF_Int *n, JXF_Int *nrhs, JXF_Real *a, JXF_Int *lda, JXF_Int *ipiv, JXF_Real *b, JXF_Int *ldb, JXF_Int *info );

/* dlasq1.c */
JXF_Int
jxf_dlasq1( JXF_Int *n, JXF_Real *d__, JXF_Real *e, JXF_Real *work, JXF_Int *info );

/* dlasq2.c */
JXF_Int
jxf_dlasq2( JXF_Int *n, JXF_Real *z__, JXF_Int *info );

/* dlasrt.c */
JXF_Int
jxf_dlasrt(const char *id, JXF_Int *n, JXF_Real *d__, JXF_Int *info );

/* dorg2l.c */
JXF_Int
jxf_dorg2l( JXF_Int *m, JXF_Int *n, JXF_Int *k, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *info );

/* dorg2r.c */
JXF_Int
jxf_dorg2r( JXF_Int *m, JXF_Int *n, JXF_Int *k, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *info );

/* dorgbr.c */
JXF_Int
jxf_dorgbr(const char *vect, JXF_Int *m, JXF_Int *n, JXF_Int *k, JXF_Real *a,
JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dorgl2.c */
JXF_Int
jxf_dorgl2( JXF_Int *m, JXF_Int *n, JXF_Int *k, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *info );

/* dorglq.c */
JXF_Int
jxf_dorglq( JXF_Int *m, JXF_Int *n, JXF_Int *k, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dorgql.c */
JXF_Int
jxf_dorgql( JXF_Int *m, JXF_Int *n, JXF_Int *k, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dorgqr.c */
JXF_Int
jxf_dorgqr( JXF_Int *m, JXF_Int *n, JXF_Int *k, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dorgtr.c */
JXF_Int
jxf_dorgtr(const char *uplo, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dorm2r.c */
JXF_Int
jxf_dorm2r(const char *side,const char *trans, JXF_Int *m, JXF_Int *n, JXF_Int *k,
JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *c__, JXF_Int *ldc, JXF_Real *work, JXF_Int *info );

/* dormbr.c */
JXF_Int
jxf_dormbr(const char *vect,const char *side,const char *trans, JXF_Int *m,
JXF_Int *n, JXF_Int *k, JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *c__, JXF_Int *ldc, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dorml2.c */
JXF_Int
jxf_dorml2(const char *side,const char *trans, JXF_Int *m, JXF_Int *n, JXF_Int *k,
JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *c__, JXF_Int *ldc, JXF_Real *work, JXF_Int *info );

/* dormlq.c */
JXF_Int
jxf_dormlq(const char *side,const char *trans, JXF_Int *m, JXF_Int *n, JXF_Int *k,
JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *c__, JXF_Int *ldc, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dormqr.c */
JXF_Int
jxf_dormqr(const char *side,const char *trans, JXF_Int *m, JXF_Int *n, JXF_Int *k,
JXF_Real *a, JXF_Int *lda, JXF_Real *tau, JXF_Real *c__, JXF_Int *ldc, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dpotf2.c */
JXF_Int
jxf_dpotf2(const char *uplo, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Int *info );

/* dpotrf.c */
JXF_Int
jxf_dpotrf(const char *uplo, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Int *info );

/* dpotrs.c */
JXF_Int
jxf_dpotrs( char *uplo, JXF_Int *n, JXF_Int *nrhs, JXF_Real *a, JXF_Int *lda, JXF_Real *b, JXF_Int *ldb, JXF_Int *info );

/* dsteqr.c */
JXF_Int
jxf_dsteqr(const char *compz, JXF_Int *n, JXF_Real *d__, JXF_Real *e, JXF_Real *z__, JXF_Int *ldz, JXF_Real *work, JXF_Int *info );

/* dsterf.c */
JXF_Int
jxf_dsterf( JXF_Int *n, JXF_Real *d__, JXF_Real *e, JXF_Int *info );

/* dsyev.c */
JXF_Int
jxf_dsyev(const char *jobz,const char *uplo, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *w, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dsygs2.c */
JXF_Int
jxf_dsygs2( JXF_Int *itype,const char *uplo, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *b, JXF_Int *ldb, JXF_Int *info );

/* dsygst.c */
JXF_Int
jxf_dsygst( JXF_Int *itype,const char *uplo, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *b, JXF_Int *ldb, JXF_Int *info );

/* dsygv.c */
JXF_Int
jxf_dsygv( JXF_Int *itype, char *jobz, char *uplo, JXF_Int *n, JXF_Real *a,
JXF_Int *lda, JXF_Real *b, JXF_Int *ldb, JXF_Real *w, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dsytd2.c */
JXF_Int
jxf_dsytd2(const char *uplo, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *d__, JXF_Real *e, JXF_Real *tau, JXF_Int *info );

/* dsytrd.c */
JXF_Int
jxf_dsytrd(const char *uplo, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Real *d__,
JXF_Real *e, JXF_Real *tau, JXF_Real *work, JXF_Int *lwork, JXF_Int *info );

/* dtrti2.c */
JXF_Int
jxf_dtrti2(const char *uplo, const char *diag, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Int *info);

/* dtrtri.c */
JXF_Int
jxf_dtrtri(const char *uplo, const char *diag, JXF_Int *n, JXF_Real *a, JXF_Int *lda, JXF_Int *info);

#ifdef __cplusplus
}
#endif

#endif
