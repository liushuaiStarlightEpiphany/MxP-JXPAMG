/***** DO NOT use this file outside of the BLAS directory *****/

/*--------------------------------------------------------------------------
 * This header renames the functions in BLAS to avoid conflicts
 *--------------------------------------------------------------------------*/

/* blas */
#define dasum_   jxf_dasum
#define daxpy_   jxf_daxpy
#define dcopy_   jxf_dcopy
#define ddot_    jxf_ddot
#define dgemm_   jxf_dgemm
#define dgemv_   jxf_dgemv
#define dger_    jxf_dger
#define dnrm2_   jxf_dnrm2
#define drot_    jxf_drot
#define dscal_   jxf_dscal
#define dswap_   jxf_dswap
#define dsymm_   jxf_dsymm
#define dsymv_   jxf_dsymv
#define dsyr2_   jxf_dsyr2
#define dsyr2k_  jxf_dsyr2k
#define dsyrk_   jxf_dsyrk
#define dtrmm_   jxf_dtrmm
#define dtrmv_   jxf_dtrmv
#define dtrsm_   jxf_dtrsm
#define dtrsv_   jxf_dtrsv
#define idamax_  jxf_idamax

/* f2c library routines */
#define s_cmp    jxf_s_cmp
#define s_copy   jxf_s_copy
#define s_cat    jxf_s_cat
#define d_lg10   jxf_d_lg10
#define d_sign   jxf_d_sign
#define pow_dd   jxf_pow_dd
#define pow_di   jxf_pow_di

/* these auxiliary routines have a different definition in LAPACK */
#define lsame_   jxf_blas_lsame
#define xerbla_  jxf_blas_xerbla
