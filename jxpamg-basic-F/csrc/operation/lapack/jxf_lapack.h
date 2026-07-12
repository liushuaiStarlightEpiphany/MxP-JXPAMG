/***** DO NOT use this file outside of the LAPACK directory *****/

/*--------------------------------------------------------------------------
 * This header renames the functions in LAPACK to avoid conflicts
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

/* lapack */
#define dbdsqr_  jxf_dbdsqr
#define dgebd2_  jxf_dgebd2
#define dgebrd_  jxf_dgebrd
#define dgelq2_  jxf_dgelq2
#define dgelqf_  jxf_dgelqf
#define dgels_   jxf_dgels
#define dgeqr2_  jxf_dgeqr2
#define dgeqrf_  jxf_dgeqrf
#define dgesvd_  jxf_dgesvd
#define dgetf2_  jxf_dgetf2
#define dgetrf_  jxf_dgetrf
#define dgetri_  jxf_dgetri
#define dgetrs_  jxf_dgetrs
#define dlasq1_  jxf_dlasq1
#define dlasq2_  jxf_dlasq2
#define dlasrt_  jxf_dlasrt
#define dorg2l_  jxf_dorg2l
#define dorg2r_  jxf_dorg2r
#define dorgbr_  jxf_dorgbr
#define dorgl2_  jxf_dorgl2
#define dorglq_  jxf_dorglq
#define dorgql_  jxf_dorgql
#define dorgqr_  jxf_dorgqr
#define dorgtr_  jxf_dorgtr
#define dorm2r_  jxf_dorm2r
#define dormbr_  jxf_dormbr
#define dorml2_  jxf_dorml2
#define dormlq_  jxf_dormlq
#define dormqr_  jxf_dormqr
#define dpotf2_  jxf_dpotf2
#define dpotrf_  jxf_dpotrf
#define dpotrs_  jxf_dpotrs
#define dsteqr_  jxf_dsteqr
#define dsterf_  jxf_dsterf
#define dsyev_   jxf_dsyev
#define dsygs2_  jxf_dsygs2
#define dsygst_  jxf_dsygst
#define dsygv_   jxf_dsygv
#define dsytd2_  jxf_dsytd2
#define dsytrd_  jxf_dsytrd
#define dtrti2_  jxf_dtrti2
#define dtrtri_  jxf_dtrtri

/* lapack auxiliary routines */
#define dlabad_  jxf_dlabad
#define dlabrd_  jxf_dlabrd
#define dlacpy_  jxf_dlacpy
#define dlae2_   jxf_dlae2
#define dlaev2_  jxf_dlaev2
#define dlamch_  jxf_dlamch
#define dlamc1_  jxf_dlamc1
#define dlamc2_  jxf_dlamc2
#define dlamc3_  jxf_dlamc3
#define dlamc4_  jxf_dlamc4
#define dlamc5_  jxf_dlamc5
#define dlange_  jxf_dlange
#define dlanst_  jxf_dlanst
#define dlansy_  jxf_dlansy
#define dlapy2_  jxf_dlapy2
#define dlarf_   jxf_dlarf
#define dlarfb_  jxf_dlarfb
#define dlarfg_  jxf_dlarfg
#define dlarft_  jxf_dlarft
#define dlartg_  jxf_dlartg
#define dlas2_   jxf_dlas2
#define dlascl_  jxf_dlascl
#define dlaset_  jxf_dlaset
#define dlasq3_  jxf_dlasq3
#define dlasq4_  jxf_dlasq4
#define dlasq5_  jxf_dlasq5
#define dlasq6_  jxf_dlasq6
#define dlasr_   jxf_dlasr
#define dlassq_  jxf_dlassq
#define dlasv2_  jxf_dlasv2
#define dlaswp_  jxf_dlaswp
#define dlatrd_  jxf_dlatrd
#define ieeeck_  jxf_ieeeck
#define ilaenv_  jxf_ilaenv

/* these auxiliary routines have a different definition in BLAS */
#define lsame_   jxf_lapack_lsame
#define xerbla_  jxf_lapack_xerbla

/* this is needed so that lapack can call external BLAS */
#include "jxf_blas2.h"
