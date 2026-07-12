//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jx_blockprec.h -- head files for block preconditioner
 *  Date: 2021/12/10
 * 
 *  Created by dyt
 */ 

#ifndef JX_BLOCKPREC_HEADER
#define JX_BLOCKPREC_HEADER

#ifndef JX_UTIL_HEADER 
#include "jx_util.h"
#endif

#ifndef JX_MV_HEADER 
#include "jx_mv.h"
#endif

#ifndef JX_HPCSRMV_HEADER 
#include "jx_hpcsr.h"
#endif

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jx_BlockPrecData
 */  
typedef struct
{
    /* precond params */
    JX_Int max_iter;    // m
    JX_Int inv_iter;    // s
    JX_Int level;       // l

    /* data generated in the setup phase */
    jx_ParVector *Df;
    jx_ParVector *DBf;
    jx_ParVector *Vtemp;
    jx_ParVector *Utemp;

} jx_BlockPrecData;


#define jx_BlockPrecDataMaxIter(bds_data)        ((bds_data) -> max_iter)
#define jx_BlockPrecDataInvIter(bds_data)        ((bds_data) -> inv_iter)
#define jx_BlockPrecDataLevel(bds_data)          ((bds_data) -> level)
#define jx_BlockPrecDataDf(bds_data)             ((bds_data) -> Df)
#define jx_BlockPrecDataDBf(bds_data)            ((bds_data) -> DBf)
#define jx_BlockPrecDataVtemp(bds_data)          ((bds_data) -> Vtemp)
#define jx_BlockPrecDataUtemp(bds_data)          ((bds_data) -> Utemp)

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/blockprec/blockprec.c */
JX_Int JX_BlockPrecCreate( MPI_Comm comm, JX_Solver *solver );

JX_Int JX_BlockPrecSetMaxIter( JX_Solver solver, JX_Int max_iter);
JX_Int JX_BlockPrecGetMaxIter( JX_Solver solver, JX_Int *max_iter);
JX_Int JX_BlockPrecSetInvIter( JX_Solver solver, JX_Int inv_iter);
JX_Int JX_BlockPrecGetInvIter( JX_Solver solver, JX_Int *inv_iter);
JX_Int JX_BlockPrecSetLevel( JX_Solver solver, JX_Int level);
JX_Int JX_BlockPrecGetLevel( JX_Solver solver, JX_Int *level);

JX_Int JX_BlockPrecDestroy( JX_Solver solver );

JX_Int 
JX_BlockPrecSetup( JX_Solver       solver, 
                   JX_hpCSRMatrix  hp_matrix );

JX_Int 
JX_BlockPrec_JAC( JX_Solver       solver,
                 JX_hpCSRMatrix  hp_matrix,
                 JX_ParVector    par_rhs,
                 JX_ParVector    par_app  );

JX_Int 
JX_BlockPrec_GS( JX_Solver       solver,
                    JX_hpCSRMatrix  hp_matrix,
                    JX_ParVector    par_rhs,
                    JX_ParVector    par_app  );

JX_Int 
JX_BlockPrec_SGS(JX_Solver       solver,
                    JX_hpCSRMatrix  hp_matrix,
                    JX_ParVector    par_rhs,
                    JX_ParVector    par_app  );

void  *jx_BlockPrecCreate(MPI_Comm comm);
JX_Int jx_BlockPrecSetMaxIter( void *solver, JX_Int max_iter );
JX_Int jx_BlockPrecGetMaxIter( void *solver, JX_Int *max_iter );
JX_Int jx_BlockPrecSetInvIter( void *solver, JX_Int inv_iter );
JX_Int jx_BlockPrecGetInvIter( void *solver, JX_Int *max_iter );
JX_Int jx_BlockPrecSetLevel( void *solver, JX_Int level );
JX_Int jx_BlockPrecGetLevel( void *solver, JX_Int *level );

JX_Int jx_BlockPrecDestroy( void *solver );

JX_Int 
jx_BlockPrecSetup( void            *solver, 
                   jx_hpCSRMatrix  *hp_matrix );

JX_Int 
jx_BlockPrec_JAC( void        *solver,
                 jx_hpCSRMatrix  *hp_matrix,
                 jx_ParVector    *par_rhs,
                 jx_ParVector    *par_app  );

JX_Int 
jx_BlockPrec_GS( void         *solver,
                 jx_hpCSRMatrix  *hp_matrix,
                 jx_ParVector    *par_rhs,
                 jx_ParVector    *par_app  );

JX_Int 
jx_BlockPrec_SGS( void        *solver,
                 jx_hpCSRMatrix  *hp_matrix,
                 jx_ParVector    *par_rhs,
                 jx_ParVector    *par_app  );
#endif