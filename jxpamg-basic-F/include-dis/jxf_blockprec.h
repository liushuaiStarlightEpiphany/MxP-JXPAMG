//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_blockprec.h -- head files for block preconditioner
 *  Date: 2021/12/10
 * 
 *  Created by dyt
 */ 

#ifndef JXF_BLOCKPREC_HEADER
#define JXF_BLOCKPREC_HEADER

#ifndef JXF_UTIL_HEADER 
#include "jxf_util.h"
#endif

#ifndef JXF_MV_HEADER 
#include "jxf_mv.h"
#endif

#ifndef JXF_HPCSRMV_HEADER 
#include "jxf_hpcsr.h"
#endif

/*----------------------------------------------------------------*
 *                   Struct Declaration                           *
 *----------------------------------------------------------------*/

/*!
 * \struct jxf_BlockPrecData
 */  
typedef struct
{
    /* precond params */
    JXF_Int max_iter;    // m
    JXF_Int inv_iter;    // s
    JXF_Int level;       // l

    /* data generated in the setup phase */
    jxf_ParVector *Df;
    jxf_ParVector *DBf;
    jxf_ParVector *Vtemp;
    jxf_ParVector *Utemp;

} jxf_BlockPrecData;


#define jxf_BlockPrecDataMaxIter(bds_data)        ((bds_data) -> max_iter)
#define jxf_BlockPrecDataInvIter(bds_data)        ((bds_data) -> inv_iter)
#define jxf_BlockPrecDataLevel(bds_data)          ((bds_data) -> level)
#define jxf_BlockPrecDataDf(bds_data)             ((bds_data) -> Df)
#define jxf_BlockPrecDataDBf(bds_data)            ((bds_data) -> DBf)
#define jxf_BlockPrecDataVtemp(bds_data)          ((bds_data) -> Vtemp)
#define jxf_BlockPrecDataUtemp(bds_data)          ((bds_data) -> Utemp)

/*----------------------------------------------------------------*
 *                   Function Declaration                         *
 *----------------------------------------------------------------*/

/* csrc/blockprec/blockprec.c */
JXF_Int JXF_BlockPrecCreate( MPI_Comm comm, JXF_Solver *solver );

JXF_Int JXF_BlockPrecSetMaxIter( JXF_Solver solver, JXF_Int max_iter);
JXF_Int JXF_BlockPrecGetMaxIter( JXF_Solver solver, JXF_Int *max_iter);
JXF_Int JXF_BlockPrecSetInvIter( JXF_Solver solver, JXF_Int inv_iter);
JXF_Int JXF_BlockPrecGetInvIter( JXF_Solver solver, JXF_Int *inv_iter);
JXF_Int JXF_BlockPrecSetLevel( JXF_Solver solver, JXF_Int level);
JXF_Int JXF_BlockPrecGetLevel( JXF_Solver solver, JXF_Int *level);

JXF_Int JXF_BlockPrecDestroy( JXF_Solver solver );

JXF_Int 
JXF_BlockPrecSetup( JXF_Solver       solver, 
                   JXF_hpCSRMatrix  hp_matrix );

JXF_Int 
JXF_BlockPrec_JAC( JXF_Solver       solver,
                 JXF_hpCSRMatrix  hp_matrix,
                 JXF_ParVector    par_rhs,
                 JXF_ParVector    par_app  );

JXF_Int 
JXF_BlockPrec_GS( JXF_Solver       solver,
                    JXF_hpCSRMatrix  hp_matrix,
                    JXF_ParVector    par_rhs,
                    JXF_ParVector    par_app  );

JXF_Int 
JXF_BlockPrec_SGS(JXF_Solver       solver,
                    JXF_hpCSRMatrix  hp_matrix,
                    JXF_ParVector    par_rhs,
                    JXF_ParVector    par_app  );

void  *jxf_BlockPrecCreate(MPI_Comm comm);
JXF_Int jxf_BlockPrecSetMaxIter( void *solver, JXF_Int max_iter );
JXF_Int jxf_BlockPrecGetMaxIter( void *solver, JXF_Int *max_iter );
JXF_Int jxf_BlockPrecSetInvIter( void *solver, JXF_Int inv_iter );
JXF_Int jxf_BlockPrecGetInvIter( void *solver, JXF_Int *max_iter );
JXF_Int jxf_BlockPrecSetLevel( void *solver, JXF_Int level );
JXF_Int jxf_BlockPrecGetLevel( void *solver, JXF_Int *level );

JXF_Int jxf_BlockPrecDestroy( void *solver );

JXF_Int 
jxf_BlockPrecSetup( void            *solver, 
                   jxf_hpCSRMatrix  *hp_matrix );

JXF_Int 
jxf_BlockPrec_JAC( void        *solver,
                 jxf_hpCSRMatrix  *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 jxf_ParVector    *par_app  );

JXF_Int 
jxf_BlockPrec_GS( void         *solver,
                 jxf_hpCSRMatrix  *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 jxf_ParVector    *par_app  );

JXF_Int 
jxf_BlockPrec_SGS( void        *solver,
                 jxf_hpCSRMatrix  *hp_matrix,
                 jxf_ParVector    *par_rhs,
                 jxf_ParVector    *par_app  );
#endif