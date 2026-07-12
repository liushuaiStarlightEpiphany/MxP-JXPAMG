#include "jx_pamg.h"
#if JX_USING_MUMPS
#include "dmumps_c.h"
#endif

JX_Int
jx_Mumps( jx_hpCSRMatrix *hp_A, jx_ParVector *pb, jx_ParVector *px )
{
    JX_Int relax_error = 0;
#if JX_USING_MUMPS
    jx_ParCSRMatrix *pA = jx_hpCSRMatrixPar(hp_A);
    MPI_Comm comm = jx_ParCSRMatrixComm(pA);
    jx_CSRMatrix *A = jx_ParCSRMatrixDiag(pA);
    JX_Real *ADATA = jx_CSRMatrixData(A);
    JX_Int *AIA = jx_CSRMatrixI(A);
    JX_Int *AJA = jx_CSRMatrixJ(A);
    JX_Int num_rows = jx_CSRMatrixNumRows(A);
    JX_Int num_nonzeros = jx_CSRMatrixNumNonzeros(A);
    jx_Vector *b = jx_ParVectorLocalVector(pb);
    JX_Real *bDATA = jx_VectorData(b);
    jx_Vector *x = jx_ParVectorLocalVector(px);
    JX_Real *xDATA = jx_VectorData(x);
    DMUMPS_STRUC_C id;
	JX_Int myid, j, k;

    JX_Int *irn = jx_CTAlloc(JX_Int, num_nonzeros);
    JX_Int *jcn = jx_CTAlloc(JX_Int, num_nonzeros);
    JX_Real *adn = jx_CTAlloc(JX_Real, num_nonzeros);
    JX_Real *rhs = jx_CTAlloc(JX_Real, num_rows);

    jx_MPI_Comm_rank(comm, &myid);

    for (j = 0; j < num_rows; j ++)
    {
        for (k = AIA[j]; k < AIA[j+1]; k ++)
        {
            irn[k] = j + 1;
            jcn[k] = AJA[k] + 1;
            adn[k] = ADATA[k];
        }
    }
    for (j = 0; j < num_rows; j ++)
    {
        rhs[j] = bDATA[j];
    }
    /* Initialize a MUMPS instance */
    id.job = -1;
    id.par = 1;
    id.sym = 0;
    id.comm_fortran = -987654;
    dmumps_c(&id);
    /* Define the problem on the host */
    if (myid == 0)
    {
        id.n = num_rows;
        id.nz = num_nonzeros;
        id.irn = irn;
        id.jcn = jcn;
        id.a = adn;
        id.rhs = rhs;
    }
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
    id.ICNTL(1) = -1;
    id.ICNTL(2) = -1;
    id.ICNTL(3) = -1;
    id.ICNTL(4) = 0;
    /* Call the MUMPS package. */
    id.job = 6;
    dmumps_c(&id);
    for (j = 0; j < num_rows; j ++)
    {
        xDATA[j] = rhs[j];
    }
    id.job = -2;
    dmumps_c(&id); /* Terminate instance */
    jx_TFree(irn);
    jx_TFree(jcn);
    jx_TFree(adn);
    jx_TFree(rhs);
#else
    jx_printf(" >>> Without linking the MUMPS library but invoking jx_Mumps!!\n");
#endif
    return(relax_error);
}
