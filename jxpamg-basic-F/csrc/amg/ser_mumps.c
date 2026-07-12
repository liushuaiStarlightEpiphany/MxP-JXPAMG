#include "jxf_pamg.h"
#if JXF_USING_MUMPS
#include "dmumps_c.h"
#endif

JXF_Int
jxf_Mumps( jxf_hpCSRMatrix *hp_A, jxf_ParVector *pb, jxf_ParVector *px )
{
    JXF_Int relax_error = 0;
#if JXF_USING_MUMPS
    jxf_ParCSRMatrix *pA = jxf_hpCSRMatrixPar(hp_A);
    MPI_Comm comm = jxf_ParCSRMatrixComm(pA);
    jxf_CSRMatrix *A = jxf_ParCSRMatrixDiag(pA);
    JXF_Real *ADATA = jxf_CSRMatrixData(A);
    JXF_Int *AIA = jxf_CSRMatrixI(A);
    JXF_Int *AJA = jxf_CSRMatrixJ(A);
    JXF_Int num_rows = jxf_CSRMatrixNumRows(A);
    JXF_Int num_nonzeros = jxf_CSRMatrixNumNonzeros(A);
    jxf_Vector *b = jxf_ParVectorLocalVector(pb);
    JXF_Real *bDATA = jxf_VectorData(b);
    jxf_Vector *x = jxf_ParVectorLocalVector(px);
    JXF_Real *xDATA = jxf_VectorData(x);
    DMUMPS_STRUC_C id;
	JXF_Int myid, j, k;

    JXF_Int *irn = jxf_CTAlloc(JXF_Int, num_nonzeros);
    JXF_Int *jcn = jxf_CTAlloc(JXF_Int, num_nonzeros);
    JXF_Real *adn = jxf_CTAlloc(JXF_Real, num_nonzeros);
    JXF_Real *rhs = jxf_CTAlloc(JXF_Real, num_rows);

    jxf_MPI_Comm_rank(comm, &myid);

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
    jxf_TFree(irn);
    jxf_TFree(jcn);
    jxf_TFree(adn);
    jxf_TFree(rhs);
#else
    jxf_printf(" >>> Without linking the MUMPS library but invoking jxf_Mumps!!\n");
#endif
    return(relax_error);
}
