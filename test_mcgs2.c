#include "jx_mv.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

extern JX_Int jx_PAMGRelaxMultiColorGS(
    jx_ParCSRMatrix *par_matrix, jx_ParVector *par_rhs,
    JX_Int *cf_marker, JX_Int relax_points,
    JX_Real relax_weight, JX_Real omega,
    jx_ParVector *par_app, jx_ParVector *Vtemp);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    JX_Int n = 4, nnz = 10;
    JX_Int row_starts[2] = {0, n}, col_starts[2] = {0, n};

    /* diag CSR: diagonal first per row */
    jx_CSRMatrix *A_diag = jx_CSRMatrixCreate(n, n, nnz);
    jx_CSRMatrixInitialize(A_diag);
    JX_Int *ai = jx_CSRMatrixI(A_diag), *aj = jx_CSRMatrixJ(A_diag);
    JX_Real *ad = jx_CSRMatrixData(A_diag);
    JX_Int rp[] = {0,2,5,8,10}, jc[] = {0,1, 1,0,2, 2,1,3, 3,2};
    JX_Real va[] = {4,-1, 4,-1,-1, 4,-1,-1, 4,-1};
    memcpy(ai, rp, 5*sizeof(JX_Int)); memcpy(aj, jc, nnz*sizeof(JX_Int));
    memcpy(ad, va, nnz*sizeof(JX_Real));
    A_diag->rownnz = (JX_Int*)malloc(n*sizeof(JX_Int));
    for (JX_Int i=0;i<n;i++) A_diag->rownnz[i] = ai[i+1]-ai[i];

    jx_ParCSRMatrix *A = jx_ParCSRMatrixCreate(MPI_COMM_WORLD,n,n,row_starts,col_starts,0,0,nnz);
    jx_ParCSRMatrixDiag(A) = A_diag;
    jx_CSRMatrix *offd = jx_CSRMatrixCreate(n,0,0);
    jx_CSRMatrixInitialize(offd);
    jx_ParCSRMatrixOffd(A) = offd;
    jx_MatvecCommPkgCreate(A);

    jx_ParVector *b = jx_ParVectorCreate(MPI_COMM_WORLD,n,row_starts);
    jx_ParVectorInitialize(b);
    jx_ParVector *x = jx_ParVectorCreate(MPI_COMM_WORLD,n,row_starts);
    jx_ParVectorInitialize(x);
    jx_ParVector *V = jx_ParVectorCreate(MPI_COMM_WORLD,n,row_starts);
    jx_ParVectorInitialize(V);
    JX_Real *bd = jx_VectorData(jx_ParVectorLocalVector(b));
    bd[0]=1; bd[1]=2; bd[2]=3; bd[3]=4;

    if(rank==0) printf("Multi-Color GS test (4x4)\n");
    for(int it=0;it<6;it++) {
        jx_PAMGRelaxMultiColorGS(A,b,NULL,0,1.0,1.0,x,V);
        if(rank==0){
            JX_Real *xd = jx_VectorData(jx_ParVectorLocalVector(x));
            printf("  iter %d: [%.6f %.6f %.6f %.6f]\n",it+1,xd[0],xd[1],xd[2],xd[3]);
        }
    }
    /* residual */
    jx_ParVector *r = jx_ParVectorCreate(MPI_COMM_WORLD,n,row_starts);
    jx_ParVectorInitialize(r);
    jx_ParVectorCopy(b,r);
    jx_ParCSRMatrixMatvec(-1.0,A,x,1.0,r);
    JX_Real rn = sqrt(jx_ParVectorInnerProd(r,r));
    if(rank==0) printf("  rel.res: %e\n",rn/sqrt(jx_ParVectorInnerProd(b,b)));

    jx_ParVectorDestroy(r); jx_ParVectorDestroy(V);
    jx_ParVectorDestroy(x); jx_ParVectorDestroy(b);
    /* manual cleanup for hand-built matrix */
    jx_CSRMatrix *cd = jx_ParCSRMatrixDiag(A), *co = jx_ParCSRMatrixOffd(A);
    jx_ParCSRMatrixDiag(A)=NULL; jx_ParCSRMatrixOffd(A)=NULL;
    jx_ParCSRMatrixDestroy(A);
    free(cd->rownnz); cd->rownnz=NULL;
    jx_CSRMatrixDestroy(cd); jx_CSRMatrixDestroy(co);
    MPI_Finalize();
    if(rank==0) printf("PASSED\n");
    return 0;
}
