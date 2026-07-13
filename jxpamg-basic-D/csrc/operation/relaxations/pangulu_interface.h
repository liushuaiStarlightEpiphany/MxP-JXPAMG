#ifndef PANGULU_INTERFACE_H_
#define PANGULU_INTERFACE_H_

#ifndef CALCULATE_TYPE_R64
#define CALCULATE_TYPE_R64
#endif

#include <stddef.h>
#include <stdint.h>

typedef unsigned long long sparse_pointer_t; /* CSR rowptr */
typedef unsigned int       sparse_index_t;   /* CSR colidx */
typedef double             sparse_value_t;   /* CSR values / vectors */

#ifdef __cplusplus
extern "C" {
#endif

int jx_pglu_ensure_factorized(
    int                       n,
    long long                 nnz,
    const sparse_pointer_t*   rowptr,
    const sparse_index_t*     colidx,
    const sparse_value_t*     aval,
    int                       nb_hint,
    int                       nthread_hint,
    int*                      n_pad_out);

int jx_pglu_solve(const double* rhs_pad, double* sol_pad);

void jx_pglu_release(void);

#ifdef __cplusplus
}
#endif

#endif /* PANGULU_INTERFACE_H_ */

