#include "spgemm_adapter.h"

int spgemm_adapter_multiply(jx_CSRMatrix *A, jx_CSRMatrix *B, jx_CSRMatrix **C)
{
    return spgemm_adapter_multiply_ex(A, B, C, spgemm_algo_from_env());
}
