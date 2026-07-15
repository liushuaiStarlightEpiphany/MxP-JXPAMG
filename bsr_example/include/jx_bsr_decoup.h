#ifndef JX_BSR_DECOUP_H
#define JX_BSR_DECOUP_H

#include "jx_mv.h"
#include "jx_bsr_mv.h"
#include "jx_parbsr_mv.h"

#ifdef __cplusplus
extern "C" {
#endif

// 解耦方法枚举
typedef enum {
    JX_DECOUP_NONE = 0,    // 不解耦
    JX_DECOUP_ABF = 1,     // ABF解耦（推荐用于大多数情况）
    JX_DECOUP_ANL = 2,     // 解析解耦
    JX_DECOUP_SEM = 3,     // 半解析解耦
    JX_DECOUP_QI = 4,      // Quasi-IMPES解耦
    JX_DECOUP_TIMPES = 5,  // True-IMPES解耦
    JX_DECOUP_TIMPES2 = 6  // True-IMPES解耦（版本2）
} JX_DecoupType;

/**
 * @brief 串行BSR矩阵解耦（同时处理矩阵和右端项）
 * 
 * @param A_bsr 输入的BSR矩阵，解耦会直接修改矩阵数据
 * @param rhs 右端项向量，解耦会直接修改向量数据
 * @param decoup_type 解耦类型
 * @param is_thermal 是否为热力模型（影响解耦算法）
 * @return JX_Int 成功返回JX_SUCCESS，失败返回错误码
 */
JX_Int jx_BSRMatrixDecouple(jx_BSRMatrix* A_bsr, 
                           jx_Vector* rhs,
                           JX_DecoupType decoup_type,
                           JX_Int is_thermal);

/**
 * @brief 并行BSR矩阵解耦（同时处理矩阵和右端项）
 * 
 * @param A_parbsr 输入的并行BSR矩阵，解耦会直接修改矩阵数据
 * @param rhs 并行的右端项向量，解耦会直接修改向量数据
 * @param decoup_type 解耦类型
 * @param is_thermal 是否为热力模型
 * @return JX_Int 成功返回JX_SUCCESS，失败返回错误码
 */
JX_Int jx_ParBSRMatrixDecouple(jx_ParBSRMatrix* A_parbsr,
                              jx_ParVector* rhs,
                              JX_DecoupType decoup_type,
                              JX_Int is_thermal);

/**
 * @brief 获取解耦方法名称
 * 
 * @param decoup_type 解耦类型
 * @return const char* 解耦方法名称字符串
 */
const char* jx_DecoupTypeToString(JX_DecoupType decoup_type);

#ifdef __cplusplus
}
#endif

#endif // JX_BSR_DECOUP_H