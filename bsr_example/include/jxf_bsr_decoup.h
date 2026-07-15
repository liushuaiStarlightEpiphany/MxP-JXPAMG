#ifndef JXF_BSR_DECOUP_H
#define JXF_BSR_DECOUP_H

#include "jxf_mv.h"
#include "jxf_bsr_mv.h"
#include "jxf_parbsr_mv.h"

#ifdef __cplusplus
extern "C" {
#endif

// 解耦方法枚举
typedef enum {
    JXF_DECOUP_NONE = 0,
    JXF_DECOUP_ABF,
    JXF_DECOUP_ANL,
    JXF_DECOUP_SEM,
    JXF_DECOUP_QI,
    JXF_DECOUP_TIMPES,   // 新增
    JXF_DECOUP_TIMPES2   // 新增
} JXF_DecoupType;
/**
 * @brief 串行BSR矩阵解耦（同时处理矩阵和右端项）
 * 
 * @param A_bsr 输入的BSR矩阵，解耦会直接修改矩阵数据
 * @param rhs 右端项向量，解耦会直接修改向量数据
 * @param decoup_type 解耦类型
 * @param is_thermal 是否为热力模型（影响解耦算法）
 * @return JXF_Int 成功返回JXF_SUCCESS，失败返回错误码
 */
JXF_Int jxf_BSRMatrixDecouple(jxf_BSRMatrix* A_bsr, 
                           jxf_Vector* rhs,
                           JXF_DecoupType decoup_type,
                           JXF_Int is_thermal);

/**
 * @brief 并行BSR矩阵解耦（同时处理矩阵和右端项）
 * 
 * @param A_parbsr 输入的并行BSR矩阵，解耦会直接修改矩阵数据
 * @param rhs 并行的右端项向量，解耦会直接修改向量数据
 * @param decoup_type 解耦类型
 * @param is_thermal 是否为热力模型
 * @return JXF_Int 成功返回JXF_SUCCESS，失败返回错误码
 */
JXF_Int jxf_ParBSRMatrixDecouple(jxf_ParBSRMatrix* A_parbsr,
                              jxf_ParVector* rhs,
                              JXF_DecoupType decoup_type,
                              JXF_Int is_thermal);

/**
 * @brief 获取解耦方法名称
 * 
 * @param decoup_type 解耦类型
 * @return const char* 解耦方法名称字符串
 */
const char* jxf_DecoupTypeToString(JXF_DecoupType decoup_type);

#ifdef __cplusplus
}
#endif

#endif // JXF_BSR_DECOUP_H