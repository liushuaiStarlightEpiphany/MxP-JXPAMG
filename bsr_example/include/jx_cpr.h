//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2024        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jx_cpr.h --  JXPAMG CPR预条件器头文件 
 *  Date: 2025/10/08
 */ 


#ifndef JX_CPR_HEADER
#define JX_CPR_HEADER

#include "jx_mv.h"
#include "jx_pamg.h"
#include "jx_parbsr_mv.h"


#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------*
 *             CPR预条件器结构体                   *
 *----------------------------------------------------------------*/

typedef struct
{
    MPI_Comm         comm;                /* MPI通信域 */
    
    /* 输入矩阵 */
    jx_ParBSRMatrix *A_bsr;               /* 原始BSR矩阵 */
    jx_ParCSRMatrix *A_pressure;          /* 压力矩阵（CSR格式） */
    
    /* 子求解器 */
    JX_Solver            *stage1_solver;       /* 压力求解器（AMG） */
    JX_Solver            *stage2_solver;       /* 全局求解器（松弛法） */
    
    /* 工作向量 */
    jx_ParVector    *rp;                  /* 压力残差向量 */
    jx_ParVector    *xp;                  /* 压力解向量 */
    jx_ParVector    *work;                /* 临时工作向量 */      // dont need currently?
    
    /* 配置参数 */
    JX_Int           pressure_index;      /* 压力变量索引 */
    JX_Int           block_size;          /* 块大小 */
    JX_Real          threshold;           /* 提取阈值 */
    JX_Int           stage1_maxit;        /* 阶段1最大迭代次数 */
    JX_Int           stage2_maxit;        /* 阶段2最大迭代次数 */
    JX_Int           stage1_solver_type;  /* 阶段1求解器类型 */
    JX_Int           stage2_solver_type;  /* 阶段2求解器类型 */
    
    /* 性能统计 */
    JX_Real          stage1_setup_time;
    JX_Real          stage1_solve_time;
    JX_Real          stage2_setup_time;
    JX_Real          stage2_solve_time;
    
    JX_Int           is_initialized;      /* 是否已初始化 */
    JX_Int           print_level;      /* 打印级别 */
    
} jx_CPRPrecond;

/*----------------------------------------------------------------*
 *             创建和销毁函数                                      *
 *----------------------------------------------------------------*/

/**
 * \brief 创建CPR预条件器
 * 
 * \param comm MPI通信域
 * \return jx_CPRPrecond* CPR预条件器指针
 */
jx_CPRPrecond* JX_CPRCreate(MPI_Comm comm);

/**
 * \brief 销毁CPR预条件器
 * 
 * \param cpr CPR预条件器指针的指针
 * \return JX_Int 错误代码
 */
JX_Int JX_CPRDestroy(jx_CPRPrecond **cpr);

/**
 * \brief 设置CPR参数（与FASPxx风格一致）
 * 
 * \param cpr CPR预条件器
 * \param param_name 参数名
 * \param value 参数值
 * \return JX_Int 错误代码
 */
JX_Int JX_CPRSetParameter(jx_CPRPrecond *cpr, const char *param_name, JX_Int value);

/**
 * \brief 设置CPR浮点参数
 * 
 * \param cpr CPR预条件器
 * \param param_name 参数名
 * \param value 参数值
 * \return JX_Int 错误代码
 */
JX_Int JX_CPRSetRealParameter(jx_CPRPrecond *cpr, const char *param_name, JX_Real value);

/*----------------------------------------------------------------*
 *             设置和求解函数（与PAMG接口一致）                    *
 *----------------------------------------------------------------*/

/**
 * \brief 设置CPR预条件器（与JX_PAMGSetup类似）
 * 
 * \param cpr CPR预条件器
 * \param par_matrix 输入矩阵（BSR格式）
 * \return JX_Int 错误代码
 */
JX_Int JX_CPRSetup(jx_CPRPrecond *cpr, jx_ParBSRMatrix* par_matrix);

/**
 * \brief 应用CPR预条件器（与JX_PAMGPrecond类似）
 * 
 * \param cpr CPR预条件器
 * \param par_matrix 输入矩阵（BSR格式）
 * \param par_rhs 右端项向量
 * \param par_sol 解向量（输出）
 * \return JX_Int 错误代码
 */
// JX_Int JX_CPRPrecond(jx_CPRPrecond *cpr, jx_ParBSRMatrix par_matrix,
//                      JX_Vector par_rhs, JX_Vector par_sol);

JX_Int JX_CPRPrecond(jx_CPRPrecond *cpr, jx_ParBSRMatrix* par_matrix,
                     jx_ParVector *par_rhs, jx_ParVector *par_sol);
/**
 * \brief 获取CPR性能统计
 * 
 * \param cpr CPR预条件器
 * \param stats 统计数组（至少4个元素）
 * \return JX_Int 错误代码
 */
JX_Int JX_CPRGetStatistics(jx_CPRPrecond *cpr, JX_Real *stats);

/**
 * \brief 打印CPR信息
 * 
 * \param cpr CPR预条件器
 * \param level 打印级别
 * \return JX_Int 错误代码
 */
JX_Int JX_CPRPrint(jx_CPRPrecond *cpr, JX_Int level);

/*----------------------------------------------------------------*
 *             工具函数                                           *
 *----------------------------------------------------------------*/

/**
 * \brief 从BSR矩阵提取压力子系统
 * 
 * \param A_bsr 输入BSR矩阵
 * \param pressure_index 压力变量索引
 * \param threshold 提取阈值
 * \return jx_ParCSRMatrix* 压力矩阵（CSR格式）
 */
jx_ParCSRMatrix* JX_ExtractPressureFromBSR(jx_ParBSRMatrix *A_bsr,
                                          JX_Int pressure_index,
                                          JX_Real threshold);

#ifdef __cplusplus
}
#endif

#endif /* JX_CPR_GMRES_HEADER */