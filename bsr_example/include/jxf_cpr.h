//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2024        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_cpr.h --  JXFPAMG CPR预条件器头文件 
 *  Date: 2025/10/08
 */ 


#ifndef JXF_CPR_HEADER
#define JXF_CPR_HEADER

#include "jxf_mv.h"
#include "jxf_pamg.h"
#include "jxf_parbsr_mv.h"


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
    jxf_ParBSRMatrix *A_bsr;               /* 原始BSR矩阵 */
    jxf_ParCSRMatrix *A_pressure;          /* 压力矩阵（CSR格式） */
    
    /* 子求解器 */
    JXF_Solver            *stage1_solver;       /* 压力求解器（AMG） */
    JXF_Solver            *stage2_solver;       /* 全局求解器（松弛法） */
    
    /* 工作向量 */
    jxf_ParVector    *rp;                  /* 压力残差向量 */
    jxf_ParVector    *xp;                  /* 压力解向量 */
    jxf_ParVector    *work;                /* 临时工作向量 */      // dont need currently?
    
    /* 配置参数 */
    JXF_Int           pressure_index;      /* 压力变量索引 */
    JXF_Int           block_size;          /* 块大小 */
    JXF_Real          threshold;           /* 提取阈值 */
    JXF_Int           stage1_maxit;        /* 阶段1最大迭代次数 */
    JXF_Int           stage2_maxit;        /* 阶段2最大迭代次数 */
    JXF_Int           stage1_solver_type;  /* 阶段1求解器类型 */
    JXF_Int           stage2_solver_type;  /* 阶段2求解器类型 */
    
    /* 性能统计 */
    JXF_Real          stage1_setup_time;
    JXF_Real          stage1_solve_time;
    JXF_Real          stage2_setup_time;
    JXF_Real          stage2_solve_time;
    
    JXF_Int           is_initialized;      /* 是否已初始化 */
    JXF_Int           print_level;      /* 打印级别 */
    
} jxf_CPRPrecond;

/*----------------------------------------------------------------*
 *             创建和销毁函数                                      *
 *----------------------------------------------------------------*/

/*! 
 * \brief 创建CPR预条件器求解器对象
 * \param solver 输出参数，指向创建的求解器对象
 * \return 错误码
 */
JXF_Int JXF_CPRCreateSolver(JXF_Solver *solver);

/**
 * \brief 创建CPR预条件器
 * 
 * \param comm MPI通信域
 * \return jxf_CPRPrecond* CPR预条件器指针
 */
jxf_CPRPrecond* JXF_CPRCreate(MPI_Comm comm);

/**
 * \brief 销毁CPR预条件器
 * 
 * \param cpr CPR预条件器指针的指针
 * \return JXF_Int 错误代码
 */
JXF_Int JXF_CPRDestroy(jxf_CPRPrecond **cpr);

/**
 * \brief 设置CPR参数（与FASPxx风格一致）
 * 
 * \param cpr CPR预条件器
 * \param param_name 参数名
 * \param value 参数值
 * \return JXF_Int 错误代码
 */
JXF_Int JXF_CPRSetParameter(jxf_CPRPrecond *cpr, const char *param_name, JXF_Int value);

/**
 * \brief 设置CPR浮点参数
 * 
 * \param cpr CPR预条件器
 * \param param_name 参数名
 * \param value 参数值
 * \return JXF_Int 错误代码
 */
JXF_Int JXF_CPRSetRealParameter(jxf_CPRPrecond *cpr, const char *param_name, JXF_Real value);

/*----------------------------------------------------------------*
 *             设置和求解函数（与PAMG接口一致）                    *
 *----------------------------------------------------------------*/

/**
 * \brief 设置CPR预条件器（与JXF_PAMGSetup类似）
 * 
 * \param cpr CPR预条件器
 * \param par_matrix 输入矩阵（BSR格式）
 * \return JXF_Int 错误代码
 */
JXF_Int JXF_CPRSetup(jxf_CPRPrecond *cpr, jxf_ParBSRMatrix* par_matrix);

/**
 * \brief 应用CPR预条件器（与JXF_PAMGPrecond类似）
 * 
 * \param cpr CPR预条件器
 * \param par_matrix 输入矩阵（BSR格式）
 * \param par_rhs 右端项向量
 * \param par_sol 解向量（输出）
 * \return JXF_Int 错误代码
 */
// JXF_Int JXF_CPRPrecond(jxf_CPRPrecond *cpr, jxf_ParBSRMatrix par_matrix,
//                      JXF_Vector par_rhs, JXF_Vector par_sol);

JXF_Int JXF_CPRPrecond(jxf_CPRPrecond *cpr, jxf_ParBSRMatrix* par_matrix,
                     jxf_ParVector *par_rhs, jxf_ParVector *par_sol);
/**
 * \brief 获取CPR性能统计
 * 
 * \param cpr CPR预条件器
 * \param stats 统计数组（至少4个元素）
 * \return JXF_Int 错误代码
 */
JXF_Int JXF_CPRGetStatistics(jxf_CPRPrecond *cpr, JXF_Real *stats);

/**
 * \brief 打印CPR信息
 * 
 * \param cpr CPR预条件器
 * \param level 打印级别
 * \return JXF_Int 错误代码
 */
JXF_Int JXF_CPRPrint(jxf_CPRPrecond *cpr, JXF_Int level);

/*----------------------------------------------------------------*
 *             工具函数                                           *
 *----------------------------------------------------------------*/

/**
 * \brief 从BSR矩阵提取压力子系统
 * 
 * \param A_bsr 输入BSR矩阵
 * \param pressure_index 压力变量索引
 * \param threshold 提取阈值
 * \return jxf_ParCSRMatrix* 压力矩阵（CSR格式）
 */
jxf_ParCSRMatrix* JXF_ExtractPressureFromBSR(jxf_ParBSRMatrix *A_bsr,
                                          JXF_Int pressure_index,
                                          JXF_Real threshold);

jxf_ParVector* jxf_CreateGlobalVector(jxf_ParBSRMatrix *A_bsr);

#ifdef __cplusplus
}
#endif

#endif /* JXF_CPR_GMRES_HEADER */