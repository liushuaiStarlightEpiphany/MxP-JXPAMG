#include <compiler/m3000.h>
#include "hthread_device.h"
#include <math.h>
#include <stdint.h>
#include "../vec_sum.h"

// 考虑并发累加？
// 16 次累加 or 分批次并行累加？
double reduce_sum_fp64(lvector double a)
{
  mov_to_svr_v16df(a); // 从向量搬移到SVR
  return mov_from_svr0_df() + mov_from_svr1_df() +
         mov_from_svr2_df() + mov_from_svr3_df() +
         mov_from_svr4_df() + mov_from_svr5_df() +
         mov_from_svr6_df() + mov_from_svr7_df() +
         mov_from_svr8_df() + mov_from_svr9_df() +
         mov_from_svr10_df() + mov_from_svr11_df() +
         mov_from_svr12_df() + mov_from_svr13_df() +
         mov_from_svr14_df() + mov_from_svr15_df();
}

// double reduce_sum_fp64_parallel(lvector double a)
// double reduce_sum_fp64(lvector double a)
// {
//   mov_to_svr_v16df(a); // 从向量搬移到SVR

//   // 第一级：两两相加，得到8个中间结果
//   double s0 = mov_from_svr0_df() + mov_from_svr1_df();
//   double s1 = mov_from_svr2_df() + mov_from_svr3_df();
//   double s2 = mov_from_svr4_df() + mov_from_svr5_df();
//   double s3 = mov_from_svr6_df() + mov_from_svr7_df();
//   double s4 = mov_from_svr8_df() + mov_from_svr9_df();
//   double s5 = mov_from_svr10_df() + mov_from_svr11_df();
//   double s6 = mov_from_svr12_df() + mov_from_svr13_df();
//   double s7 = mov_from_svr14_df() + mov_from_svr15_df();

//   // 第二级：两两相加，得到4个中间结果
//   double r0 = s0 + s1;
//   double r1 = s2 + s3;
//   double r2 = s4 + s5;
//   double r3 = s6 + s7;

//   // 第三级：两两相加，得到2个中间结果
//   double result_0 = r0 + r1;
//   double result_1 = r2 + r3;

//   // 第四级：最终相加，得到结果
//   return result_0 + result_1;
// }

double reduce_sum_fp32(lvector float a)
{
  lvector double l = vec_fstdl(a);
  lvector double h = vec_fstdh(a);
  return reduce_sum_fp64(l) + reduce_sum_fp64(h);
}

float reduce_sum_fp16(lvector_float16 a)
{
  lvector float l = vec_fstdl(a);
  lvector float h = vec_fstdh(a);
  return reduce_sum_fp32(l) + reduce_sum_fp32(h);
}

union FloatToInt
{
  float floatValue;
  int intValue;
};

float reduce_sum_float(lvector float a)
{
  lvector signed long int b = (lvector signed long int)a; // 64位有符号整数
                                                          // lvector signed int b = (lvector signed int)a;
  mov_to_svr_v16di(b);                                    // 向量搬移到SVR
  float sum = 0;
  union FloatToInt converter;

  converter.intValue = (int)(mov_from_svr0() & 0xFFFFFFFF); // SVR到目的地址
  // 位掩码：提取64位值的低32位，将64位值的高32位清零，存入intValue
  sum += converter.floatValue; // 通过floatValue解释为浮点数并累加
  converter.intValue = (int)(mov_from_svr0() >> 32);
  // 移位：提取64位值的高32位，存入intValue
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr1() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr1() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr2() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr2() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr3() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr3() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr4() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr4() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr5() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr5() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr6() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr6() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr7() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr7() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr8() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr8() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr9() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr9() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr10() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr10() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr11() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr11() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr12() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr12() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr13() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr13() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr14() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr14() >> 32);
  sum += converter.floatValue;

  converter.intValue = (int)(mov_from_svr15() & 0xFFFFFFFF);
  sum += converter.floatValue;
  converter.intValue = (int)(mov_from_svr15() >> 32);
  sum += converter.floatValue;

  return sum;
}

union HalfToInt
{
  __fp16 halfValue;  // 16-bit half-precision float
  uint16_t intValue; // 16-bit unsigned integer
};

__fp16 reduce_sum_half(lvector __fp16 a)
{
  lvector signed long int b = (lvector signed long int)a; // 转换为 64 位整数向量
  mov_to_svr_v16di(b);                                    // 向量搬移到 SVR
  float sum = 0.0f;                                       // 使用 float 累加
  union HalfToInt converter;

  // 处理 svr0 到 svr15
  long long value;

  // svr0
  value = mov_from_svr0();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr1
  value = mov_from_svr1();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr2
  value = mov_from_svr2();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr3
  value = mov_from_svr3();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr4
  value = mov_from_svr4();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr5
  value = mov_from_svr5();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr6
  value = mov_from_svr6();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr7
  value = mov_from_svr7();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr8
  value = mov_from_svr8();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr9
  value = mov_from_svr9();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr10
  value = mov_from_svr10();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr11
  value = mov_from_svr11();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr12
  value = mov_from_svr12();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr13
  value = mov_from_svr13();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr14
  value = mov_from_svr14();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  // svr15
  value = mov_from_svr15();
  converter.intValue = (uint16_t)(value & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 16) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 32) & 0xFFFF);
  sum += converter.halfValue;
  converter.intValue = (uint16_t)((value >> 48) & 0xFFFF);
  sum += converter.halfValue;

  return sum; // 返回 float 类型
}

// ================== mBSR =================================
lvector double Reduction_Vector4(lvector double vec)
{
  double r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
  lvector double ret;
  asm(
      "VMVCGC		%16,SVR\n\t"
      "SNOP		2\n\t"
      "SMVCCG		SVR0,%0\n\t"
      "SMVCCG		SVR1,%1\n\t"
      "SMVCCG		SVR2,%2\n\t"
      "SMVCCG		SVR3,%3\n\t"
      "SMVCCG		SVR4,%4\n\t"
      "SMVCCG		SVR5,%5\n\t"
      "SMVCCG		SVR6,%6\n\t"
      "SMVCCG		SVR7,%7\n\t"
      "SMVCCG		SVR8,%8\n\t"
      "SMVCCG		SVR9,%9\n\t"
      "SMVCCG		SVR10,%10\n\t"
      "SMVCCG		SVR11,%11\n\t"
      "SMVCCG		SVR12,%12\n\t"
      "SMVCCG		SVR13,%13\n\t"
      "SMVCCG		SVR14,%14\n\t"
      "SMVCCG		SVR15,%15\n\t"
      "SNOP		1\n\t"
      : "=r"(r0), "=r"(r1), "=r"(r2), "=r"(r3), "=r"(r4), "=r"(r5), "=r"(r6), "=r"(r7), "=r"(r8), "=r"(r9), "=r"(r10), "=r"(r11), "=r"(r12), "=r"(r13), "=r"(r14), "=r"(r15)
      : "v"(vec) // 16
  );
  r0 = r0 + r1 + r2 + r3;
  r1 = r4 + r5 + r6 + r7;
  r2 = r8 + r9 + r10 + r11;
  r3 = r12 + r13 + r14 + r15;
  r4 = 0;
  r5 = 0;
  r6 = 0;
  r7 = 0;
  r8 = 0;
  r9 = 0;
  r10 = 0;
  r11 = 0;
  r12 = 0;
  r13 = 0;
  r14 = 0;
  r15 = 0;
  asm(
      "SMVCGC		%1, SVR0\n\t"
      "SMVCGC		%2, SVR1\n\t"
      "SMVCGC		%3, SVR2\n\t"
      "SMVCGC		%4,  SVR3\n\t"
      "SMVCGC		%5,  SVR4\n\t"
      "SMVCGC		%6,  SVR5\n\t"
      "SMVCGC		%7,  SVR6\n\t"
      "SMVCGC		%8,  SVR7\n\t"
      "SMVCGC		%9,  SVR8\n\t"
      "SMVCGC		%10, SVR9\n\t"
      "SMVCGC		%11, SVR10\n\t"
      "SMVCGC		%12, SVR11\n\t"
      "SMVCGC		%13, SVR12\n\t"
      "SMVCGC		%14, SVR13\n\t"
      "SMVCGC		%15, SVR14\n\t"
      "SMVCGC		%16, SVR15\n\t"
      "SNOP		2\n\t"
      "VMVCCG		SVR,%0\n\t"
      "SNOP		1\n\t"
      : "=v"(ret) //
      : "r"(r0), "r"(r1), "r"(r2), "r"(r3), "r"(r4), "r"(r5), "r"(r6), "r"(r7), "r"(r8), "r"(r9), "r"(r10), "r"(r11), "r"(r12), "r"(r13), "r"(r14), "r"(r15));
  return ret;
}