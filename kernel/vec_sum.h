#include <compiler/m3000.h>
#include "hthread_device.h"
#include <math.h>

double reduce_sum_fp64(lvector double a);

double reduce_sum_fp32(lvector float a);

float reduce_sum_fp16(lvector_float16 a);

float reduce_sum_float(lvector float a);

__fp16 reduce_sum_half(lvector __fp16 a);

lvector double Reduction_Vector4(lvector double vec);