# 设置环境变量
export LD_LIBRARY_PATH=/vol8/appsoftware/mt3000_programming_env/third-party-lib:$LD_LIBRARY_PATH
GCCROOT =/vol8/appsoftware/mt3000_programming_env/dsp_compiler
ENV_ROOT=/vol8/appsoftware/mt3000_programming_env/hthreads

# 指定编译器/归档工具/链接器
CC=MT-3000-gcc
AR=MT-3000-ar
LD=MT-3000-ld

# 编译器标志:
CFLAGS = -c -O3 -fenable-m3000 -ffunction-sections -flax-vector-conversions -I./ \
		 -I$(ENV_ROOT)/include \
		 -I$(GCCROOT)/include/ \
		 -I$(GCCROOT)/include/compiler

# 链接器标志:
LDFLAGS=-L$(ENV_ROOT)/lib --gc-sections -Tdsp.lds

# 定义源文件/目标文件/可执行文件
# SRCS = kernel_spmv_FP64.c kernel_spmv_FP32.c kernel_spmv_FP16.c  \
#        kernel_dot_FP64.c  kernel_dot_FP32.c  kernel_dot_FP16.c vec_sum.c 
# OBJS = kernel_spmv_FP64.o kernel_spmv_FP32.o kernel_spmv_FP16.o  \
#        kernel_dot_FP64.o  kernel_dot_FP32.o  kernel_dot_FP16.o vec_sum.o 
SRCS = kernel_spmv_FP64.c kernel_spmv_FP32.c kernel_spmv_FP16.c vec_sum.c
OBJS = kernel_spmv_FP64.o kernel_spmv_FP32.o kernel_spmv_FP16.o vec_sum.o
EXE = kernel.out
DAT = kernel.dat

# 编译生成可执行文件/使用MT-3000-makedat工具生成数据文件/将生成的数据文件复制到../bin目录
# ALL: $(EXE)
# 	$(GCCROOT)/bin/MT-3000-makedat -J $(EXE)
# 	cp ${DAT} ../bin
# 	cp ${DAT} ../../

# 默认目标
ALL: $(DAT)
	cp $(DAT) ../bin || true
	cp $(DAT) ../../ || true

# 生成 .dat 文件
$(DAT): $(EXE)
	$(GCCROOT)/bin/MT-3000-makedat -J $(EXE)

# 链接生成单一的可执行文件
$(EXE): $(OBJS)
	$(GCCROOT)/bin/$(LD) $(LDFLAGS) $(OBJS) $(ENV_ROOT)/lib/libhthread_device.a \
	$(GCCROOT)/lib/vlib3000.a $(GCCROOT)/lib/slib3000.a -o $(EXE)

# 编译源文件生成目标文件
%.o: %.c
	$(GCCROOT)/bin/$(CC) $(CFLAGS) $< -o $@

# 清理生成的文件
clean:
	rm -f $(EXE) $(OBJS) $(DAT)
	rm -f ../bin/$(DAT) ../../$(DAT)

# 声明伪目标
.PHONY: ALL clean