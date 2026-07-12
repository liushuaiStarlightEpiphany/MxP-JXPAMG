1. 概述
   JXPAMG作为JASMIN框架的依赖库, 以往的基于makefile的编译系统
无法实现编译/安装/测试/打包的自动化, 既不利于JXPAMG库本身的研
发和质量管控, 也不利于JASMIN的自动化部署.
   为此, 特针对JXPAMG的编译系统进行了重构, 采用cmake替代makefile,
并提供了一键编译脚本. 同时为了照顾研发人员的科研习惯, 继续保留
了原来的makefile编译系统.

   新增文件及功能说明:
   $TOPDIR/
   |-- CMakeLists.txt    # cmake描述文件, 研发人员请认真阅读此文件
   `-- scripts/          # 编译相关的脚本
       |-- autobuild.sh  # 自动化编译脚本, 完成编译/安装/测试/打包的功能
       |-- cmake-modules # cmake资源模块
       |   `-- Findmpi.cmake
       |-- CPack.STGZ_Header.sh.in # cmake 打包时的脚本文件
       |-- get_svn_reversion.sh    # 获取svn版本号
       `-- get_system_name.sh      # 获取系统平台信息
   
   研发人员请认真阅读理解CMakeLists.txt, 掌握cmake的常用命令.
   同时, 还请掌握autobuild.sh脚本的用法.

2. 依赖关系库
   JXPAMG依赖于MPI和openmp. 在CMakeLists.txt中对MPI和OpenMP分别进行了处理.
   其中缺省的mpi_ROOT为/usr/local/jasmin/gnu/thirdparty/mpich2-1.5

3. 编译/安装/测试/打包

3.1 如何编译
   $ cd $TOPDIR
   $ sh scripts/autobuild.sh -T make
   如果需要指定mpi, 则:
   $ sh scripts/autobuild.sh -T make --with-mpi <mpi_root>
   如果需要指定编译器, 可以:
   $ sh scripts/autobuild.sh -T make --with-CC=gcc|icc --with-CXX=g++|icpc --with-FC=gfortran|ifort
   ATTENTION: 如果指定的编译器为mpicc,mpicxx,mpif90,则无需再次指定--with-mpi.

3.2 如何安装
   $ cd $TOPDIR
   # 缺省情况下安装到$TOPDIR/build/jxpamg
   $ sh scripts/autobuild.sh -T install --with-mpi <mpi_root>
   # 如需指定安装路径, 可执行:
   $ sh scripts/autobuild.sh -T install --with-mpi <mpi_root> --install <prefix>

3.3 如何测试
   $ cd $TOPDIR
   $ sh scripts/autobuild.sh -T test --with-mpi <mpi_root>

3.4 如何制作安装包
   $ cd $TOPDIR
   $ sh scripts/autobuild.sh -T pack --with-mpi <mpi_root>

3.5 如何自定义安装包的名字
   $ cd $TOPDIR
   $ sh scripts/autobuild.sh -T pack -e "-DCPACK_PACKAGE_NAME=<XXX> -DCPACK_PACKAGE_VERSION=<XXX> -DCPACK_SYSTEM_NAME=<XXX>"

3.6 如何添加新的测试用例
   请参阅CMakeLists.txt中添加测试程序和添加测试用例.

3.7 如何修改软件发布版本
   请参阅CMakeLists.txt中JXPAMG_MAJOR_VERSION, JXPAMG_MINOR_VERSION, JXPAMG_PATCH_VERSION


