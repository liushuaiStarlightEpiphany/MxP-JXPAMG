cmake_minimum_required(VERSION 2.8)
project(jxpamg)

###########################################################
# 1. 设置版本号, 每次发布前请更新此版本号
set(JXPAMG_MAJOR_VERSION 1)
set(JXPAMG_MINOR_VERSION 2)
set(JXPAMG_PATCH_VERSION 20)
set(JXPAMG_VERSION "${JXPAMG_MAJOR_VERSION}.${JXPAMG_MINOR_VERSION}.${JXPAMG_PATCH_VERSION}")

###########################################################
# 2. 处理依赖关系
# 依赖库1: mpi
# 设置mpi的根目录
if (NOT mpi_ROOT)
  set(mpi_ROOT "/usr/local/jasmin/thirdparty/mpich2-1.5")
  message(STATUS "set mpi_ROOT as ${mpi_ROOT} by default!"
    "To change CMAKE_CXX_COMPILER, please run "
    "$ cmake -Dmpi_ROOT=<dir> <project_src_dir>")
endif ()
# 调用scripts/cmake-modules/Findmpi.cmake查找mpi库
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts/cmake-modules)
find_package(mpi)
if (mpi_FOUND)
  include_directories(${mpi_INCLUDE_DIRS})
endif ()

# 依赖库2: openmp
find_package(OpenMP)
if (OPENMP_FOUND)
  add_definitions(-DJX_USING_OPENMP=1)
  list(APPEND CMAKE_C_FLAGS "${OpenMP_C_FLAGS}")
  list(APPEND CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif ()

###########################################################
# 3. 编译JXPAMG库
# 添加头文件搜索路径
include_directories(include)
add_definitions(-DJX_USING_BIG_INT=0)
add_definitions(-DJX_USING_BIG_DOUBLE=0)
# 检索源文件, 添加库目标
file(GLOB_RECURSE C_src_files csrc/amg/*.c csrc/apctl/*.c csrc/combined/*.c csrc/communication/*.c csrc/diagscale/*.c csrc/euclid/*.c csrc/ilu/*.c csrc/krylov/*.c csrc/mat_vec_data/*.c csrc/multils/*.c csrc/utilities/*.c csrc/operation/matops/*.c csrc/operation/mvops/*.c csrc/operation/relaxations/*.c csrc/operation/vecops/*.c csrc/operation/blas/*.c csrc/operation/lapack/*.c)
file(GLOB_RECURSE Fortran_src_files fsrc/*.for fsrc/*.f fsrc/*.f90)
add_library(JXPAMG ${C_src_files} ${Fortran_src_files})
target_link_libraries(JXPAMG ${mpi_LIBRARIES} pthread rt m)

###########################################################
# 4. 安装
install(TARGETS JXPAMG DESTINATION lib)
install(DIRECTORY include DESTINATION .)
install(DIRECTORY example DESTINATION .)

###########################################################
# 5. 打包
# 使cpack使用我们自定义的脚本CPack.STGZ_Header.sh.in制作stgz包
# 可以把pre-install和post-install都放在这个脚本里面
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts)

# 安装包的命名规则
# ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}
if (NOT CPACK_PACKAGE_NAME)
  set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
endif ()
if (NOT CPACK_PACKAGE_VERSION)
  set(CPACK_PACKAGE_VERSION "${JXPAMG_VERSION}")
endif ()
if (NOT CPACK_SYSTEM_NAME)
  execute_process(
    COMMAND bash ${PROJECT_SOURCE_DIR}/scripts/get_system_name.sh
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    OUTPUT_VARIABLE system_name
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CPACK_SYSTEM_NAME "${system_name}")
endif ()
include(InstallRequiredSystemLibraries)
include(CPack)

# 制作自解压格式的STGZ包
add_custom_target(stgz COMMAND cpack -G STGZ)

###########################################################
# 6. 测试
include(CTest)

# 添加测试程序
#add_executable(zzy_apctl example/apctl-jasminIF/driver.c)
#target_link_libraries(zzy_apctl JXPAMG)
#
#add_executable(solver_weak example/solver_weak.c)
#target_link_libraries(solver_weak JXPAMG)
#
#add_executable(solver_strong example/solver_strong.c)
#target_link_libraries(solver_strong JXPAMG)
#
add_executable(solver_laplace example/par_laplace.c example/solver_laplace.c)
target_link_libraries(solver_laplace JXPAMG)

# 添加测试用例
set(nplist 1 2 4 8) # mpi procs
set(ntlist 1 2 4)   # pthreads
#foreach (np ${nplist})
#  add_test(NAME test_zzy_apctl-np${np}
#    COMMAND ${mpi_ROOT}/bin/mpiexec -np ${np} ${PROJECT_BINARY_DIR}/zzy_apctl)
#endforeach ()
#
#foreach (np ${nplist})
#  add_test(NAME test_solver_weak-np${np}
#    COMMAND ${mpi_ROOT}/bin/mpiexec -np ${np} ${PROJECT_BINARY_DIR}/solver_weak)
#endforeach ()
#
#foreach (np ${nplist})
#  foreach (nt ${ntlist})
#    add_test(NAME test_solver_strong-np${np}-nts${nt}
#      COMMAND ${mpi_ROOT}/bin/mpiexec -np ${np} ${PROJECT_BINARY_DIR}/solver_strong -nts ${nt})
#  endforeach ()
#endforeach ()
#
foreach (np ${nplist})
  foreach (nt ${ntlist})
    add_test(NAME test_solver_laplace-np${np}-nts${nt}
      COMMAND ${mpi_ROOT}/bin/mpiexec -np ${np} ${PROJECT_BINARY_DIR}/solver_laplace -nts ${nt})
  endforeach ()
endforeach ()

