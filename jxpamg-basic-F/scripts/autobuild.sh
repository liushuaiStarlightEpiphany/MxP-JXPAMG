#!/bin/bash

NP=$(cat /proc/cpuinfo | grep processor | wc -l)
SCRIPT=$(basename ${BASH_SOURCE})
SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE}); pwd)
SOURCE_ROOT=$(cd $(dirname $SCRIPT_DIR); pwd)
BUILD_ROOT=$SOURCE_ROOT/build     # 缺省编译目录
INSTALL_PREFIX=$BUILD_ROOT/jxfpamg # 缺省安装目录
mpi_ROOT=/usr/local/jasmin/thirdparty/mpich2-1.5
TARGET=pack    # 默认情况下, 执行cmake, make, make package
BUILD_TYPE=rel # 默认情况下, 编译类型设为Release

function usage()
{
  cat << EOF
USAGE
  $SCRIPT <-T|--target <cmake|make|test|install|pack>> [OPTIONS]
      cmake            只执行cmake
      make             执行cmake, make
      test             执行cmake, make, make test
      install          执行cmake, make, make install
      pack             执行cmake, make, make package

DEPENDS
  --with-CC  <XXX>     设置C编译器, 缺省值gcc
  --with-CXX <XXX>     设置CXX编译器, 缺省值: g++
  --with-FC  <XXX>     设置Fortran编译器, 缺省值: gfortran
  --with-mpi <DIR>     设置mpi库的根目录, 缺省值: $mpi_ROOT

OPTIONS
  -s|--source  <dir>   设置源码目录, 缺省值: ${SOURCE_ROOT}
  -b|--build   <dir>   设置编译目录, 缺省值: ${BUILD_ROOT}
  -i|--install <dir>   设置安装目录, 缺省值: ${INSTALL_PREFIX}
  -bt|--buildtype <rel|dbg|reldbg|minrel> 设置CMAKE_BUILD_TYPE, 缺省值: rel
      rel              CMAKE_BUILD_TYPE=RELEASE
      dbg              CMAKE_BUILD_TYPE=DEBUG
      reldbg           CMAKE_BUILD_TYPE=RELWITHDEBINFO
      minrel           CMAKE_BUILD_TYPE=MINSIZEREL
  -e|--extra_cmake_vars <str> 传给cmake的其它命令行参数, 多个参数用空格
                       分开并用引号引起来, 例如: "-DCPACK_PACKAGE_NAME=JXFPAMG"
  -h|--help            打印帮助信息并退出
EOF
}


function check_source_root()
{
  if [ ! -d ${SOURCE_ROOT} ]; then
    cat << EOF

Fatal error: SOURCE_ROOT=${SOURCE_ROOT} not exist!
Please check -s|--source set properly or not!

EOF
    usage; exit
  fi

  if [ ! -e ${SOURCE_ROOT}/CMakeLists.txt ]; then
    cat << EOF

Fatal error: invalid SOURCE_ROOT - ${SOURCE_ROOT}
Please check -s|--source set properly or not!

EOF
    usage; exit
  fi
}

function create_build_root_if_not_exists()
{
  if [ ! -d ${BUILD_ROOT} ]; then
    mkdir -p ${BUILD_ROOT} || exit
  else
    \rm -rf ${BUILD_ROOT}/*
  fi
}

function check_autobuild_target()
{
  case ${TARGET} in
    cmake|make|test|install|pack)
      echo "-- do ${TARGET}"
      ;;
    *)
      usage; exit
      ;;
  esac
}

function check_cmake_build_type()
{
  case ${BUILD_TYPE} in
    rel)
      CMAKE_BUILD_TYPE="RELEASE"
      ;;
    dbg)
      CMAKE_BUILD_TYPE="DEBUG"
      ;;
    reldbg)
      CMAKE_BUILD_TYPE="RELWITHDEBINFO"
      ;;
    minrel)
      CMAKE_BUILD_TYPE="MINSIZEREL"
      ;;
    *)
      usage; exit
      ;;
  esac
}

function execute_cmake_command()
{
  # compiler
  COMPILE_ENV=""
  if [ "x${C_COMPILER}x" != "xx" ]; then COMPILE_ENV="CC=${C_COMPILER}"; fi
  if [ "x${CXX_COMPILER}x" != "xx" ]; then COMPILE_ENV="${COMPILE_ENV} CXX=${CXX_COMPILER}"; fi
  if [ "x${Fortran_COMPILER}x" != "xx" ]; then COMPILE_ENV="${COMPILE_ENV} FC=${Fortran_COMPILER}"; fi

  # cmake vars
  CMAKE_VARS="-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} -Dmpi_ROOT=${mpi_ROOT}"

  export ${COMPILE_ENV}; cmake ${CMAKE_VARS} ${EXTRA_CMAKE_VARS} ${SOURCE_ROOT}
}

function run_autobuild_target()
{
  pushd $BUILD_ROOT &> /dev/null
    
  execute_cmake_command

  case ${TARGET} in
    make)
      make -j $NP
      ;;
    test)
      make -j $NP
      make test
      ;;
    install)
      make -j $NP
      make install
      ;;
    pack)
      make -j $NP
      make stgz
      ;;
  esac

  popd &> /dev/null
}


# script entry
if [ $# -eq 0 ]; then
  usage; exit
fi
while [ $# -ge 1 ]
do
  case $1 in
    -T|--target)
      TARGET=$2
      shift 2
      ;;
    --with-CC)
      C_COMPILER=$2
      shift 2
      ;;
    --with-CXX)
      CXX_COMPILER=$2
      shift 2
      ;;
    --with-FC)
      Fortran_COMPILER=$2
      shift 2
      ;;
    --with-mpi)
      mpi_ROOT=$2
      shift 2
      ;;
    -s|--source)
      SOURCE_ROOT=$2
      shift 2
      ;;
    -b|--build)
      BUILD_ROOT=$2
      shift 2
      ;;
    -i|--install)
      INSTALL_PREFIX=$2
      shift 2
      ;;
    -bt|--buildtype)
      BUILD_TYPE=$2
      shift 2
      ;;
    -e|--extra_cmake_vars)
      EXTRA_CMAKE_VARS="$2"
      shift 2
      ;;
    -h|--help|*)
      usage
      exit
      ;;
  esac
done

# 检查SOURCE_ROOT是否正确设置
check_source_root

# 检查BUILD_ROOT
create_build_root_if_not_exists

# 检查autobuild目标
check_autobuild_target

# 设置CMAKE_BUILD_TYPE
check_cmake_build_type

# 执行编译构建
run_autobuild_target

