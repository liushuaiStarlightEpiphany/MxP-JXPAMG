# 输入变量:
# - mpi_ROOT         依据这输入路径搜索需要的头文件和库文件,
#
# 输出变量:
# - mpi_FOUND        是否查找到mpi头文件及库
# - mpi_INCLUDE_DIRS mpi头文件路径
# - mpi_LIBRARIES    mpi库
# 模块说明:
#   1. 如果用户设置了mpicc或者mpicxx编译器, 则直接返回
#      用户可以通过环境变量CC,CXX,FC影响编译器的设置.
#      用户还可以设置CMAKE_{C|CXX|Fortran}_COMPILER来影响编译器
#   2. 如果用户仅设置了mpi_ROOT, 首先查找mpicc, mpicxx, mpif90
#      并将它们设为MPI_{C|CXX|Fortran}_COMPILER, 优先调用系统自带
#      的FindMPI.cmake来查找MPI库.
#   3. 如果cmake官方提供的FindMPI.cmake模块查找失败, 则尝试查找穷举
#      列表中的库.

# 1. 判断用户有没有自行设置mpi编译器, 如果是则返回
set(_IS_MPI_COMPILER FALSE)
foreach(lang C CXX)
  message(STATUS "CMAKE_${lang}_COMPILER : ${CMAKE_${lang}_COMPILER}")
  get_filename_component(_${lang}_Compiler_Name ${CMAKE_${lang}_COMPILER} NAME)
  string(REGEX MATCH "^mpi.+" _Is_MPI_${lang}_Compiler ${_${lang}_Compiler_Name})
  if (_Is_MPI_${lang}_Compiler)
    message(STATUS "User set ${CMAKE_${lang}_COMPILER} as MPI Compiler.")
    set(_IS_MPI_COMPILER TRUE)
  endif ()
endforeach ()
if (_IS_MPI_COMPILER)
  return()
endif ()

# 设置搜索路径
set(mpi_search_PATHS ${mpi_ROOT})
set(_OLD_PREFIX_PATH ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH ${mpi_search_PATHS})

# 首先搜索需要的头文件
# 不使用FindMPI.cmake返回的MPI_{lang}_INCLUDE_PATH是因为有些机器上的mpi
# 是从别的机器上直接拷贝安装的, mpicc, mpicxx, mpif90 -show给出的头文件
# 搜索路径是原机器的, 与目标机器的路径不符.
find_file(mpi_INCLUDE_DIRS NAMES include)
if (NOT mpi_INCLUDE_DIRS AND mpi_FIND_REQUIRED)
  message(FATAL_ERROR "Not found 'include' directory in ${mpi_search_PATHS}")
endif ()

# FIXME
# 2. 搜索mpicc, mpicxx, mpif90作为MPI_COMPILER, 调用cmake自带的FindMPI模块查找MPI库
message(STATUS "Try to find mpicc, mpicxx and mpif90 in ${CMAKE_PREFIX_PATH}")
find_program(MPI_C_COMPILER NAMES mpicc PATH_SUFFIXES bin)
find_program(MPI_CXX_COMPILER NAMES mpicxx PATH_SUFFIXES bin)
find_program(MPI_Fortran_COMPILER NAMES mpif90 PATH_SUFFIXES bin)
foreach(lang C CXX Fortran)
  message(STATUS "Found MPI_${lang}_COMPILER : ${MPI_${lang}_COMPILER}")
endforeach ()

message(STATUS "Try to find mpi libraries with FindMPI.cmake")
find_package(MPI)
foreach(lang C CXX Fortran)
  message(STATUS "MPI_${lang}_FOUND : ${MPI_${lang}_FOUND}")
endforeach ()

if (MPI_C_FOUND OR MPI_CXX_FOUND OR MPI_Fortran_FOUND)
  set(mpi_LIBRARIES ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
  list(REMOVE_DUPLICATES mpi_LIBRARIES)
  message(STATUS "Found mpi_LIBRARIES: ${mpi_LIBRARIES}")
else ()
  unset(mpi_LIBRARIES)
endif ()


# 3. 如果使用cmake自带的FindMPI.cmake查找失败, 再按照自己定义的去找
if (NOT mpi_LIBRARIES)
  # 设置需要搜索的mpi库
  if (NOT mpi_FIND_COMPONENTS)
    set(mpi_FIND_COMPONENTS
      mpichf90 mpichcxx mpich opa mpl dl rt pthread # mpich2
      fmpich mpichf90 mpich opa #openmpi
      mpi mpl         # mpi3
      mpigf mpi mpigi # intel impi
      )
  endif ()
  
  message(STATUS "Try to find mpi libraries with exhaustivity method.")
  message(STATUS "mpi_FIND_COMPONENTS : ${mpi_FIND_COMPONENTS}")
  foreach (module ${mpi_FIND_COMPONENTS})
    find_library(mpi_LIBRARIES_${module} NAMES ${module})
    if (mpi_LIBRARIES_${module})
      list(APPEND mpi_LIBRARIES ${mpi_LIBRARIES_${module}})
    else ()
      list(APPEND mpi_LIBRARIES_NOT_FOUND ${mpi_LIBRARIES_${module}})
    endif ()
  endforeach ()
  list(REMOVE_DUPLICATES mpi_LIBRARIES)
  list(REMOVE_DUPLICATES mpi_LIBRARIES_NOT_FOUND)
endif()

# 恢复CMAKE_PREFIX_PATH
set(CMAKE_PREFIX_PATH ${_OLD_PREFIX_PATH})

# 输出信息
if (NOT mpi_FIND_QUIETLY)
  message(STATUS "Found mpi_INLUDE_DIRS: ${mpi_INCLUDE_DIRS} in ${mpi_search_PATHS}")
  message(STATUS "Found mpi_LIBRARIES : ${mpi_LIBRARIES}")
  message(STATUS "NOT Found mpi_LIBRARIES : ${mpi_LIBRARIES_NOT_FOUND}")
endif ()

# 定义输出变量
if (mpi_INCLUDE_DIRS AND mpi_LIBRARIES)
  set(mpi_FOUND TRUE)
else ()
  set(mpi_FOUND FALSE)
endif ()

