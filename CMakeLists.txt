cmake_minimum_required (VERSION 2.6)
project (NYVONNA_CPU_ABM)
include(CTest)
#project(ProjectName)
#SET(CMAKE_C_COMPILER gcc-7)
#SET(CMAKE_CXX_COMPILER g++-7)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/cuda" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

message("C_COMPILER = ${CMAKE_C_COMPILER}")
message("CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
#include(CMakeForceCompiler)
#CMAKE_FORCE_C_COMPILER(icc "Intel C Compiler")
#CMAKE_FORCE_CXX_COMPILER(icpc "Intel C++ Compiler")
#project (NYVONNA_CPU_ABM)

message(" (1) cxx compiler = ${CMAKE_CXX_COMPILER}")

#added by Yvonna June 2 2013 to enable nvcc
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/cuda" ${CMAKE_MODULE_PATH})
# message("CMAKE_SIZEOF_VOID_P = ${CMAKE_SIZEOF_VOID_P}")
# message("CUDA_64_BIT_DEVICE_CODE = ${CUDA_64_BIT_DEVICE_CODE}")

message(" (2) cxx compiler = ${CMAKE_CXX_COMPILER}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

# find_package(CUDA)
find_package(CUDA QUIET REQUIRED)
# find_package(CUDA 2.1)
# find_package(CUDA 2.2)
# find_package(CUDA 2.3)
# find_package(CUDA 2.1 REQUIRED)
# find_package(CUDA 2.2 REQUIRED)
# find_package(CUDA 2.3 REQUIRED)
# find_package(CUDA 2.1 EXACT)
# find_package(CUDA 2.2 EXACT)
# find_package(CUDA 2.3 EXACT)
# find_package(CUDA 2.1 EXACT REQUIRED)
# find_package(CUDA 2.2 EXACT REQUIRED)
# find_package(CUDA 2.3 EXACT REQUIRED)
# find_package(CUDA 2.2 EXACT QUIET REQUIRED)

SET(CUDA_PROPAGATE_HOST_FLAGS OFF)

 if (CUDA_FOUND)
   message("CUDA found")
 else()
   message("CUDA not found")
 endif()


#added by Yvonna June 2 2013 to enable nvcc 
message(" (3) cxx compiler = ${CMAKE_CXX_COMPILER}")

#SET(CMAKE_BUILD_TYPE DEBUG) #to enable debugging. Change to release for production
SET(CMAKE_BUILD_TYPE RELEASE)


set(NYVONNA_LIBS World Agent FieldVariable Patch ECM ArrayChain Utilities Diffusion)


#added "Patch" as part of library - Yvonna, Jun 6 2013

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH}; src)

#####################
#
#  MPI
#####################
#find_package ( MPI )
message(" --->")
message(" CXX is set to ${CXX}")
message(" cxx_compiler is set to ${CMAKE_CXX_COMPILER}")


#set(CUDA_SAMPLES_INC ${CUDA_TOOLKIT_ROOT_DIR}/samples/common/inc)
set(CUDA_SAMPLES_INC /home/meghana7/projects/def-nicoleli/meghana7/cuda-samples/Common)
message(STATUS "CUDA_SAMPLES_INC = ${CUDA_SAMPLES_INC}")
cuda_include_directories(${CUDA_SAMPLES_INC})
include_directories(${CUDA_SAMPLES_INC})

# Add your custom CUDA samples include path
include_directories(/home/meghana7/projects/def-nicoleli/meghana7/cuda-samples)
include_directories(/home/meghana7/projects/def-nicoleli/meghana7/cuda-samples/common/inc)


# Set CUDA NVCC flags to support compute capability 9.0 (for NVIDIA H100 on Nibi - Aug 2025)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11; -gencode arch=compute_90,code=sm_90 "-I/home/meghana7/projects/def-nicoleli/meghana7/cuda-samples" "-I/home/meghana7/projects/def-nicoleli/meghana7/cuda-samples/common/inc")

# Set CUDA NVCC flags to support compute capability 6.0 (for Tesla P100)
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11; -gencode arch=compute_60,code=sm_60 "-I/home/meghana7/projects/def-nicoleli/meghana7/cuda-samples" "-I/home/meghana7/projects/def-nicoleli/meghana7/cuda-samples/common/inc")
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11;-gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61) 
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11;-gencode arch=compute_70,code=sm_70 -gencode arch=compute_72,code=sm_72) 
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11; -gencode arch=compute_70,code=sm_70 -gencode arch=compute_72,code=sm_72 "-I/home/meghana7/projects/def-nicoleli/meghana7/cuda-samples" "-I/home/meghana7/projects/def-nicoleli/meghana7/cuda-samples/common/inc")
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11;-gencode arch=compute_37,code=sm_37 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_61,code=sm_61)
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35)
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11;-gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61)


set (CUDA_TARGET_LIB ${CUDA_TOOLKIT_ROOT_DIR}/include)
message(STATUS "CUDA_TARGET_LIB = ${CUDA_TARGET_LIB}")
include_directories(${CUDA_TARGET_LIB})

add_subdirectory(src)

add_executable(Nyvonna main.cpp)
add_executable(ex_diffuse_temperature examples/ex_diffuse_temperature.cpp)
add_executable(testRun examples/test_not_main.cpp)

message("here 1")

set( allLIBS ${NYVONNA_LIBS} cufft cudart)


target_link_libraries(Nyvonna ${allLIBS})
target_link_libraries(ex_diffuse_temperature ${allLIBS})
target_link_libraries(testRun ${allLIBS})

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lGL -lGLU -lglut")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lGL -lGLU -lglut")

## Added by Nuttiiya for OMP
#find_package(OpenMP)
#if(OPENMP_FOUND)
#  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

## Added by Nuttiiya to run profiler
## Compiler Flags
if(CMAKE_COMPILER_IS_GNUCXX)
  #message("Setting -lpthread flag")
  #set(CMAKE_CXX_FLAGS "-lpthread")
endif()


