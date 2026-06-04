cmake_minimum_required (VERSION 2.6)
project (NYVONNA_CPU_ABM)
include(CTest)
#project(ProjectName)
#SET(CMAKE_C_COMPILER gcc-7)
#SET(CMAKE_CXX_COMPILER g++-7)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

message("C_COMPILER = ${CMAKE_C_COMPILER}")
message("CXX_COMPILER = ${CMAKE_CXX_COMPILER}")

#SET(CMAKE_BUILD_TYPE DEBUG) #to enable debugging. Change to release for production
SET(CMAKE_BUILD_TYPE RELEASE)

set(NYVONNA_LIBS World Agent FieldVariable Patch ECM ArrayChain Utilities)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH}; src)

message(" --->")
message(" CXX is set to ${CXX}")
message(" cxx_compiler is set to ${CMAKE_CXX_COMPILER}")

add_subdirectory(src)

add_executable(Nyvonna main.cpp)
if(EXISTS "${CMAKE_SOURCE_DIR}/examples/ex_diffuse_temperature.cpp")
  add_executable(ex_diffuse_temperature examples/ex_diffuse_temperature.cpp)
endif()
if(EXISTS "${CMAKE_SOURCE_DIR}/examples/test_not_main.cpp")
  add_executable(testRun examples/test_not_main.cpp)
endif()
if(EXISTS "${CMAKE_SOURCE_DIR}/examples/ex_chem_random_spikes.cpp")
  add_executable(ex_chem_random_spikes examples/ex_chem_random_spikes.cpp)
  target_compile_definitions(ex_chem_random_spikes PRIVATE
    IVDBM_CHEM_CONFIG_DIR="${CMAKE_SOURCE_DIR}/configFiles")
  set_target_properties(ex_chem_random_spikes PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON)
endif()

message("here 1")

set(allLIBS ${NYVONNA_LIBS})

target_link_libraries(Nyvonna ${allLIBS})
if(TARGET ex_diffuse_temperature)
  target_link_libraries(ex_diffuse_temperature ${allLIBS})
endif()
if(TARGET testRun)
  target_link_libraries(testRun ${allLIBS})
endif()
if(TARGET ex_chem_random_spikes)
  target_link_libraries(ex_chem_random_spikes Chemistry FieldVariable)
endif()

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


