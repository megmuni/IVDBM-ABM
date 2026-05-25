set(WORLD_SRC World.cpp PatchFieldDiffusion.cpp)
file( GLOB USR_WORLD_SRC Usr_World/*.cpp)

set(world_lib_type SHARED)
add_library(World ${world_lib_type} ${WORLD_SRC} ${USR_WORLD_SRC})

target_include_directories(World PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(World FieldVariable Chemistry)

if(BUILD_DIFFUSION3D)
    target_link_libraries(World diffusion3d_core_lib)
    set_source_files_properties(PatchFieldDiffusion.cpp PROPERTIES COMPILE_FLAGS "-std=c++17")
endif()

## Added by Nuttiiya for atomic flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# Added by Nuttiiya
#find_package(OpenMP)
#if (OPENMP_FOUND)
#  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

# Added by Nuttiiya
## Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
#  message("        Setting -pg flag in src/World/")
#  set(CMAKE_CXX_FLAGS "-pg")
endif()
