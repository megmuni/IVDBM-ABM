cmake_minimum_required (VERSION 2.6)
project (Nyvonna_CPU_ABM)


add_subdirectory(Agent)
add_subdirectory(FieldVariable)
add_subdirectory(World)
add_subdirectory(Patch)
add_subdirectory(ECM)
add_subdirectory(ArrayChain)
add_subdirectory(Utilities)
add_subdirectory(Diffusion)

# Note: Will not hide this behind a flag until Diffusion3D is fully complete
option(BUILD_DIFFUSION3D "Build Diffusion3D module and Catch2 tests" ON)
if(BUILD_DIFFUSION3D)
    add_subdirectory(Diffusion3D)
endif()

#added the subdir - Yvonna Jun 6 2013
