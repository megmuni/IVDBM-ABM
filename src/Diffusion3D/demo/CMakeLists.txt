# ParaView demo executable (all three diffusion algorithms when CUDA is enabled)

add_executable(diffusion3d_demo main.cpp)

target_include_directories(diffusion3d_demo PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../core
    ${CMAKE_CURRENT_SOURCE_DIR}/../core/gpu
)

target_link_libraries(diffusion3d_demo PRIVATE diffusion3d_core_lib)

set_target_properties(diffusion3d_demo PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

if(DIFFUSION3D_CUDA)
    target_compile_definitions(diffusion3d_demo PRIVATE DIFFUSION3D_CUDA=1)
    set_target_properties(diffusion3d_demo PROPERTIES
        CUDA_ARCHITECTURES "${DIFFUSION3D_CUDA_ARCHITECTURES}"
    )
endif()
