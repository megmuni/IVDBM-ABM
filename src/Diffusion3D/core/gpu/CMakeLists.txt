# GPU / FFT numerical kernels (optional CUDA)

set(DIFFUSION3D_GPU_INCLUDES
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

if(DIFFUSION3D_CUDA)
    enable_language(CUDA)
    if(NOT TARGET CUDA::cudart)
        find_package(CUDAToolkit REQUIRED)
    endif()
    add_library(diffusion3d_core_gpu OBJECT
        gpu_diffusion3d.cu
        convolutionFFT3D.cu
        diffusion3d_fft_scratch.cu
    )
    target_compile_definitions(diffusion3d_core_gpu PUBLIC DIFFUSION3D_CUDA=1)
    target_link_libraries(diffusion3d_core_gpu PUBLIC CUDA::cudart CUDA::cufft)
    set_target_properties(diffusion3d_core_gpu PROPERTIES
        CUDA_ARCHITECTURES "${DIFFUSION3D_CUDA_ARCHITECTURES}"
        CXX_STANDARD 17
        CUDA_STANDARD 17
    )
else()
    add_library(diffusion3d_core_gpu OBJECT
        gpu_diffusion3d_stub.cpp
        diffusion3d_fft_scratch_stub.cpp
    )
endif()

target_include_directories(diffusion3d_core_gpu PUBLIC ${DIFFUSION3D_GPU_INCLUDES})
target_link_libraries(diffusion3d_core_gpu PUBLIC diffusion3d_core)
# The GPU code is not thread-safe, so we need to mark it as position-independent code.
# This is a hack to get around the fact that the GPU code is not thread-safe.
# It is not a good solution, but it is the best we can do for now.
set_target_properties(diffusion3d_core_gpu PROPERTIES POSITION_INDEPENDENT_CODE ON)
