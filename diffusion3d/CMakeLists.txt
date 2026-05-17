cmake_minimum_required(VERSION 3.18)

include(FetchContent)

option(GPU_DIFFUSE "Enable CUDA GPU diffusion code path" OFF)
option(LOGGING "Enable detailed logging output" OFF)
option(DIFFUSION3D_TESTS "Build diffusion3d tests (CTest)" ON)

if(GPU_DIFFUSE)
    project(diffusion3d LANGUAGES CXX CUDA)
else()
    project(diffusion3d LANGUAGES CXX)
endif()

message(STATUS "Configuring diffusion3d project")

if(GPU_DIFFUSE)
    find_package(CUDAToolkit REQUIRED)
endif()

set(DIFFUSION3D_SOURCES
    src/main.cpp
    src/cpu_diffusion3d.cpp
    src/analytic_diffusion3d.cpp
    src/diffusion3d_timestep.cpp
    src/diffusion3d_solver.cpp
    src/diffusion3d_common.h
    src/diffusion3d_cuda_check.h
    src/gpu_diffusion3d.cuh
    src/logger.h
    src/vtk_export.h
    src/convolutionFFT3D_common.h
    src/analytic_diffusion3d.h
)

if(GPU_DIFFUSE)
    list(APPEND DIFFUSION3D_SOURCES
        src/gpu_diffusion3d.cu
        src/convolutionFFT3D.cu
        src/diffusion3d_fft_scratch.cu
    )
else()
    list(APPEND DIFFUSION3D_SOURCES
        src/gpu_diffusion3d_stub.cpp
        src/diffusion3d_fft_scratch_stub.cpp
    )
endif()

add_executable(diffusion3d ${DIFFUSION3D_SOURCES})
target_include_directories(diffusion3d PRIVATE src)
set_target_properties(diffusion3d PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

if(GPU_DIFFUSE)
    target_compile_definitions(diffusion3d PRIVATE GPU_DIFFUSE=1)
    target_link_libraries(diffusion3d PRIVATE CUDA::cudart CUDA::cufft)
    set_target_properties(diffusion3d PROPERTIES CUDA_ARCHITECTURES native)
endif()

if(LOGGING)
    message(STATUS "Logging is enabled")
    target_compile_definitions(diffusion3d PRIVATE LOGGING=1)
endif()

if(DIFFUSION3D_TESTS)
    enable_testing()

    # Catch2 v2.x, header-only style include: <catch2/catch.hpp>
    FetchContent_Declare(
        Catch2
        URL https://github.com/catchorg/Catch2/archive/refs/tags/v2.13.10.tar.gz
    )
    FetchContent_MakeAvailable(Catch2)

    add_executable(diffusion3d_tests
        tests/catch2_main.cpp
        tests/convolutionFFT3D_test.cpp
        tests/diffusion3d_step_euler_test.cpp
        tests/diffusion3d_timestep_test.cpp
        tests/diffusion3d_solver_test.cpp
        src/convolutionFFT3D_common.h
        src/cpu_diffusion3d.cpp
        src/diffusion3d_timestep.cpp
        src/diffusion3d_solver.cpp
        src/diffusion3d_solver.h
        src/diffusion3d_common.h
        src/diffusion3d_cuda_check.h
        src/analytic_diffusion3d.cpp
        src/analytic_diffusion3d.h
    )
    target_include_directories(diffusion3d_tests PRIVATE src)
    target_link_libraries(diffusion3d_tests PRIVATE Catch2::Catch2)

    if(GPU_DIFFUSE)
        target_sources(diffusion3d_tests PRIVATE
            src/convolutionFFT3D.cu
            src/gpu_diffusion3d.cu
            src/diffusion3d_fft_scratch.cu
        )
        target_compile_definitions(diffusion3d_tests PRIVATE GPU_DIFFUSE=1)
        target_link_libraries(diffusion3d_tests PRIVATE CUDA::cudart CUDA::cufft)
        set_target_properties(diffusion3d_tests PROPERTIES CUDA_ARCHITECTURES native)
    else()
        target_sources(diffusion3d_tests PRIVATE
            src/gpu_diffusion3d_stub.cpp
            src/diffusion3d_fft_scratch_stub.cpp
        )
    endif()

    # Register individual Catch2 test cases as separate CTest entries.
    add_test(NAME diffusion3d_modulateAndNormalize3D
             COMMAND diffusion3d_tests -c "modulateAndNormalize3D matches expected output")
    add_test(NAME diffusion3d_spectrumRaisePowInt3D
             COMMAND diffusion3d_tests -c "spectrumRaisePowInt3D matches host complex powers")
    add_test(NAME diffusion3d_self_convolution_matches_spectrum_power
             COMMAND diffusion3d_tests -c "Self-convolution: repeated one-step FFT matches spectrum power")
    add_test(NAME diffusion3d_cpuDiffusion_vs_analyticGaussian
             COMMAND diffusion3d_tests -c "CPU diffusion is qualitatively close to analytic Gaussian (early time)")
    add_test(NAME diffusion3d_unpadResult3D
             COMMAND diffusion3d_tests -c "unpadResult3D extracts the top-left-front volume")
    add_test(NAME diffusion3d_padKernel3D
             COMMAND diffusion3d_tests -c "padKernel3D performs ifftshift-on-write")
    add_test(NAME diffusion3d_padDataClampToBorder3D
             COMMAND diffusion3d_tests -c "padDataClampToBorder3D matches CUDA-sample border semantics")
    add_test(NAME diffusion3d_fft3d_convolution_matches_cpu
             COMMAND diffusion3d_tests -c "3D FFT convolution matches CPU circular convolution (padded domain)")
    add_test(NAME diffusion3d_fft3d_matches_cpu_laplacian
             COMMAND diffusion3d_tests -c "3D FFT convolution matches CPU 6-point Laplacian update (interior voxels)")
    add_test(NAME diffusion3d_step_euler_cpu_determinism
             COMMAND diffusion3d_tests -c "diffusion3d_step_euler_cpu identical ICs stay aligned over multiple steps")
    add_test(NAME diffusion3d_timestep_suite
             COMMAND diffusion3d_tests "[timestep]")
    add_test(NAME diffusion3d_solver_suite
             COMMAND diffusion3d_tests "[solver]")
    if(GPU_DIFFUSE)
        add_test(NAME diffusion3d_step_euler_gpu_vs_cpu_chain
                 COMMAND diffusion3d_tests -c "diffusion3d_step_euler_gpu chained matches CPU euler chain")
    endif()
endif()