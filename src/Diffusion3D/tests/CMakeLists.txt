# Catch2 tests for Diffusion3D (legacy numerical kernels + Phase I core placeholders)
# Currently testing the demo code! We will move the LEGACY_DIFFUSION3D_SRC to the new src soon.

set(DIFFUSION3D_TEST_SOURCES
    catch2_main.cpp
    convolutionFFT3D_test.cpp
    diffusion3d_step_euler_test.cpp
    diffusion3d_timestep_test.cpp
    diffusion3d_solver_test.cpp
    test_diffusion3d.cpp
    ${LEGACY_DIFFUSION3D_SRC}/cpu_diffusion3d.cpp
    ${LEGACY_DIFFUSION3D_SRC}/analytic_diffusion3d.cpp
    ${LEGACY_DIFFUSION3D_SRC}/diffusion3d_timestep.cpp
    ${LEGACY_DIFFUSION3D_SRC}/diffusion3d_solver.cpp
)

if(DIFFUSION3D_CUDA)
    find_package(CUDAToolkit REQUIRED)
    list(APPEND DIFFUSION3D_TEST_SOURCES
        ${LEGACY_DIFFUSION3D_SRC}/gpu_diffusion3d.cu
        ${LEGACY_DIFFUSION3D_SRC}/convolutionFFT3D.cu
        ${LEGACY_DIFFUSION3D_SRC}/diffusion3d_fft_scratch.cu
    )
else()
    list(APPEND DIFFUSION3D_TEST_SOURCES
        ${LEGACY_DIFFUSION3D_SRC}/gpu_diffusion3d_stub.cpp
        ${LEGACY_DIFFUSION3D_SRC}/diffusion3d_fft_scratch_stub.cpp
    )
endif()

add_executable(diffusion3d_tests ${DIFFUSION3D_TEST_SOURCES})

target_include_directories(diffusion3d_tests PRIVATE
    ${LEGACY_DIFFUSION3D_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/../core
)

target_link_libraries(diffusion3d_tests PRIVATE
    Catch2::Catch2
    diffusion3d_core_objs
)

set_target_properties(diffusion3d_tests PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

if(DIFFUSION3D_CUDA)
    target_compile_definitions(diffusion3d_tests PRIVATE GPU_DIFFUSE=1)
    target_link_libraries(diffusion3d_tests PRIVATE CUDA::cudart CUDA::cufft)
    set_target_properties(diffusion3d_tests PROPERTIES CUDA_ARCHITECTURES native)
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
add_test(NAME diffusion3d_core_suite
         COMMAND diffusion3d_tests "[core]")

if(DIFFUSION3D_CUDA)
    add_test(NAME diffusion3d_step_euler_gpu_vs_cpu_chain
             COMMAND diffusion3d_tests -c "diffusion3d_step_euler_gpu chained matches CPU euler chain")
endif()
