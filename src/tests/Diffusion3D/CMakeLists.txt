# Catch2 tests for Diffusion3D (lives under src/tests/, links diffusion3d_core_lib)

set(DIFFUSION3D_TEST_SOURCES
    catch2_main.cpp
    convolutionFFT3D_test.cpp
    diffusion3d_step_euler_test.cpp
    diffusion3d_timestep_test.cpp
    engine_stepper_test.cpp
    test_diffusion3d.cpp
    support/analytic_diffusion3d.cpp
)

add_executable(diffusion3d_tests ${DIFFUSION3D_TEST_SOURCES})

target_include_directories(diffusion3d_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Diffusion3D/core
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Diffusion3D/core/gpu
    ${CMAKE_CURRENT_SOURCE_DIR}/support
)

target_link_libraries(diffusion3d_tests PRIVATE
    Catch2::Catch2
    diffusion3d_core_lib
)

set_target_properties(diffusion3d_tests PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

if(DIFFUSION3D_CUDA)
    target_compile_definitions(diffusion3d_tests PRIVATE DIFFUSION3D_CUDA=1)
    set_target_properties(diffusion3d_tests PROPERTIES
        CUDA_ARCHITECTURES "${DIFFUSION3D_CUDA_ARCHITECTURES}"
    )
endif()

add_test(NAME diffusion3d_modulateAndNormalize3D
         COMMAND diffusion3d_tests "modulateAndNormalize3D matches expected output")
add_test(NAME diffusion3d_spectrumRaisePowInt3D
         COMMAND diffusion3d_tests "spectrumRaisePowInt3D matches host complex powers")
add_test(NAME diffusion3d_self_convolution_matches_spectrum_power
         COMMAND diffusion3d_tests "Self-convolution: repeated one-step FFT matches spectrum power")
add_test(NAME diffusion3d_cpuDiffusion_vs_analyticGaussian
         COMMAND diffusion3d_tests "CPU diffusion is qualitatively close to analytic Gaussian (early time)")
add_test(NAME diffusion3d_unpadResult3D
         COMMAND diffusion3d_tests "unpadResult3D extracts the top-left-front volume")
add_test(NAME diffusion3d_padKernel3D
         COMMAND diffusion3d_tests "padKernel3D performs ifftshift-on-write")
add_test(NAME diffusion3d_padDataClampToBorder3D
         COMMAND diffusion3d_tests "padDataClampToBorder3D matches CUDA-sample border semantics")
add_test(NAME diffusion3d_fft3d_convolution_matches_cpu
         COMMAND diffusion3d_tests "3D FFT convolution matches CPU circular convolution (padded domain)")
add_test(NAME diffusion3d_fft3d_matches_cpu_laplacian
         COMMAND diffusion3d_tests "3D FFT convolution matches CPU 6-point Laplacian update (interior voxels)")
add_test(NAME diffusion3d_step_euler_cpu_determinism
         COMMAND diffusion3d_tests "diffusion3d_step_euler_scalar identical ICs stay aligned over multiple steps")
add_test(NAME diffusion3d_timestep_suite
         COMMAND diffusion3d_tests "[timestep]")
add_test(NAME diffusion3d_engine_suite
         COMMAND diffusion3d_tests "[engine]")
add_test(NAME diffusion3d_core_suite
         COMMAND diffusion3d_tests "[core]")

if(DIFFUSION3D_CUDA)
    add_test(NAME diffusion3d_step_euler_gpu_vs_cpu_chain
             COMMAND diffusion3d_tests "diffusion3d_step_euler_gpu chained matches CPU euler chain")
endif()
