set(DIFFUSION_SRC "convolutionFFT2D.cu")
file( GLOB USR_DIFFUSION_SRC "Usr_Diffusion/*.cu")

set(lib_type SHARED)
#add_library(Diffusion ${lib_type} ${DIFFUSION_SRC} ${USR_DIFFUSION_SRC})

#target_link_libraries(Diffusion World)

#include_directories( ${PROJECT_SOURCE_DIR}/src)


#set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE ON)

#set(BUILD_SHARED_LIBS OFF)

#set(CUDA_SEPARABLE_COMPILATION ON)
#list(APPEND CUDA_NVCC_FLAGS -arch=sm_20)

set(CUDA_SAMPLES_INC "/usr/local/cuda/cuda-7.0/samples/common/inc")
cuda_include_directories(${CUDA_SAMPLES_INC})

set(LIB_NAME "Diffusion")

cuda_add_library(${LIB_NAME}
  ${lib_type}
  ${DIFFUSION_SRC}
  ${USR_DIFFUSION_SRC}
  #OPTIONS -DSTUFF="blah blah"
  #RELEASE -DNDEBUG
  #DEBUG -g -DDEBUG
)

target_link_libraries(${LIB_NAME} World)

#set(PUBLIC_HEADERS "test.cuh")

#INSTALL(FILES ${PUBLIC_HEADERS} DESTINATION include)
#INSTALL(FILES "${CMAKE_BINARY_DIR}/lib${LIB_NAME}.a" DESTINATION lib)


