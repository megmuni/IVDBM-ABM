set(ARRAY_CHAIN_SRC ArrayChain.cpp)
file( GLOB USR_ARRAY_CHAIN_SRC Usr_ArrayChain/*.cpp)

set(lib_type SHARED)
add_library(ArrayChain ${lib_type} ${ARRAY_CHAIN_SRC} ${USR_ARRAY_CHAIN_SRC})

target_link_libraries(ArrayChain World)

include_directories( ${PROJECT_SOURCE_DIR}/src)

## Added by Nuttiiya for atomic flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
