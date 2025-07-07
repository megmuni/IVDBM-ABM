set(ECM_SRC ECM.cpp)
file( GLOB USR_ECM_SRC Usr_ECM/*.cpp)

set(lib_type SHARED)
add_library(ECM ${lib_type} ${ECM_SRC} ${USR_ECM_SRC})

target_link_libraries(ECM World)

include_directories( ${PROJECT_SOURCE_DIR}/src)

## Added by Nuttiiya for atomic flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
