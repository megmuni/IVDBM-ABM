set(PATCH_SRC Patch.cpp)
file( GLOB USR_PATCH_SRC Usr_Patch/*.cpp)

add_library(Patch ${PATCH_SRC} ${USR_PATCH_SRC})

## Added by Nuttiiya for atomic flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
