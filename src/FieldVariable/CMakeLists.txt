set(FV_SRC FieldVariable.cpp)
file( GLOB USR_FV_SRC Usr_FieldVariables/*.cpp)

set(fv_lib_type SHARED)
message("glob of FV is ${USR_FV_SRC}")
add_library(FieldVariable ${fv_lib_type} ${FV_SRC} ${USR_FV_SRC})

## Added by Nuttiiya for atomic flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

