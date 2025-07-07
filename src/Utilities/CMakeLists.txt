set(UTILITIES_SRC timer.cpp)
file( GLOB USR_UTILITIES_SRC Usr_Utilities/*.cpp)

set(lib_type SHARED)
add_library(Utilities ${lib_type} ${UTILITIES_SRC} ${UTILITIES_SRC})

target_link_libraries(Utilities World)

include_directories( ${PROJECT_SOURCE_DIR}/src)

