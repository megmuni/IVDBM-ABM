set(AGENT_SRC Agent.cpp)
file( GLOB USR_AGENT_SRC Usr_Agents/*.cpp)

set(lib_type SHARED)
add_library(Agent ${lib_type} ${AGENT_SRC} ${USR_AGENT_SRC})

target_link_libraries(Agent World)
#World is shared and is referenced by Agent! Without the last command, the linking process wouldnt work!

include_directories( ${PROJECT_SOURCE_DIR}/src)

## Added by Nuttiiya for atomic flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#added by Yvonna Jun 19th 2013
