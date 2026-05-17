# Adapters for memory mapping and type conversion

set(DIFFUSION3D_ADAPTER_SOURCES
    memory_mapped_adapter.h
    memory_mapped_adapter.cpp
)

add_library(diffusion3d_adapters_objs OBJECT
    ${DIFFUSION3D_ADAPTER_SOURCES}
)

target_include_directories(diffusion3d_adapters_objs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
