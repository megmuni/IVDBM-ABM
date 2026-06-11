# src-level Catch2 tests (module tests live in subdirectories)

if(NOT TARGET Catch2::Catch2)
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        URL https://github.com/catchorg/Catch2/archive/refs/tags/v2.13.10.tar.gz
    )
    FetchContent_MakeAvailable(Catch2)
endif()

if(TARGET diffusion3d_core_lib)
    add_subdirectory(Diffusion3D)
endif()

if(TARGET Chemistry)
    add_subdirectory(Chemistry)
endif()

if(TARGET diffusion3d_core_lib AND TARGET Chemistry)
    add_subdirectory(World)
endif()
