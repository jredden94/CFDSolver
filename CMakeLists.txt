cmake_minimum_required(VERSION 3.10)

project(CFD VERSION 1.0 LANGUAGES CXX)

include(FetchContent)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

FetchContent_Declare(XAD
    GIT_REPOSITORY https://github.com/auto-differentiation/xad.git
    GIT_TAG v1.6.0
    )
FetchContent_Declare(CODI
    GIT_REPOSITORY https://github.com/SciCompKL/CoDiPack.git
    GIT_TAG v2.2.0
    )

FetchContent_MakeAvailable(XAD)
FetchContent_MakeAvailable(CODI)


add_subdirectory(src)
