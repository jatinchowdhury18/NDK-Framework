cmake_minimum_required(VERSION 3.15)
project(NDK-Framework-tests)
set(CMAKE_CXX_STANDARD 20)

include(cmake/CPM.cmake)

CPMAddPackage(
    NAME eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen
    GIT_TAG 969c31eefcdfaab11da763bea3f7502086673ab0
)

CPMAddPackage(
    NAME matplotlib-cpp
    GIT_REPOSITORY https://github.com/jatinchowdhury18/matplotlib-cpp
    GIT_TAG main
)

CPMAddPackage(
    NAME kissfft
    GIT_REPOSITORY https://github.com/mborgerding/kissfft
    GIT_TAG master
    OPTIONS "KISSFFT_DATATYPE double" "KISSFFT_STATIC ON" "KISSFFT_TEST OFF" "KISSFFT_TOOLS OFF" "KISSFFT_USE_ALLOCA ON"
)

add_subdirectory(cry_baby_test)
add_subdirectory(big_muff_test)
add_subdirectory(rc_lowpass_test)
