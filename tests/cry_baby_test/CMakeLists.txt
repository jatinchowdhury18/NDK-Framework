message(STATUS "Setting up Cry Baby test...")

add_executable(cry_baby_test main.cpp CryBabyNDK.cpp)

CPMAddPackage(
    NAME kissfft
    GIT_REPOSITORY https://github.com/mborgerding/kissfft
    GIT_TAG master
    OPTIONS "KISSFFT_DATATYPE double" "KISSFFT_STATIC ON" "KISSFFT_TEST OFF" "KISSFFT_TOOLS OFF" "KISSFFT_USE_ALLOCA ON"
)

target_link_libraries(cry_baby_test
    PRIVATE
        Eigen3::Eigen
        matplotlib-cpp
        kissfft
)
