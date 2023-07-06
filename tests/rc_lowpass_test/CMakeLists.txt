message(STATUS "Setting up RC-Lowpass test...")

add_executable(rc_lowpass_test main.cpp RCLowpassNDK.cpp)

target_link_libraries(rc_lowpass_test
    PRIVATE
        Eigen3::Eigen
        matplotlib-cpp
        kissfft
)
