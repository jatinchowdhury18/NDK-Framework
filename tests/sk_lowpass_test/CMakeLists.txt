message(STATUS "Setting up Sallen-Key Lowpass test...")

add_executable(sk_lowpass_test main.cpp SKLowpassNDK.cpp)

target_link_libraries(sk_lowpass_test
    PRIVATE
        Eigen3::Eigen
        matplotlib-cpp
        kissfft
)
