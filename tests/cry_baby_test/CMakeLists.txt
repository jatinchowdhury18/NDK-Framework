message(STATUS "Setting up Cry Baby test...")

add_executable(cry_baby_test main.cpp CryBabyNDK.cpp)

target_link_libraries(cry_baby_test
    PRIVATE
        Eigen3::Eigen
        matplotlib-cpp
        kissfft
)
