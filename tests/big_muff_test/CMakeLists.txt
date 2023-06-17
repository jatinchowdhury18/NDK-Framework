message(STATUS "Setting up Big Muff test...")

add_executable(big_muff_test main.cpp BigMuff2D.cpp BigMuffDP.cpp)

target_link_libraries(big_muff_test
    PRIVATE
        Eigen3::Eigen
        matplotlib-cpp
)
