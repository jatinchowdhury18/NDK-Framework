message(STATUS "Setting up Common Cathode test...")

add_executable(common_cathode_test main.cpp CommonCathodeNDK.cpp)

target_link_libraries(common_cathode_test
    PRIVATE
        Eigen3::Eigen
        matplotlib-cpp
)
