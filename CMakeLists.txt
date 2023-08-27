cmake_minimum_required(VERSION 3.10)
project(Matrix)

find_package(GTest REQUIRED)

add_library(matrixLib INTERFACE)
target_include_directories(matrixLib INTERFACE matrix)

enable_testing()

add_executable(tests tests.cpp)

target_include_directories(tests PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(tests PRIVATE ${GTEST_BOTH_LIBRARIES})
target_link_libraries(tests PRIVATE matrixLib)

gtest_discover_tests(tests)
