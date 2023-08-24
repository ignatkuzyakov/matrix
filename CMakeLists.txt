cmake_minimum_required(VERSION 3.10)
project(Matrix)

enable_testing()
find_package(GTest REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main GTest::gtest GTest::gtest_main)


gtest_discover_tests(main)
