cmake_minimum_required(VERSION 3.11)
project(Matrix)

find_package(GTest REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(SRC main.cpp)

enable_testing()

set(UNIT_TARGET ${PROJECT_NAME}_unit)
set(TEST_TARGET ${PROJECT_NAME}_test)

add_executable(${UNIT_TARGET} ${SRC})
target_include_directories(${UNIT_TARGET} PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(${UNIT_TARGET} PRIVATE ${GTEST_BOTH_LIBRARIES})
target_link_libraries(${UNIT_TARGET} PRIVATE Threads::Threads)
target_compile_features(${UNIT_TARGET} PRIVATE cxx_std_11)

add_test(
    NAME ${TEST_TARGET}
    COMMAND ${UNIT_TARGET}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

set_tests_properties(${TEST_TARGET} PROPERTIES DEPENDS ${UNIT_TARGET})