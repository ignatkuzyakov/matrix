cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)

project(MatrixTests LANGUAGES CXX)

include_directories(../googletest/googletest/include)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/../googletest/build/lib/libgtest.a ${CMAKE_SOURCE_DIR}/../googletest/build/lib/libgtest_main.a)

set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)