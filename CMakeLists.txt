cmake_minimum_required(VERSION 3.14)
project(glmx)

set(CMAKE_CXX_STANDARD 17)

add_library(glmx INTERFACE)
target_include_directories(glmx INTERFACE include)

add_executable(glmx_test tests/test_main.cpp tests/test_spatial.cpp tests/test_transform.cpp)
target_link_libraries(glmx_test PRIVATE glmx)


