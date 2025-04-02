cmake_minimum_required(VERSION 3.10)
project(Matmul)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include include/external)

add_subdirectory(src)
add_subdirectory(tests)

enable_testing()