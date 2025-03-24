cmake_minimum_required(VERSION 3.10)
project(Matmul)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Matmul main.cpp MurmurHash3.cpp)
