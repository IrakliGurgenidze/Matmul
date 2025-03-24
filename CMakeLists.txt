cmake_minimum_required(VERSION 3.10)
project(Matmul)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include include/third_party)

# Add source and test directories
add_subdirectory(src)
add_subdirectory(tests)

# Enable testing support
enable_testing()

# Main executable
add_executable(Matmul
        src/main.cpp
        include/third_party/MurmurHash3.cpp
)