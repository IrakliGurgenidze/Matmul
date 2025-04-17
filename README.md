# Better Size Estimation for Sparse Matrix Products

This project implements a fast, probabilistic algorithm for estimating the number of non-zero entries in the product of two sparse boolean matrices. The algorithm is based on the 2011 paper [*Better size estimation for sparse matrix products*](https://arxiv.org/abs/1006.4173) by Amossen, Campagna, and Pagh.

We have also provided a series of sparse matrix formats and multiplication functions to take leverage this estimate for greater efficiency. Sample usage of all matrix classes and functions can be found in `src/main.cpp`.

## Project Layout
Matmul/ 
- `data/` – Optional .mtx files
- `include/` – Header files
  - `external/` – Third-party headers (e.g., hashing, Catch2)
- `src/` – Source files
- `tests/` – Unit tests and benchmark drivers
- `CMakeLists.txt` – Root build configuration
- `README.md` – Project documentation

## Build and Run Instructions

#### Requirements
- C++20 compiler (e.g. `g++ >= 10`, `clang++ >= 11`
- CMake 3.10 or higher

#### With CLion
Using CLion, simply cloning the repository and pressing `build`. CLion will automatically detect program structure and download required dependencies.

From there, run the built program from the GUI.

Otherwise,

#### Compilation (with CMake)
```
mkdir -p build
cd build
cmake ..
make
```
This will compile the project.

#### Running the Program (will execute `src/main.cpp`)
From the `build/src` directory, run `./Matmul`

## Usage
Please refer to `src/main.cpp` for basic usage. 
