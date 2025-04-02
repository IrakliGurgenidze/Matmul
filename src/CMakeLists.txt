# src/CMakeLists.txt

set(SOURCES
        main.cpp
        Estimator.cpp
        HashUtils.cpp
        CoordListMatrix.cpp
        Types.cpp
        ${PROJECT_SOURCE_DIR}/include/external/MurmurHash3.cpp
)

add_executable(Matmul ${SOURCES})
