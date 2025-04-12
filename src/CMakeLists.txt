# src/CMakeLists.txt

set(SOURCES
        main.cpp
        Estimator.cpp
        HashUtils.cpp
        CoordListMatrix.cpp
        CSRMatrix.cpp
        Types.cpp
        ${PROJECT_SOURCE_DIR}/include/external/MurmurHash3.cpp
        CSRMatrix.cpp
        ../include/CSRMatrix.h
        ../include/HashContext.h
)

add_executable(Matmul ${SOURCES})
