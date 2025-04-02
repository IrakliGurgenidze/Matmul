include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.3.0 # Latest stable version as of now
)

# Make Catch2 available
FetchContent_MakeAvailable(Catch2)

# Add your test executable
add_executable(tests
        test_main.cpp
        # test_cardinality.cpp  # Add your test source files here
)

# Link Catch2 to your test executable
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

# Register tests with CTest
include(CTest)
include(Catch)
catch_discover_tests(tests)