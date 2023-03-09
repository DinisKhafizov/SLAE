cmake_minimum_required(VERSION 3.0.0)
project(tridiagonalmatrix VERSION 0.1.0)

add_subdirectory(src)


OPTION(TESTS "Enable testing of the project" ON)
if (TESTS)
    add_compile_definitions(COMPILE_WITH_TESTS)
    add_subdirectory(tests)
    message("TESTING ENABLED")
else()
    MESSAGE("TESTING DISABLED")
endif()

