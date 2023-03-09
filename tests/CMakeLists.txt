# Add Google Test
add_subdirectory(googletest)

# Link to Google Test Library
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Add test executable
#add_executable(run_tests tests_Tridiag.cpp)
#add_executable(run_tests tests_CSR.cpp Problems_3_4.cpp)
add_executable(run_tests Problem3.cpp)


# Link test executable against Google Test and our library
#target_link_libraries(run_tests gtest gtest_main SLAE_methods)
target_link_libraries(run_tests gtest gtest_main SLAE_methods)