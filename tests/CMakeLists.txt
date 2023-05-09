# Add Google Test
add_subdirectory(googletest)

# Link to Google Test Library
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Add test executables
add_executable(run_tests 
test_GaussRev.cpp
#ProblemsKR2.cpp
test_IterMethods.cpp
#tests_CSR.cpp
test_Dense.cpp
tests_Tridiag.cpp
tests_vect_op.cpp 
test_GMRES.cpp
test_cheb.cpp
test_QR.cpp)


# Link test executable against Google Test and our library
#target_link_libraries(run_tests gtest gtest_main SLAE_methods)
target_link_libraries(run_tests gtest gtest_main SLAE_methods)