add_executable(matrix_read matrix_read.cpp)   
target_link_libraries(matrix_read PUBLIC Eigen3::Eigen)
target_link_libraries(matrix_read PUBLIC FFTWpp)
# find_package(OpenMP)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(matrix_read PUBLIC OpenMP::OpenMP_CXX)
# endif()
# message(STATUS "BLAS VENDOR:    " ${BLA_VENDOR})
# message(STATUS "BLAS LIBRARIES: " ${BLAS_LIBRARIES})
# message(STATUS "FFTW LIBRARIES: " ${FFTW_LIBRARIES})
# target_link_libraries(matrix_read PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(test test.cpp)   
target_link_libraries(test PUBLIC Eigen3::Eigen)
target_link_libraries(test PUBLIC FFTWpp)
target_link_libraries(test PUBLIC OpenMP::OpenMP_CXX)

message(STATUS "BLAS VENDOR:    " ${BLA_VENDOR})
message(STATUS "BLAS LIBRARIES: " ${BLAS_LIBRARIES})
message(STATUS "FFTW LIBRARIES: " ${FFTW_LIBRARIES})
target_link_libraries(test PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
# target_link_libraries(test PUBLIC lapack)
# target_link_libraries(test PUBLIC ${LAPACKE_LIBRARIES})

add_executable(testeig testeigsolve.cpp)   
target_link_libraries(testeig PUBLIC Eigen3::Eigen)
target_link_libraries(testeig PUBLIC FFTWpp)
target_link_libraries(testeig PUBLIC OpenMP::OpenMP_CXX)


target_link_libraries(testeig PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(testRK testRK.cpp)   
target_link_libraries(testRK PUBLIC Eigen3::Eigen)
target_link_libraries(testRK PUBLIC FFTWpp)
target_link_libraries(testRK PUBLIC OpenMP::OpenMP_CXX)


target_link_libraries(testeig PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})