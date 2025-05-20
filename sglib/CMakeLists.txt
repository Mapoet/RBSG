cmake_minimum_required(VERSION 3.12)
project(SGLib VERSION 0.2.0 LANGUAGES Fortran)

# Options
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_TESTING "Build the testing tree" ON)

# Find required packages
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Set Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -g -Wall -Wextra -pedantic")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g3 -Wall -Wextra -pedantic -fcheck=all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-loops")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -check all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo")
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Library source files
set(LIB_SOURCES
    src/sg_coeff.f90
    src/sg_filter.f90
)

# Test program source files
set(TEST_SOURCES
    src/main_test.f90
)

# Custom data test program
set(CUSTOM_TEST_SOURCES
    src/test_sg_filter_data.f90
)

# Add library target
add_library(sglib ${LIB_SOURCES})
target_link_libraries(sglib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Set public headers
set_target_properties(sglib PROPERTIES
    PUBLIC_HEADER "include/sglib_mod.f90"
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Add test executable
add_executable(test_sg_filter ${TEST_SOURCES})
target_link_libraries(test_sg_filter sglib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Add custom data test executable
add_executable(test_sg_filter_data ${CUSTOM_TEST_SOURCES})
target_link_libraries(test_sg_filter_data sglib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Installation targets
install(TARGETS sglib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include
)

install(TARGETS test_sg_filter test_sg_filter_data
    RUNTIME DESTINATION bin
)

# Add testing
if(BUILD_TESTING)
    enable_testing()
    add_test(NAME RunBasicTest 
             COMMAND test_sg_filter)
    add_test(NAME GenerateData
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_test_data.py -n 200 --noise 0.1 --outlier-prob 0.05)
    add_test(NAME RunDataTest
             COMMAND test_sg_filter_data -i test_data.dat -o test_results.dat -w 7 -p 1 -r 3.0)
endif()

# Print summary
message(STATUS "")
message(STATUS "SGLib configuration summary:")
message(STATUS "  Fortran compiler:      ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  Build type:            ${CMAKE_BUILD_TYPE}")
message(STATUS "  Build shared libs:     ${BUILD_SHARED_LIBS}")
message(STATUS "  Build testing:         ${BUILD_TESTING}")
message(STATUS "  Install prefix:        ${CMAKE_INSTALL_PREFIX}")
message(STATUS "") 