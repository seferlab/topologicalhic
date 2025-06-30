cmake_minimum_required(VERSION 3.15)
project(HiCAnalysis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# Find required packages
find_package(Boost REQUIRED COMPONENTS system filesystem thread)
find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(TBB)

# Create executable
add_executable(hic_analysis main.cpp)

# Include directories for GUDHI - Updated for Cech Complex
target_include_directories(hic_analysis PRIVATE
    "${CMAKE_SOURCE_DIR}/gudhi-devel/src/Cech_complex/include"
    "${CMAKE_SOURCE_DIR}/gudhi-devel/src/Simplex_tree/include"
    "${CMAKE_SOURCE_DIR}/gudhi-devel/src/common/include"
    "${CMAKE_SOURCE_DIR}/gudhi-devel/src/Persistent_cohomology/include"
    "${CMAKE_SOURCE_DIR}/gudhi-devel/src/Bottleneck_distance/include"
    "${CMAKE_SOURCE_DIR}/gudhi-devel/src/Alpha_complex/include"
    "${CMAKE_SOURCE_DIR}/gudhi-devel/src/Spatial_searching/include"
    "${CMAKE_SOURCE_DIR}/gudhi-devel/src"
    ${Boost_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    ${CGAL_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(hic_analysis
    ${Boost_LIBRARIES}
    CGAL::CGAL
    Eigen3::Eigen
)

# Link TBB if available
if(TBB_FOUND)
    target_link_libraries(hic_analysis TBB::tbb)
endif()

# Compiler options
target_compile_options(hic_analysis PRIVATE
    -Wall -Wextra -O3 -DNDEBUG
)

# Preprocessor definitions
target_compile_definitions(hic_analysis PRIVATE
    GUDHI_USE_TBB
    WITH_GUDHI_USE_TBB
    CGAL_EIGEN3_ENABLED
)

# Add CGAL compile definitions
target_compile_definitions(hic_analysis PRIVATE ${CGAL_DEFINITIONS})

# Set additional compile flags for CGAL
if(CGAL_FOUND)
    target_compile_options(hic_analysis PRIVATE ${CGAL_CXX_FLAGS_INIT})
endif()