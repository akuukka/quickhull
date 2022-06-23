cmake_minimum_required(VERSION 3.11)
project(QuickHull)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

cmake_policy(SET CMP0076 NEW)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

add_library(QuickHull SHARED)
target_include_directories(QuickHull
    PUBLIC
        inc
)
target_sources(QuickHull
	PUBLIC
        inc/QuickHull/ConvexHull.hpp
        inc/QuickHull/HalfEdgeMesh.hpp
        inc/QuickHull/MathUtils.hpp
        inc/QuickHull/QuickHull.hpp
        inc/QuickHull/Structs/Mesh.hpp
        inc/QuickHull/Structs/Plane.hpp
        inc/QuickHull/Structs/Pool.hpp
        inc/QuickHull/Structs/Ray.hpp
        inc/QuickHull/Structs/Vector3.hpp
        inc/QuickHull/Structs/VertexDataSource.hpp
        src/QuickHull.cpp
)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
