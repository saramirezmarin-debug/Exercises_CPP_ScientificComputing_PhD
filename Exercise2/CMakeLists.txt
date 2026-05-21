cmake_minimum_required(VERSION 3.16)

project(Exercise2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Eigen
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen-5.0.1")

# {fmt}
set(LOCAL_FMT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/fmt-12.1.0")
set(LOCAL_FMT_BUILD_DIR "${LOCAL_FMT_DIR}/build")

if(EXISTS "${LOCAL_FMT_BUILD_DIR}/fmt-config.cmake")
  find_package(fmt CONFIG REQUIRED
    PATHS "${LOCAL_FMT_BUILD_DIR}"
    NO_DEFAULT_PATH
  )
else()
  add_subdirectory("${LOCAL_FMT_DIR}" "${CMAKE_BINARY_DIR}/fmt-build")
endif()

add_executable(BVP
    src/BVP.cpp
    include/BVP.hh
)

target_include_directories(BVP PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
)
target_link_libraries(BVP PRIVATE fmt::fmt)


if(UNIX)
  target_compile_options(BVP PRIVATE -Wall -Wextra -Wpedantic)
endif()

message(STATUS "Eigen:     ${EIGEN3_INCLUDE_DIR}")
message(STATUS "{fmt}:      ${LOCAL_FMT_DIR}")

