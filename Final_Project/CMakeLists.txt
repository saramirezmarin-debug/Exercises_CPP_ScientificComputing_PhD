cmake_minimum_required(VERSION 3.16)

set( TARGET CSC)

project(${TARGET} VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Eigen
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen-5.0.1")

add_executable(${TARGET}
    src/ODE4.cc
    src/main.cpp
    src/CSC_RL.cc
    src/Equilibrium.cc
    src/CSC_RL_Case.cc

    include/ODE4.hh
    include/CSC_RL.hh
    include/CSC_RL_Case.hh
    include/Equilibrium.hh
    include/Newton.hh
    include/dual.hh
    
)

target_include_directories(${TARGET} PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
)

if(UNIX)
  target_compile_options(${TARGET} PRIVATE -Wall -Wextra -Wpedantic)
endif()

# ----------------------------------------------------------------------------
# Doxygen documentation
# ----------------------------------------------------------------------------
find_package(Doxygen)

if(DOXYGEN_FOUND)

  set(DOXYGEN_PROJECT_NAME           "${PROJECT_NAME}")
  set(DOXYGEN_PROJECT_BRIEF          "Current source converter dynamic simulation")
  set(DOXYGEN_PROJECT_NUMBER         "${PROJECT_VERSION}")

  # Output folder for generated documentation.
  # It is better to keep generated files outside the source docs/ folder.
  set(DOXYGEN_OUTPUT_DIRECTORY       "${CMAKE_CURRENT_BINARY_DIR}/docs")

  set(DOXYGEN_CREATE_SUBDIRS         NO)
  set(DOXYGEN_GENERATE_HTML          YES)
  set(DOXYGEN_GENERATE_LATEX         NO)

  # Mathematical equations.
  set(DOXYGEN_USE_MATHJAX            YES)
  set(DOXYGEN_HTML_FORMULA_FORMAT    svg)

  # General extraction options.
  set(DOXYGEN_RECURSIVE              YES)
  set(DOXYGEN_EXTRACT_ALL            YES)
  set(DOXYGEN_EXTRACT_PRIVATE        NO)
  set(DOXYGEN_EXTRACT_STATIC         YES)

  # Path formatting.
  set(DOXYGEN_FULL_PATH_NAMES        YES)
  set(DOXYGEN_STRIP_FROM_PATH        "${CMAKE_CURRENT_SOURCE_DIR}")

  # Warnings.
  set(DOXYGEN_QUIET                  NO)
  set(DOXYGEN_WARN_IF_UNDOCUMENTED   NO)
  set(DOXYGEN_WARN_IF_DOC_ERROR      YES)

  doxygen_add_docs(
    doc
    "${CMAKE_CURRENT_SOURCE_DIR}/docs/mainpage.dox"
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    COMMENT "Generating API documentation with Doxygen"
  )

else()
  message(STATUS "Doxygen not found: target 'doc' will not be available")
endif()

# post build command, copy the executable to ${CMAKE_SOURCE_DIR}/exe
add_custom_command(TARGET ${TARGET} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory    "${CMAKE_SOURCE_DIR}/exe"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:${TARGET}>" "${CMAKE_SOURCE_DIR}/exe/$<TARGET_FILE_NAME:${TARGET}>"
  COMMENT "Copying ${TARGET} executable to ${CMAKE_SOURCE_DIR}/exe"
)

message(STATUS "Project: ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

