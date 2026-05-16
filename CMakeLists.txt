cmake_minimum_required(VERSION 3.14)

project(Zipf VERSION 1.0.0 LANGUAGES CXX)  

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON) 
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(Zipf src/Zipf.cpp src/main.cpp include/Zipf.hh)  
target_include_directories(Zipf PRIVATE include)

if (UNIX)
  target_compile_options(Zipf PRIVATE -Wall -Wextra -Wpedantic )
endif()



add_custom_command(TARGET Zipf POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory    "${CMAKE_SOURCE_DIR}/exe"  
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:Zipf>" "${CMAKE_SOURCE_DIR}/exe/$<TARGET_FILE_NAME:Zipf>"
  COMMENT "Copying Zipf executable to ${CMAKE_SOURCE_DIR}/exe"
)


message(STATUS "Project: ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")