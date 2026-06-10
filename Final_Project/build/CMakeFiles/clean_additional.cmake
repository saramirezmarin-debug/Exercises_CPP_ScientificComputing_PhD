# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "")
  file(REMOVE_RECURSE
  "C:\\CPP_Projects\\CPP_Course\\Exercises\\Final_Project\\docs\\html"
  )
endif()
