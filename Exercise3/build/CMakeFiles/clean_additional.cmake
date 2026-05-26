# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "")
  file(REMOVE_RECURSE
  "C:\\CPP_Projects\\CPP_Course\\Exercises\\Exercise3\\docs\\html"
  )
endif()
