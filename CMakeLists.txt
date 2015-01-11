cmake_minimum_required(VERSION 2.6.4)

include_directories (.)

set (PROJECT cylinder)

project (${PROJECT})

include_directories (./math)
include_directories (./gauss)
add_subdirectory(gauss)

if(CMAKE_COMPILER_IS_GNUCC)
  add_definitions (-Wall -pedantic -std=gnu++11 -O2)
endif()


set (HEADERS 
  ./given.h
  ./enter_data.h

  ./math/hankel.h
  ./math/chebyshev.h
  ./math/simple_math.h
  ./math/point.h
)

set (SOURCES
  ./main.cpp
  ./enter_data.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

add_executable(${PROJECT} ${HEADERS} ${SOURCES})

if(CMAKE_COMPILER_IS_GNUCC)
  target_link_libraries(${PROJECT} pthread)
endif()
