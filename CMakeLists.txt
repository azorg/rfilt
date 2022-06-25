cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(PROJECT_NAME "rfilt_test")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Os")

add_definitions(-DRFILT_DEBUG)

set(SOURCES
  rfilt.c
  rfilt_test.c
)

include_directories(
 ./
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} m pthread)

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE C)

