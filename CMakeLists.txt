cmake_minimum_required(VERSION 3.10)
project(circular_project)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR}/../sdsl/include)
link_directories(${CMAKE_SOURCE_DIR}/../sdsl/lib)

add_executable(main "delta_core_rlz_test.cpp")
target_link_libraries(main sdsl divsufsort divsufsort64)