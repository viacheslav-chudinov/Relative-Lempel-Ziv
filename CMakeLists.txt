cmake_minimum_required(VERSION 3.10)
project(circular_project)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

include_directories(${CMAKE_SOURCE_DIR}/../sdsl/include)
link_directories(${CMAKE_SOURCE_DIR}/../sdsl/lib)

add_executable(main "delta_core_rlz_test.cpp")
target_link_libraries(main sdsl divsufsort divsufsort64)


add_custom_target(profile
    COMMAND perf record -g -- $<TARGET_FILE:main>
    COMMAND perf report
    DEPENDS main
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running perf record + report on main"
)