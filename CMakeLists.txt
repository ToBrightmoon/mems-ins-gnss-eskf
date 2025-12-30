cmake_minimum_required(VERSION 3.10.0)

project(mems-ins-gnss-eskf)

set(CMAKE_CXX_STANDARD 17)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/third_party)

file(GLOB_RECURSE SRC_SOURCE main.cpp ins/ins_mech.cpp)

add_executable(${PROJECT_NAME} ${SRC_SOURCE})