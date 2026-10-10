cmake_minimum_required(VERSION 3.15)
project(nbody LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(HDF5 REQUIRED COMPONENTS CXX)

set(SOURCES
    src/main.cpp
    src/barneshut.cpp
    src/integrator.cpp
    src/morton_code.cpp
    src/io.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE 
    ${HDF5_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(${PROJECT_NAME} PRIVATE ${HDF5_CXX_LIBRARIES})