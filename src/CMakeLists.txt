cmake_minimum_required(VERSION 3.25)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Matplot++
        PATHS "D:/Matplotplusplus1_1/lib/cmake/Matplot++"
        NO_DEFAULT_PATH
        REQUIRED)
add_subdirectory(tool)

include_directories(tool)
add_executable(mainApp main.cpp)

target_link_libraries(mainApp PRIVATE Matplot++::matplot Tool)