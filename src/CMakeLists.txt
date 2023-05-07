cmake_minimum_required(VERSION 3.24)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Matplot++
        PATHS "D:/Matplotplusplus1_1/lib/cmake/Matplot++"
        NO_DEFAULT_PATH
        REQUIRED)
add_subdirectory(tool)
find_path(DBG_MACRO_INCLUDE_DIRS "dbg.h" PATHS D:\\vcpkg\\installed\\x64-windows\\include)
set(spdlog_DIR "D:\\vcpkg\\installed\\x64-windows\\share\\spdlog")
set(fmt_DIR "D:\\vcpkg\\installed\\x64-windows\\share\\fmt")
find_package(spdlog CONFIG REQUIRED)

include_directories(tool)
add_executable(mainApp main.cpp)

target_link_libraries(mainApp PRIVATE Matplot++::matplot Tool)
target_include_directories(mainApp PRIVATE ${DBG_MACRO_INCLUDE_DIRS})
target_link_libraries(mainApp PRIVATE spdlog::spdlog spdlog::spdlog_header_only)
target_link_libraries(mainApp PRIVATE fmt::fmt fmt::fmt-header-only)