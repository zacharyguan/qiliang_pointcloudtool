cmake_minimum_required(VERSION 3.22)
#message(STATUS "当前源代码所在的路径：${CMAKE_CURRENT_SOURCE_DIR}")
#set(CMAKE_TOOLCHAIN_FILE "E:/Zachary/Git/vcpkg/scripts/buildsystems/vcpkg.cmake")
#message(STATUS "当前CMAKE_TOOLCHAIN_FILE路径：${CMAKE_TOOLCHAIN_FILE}")
project(PlyFastLoader)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

file(GLOB_RECURSE src_files src/*.cpp)
file(GLOB_RECURSE include_files include/*.h)
add_executable(PlyFastLoader ${src_files} ${include_files})

find_package(jsoncpp CONFIG REQUIRED)
target_link_libraries(PlyFastLoader PRIVATE JsonCpp::JsonCpp)

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(PlyFastLoader PRIVATE Eigen3::Eigen)