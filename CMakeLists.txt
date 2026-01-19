cmake_minimum_required(VERSION 3.0)

project(a4)
set(CMAKE_CXX_STANDARD 11)

find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)

add_library(col781 src/hw.cpp src/camera.cpp deps/src/gl.c)
target_include_directories(col781 PUBLIC deps/include)
target_link_libraries(col781 glm::glm OpenGL::GL SDL2::SDL2)

add_executable(part1 executables/part1.cpp)
target_include_directories(part1 PUBLIC src)
target_link_libraries(part1 col781)

add_executable(part2 executables/part2.cpp)
target_include_directories(part2 PUBLIC src)
target_link_libraries(part2 col781)


add_executable(part3 executables/part3.cpp)
target_include_directories(part3 PUBLIC src)
target_link_libraries(part3 col781)