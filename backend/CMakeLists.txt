cmake_minimum_required(VERSION 3.10)
project(ShockSimAPI CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 march=native -DNDEBUG")

find_package(Threads REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURE_DIR}/include/cpp-httplib
  ${CMAKE_SOURCE_DIR}/include/json/include
)

add_executable(shocktube
  src/main.cpp
  src/solver.cpp
)

target_include_directories(shocktube PRIVATE
  ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(shocktube PRIVATE
  Threads::Threads
)