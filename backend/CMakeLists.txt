cmake_minimum_required(VERSION 3.10)
project(ShockSimAPI CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 march=native -DNDEBUG")

if(APPLE)
  if(EXISTS "/opt/homebrew/opt/libomp")
    set(LIBOMP_PREFIX "/opt/homebrew/opt/libomp")
  elseif(EXISTS "user/local/opt/libomp")
    set(LIBOMP_PREFIX "/usr/local/opt/libomp")
  else()
    message(FATAL_ERROR "libomp not found. Installation required")
  endif()

  message(STATUS "Found libomp at ${LIBOMP_PREFIX}")

  set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
  set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_CXX_LIB_NAMES "omp")
  set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
endif()


find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURE_DIR}/include/cpp-httplib
  ${CMAKE_SOURCE_DIR}/include/json/include
)

add_executable(shocktube
  src/main.cpp
  src/solver.cpp
)

add_executable(solver_test
  src/solver_test.cpp
  src/solver.cpp
)

add_executable(toro_test
  src/toro_test.cpp
  src/solver.cpp
)

add_executable(mc_test
  src/montecarlo_test.cpp
)

target_include_directories(toro_test PRIVATE ${CMAKE_SOURCE_DIR}/include)

target_include_directories(solver_test PRIVATE ${CMAKE_SOURCE_DIR}/include)

target_include_directories(shocktube PRIVATE
  ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(shocktube PRIVATE
  Threads::Threads
  OpenMP::OpenMP_CXX
)

target_compile_options(shocktube PRIVATE -O3)

target_link_libraries(mc_test PRIVATE
  OpenMP::OpenMP_CXX
)