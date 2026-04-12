cmake_minimum_required(VERSION 3.14)
project(BlackScholesSolver)

# GoogleTest requires at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Create a library with all source files
add_library(solver_lib src/linear_algebra.cpp src/black_scholes.cpp)
target_include_directories(solver_lib PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Create executable for main.cpp
add_executable(solver src/main.cpp)
target_link_libraries(solver solver_lib)
# optional: ensure include dir if other targets don't transitively see it
target_include_directories(solver PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Create a single unified test executable with all test files
add_executable(
  all_tests
  tests/decomp_tests.cpp
  tests/forward_sub_tests.cpp
  tests/back_sub_tests.cpp
)
target_link_libraries(
  all_tests
  solver_lib
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(all_tests)


# Fetch pybind11
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.11.1
)
FetchContent_MakeAvailable(pybind11)

# Tell pybind11 to build a Python extension module named 'black_scholes_solver'
# We include the binding file AND the core library files it depends on
pybind11_add_module(black_scholes_solver 
    src/bindings.cpp 
    src/linear_algebra.cpp 
    src/black_scholes.cpp
)

# Install the compiled C++ library 
install(TARGETS black_scholes_solver DESTINATION .)

# Install the Python type hints stub file alongside it
install(FILES src/black_scholes_solver.pyi DESTINATION .)

# Ensure the module can find your headers
target_include_directories(black_scholes_solver PRIVATE ${CMAKE_SOURCE_DIR}/include)