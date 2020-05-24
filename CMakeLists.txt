cmake_minimum_required(VERSION 3.16)
project(untitled2)

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall -Wextra -pedantic -O3)

add_library(big BigInt/BigUInt.cpp BigInt/BigUInt.h timer.cpp timer.h BigInt/DigitVector.cpp BigInt/DigitVector.h BigInt/BigInt.cpp BigInt/BigInt.h)
link_libraries(big gmp gmpxx)

add_executable(untitled2 main.cpp)

add_subdirectory(Catch_tests)