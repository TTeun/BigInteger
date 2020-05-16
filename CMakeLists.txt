cmake_minimum_required(VERSION 3.16)
project(untitled2)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wall -Wextra -pedantic -Werror -O3)

add_executable(untitled2 main.cpp BigInt/BigUInt.cpp BigInt/BigUInt.h timer.cpp timer.h BigInt/DigitVector.cpp BigInt/DigitVector.h BigInt/BigUIntBase.cpp BigInt/BigUIntBase.h BigInt/BigInt.cpp BigInt/BigInt.h)