cmake_minimum_required(VERSION 3.16)
project(untitled2)

set(CMAKE_CXX_STANDARD 14)

add_executable(untitled2 main.cpp BigInt/BigUnsignedInt.cpp BigInt/BigUnsignedInt.h timer.cpp timer.h BigInt/DigitVector.cpp BigInt/DigitVector.h)