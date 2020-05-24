#define NDEBUG

#include "BigInt/BigUInt.h"
#include "Time/BenchMark.h"

#include <iostream>

int main() {
        BenchMark::run({{10, 100}, {100, 100}, {1000, 100}, {10000, 50}, {100000, 10}, {1000000, 5}});
//    big::BigUInt b("1111111111111111111111111111");
//    std::cout << b << '\n';
//    b *= big::DigitVector::s_base;
//
//    std::cout << b << '\n';
}
