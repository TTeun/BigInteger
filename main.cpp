#include "BigInt/BigUInt.h"
#include "Time/BenchMark.h"

#include <iostream>

int main() {
//    big::BigUInt b("111111111222222222333333333");
//    std::cout << b << '\n';
//    std::cout << big::BigUInt::longDivision(b, big::BigUInt("3333333333")) << '\n';
//    std::cout << b << '\n';
//
//
//    std::cout << b.squareRootRemainder().first << '\n';

        BenchMark::run({{10, 100}, {100, 100}, {1000, 100}, {10000, 50}, {100000, 10}, {1000000, 5}});
}
