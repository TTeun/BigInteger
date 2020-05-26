#include "BigInt/BigUInt.h"
#include "Time/BenchMark.h"

#include <iostream>

int main() {
//    big::BigUInt b{{1, 1, 1, 1}, true};
//    big::BigUInt c{{1, 1, 1, 1}, true};
//    big::BigUInt d{{0, 0,0,0,0,0, 0, 0}, true};
//
//    big::BigUInt::toomCook_4(d.rlBegin(), d.rlEnd(), b.rlcBegin(), b.rlcEnd(), c.rlcBegin(), c.rlcEnd());
//    d.resizeToFit();
//    std::cout << d << '\n';
        BenchMark::run({{10, 100}, {100, 100}, {1000, 100}, {10000, 50}, {100000, 10}, {1000000, 5}});
//                BenchMark::run({{10, 100}, {100, 100}, {1000, 100}, {10000, 50}});
}
