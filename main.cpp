//#define NDEBUG

#include "BigInt/BigInt.h"
#include "BigInt/BigUInt.h"
#include "timer.h"

#include <ctime>
#include <gmpxx.h>
#include <iostream>

int main() {
    mpz_t n1;
    mpz_init(n1);
    mpz_set_str(n1, "21983712987319237912389139871982398127398712312312312312", 10);

    Timer t1;
    for (size_t i = 0; i != 7; ++i) {
        mpz_mul(n1, n1, n1);
    }
    const double gmpTime = t1.elapsed();
    BigUInt      b       = BigUInt{"21983712987319237912389139871982398127398712312312312312"};
    Timer        t;
    for (size_t i = 0; i != 7; ++i) {
        b *= b;
    }
    const double bTime = t.elapsed();
    std::cout << bTime / gmpTime << '\n';
}
