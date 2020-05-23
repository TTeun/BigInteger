#include "BigInt/BigInt.h"
#include "BigInt/BigUInt.h"
#include "timer.h"

#include <ctime>
#include <gmpxx.h>
#include <iostream>

int main() {

    //    std::vector<size_t> vec = {0, DigitVector::s_base - 1ul,DigitVector::s_base - 1ul};
    //    for (auto it : vec){
    //        std::cout << it << " " ;
    //    }
    //    std::cout << '\n';
    //
    //    BigUIntBase::multiplyBySmallNumberViaIterators(vec.rbegin(), vec.rend(), 6);
    //
    //    for (auto it : vec){
    //        std::cout << it << " " ;
    //    }
    //    std::cout << '\n';
    //
    //    return 0;

    mpz_t n1;
    mpz_init(n1);
    mpz_set_str(n1, "21983712987319237912389139871982398127398712312312312312", 10);

    Timer t1;
    for (size_t i = 0; i != 5; ++i) {
        mpz_mul(n1, n1, n1);
    }
    const double gmpTime = t1.elapsed();
    BigUInt      b       = BigUInt{"21983712987319237912389139871982398127398712312312312312"};
    Timer        t;
    for (size_t i = 0; i != 5; ++i) {
        b *= b;
    }
    const double bTime = t.elapsed();
    std::cout << bTime / gmpTime << '\n';
}
