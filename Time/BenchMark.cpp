#include "BenchMark.h"

#include "../BigInt/BigUInt.h"
#include "timer.h"

#include <gmpxx.h>
#include <iomanip>

double BenchMark::multiply(size_t powerOfTen, size_t numberOfRepetitions) {
    double gmpTime = 0.0;
    double bigTime = 0.0;
    for (size_t i = 0; i != numberOfRepetitions; ++i) {
        const BigUInt a = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        BigUInt       b = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        mpz_t         n1;
        mpz_init(n1);
        mpz_set_str(n1, a.toString().c_str(), 10);
        mpz_t n2;
        mpz_init(n2);
        mpz_set_str(n2, b.toString().c_str(), 10);
        {
            Timer t;
            mpz_mul(n1, n1, n2);
            gmpTime += t.elapsed();
        }
        {
            Timer t;
            b *= a;
            bigTime += t.elapsed();
        }
    }
    return bigTime / gmpTime;
}

double BenchMark::add(size_t powerOfTen, size_t numberOfRepetitions) {
    double gmpTime = 0.0;
    double bigTime = 0.0;
    for (size_t i = 0; i != numberOfRepetitions; ++i) {
        const BigUInt a = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        BigUInt       b = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        mpz_t         n1;
        mpz_init(n1);
        mpz_set_str(n1, a.toString().c_str(), 10);
        mpz_t n2;
        mpz_init(n2);
        mpz_set_str(n2, b.toString().c_str(), 10);
        {
            Timer t;
            mpz_add(n1, n1, n2);
            gmpTime += t.elapsed();
        }
        {
            Timer t;
            b += a;
            bigTime += t.elapsed();
        }
    }
    return bigTime / gmpTime;
}

double BenchMark::modulo(size_t powerOfTen, size_t numberOfRepetitions) {
    double gmpTime = 0.0;
    double bigTime = 0.0;
    for (size_t i = 0; i != numberOfRepetitions; ++i) {
        const BigUInt a = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        BigUInt       b = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        mpz_t         n1;
        mpz_init(n1);
        mpz_set_str(n1, a.toString().c_str(), 10);
        mpz_t n2;
        mpz_init(n2);
        mpz_set_str(n2, b.toString().c_str(), 10);
        {
            Timer t;
            mpz_mod(n1, n1, n2);
            gmpTime += t.elapsed();
        }
        {
            Timer t;
            b %= a;
            bigTime += t.elapsed();
        }
    }
    return bigTime / gmpTime;
}

double BenchMark::divide(size_t powerOfTen, size_t numberOfRepetitions) {
    double gmpTime = 0.0;
    double bigTime = 0.0;
    for (size_t i = 0; i != numberOfRepetitions; ++i) {
        const BigUInt a = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        BigUInt       b = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        mpz_t         n1;
        mpz_init(n1);
        mpz_set_str(n1, a.toString().c_str(), 10);
        mpz_t n2;
        mpz_init(n2);
        mpz_set_str(n2, b.toString().c_str(), 10);
        {
            Timer t;
            mpz_div(n1, n1, n2);
            gmpTime += t.elapsed();
        }
        {
            Timer t;
            b /= a;
            bigTime += t.elapsed();
        }
    }
    return bigTime / gmpTime;
}

void BenchMark::run() {
    std::cout << "Decimal digits:\t";
    std::cout << std::setfill(' ') << std::setw(15) << 10;
    std::cout << std::setfill(' ') << std::setw(15) << 100;
    std::cout << std::setfill(' ') << std::setw(15) << 1000;
    std::cout << std::setfill(' ') << std::setw(15) << 10000;
    std::cout << std::setfill(' ') << std::setw(15) << 100000 << '\n';

    std::cout << "Multiply:\t\t";
    std::cout << std::setfill(' ') << std::setw(15) << multiply(10, 250);
    std::cout << std::setfill(' ') << std::setw(15) << multiply(100, 250);
    std::cout << std::setfill(' ') << std::setw(15) << multiply(1000, 250);
    std::cout << std::setfill(' ') << std::setw(15) << multiply(10000, 250);
    std::cout << std::setfill(' ') << std::setw(15) << multiply(100000, 250) << '\n';

    std::cout << "Add:\t\t\t";
    std::cout << std::setfill(' ') << std::setw(15) << add(10, 250);
    std::cout << std::setfill(' ') << std::setw(15) << add(100, 250);
    std::cout << std::setfill(' ') << std::setw(15) << add(1000, 250);
    std::cout << std::setfill(' ') << std::setw(15) << add(10000, 250);
    std::cout << std::setfill(' ') << std::setw(15) << add(100000, 250) << "\n";

    std::cout << "Divide:\t\t\t";
    std::cout << std::setfill(' ') << std::setw(15) << divide(10, 250);
    std::cout << std::setfill(' ') << std::setw(15) << divide(100, 250);
    std::cout << std::setfill(' ') << std::setw(15) << divide(1000, 250);
    std::cout << std::setfill(' ') << std::setw(15) << divide(10000, 250);
    std::cout << std::setfill(' ') << std::setw(15) << divide(100000, 250) << "\n";

    std::cout << "Modulo:\t\t\t";
    std::cout << std::setfill(' ') << std::setw(15) << modulo(10, 250);
    std::cout << std::setfill(' ') << std::setw(15) << modulo(100, 250);
    std::cout << std::setfill(' ') << std::setw(15) << modulo(1000, 250);
    std::cout << std::setfill(' ') << std::setw(15) << modulo(10000, 250);
    std::cout << std::setfill(' ') << std::setw(15) << modulo(100000, 250) << "\n";
}
