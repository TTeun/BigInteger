#ifndef TEUN_GAME_BIGUNSIGNEDINT_H
#define TEUN_GAME_BIGUNSIGNEDINT_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <vector>

class BigUnsignedInt {

public:
    static const size_t s_base = std::sqrt(std::numeric_limits<size_t>::max()) - 1ul;

public:
    friend BigUnsignedInt subRoutine(const BigUnsignedInt & A, const BigUnsignedInt & B);

    BigUnsignedInt();

    explicit BigUnsignedInt(std::vector<size_t> && digits);

    BigUnsignedInt(size_t val);

    BigUnsignedInt(const std::string & val);

    BigUnsignedInt & operator=(const BigUnsignedInt & rhs);

    BigUnsignedInt & operator=(size_t rhs);

    BigUnsignedInt & operator+=(const BigUnsignedInt & rhs);

    BigUnsignedInt & operator-=(const BigUnsignedInt & rhs);

    BigUnsignedInt operator-(const BigUnsignedInt & rhs) const;

    BigUnsignedInt & operator+=(size_t rhs);

    BigUnsignedInt & operator*=(const BigUnsignedInt & rhs);

    BigUnsignedInt & operator*=(size_t rhs);

    BigUnsignedInt operator+(size_t rhs) const;

    BigUnsignedInt operator+(const BigUnsignedInt & rhs) const;

    BigUnsignedInt operator*(size_t rhs) const;

    BigUnsignedInt operator*(const BigUnsignedInt & rhs) const;

    size_t operator%(size_t mod) const;

    BigUnsignedInt & operator%=(size_t mod);

    BigUnsignedInt operator/(const BigUnsignedInt & divisor) const;

    bool operator==(const BigUnsignedInt & rhs) const;

    bool operator!=(const BigUnsignedInt & rhs) const;

    bool operator<(const BigUnsignedInt & rhs) const;

    bool operator>(const BigUnsignedInt & rhs) const;

    bool operator<=(const BigUnsignedInt & rhs) const;

    bool operator>=(const BigUnsignedInt & rhs) const;

    size_t mostSignificantDigit() const;

    friend std::ostream & operator<<(std::ostream & os, const BigUnsignedInt & anInt);

    friend BigUnsignedInt power(const BigUnsignedInt & base, size_t exponent);

    size_t digitCount() const;
protected:
    void shiftAdd(const BigUnsignedInt & rhs, size_t shiftAmount);

    void square();

    bool isCorrectlySized() const;

    void resizeToFit();

    void init(size_t val);

    void bubble(size_t startIndex = 0ul);

    BigUnsignedInt & shift(size_t shiftAmount);

    std::vector<size_t> m_digits;
};

#endif // TEUN_GAME_BIGUNSIGNEDINT_H
