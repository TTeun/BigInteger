#ifndef __BIG_U_INT__H__
#define __BIG_U_INT__H__

#include "BigUIntBase.h"
#include "DigitVector.h"

#include <ostream>

class BigUInt : public BigUIntBase {

public:
    BigUInt() {
        init(0);
        assert(isWellFormed());
    }

    BigUInt(size_t val) {
        init(val);
        assert(isWellFormed());
    }

    BigUInt(BigUInt &&other) noexcept : BigUIntBase(std::move(other.m_digits)) {
    }

    BigUInt(const BigUInt &other)
        : BigUIntBase(std::vector<size_t>(other.m_digits.begin(), other.m_digits.end())) {
        assert(isWellFormed());
    }

    explicit BigUInt(std::vector<size_t> &&digits, bool isAlreadyCorrectlySized);

    explicit BigUInt(const std::string &val);

    static BigUInt createRandom(size_t numberOfDigits);

    static BigUInt createRandomFromDecimalDigits(size_t orderOfMagnitude);

    BigUInt &operator=(const BigUInt &rhs);

    BigUInt &operator=(size_t rhs) {
        init(rhs);
        return *this;
    }

    BigUInt &operator+=(const BigUInt &rhs);

    BigUInt &operator-=(const BigUInt &rhs);

    BigUInt operator-(const BigUInt &rhs) const;

    BigUInt &operator+=(size_t rhs);

    BigUInt &operator*=(const BigUInt &rhs);

    BigUInt &operator*=(size_t rhs);

    BigUInt operator+(size_t rhs) const;

    BigUInt operator+(const BigUInt &rhs) const;

    BigUInt operator*(size_t rhs) const;

    BigUInt operator*(const BigUInt &rhs) const;

    size_t operator%(size_t mod) const;

    BigUInt operator%(const BigUInt &mod) const;

    BigUInt &operator%=(size_t mod);

    BigUInt &operator%=(const BigUInt &mod);

    BigUInt operator/(const BigUInt &divisor) const;

    BigUInt operator/(size_t divisor) const;

    BigUInt &operator/=(const BigUInt &divisor);

    BigUInt &operator/=(size_t divisor);

    bool operator==(const BigUInt &rhs) const;

    bool operator!=(const BigUInt &rhs) const {
        return !(rhs == *this);
    }

    bool operator<(const BigUInt &rhs) const;

    bool operator>(const BigUInt &rhs) const {
        return rhs < *this;
    }

    bool operator<=(const BigUInt &rhs) const;

    void divideByLessThanBase(size_t factor);

    bool operator>=(const BigUInt &rhs) const {
        return !(*this < rhs);
    }

    friend std::ostream &operator<<(std::ostream &os, const BigUInt &anInt);

    friend BigUInt operator+(size_t lhs, const BigUInt &rhs) {
        return rhs + lhs;
    }

    friend BigUInt operator-(size_t lhs, const BigUInt &rhs) {
        return BigUInt(lhs) - rhs;
    }

    friend BigUInt operator*(size_t lhs, const BigUInt &rhs) {
        return rhs * lhs;
    }

    friend BigUInt power(const BigUInt &base, size_t exponent);

private:
    static size_t divisionSubRoutine(leftToRightConstIterator leftToRightConstIt,
                                     leftToRightConstIterator leftToRightConstEnd,
                                     rightToLeftIterator      rightToLeftIt,
                                     rightToLeftIterator      rightToLeftEnd,
                                     const BigUInt &          divisor);

    static BigUInt longDivision(BigUInt &dividend, const BigUInt &divisor);

    static BigUInt longDivisionAfterAdjustingDivisor(BigUInt &dividend, const BigUInt &divisor);

    void square();

    void init(size_t val);

    void bubble(size_t startIndex = 0ul);

public:
    static BigUInt toomCook_3(const BigUInt &m, const BigUInt &n);
};

#endif // __BIG_U_INT__H__