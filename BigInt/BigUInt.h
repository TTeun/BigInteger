#ifndef __BIG_U_INT__H__
#define __BIG_U_INT__H__

#include "DigitVector.h"

#include <ostream>

class BigUInt : public DigitVector {

public:
    static const size_t s_karatsubaLowerLimit = 200ul;
    static const size_t s_toomCookLowerLimit  = 250ul;

    BigUInt() {
        init(0);
        assert(isWellFormed());
    }

    BigUInt(size_t val) {
        init(val);
        assert(isWellFormed());
    }

    BigUInt(BigUInt &&other) noexcept : DigitVector(std::move(other.m_digits)) {
    }

    BigUInt(const BigUInt &other)
        : DigitVector(std::vector<size_t>(other.m_digits.begin(), other.m_digits.end())) {
        assert(isWellFormed());
    }

    explicit BigUInt(std::vector<size_t> &&digits, bool isAlreadyCorrectlySized);

    explicit BigUInt(const std::string &val);

    BigUInt(rightToLeftConstIterator it, rightToLeftConstIterator endIt) : DigitVector({it, endIt}) {
    }

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
    static BigUInt multiply(BigUInt smaller, BigUInt larger);

    static void karatsubaMultiplyViaIterators(rightToLeftIterator            resultIt,
                                              const rightToLeftIterator      resultEnd,
                                              rightToLeftConstIterator       rhsIt,
                                              const rightToLeftConstIterator rhsEnd,
                                              rightToLeftConstIterator       copyIt,
                                              const rightToLeftConstIterator copyEnd);

    static void splitOneMultiplicationViaIterators(rightToLeftIterator            resultIt,
                                                   const rightToLeftIterator      resultEnd,
                                                   rightToLeftConstIterator       rhsIt,
                                                   const rightToLeftConstIterator rhsEnd,
                                                   rightToLeftConstIterator       largeIt,
                                                   const rightToLeftConstIterator largeEnd);

    static void multiplyViaIterators(rightToLeftIterator            resultIt,
                                     const rightToLeftIterator      resultEnd,
                                     rightToLeftConstIterator       rhsIt,
                                     const rightToLeftConstIterator rhsEnd,
                                     rightToLeftConstIterator       copyIt,
                                     const rightToLeftConstIterator copyEnd);

    static bool lessThanShiftedRhsViaIterators(leftToRightConstIterator       thisIt,
                                               const leftToRightConstIterator thisEnd,
                                               leftToRightConstIterator       rhsIt,
                                               const leftToRightConstIterator rhsEnd,
                                               size_t                         trailingZeroesOfRhs);

    static void
    carryAdditionViaIterators(rightToLeftIterator thisIt, const rightToLeftIterator thisEnd, size_t carry);

    static void addViaIterators(rightToLeftIterator                       thisIt,
                                const rightToLeftIterator                 thisEnd,
                                std::vector<size_t>::const_iterator       rhsIt,
                                const std::vector<size_t>::const_iterator rhsEnd);
    static void addMultipleViaIterators(rightToLeftIterator                       thisIt,
                                        const rightToLeftIterator                 thisEnd,
                                        std::vector<size_t>::const_iterator       rhsIt,
                                        const std::vector<size_t>::const_iterator rhsEnd,
                                        const size_t                              multiplier);

    static void subtractViaIterators(rightToLeftIterator                       thisIt,
                                     const rightToLeftIterator                 thisEnd,
                                     std::vector<size_t>::const_iterator       rhsIt,
                                     const std::vector<size_t>::const_iterator rhsEnd);

    static bool lessThanViaIterators(const leftToRightConstIterator thisIt,
                                     const leftToRightConstIterator thisEnd,
                                     const leftToRightConstIterator rhsIt,
                                     const leftToRightConstIterator rhsEnd) {
        return lessThanShiftedRhsViaIterators(thisIt, thisEnd, rhsIt, rhsEnd, 0ul);
    }
    static bool greaterThanViaIterators(const leftToRightConstIterator thisIt,
                                        const leftToRightConstIterator thisEnd,
                                        const leftToRightConstIterator rhsIt,
                                        const leftToRightConstIterator rhsEnd) {
        return lessThanViaIterators(rhsIt, rhsEnd, thisIt, thisEnd);
    }

    void divideByLessThanBase(size_t factor);

    bool operator>=(const BigUInt &rhs) const {
        return !(*this < rhs);
    }

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
    static BigUInt toomCook_3(rightToLeftIterator      resultIt,
                              rightToLeftIterator      resultEnd,
                              rightToLeftConstIterator rhsIt,
                              rightToLeftConstIterator rhsEnd,
                              rightToLeftConstIterator copyIt,
                              rightToLeftConstIterator copyEnd);
};

#endif // __BIG_U_INT__H__