#ifndef __BIG_U_INT__H__
#define __BIG_U_INT__H__

#include "BigUIntBase.h"

#include <ostream>

namespace big {

    class BigUInt : public BigUIntBase {

    public:
        static const size_t s_karatsubaLowerLimit = 160ul;
        static const size_t s_toomCookLowerLimit = 800ul;

    public:
        /***************** Constructors *****************/
        BigUInt() {
            init(0);
        }
        BigUInt(size_t val) {
            init(val);
        }
        BigUInt(BigUInt &&other) noexcept : BigUIntBase(std::move(other.m_digits)) {
        }
        BigUInt(const BigUInt &other) : BigUIntBase(std::vector<size_t>(other.m_digits.begin(), other.m_digits.end())) {
            assert(isWellFormed());
        }
        BigUInt(std::vector<size_t> &&digits, bool isAlreadyCorrectlySized);
        explicit BigUInt(const std::string &val);
        BigUInt(rlcIterator it, rlcIterator endIt) : BigUIntBase({it, endIt}) {
            assert(isWellFormed());
        }

        /***************** Operators *****************/
        BigUInt &operator=(const BigUInt &rhs);
        BigUInt &operator=(size_t rhs) {
            init(rhs);
            return *this;
        }
        /** Addition **/
        BigUInt &operator+=(size_t rhs);
        BigUInt &operator+=(const BigUInt &rhs);
        BigUInt operator+(size_t rhs) const;
        BigUInt operator+(const BigUInt &rhs) const;

        /** Subtraction **/
        BigUInt &operator-=(const BigUInt &rhs);
        BigUInt operator-(const BigUInt &rhs) const;

        /** Multiplication **/
        BigUInt &operator*=(size_t rhs);
        BigUInt &operator*=(const BigUInt &rhs);
        BigUInt operator*(size_t rhs) const;
        BigUInt operator*(const BigUInt &rhs) const;

        /** Division **/
        BigUInt &operator/=(size_t divisor);
        BigUInt &operator/=(const BigUInt &divisor);
        BigUInt operator/(const BigUInt &divisor) const;
        BigUInt operator/(size_t divisor) const;

        /** Modulo **/
        BigUInt &operator%=(size_t mod);
        BigUInt &operator%=(const BigUInt &mod);
        size_t operator%(size_t mod) const;
        BigUInt operator%(const BigUInt &mod) const;

        /** Comparison **/
        bool operator==(const BigUInt &rhs) const;
        bool operator!=(const BigUInt &rhs) const {
            return !(rhs == *this);
        }
        bool operator<(const BigUInt &rhs) const;
        bool operator>(const BigUInt &rhs) const {
            return rhs < *this;
        }
        bool operator<=(const BigUInt &rhs) const;
        bool operator>=(const BigUInt &rhs) const {
            return !(*this < rhs);
        }

        /** Friends **/
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

        /***************** Builders *****************/
        static BigUInt createRandom(size_t numberOfDigits);
        static BigUInt createRandomFromDecimalDigits(size_t orderOfMagnitude);
        static BigUInt createWithRoom(size_t digitCount);

        /***************** Output *****************/
        std::string toString() const;
        friend std::ostream &operator<<(std::ostream &os, const BigUInt &anInt);

    private:
        /***************** Internal *****************/
        void init(size_t val);
        void bubble(size_t startIndex = 0ul);
        void square();
        void divideByLessThanBase(size_t factor);

        /***************** Static helpers *****************/
        /** Addition **/
        static void carryAdditionViaIterators(rlIterator thisIt, const rlIterator thisEnd, size_t carry);
        static void addViaIterators(rlIterator thisIt, const rlIterator thisEnd, rlcIterator rhsIt, const rlcIterator rhsEnd);
        static void addMultipleViaIterators(rlIterator resultIt,
                                            const rlIterator resultEnd,
                                            rlcIterator rhsIt,
                                            const rlcIterator rhsEnd,
                                            const size_t multiplier);
        static void
        subtractViaIterators(rlIterator thisIt, const rlIterator thisEnd, rlcIterator rhsIt, const rlcIterator rhsEnd);

        /** Multiplication **/
        static BigUInt multiply(const BigUInt &smaller, const BigUInt &larger);
        static void multiplyBySingleDigitViaIterators(rlIterator resultIt, const rlIterator resultEnd, const size_t rhs);
        static void multiplyByDoubleDigitsViaIterators(rlIterator resultIt,
                                                       const rlIterator resultEnd,
                                                       const size_t least,
                                                       const size_t most);

        static void karatsubaMultiplyViaIterators(rlIterator resultIt,
                                                  const rlIterator resultEnd,
                                                  rlcIterator smallIt,
                                                  const rlcIterator smallEnd,
                                                  rlcIterator largeIt,
                                                  const rlcIterator largeEnd);
        static void splitOneMultiplicationViaIterators(rlIterator resultIt,
                                                       const rlIterator resultEnd,
                                                       rlcIterator smallIt,
                                                       const rlcIterator smallEnd,
                                                       rlcIterator largeIt,
                                                       const rlcIterator largeEnd);
        static void multiplyViaIterators(rlIterator resultIt,
                                         const rlIterator resultEnd,
                                         rlcIterator rhsIt,
                                         const rlcIterator rhsEnd,
                                         rlcIterator copyIt,
                                         const rlcIterator copyEnd);
        static void multiplySortedViaIterators(rlIterator resultIt,
                                               const rlIterator resultEnd,
                                               rlcIterator smallIt,
                                               const rlcIterator smallEnd,
                                               rlcIterator largeIt,
                                               const rlcIterator largeEnd);
        static void toomCook_3(rlIterator resultIt,
                               rlIterator resultEnd,
                               rlcIterator rhsIt,
                               rlcIterator rhsEnd,
                               rlcIterator copyIt,
                               rlcIterator copyEnd);
        static void schoolMultiply(rlIterator resultIt,
                                   rlIterator resultEnd,
                                   rlcIterator smallIt,
                                   rlcIterator smallEnd,
                                   rlcIterator largeIt,
                                   rlcIterator largeEnd);

        /** Division **/
        static size_t divisionSubRoutine(lrcIterator leftToRightConstIt,
                                         lrcIterator leftToRightConstEnd,
                                         rlIterator rightToLeftIt,
                                         rlIterator rightToLeftEnd,
                                         const BigUInt &divisor);
        static BigUInt longDivision(BigUInt &dividend, const BigUInt &divisor);
        static BigUInt longDivisionAfterAdjustingDivisor(BigUInt &dividend, const BigUInt &divisor);

        /** Comparison **/
        static bool lessThanShiftedRhsViaIterators(lrcIterator thisIt,
                                                   const lrcIterator thisEnd,
                                                   lrcIterator rhsIt,
                                                   const lrcIterator rhsEnd,
                                                   size_t trailingZeroesOfRhs);
        static bool lessThanViaIterators(const lrcIterator thisIt,
                                         const lrcIterator thisEnd,
                                         const lrcIterator rhsIt,
                                         const lrcIterator rhsEnd) {
            return lessThanShiftedRhsViaIterators(thisIt, thisEnd, rhsIt, rhsEnd, 0ul);
        }
        static bool greaterThanViaIterators(const lrcIterator thisIt,
                                            const lrcIterator thisEnd,
                                            const lrcIterator rhsIt,
                                            const lrcIterator rhsEnd) {
            return lessThanViaIterators(rhsIt, rhsEnd, thisIt, thisEnd);
        }
    };
} // namespace big

#endif // __BIG_U_INT__H__