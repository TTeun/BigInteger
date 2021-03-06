#include "BigUInt.h"

#include "BigInt.h"

#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

namespace big {
    const size_t BigUIntBase::s_maxDigit; // So that we can use it in std::min in BigUInt::divisionSubRoutine

    /***************** Some static functions *****************/
    static size_t ceilingIntegerDivision(size_t a, size_t b) { return (a + b - 1) / b; }

    static size_t modPower(size_t base, size_t exponent, size_t mod) {
        if (exponent == 0ul) {
            return 1;
        } else if (exponent == 1ul) {
            return base % mod;
        }
        size_t aux(1);
        size_t copy(base);
        while (exponent > 1ul) {
            if (exponent % 2ul != 0ul) {
                aux *= copy;
                aux %= mod;
            }
            exponent /= 2ul;
            copy *= copy;
            copy %= mod;
        }
        copy *= aux;
        return copy % mod;
    }

    void bubbleViaIterators(rlIterator thisIt, const rlIterator &thisEnd) {
        auto next = thisIt + 1ul;
        for (; next != thisEnd; ++thisIt, ++next) {
            if (*thisIt > BigUIntBase::s_maxDigit) {
                *next += *thisIt / BigUIntBase::s_base;
                *thisIt %= BigUIntBase::s_base;
            }
        }
        assert(*next <= BigUIntBase::s_maxDigit);
    }

    /***************** Constructors *****************/
    BigUInt::BigUInt(std::vector<size_t> &&digits, bool isAlreadyCorrectlySized) : BigUIntBase(std::move(digits)) {
        if (m_digits.empty()) {
            init(0);
        }
        if (isAlreadyCorrectlySized) {
            assert(isWellFormed());
        } else {
            resizeToFit();
        }
    }

    BigUInt::BigUInt(const std::string &val) {
        static const auto   maximumDecimalDigitsInSizeType = static_cast<size_t>(std::log10(std::numeric_limits<size_t>::max()));
        static const size_t largestMultipleOfTenInSizeType =
            modPower(10ul, maximumDecimalDigitsInSizeType, std::numeric_limits<size_t>::max());

        m_digits          = {0ul};
        size_t startIndex = val.length() % maximumDecimalDigitsInSizeType;
        size_t subVal;
        if (startIndex != 0ul) {
            *this *= largestMultipleOfTenInSizeType;
            std::istringstream iss(val.substr(0, startIndex));
            iss >> subVal;
            *this += subVal;
        }
        while (startIndex + maximumDecimalDigitsInSizeType <= val.length()) {
            *this *= largestMultipleOfTenInSizeType;
            std::istringstream iss(val.substr(startIndex, maximumDecimalDigitsInSizeType));
            iss >> subVal;
            *this += subVal;
            startIndex += maximumDecimalDigitsInSizeType;
        }
        assert(isWellFormed());
    }

    /***************** Operators *****************/
    BigUInt &BigUInt::operator=(const BigUInt &rhs) {
        m_digits = rhs.m_digits;
        return *this;
    }

    /** Addition **/
    BigUInt &BigUInt::operator+=(const BigUInt &rhs) {
        assert(isWellFormed() && rhs.isWellFormed());

        resize(std::max(digitCount(), rhs.digitCount()) + 1ul);
        addViaIterators(rlBegin(), rlEnd(), rhs.rlcBegin(), rhs.rlcEnd());
        if (mostSignificantDigit() == 0ul) {
            resize(digitCount() - 1ul);
        }

        return *this;
    }

    BigUInt &BigUInt::operator+=(size_t rhs) {
        if (rhs == 0ul) {
            return *this;
        }
        static const size_t additionRoom = std::numeric_limits<size_t>::max() - s_base;
        if (rhs <= additionRoom) {
            resize(digitCount() + 1ul);
            carryAdditionViaIterators(rlBegin(), rlEnd(), rhs);
            if (mostSignificantDigit() == 0ul) {
                resize(digitCount() - 1ul);
            }
        } else {
            *this += BigUInt(rhs);
        }
        return *this;
    }

    BigUInt BigUInt::operator+(size_t rhs) const {
        auto copy = *this;
        return copy += rhs;
    }

    BigUInt BigUInt::operator+(const BigUInt &rhs) const {
        auto copy = *this;
        return copy += rhs;
    }

    /** Subtraction **/
    BigUInt &BigUInt::operator-=(const BigUInt &rhs) {
        assert(rhs <= *this);

        subtractViaIterators(rlBegin(), rlEnd(), rhs.rlcBegin(), rhs.rlcEnd());
        resizeToFit();
        return *this;
    }

    BigUInt BigUInt::operator-(const BigUInt &rhs) const {
        auto copy = *this;
        return copy -= rhs;
    }

    /** Multiplication **/
    BigUInt &BigUInt::operator*=(const size_t rhs) {
        if (rhs == 0 || *this == BigUInt(0ul)) {
            init(0);
            return *this;
        }
        if (rhs != 1ul) {
            if (rhs < s_base) {
                resize(digitCount() + 1ul);
                multiplyBySingleDigitViaIterators(rlBegin(), rlEnd(), rhs);
                resizeToFit();
            } else if (rhs < s_base * s_base) {
                resize(digitCount() + 2ul);
                multiplyByDoubleDigitsViaIterators(rlBegin(), rlEnd(), rhs % s_base, rhs / s_base);
                resizeToFit();
            } else {
                reserve(digitCount() + 3ul);
                *this *= BigUInt({rhs % s_base, (rhs / s_base) % s_base, rhs / (s_base * s_base)}, true);
            }
        }
        return *this;
    }

    BigUInt &BigUInt::operator*=(const BigUInt &rhs) {
        assert(isWellFormed());
        assert(rhs.isWellFormed());
        if (this == &rhs) {
            square();
        }

        switch (rhs.digitCount()) {
            case 1ul:
                resize(digitCount() + 1ul);
                multiplyBySingleDigitViaIterators(rlBegin(), rlEnd(), rhs.mostSignificantDigit());
                if (mostSignificantDigit() == 0ul) {
                    resize(digitCount() - 1ul);
                }
                return *this;
            case 2ul:
                resize(digitCount() + 2ul);
                multiplyByDoubleDigitsViaIterators(rlBegin(), rlEnd(), rhs.leastSignificantDigit(), rhs.mostSignificantDigit());
                resizeToFit();
                return *this;
            default: *this = rhs.digitCount() > digitCount() ? multiply(*this, rhs) : multiply(rhs, *this); return *this;
        }
    }

    BigUInt BigUInt::operator*(size_t rhs) const {
        auto copy = *this;
        return copy *= rhs;
    }

    BigUInt BigUInt::operator*(const BigUInt &rhs) const {
        if (*this < rhs) {
            return multiply(*this, rhs);
        } else {
            return multiply(rhs, *this);
        }
    }

    std::pair<BigUInt, BigUInt> BigUInt::squareRootRemainder() const {
        return recursiveSquareRoot(BigUInt(shiftedCopy((-digitCount()) % 4), true));
    }

    /** Division **/
    BigUInt &BigUInt::operator/=(size_t divisor) {
        if (divisor < s_base) {
            divideByLessThanBase(divisor);
        } else {
            *this = *this / divisor;
        }
        return *this;
    }

    BigUInt &BigUInt::operator/=(const BigUInt &divisor) {
        *this = *this / divisor;
        return *this;
    }

    BigUInt BigUInt::operator/(const BigUInt &divisor) const {
        assert(divisor != 0ul);
        if (&divisor == this) {
            return BigUInt(1);
        }
        if (divisor == 1ul) {
            return BigUInt(*this);
        }
        if (divisor.digitCount() == 1ul && divisor.mostSignificantDigit() < s_base) {
            auto copy = *this;
            copy.divideByLessThanBase(divisor.mostSignificantDigit());
            return copy;
        }
        if (*this < divisor) {
            return 0;
        }
        const size_t m = this->digitCount();
        const size_t n = divisor.digitCount();

        if (m < n) {
            return 0ul;
        }
        if (m == 1) {
            return mostSignificantDigit() / divisor.mostSignificantDigit();
        }
        BigUInt copy(*this);
        return longDivision(copy, divisor);
    }

    BigUInt BigUInt::operator/(size_t divisor) const {
        if (divisor < s_base) {
            auto copy = *this;
            copy.divideByLessThanBase(divisor);
            return copy;
        } else {
            return *this / BigUInt(divisor);
        }
    }

    /** Modulo **/
    BigUInt &BigUInt::operator%=(size_t mod) {
        *this = *this % mod;
        return *this;
    }

    BigUInt &BigUInt::operator%=(const BigUInt &mod) {
        assert(mod != 0ul);

        if (*this < mod) {
            return *this;
        }
        const size_t factor = ceilingIntegerDivision(BigUInt::s_base, 2ul * mod.mostSignificantDigit());
        if (factor != 1ul) {
            *this *= factor;
            longDivisionAfterAdjustingDivisor(*this, factor * mod);
            *this /= factor;
        } else {
            longDivisionAfterAdjustingDivisor(*this, mod);
        }
        resizeToFit();
        return *this;
    }

    size_t BigUInt::operator%(size_t mod) const {
        if (mod == 2ul) {
            return leastSignificantDigit() % 2ul;
        }
        if (mod >= std::numeric_limits<size_t>::max() / mod) {
            auto copy = *this % BigUInt(mod);
            return copy.leastSignificantDigit() + s_base * copy.digitAt(1) + s_base * s_base * copy.mostSignificantDigit();
        }
        size_t result = 0ul;
        for (size_t i = 0; i != digitCount(); ++i) {
            result += (digitAt(i) % mod) * modPower(s_base, i, mod);
            result %= mod;
        }
        return result;
    }

    BigUInt BigUInt::operator%(const BigUInt &mod) const {
        auto copy = *this;
        if (mod.digitCount() == 1ul) {
            copy %= mod.leastSignificantDigit();
        } else if (mod.digitCount() == 2ul) {
            copy %= mod.mostSignificantDigit() * s_base + mod.leastSignificantDigit();
        } else {
            copy %= mod;
        }
        return copy;
    }

    /** Comparison **/
    bool BigUInt::operator==(const BigUInt &rhs) const {
        assert(isWellFormed());
        assert(rhs.isWellFormed());

        return m_digits == rhs.m_digits;
    }

    bool BigUInt::operator<(const BigUInt &rhs) const {
        assert(isWellFormed() && rhs.isWellFormed());

        return lessThanViaIterators(lrcBegin(), lrcEnd(), rhs.lrcBegin(), rhs.lrcEnd());
    }

    bool BigUInt::operator<=(const BigUInt &rhs) const {
        assert(isWellFormed() && rhs.isCorrectlySized());

        if (digitCount() != rhs.digitCount()) {
            return digitCount() < rhs.digitCount();
        }

        auto thisIt = lrcBegin();
        auto rhsIt  = rhs.lrcBegin();
        for (; thisIt != lrcEnd(); ++thisIt, ++rhsIt) {
            if (*thisIt != *rhsIt) {
                return *thisIt <= *rhsIt;
            }
        }
        return true;
    }

    /** Friends **/
    BigUInt power(const BigUInt &base, size_t exponent) {
        if (exponent == 0ul) {
            return BigUInt(1);
        } else if (exponent == 1ul) {
            return BigUInt(base);
        }
        BigUInt aux(1);
        BigUInt copy(base);

        while (exponent > 1ul) {
            if (exponent % 2ul == 1ul) {
                aux *= copy;
            }
            exponent /= 2ul;
            copy.square();
        }
        copy *= aux;

        return copy;
    }

    /***************** Builders *****************/
    BigUInt BigUInt::createRandom(size_t numberOfDigits) {
        assert(numberOfDigits > 0ul);
        std::random_device                    rd;
        std::mt19937                          gen(rd());
        std::uniform_int_distribution<size_t> dis(0, s_maxDigit);

        BigUInt result;
        result.resize(numberOfDigits);
        for (size_t i = 0; i != numberOfDigits; ++i) {
            result.m_digits[i] = dis(gen);
        }
        result.resizeToFit();
        return result;
    }

    BigUInt BigUInt::createRandomFromDecimalDigits(size_t orderOfMagnitude) {
        assert(orderOfMagnitude > 0);
        const auto numberOfDigits = orderOfMagnitude * (std::log(10) / log(s_base)) + 1ul;
        auto       result         = createRandom((size_t)numberOfDigits);
        assert(result.isWellFormed());
        return result;
    }

    BigUInt BigUInt::createWithRoom(size_t digitCount) {
        BigUInt result;
        result.reserve(digitCount);
        return result;
    }

    /***************** Output *****************/
    std::string BigUInt::toString() const {
        std::stringstream ss;
        ss << *lrcBegin();
        for (auto it = lrcBegin() + 1ul; it != lrcEnd(); ++it) {
            ss << std::setfill('0') << std::setw(s_digitsPerLimb) << *it;
        }
        return ss.str();
    }

    std::ostream &operator<<(std::ostream &os, const BigUInt &bigUnsignedInt) {
        os << "( ";
        for (auto it = bigUnsignedInt.lrcBegin(); it != bigUnsignedInt.lrcEnd(); ++it) {
            os << *it << " ";
        }
        os << ")_" << BigUInt::s_base;
        return os;
    }

    /***************** Internal *****************/
    void BigUInt::init(size_t val) {
        static_assert(s_base <= std::numeric_limits<size_t>::max() / s_base, "s_base^2 should not exceed the maximum size_t");
        static_assert(s_base % 2 == 0, "s_base should be even");
        static_assert(s_base * s_base > std::numeric_limits<size_t>::max() / s_base, "s_base^3 should exceed the maximum size_t");
        if (val < s_base) {
            m_digits = {val};
        } else if (val < s_base * s_base) {
            m_digits = {val % s_base, val / s_base};
        } else {
            m_digits = {val % s_base, (val / s_base) % s_base, val / (s_base * s_base)};
        }
    }

    void BigUInt::bubble(size_t startIndex) {
        assert(not m_digits.empty());
        assert(startIndex < digitCount());

        auto it   = rlBegin() + startIndex;
        auto next = it + 1;

        for (; next != rlEnd(); ++it, ++next) {
            if (*it > s_maxDigit) {
                *next += *it / s_base;
                *it %= s_base;
            }
        }

        if (mostSignificantDigit() > s_maxDigit) {
            const size_t continueIndex = digitCount() - 1ul;
            resize(1ul + static_cast<size_t>(digitCount() + std::log(mostSignificantDigit()) / std::log(s_base)));
            bubble(continueIndex);
        }
        resizeToFit();
    }

    void BigUInt::divideByLessThanBase(size_t factor) {
        assert(factor < s_base);

        size_t carry  = 0ul;
        auto   thisIt = lrBegin();
        for (; thisIt != lrEnd(); ++thisIt) {
            *thisIt += carry * s_base;
            carry = *thisIt % factor;
            *thisIt /= factor;
        }
        resizeToFit();
    }

    void BigUInt::square() {
        const auto copy = *this;
        *this *= copy;
        bubble(0);
    }

    void BigUInt::append(const std::vector<size_t> &highDigits) {
        m_digits.resize(digitCount() + highDigits.size());
        std::copy(highDigits.cbegin(), highDigits.cend(), m_digits.begin() + digitCount() - highDigits.size());
        assert(isWellFormed());
    }

    std::tuple<BigUInt, BigUInt, BigUInt, BigUInt> BigUInt::splitFour(rlcIterator begin, rlcIterator end, size_t i) {
        assert(3ul * i < static_cast<size_t>(end - begin));
        return {BigUInt({begin, begin + i}, false),
                BigUInt({begin + i, begin + 2ul * i}, false),
                BigUInt({begin + 2ul * i, begin + 3ul * i}, false),
                BigUInt({begin + 3ul * i, end}, false)};
    }

    std::tuple<BigUInt, BigUInt, BigUInt> BigUInt::splitThree(rlcIterator begin, rlcIterator end, size_t i) {
        assert(2ul * i < static_cast<size_t>(end - begin));
        return {BigUInt({begin, begin + i}, false),
                BigUInt({begin + i, begin + 2ul * i}, false),
                BigUInt({begin + 2ul * i, end}, false)};
    }

    std::pair<BigUInt, BigUInt> BigUInt::recursiveSquareRoot(const BigUInt &val) {
        //        assert(val.digitCount() % 4ul == 0ul);
        //        const size_t  i = val.digitCount() / 4;
        //        const BigUInt a0({val.rlcBegin(), val.rlcBegin() + i}, false);
        //        BigUInt       a1({val.rlcBegin() + 1ul * i, val.rlcBegin() + 2ul * i}, false);
        //        const BigUInt a32({val.rlcBegin() + 2ul * i, val.rlcEnd()}, false);
        //
        //        std::pair<BigUInt, BigUInt> p = a32.squareRootRemainder();
        //        a1.append(p.second.m_digits);
        //        BigUInt quot = longDivision(a1, 2ul * p.first);
        //        BigUInt u =
        //        quot.append(p.first.m_digits);
    }

    /***************** Static helpers *****************/
    /** Addition **/
    void BigUInt::carryAdditionViaIterators(rlIterator thisIt, rlIterator thisEnd, size_t carry) {
        assert(carry != 0ul);
        assert(carry + s_base < std::numeric_limits<size_t>::max());
        for (; thisIt != thisEnd; ++thisIt) {
            *thisIt += carry;
            if (*thisIt < s_base) {
                return;
            }
            carry = (*thisIt) / s_base;
            *thisIt %= s_base;
        }
        assert(false);
    }

    void BigUInt::addViaIterators(rlIterator resultIt, rlIterator resultEnd, rlcIterator rhsIt, rlcIterator rhsEnd) {
        assert(std::distance(resultIt, resultEnd) >= std::distance(rhsIt, rhsEnd) + 1l);
        bool   carry = false;
        size_t s0;
        for (; rhsEnd - rhsIt >= 4;) {
            s0        = *resultIt + *rhsIt + carry;
            *resultIt = s0 % s_base;
            ++resultIt;
            ++rhsIt;
#define REPEAT_THREE(code) code code code
            REPEAT_THREE(s0 = *resultIt + *rhsIt + (s0 > s_maxDigit); *resultIt = s0 % s_base; ++resultIt; ++rhsIt;)
#undef REPEAT_THREE
            carry = s0 > s_maxDigit;
        }

        switch (rhsEnd - rhsIt) {
            case 0: break;
            case 1:
                s0        = *resultIt + carry + *rhsIt;
                *resultIt = s0 % s_base;
                ++resultIt;
                carry = s0 > s_maxDigit;
                break;
            case 2:
                s0        = *resultIt + carry + *rhsIt;
                *resultIt = s0 % s_base;
                ++resultIt;
                ++rhsIt;
                s0        = *resultIt + *rhsIt + (s0 > s_maxDigit);
                *resultIt = s0 % s_base;
                ++resultIt;
                carry = s0 > s_maxDigit;
                break;
            case 3:
                s0        = *resultIt + carry + *rhsIt;
                *resultIt = s0 % s_base;
                ++resultIt;
                ++rhsIt;
                s0        = *resultIt + *rhsIt + (s0 > s_maxDigit);
                *resultIt = s0 % s_base;
                ++resultIt;
                ++rhsIt;
                s0        = *resultIt + *rhsIt + (s0 > s_maxDigit);
                *resultIt = s0 % s_base;
                ++resultIt;
                carry = s0 > s_maxDigit;
                break;
            default: assert(false);
        }
        if (carry) {
            carryAdditionViaIterators(resultIt, resultEnd, 1ul);
        }
    }

    void BigUInt::addMultipleViaIterators(
        rlIterator resultIt, rlIterator resultEnd, rlcIterator rhsIt, rlcIterator rhsEnd, size_t multiplier) {
        if (multiplier == 0ul) {
            return;
        }
        assert(multiplier < s_base);
        assert(std::distance(resultIt, resultEnd) >= std::distance(rhsIt, rhsEnd) + 1l);
        size_t carry = 0ul;
        size_t s0;
        for (; rhsEnd - rhsIt >= 4;) {
            s0        = *resultIt + *rhsIt * multiplier + carry;
            *resultIt = s0 % s_base;
            ++resultIt;
            ++rhsIt;
#define REPEAT_THREE(code) code code code
            REPEAT_THREE(s0 = *resultIt + *rhsIt * multiplier + (s0 / s_base); *resultIt = s0 % s_base; ++resultIt; ++rhsIt;)
#undef REPEAT_THREE
            carry = s0 / s_base;
        }
        switch (rhsEnd - rhsIt) {
            case 0: break;
            case 1:
                s0        = *resultIt + *rhsIt * multiplier + carry;
                *resultIt = s0 % s_base;
                ++resultIt;
                ++rhsIt;
                carry = s0 / s_base;
                break;
            case 2:
                s0        = *resultIt + *rhsIt * multiplier + carry;
                *resultIt = s0 % s_base;
                ++resultIt;
                ++rhsIt;
                s0        = *resultIt + *rhsIt * multiplier + (s0 / s_base);
                *resultIt = s0 % s_base;
                ++resultIt;
                ++rhsIt;
                carry = s0 / s_base;
                break;
            case 3:
                s0        = *resultIt + *rhsIt * multiplier + carry;
                *resultIt = s0 % s_base;
                ++resultIt;
                ++rhsIt;
                s0        = *resultIt + *rhsIt * multiplier + (s0 / s_base);
                *resultIt = s0 % s_base;
                ++resultIt;
                ++rhsIt;
                s0        = *resultIt + *rhsIt * multiplier + (s0 / s_base);
                *resultIt = s0 % s_base;
                ++resultIt;
                ++rhsIt;
                carry = s0 / s_base;
                break;
            default: assert(false);
        }
        if (carry > 0UL) {
            carryAdditionViaIterators(resultIt, resultEnd, carry);
        }
    }

    void BigUInt::subtractViaIterators(rlIterator thisIt, rlIterator thisEnd, rlcIterator rhsIt, rlcIterator rhsEnd) {
        assert(std::distance(thisIt, thisEnd) >= std::distance(rhsIt, rhsEnd));
        size_t carry = 0ul;
        for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
            if (*thisIt >= *rhsIt + carry) {
                *thisIt -= *rhsIt + carry;
                carry = 0;
            } else {
                *thisIt += s_base;
                *thisIt -= *rhsIt + carry;
                carry = 1ul;
            }
        }
        if (carry != 0ul) {
            while (true) {
                assert(thisIt != thisEnd);
                if (*thisIt > 0) {
                    --*thisIt;
                    break;
                }
                *thisIt = s_maxDigit;
                ++thisIt;
            }
        }
    }

    /** Multiplication **/
    BigUInt BigUInt::multiply(const BigUInt &smaller, const BigUInt &larger) {
        BigUInt result;
        result.resize(smaller.digitCount() + larger.digitCount() + 1);

        multiplyViaIterators(
            result.rlBegin(), result.rlEnd(), smaller.rlcBegin(), smaller.rlcEnd(), larger.rlcBegin(), larger.rlcEnd());
        result.resizeToFit();
        assert(result.isWellFormed());
        return result;
    }

    void BigUInt::multiplyBySingleDigitViaIterators(rlIterator resultIt, const rlIterator resultEnd, const size_t rhs) {
        size_t carry = 0ul;
        for (; resultIt != resultEnd; ++resultIt) {
            *resultIt *= rhs;
            *resultIt += carry;
            if (*resultIt > s_maxDigit) {
                carry = *resultIt / s_base;
                *resultIt %= s_base;
            } else {
                carry = 0ul;
            }
        }
    }

    void BigUInt::multiplyByDoubleDigitsViaIterators(rlIterator       resultIt,
                                                     const rlIterator resultEnd,
                                                     const size_t     least,
                                                     const size_t     most) {
        size_t previousVal = 0ul;
        size_t currentVal;
        size_t carry = 0ul;
        for (; resultIt != resultEnd; ++resultIt) {
            currentVal = *resultIt;

            *resultIt = carry + currentVal * least + previousVal * most;
            if (*resultIt > s_maxDigit) {
                carry = *resultIt / s_base;
                *resultIt %= s_base;
            } else {
                carry = 0ul;
            }
            previousVal = currentVal;
        }
    }

    void BigUInt::karatsubaMultiplyViaIterators(rlIterator  resultIt,
                                                rlIterator  resultEnd,
                                                rlcIterator smallIt,
                                                rlcIterator smallEnd,
                                                rlcIterator largeIt,
                                                rlcIterator largeEnd) {
        assert((largeEnd - largeIt) >= (smallEnd - smallIt));
        const auto m = static_cast<size_t>(smallEnd - smallIt);
        const auto n = static_cast<size_t>(largeEnd - largeIt);
        assert(m >= s_karatsubaLowerLimit);
        const size_t splitIndex = m / 2ul;

        BigUInt high1(largeIt + splitIndex, largeEnd);
        BigUInt high2(smallIt + splitIndex, smallEnd);

        BigUInt z0;
        z0.resize(2ul * splitIndex + 1ul);
        multiplyViaIterators(z0.rlBegin(), z0.rlEnd(), smallIt, smallIt + splitIndex, largeIt, largeIt + splitIndex);
        z0.resizeToFit();

        BigUInt z2;
        z2.resize(m % 2ul + n + 1ul);
        multiplyViaIterators(z2.rlBegin(), z2.rlEnd(), smallIt + splitIndex, smallEnd, largeIt + splitIndex, largeEnd);
        z2.resizeToFit();

        high1.resize(high1.digitCount() + 1ul);
        high2.resize(high2.digitCount() + 1ul);

        addViaIterators(high1.rlBegin(), high1.rlEnd(), largeIt, largeIt + splitIndex);
        addViaIterators(high2.rlBegin(), high2.rlEnd(), smallIt, smallIt + splitIndex);

        BigUInt z1;
        z1.resize(m % 2ul + n + 3ul);

        high1.resizeToFit();
        high2.resizeToFit();

        multiplyViaIterators(z1.rlBegin(), z1.rlEnd(), high1.rlcBegin(), high1.rlcEnd(), high2.rlcBegin(), high2.rlcEnd());

        subtractViaIterators(z1.rlBegin(), z1.rlEnd(), z2.rlcBegin(), z2.rlcEnd());
        subtractViaIterators(z1.rlBegin(), z1.rlEnd(), z0.rlcBegin(), z0.rlcEnd());

        addViaIterators(resultIt, resultEnd, z0.rlcBegin(), z0.rlcEnd());
        addViaIterators(resultIt + 2ul * splitIndex, resultEnd, z2.rlcBegin(), z2.rlcEnd());
        addViaIterators(resultIt + splitIndex, resultEnd, z1.rlcBegin(), z1.rlcEnd());
    }

    void BigUInt::splitOneMultiplicationViaIterators(rlIterator  resultIt,
                                                     rlIterator  resultEnd,
                                                     rlcIterator smallIt,
                                                     rlcIterator smallEnd,
                                                     rlcIterator largeIt,
                                                     rlcIterator largeEnd) {
        const auto m = static_cast<size_t>(largeEnd - largeIt);
        assert(m >= s_karatsubaLowerLimit);
        const auto   n          = static_cast<size_t>(smallEnd - smallIt);
        const size_t splitIndex = m / 2ul;

        BigUInt z0;
        z0.resize(splitIndex + n + 1ul);
        multiplyViaIterators(z0.rlBegin(), z0.rlEnd(), smallIt, smallEnd, largeIt, largeIt + splitIndex);
        z0.resizeToFit();
        addViaIterators(resultIt, resultEnd, z0.rlcBegin(), z0.rlcEnd());

        BigUInt z1;
        z1.resize(m - splitIndex + n + 1ul);
        multiplyViaIterators(z1.rlBegin(), z1.rlEnd(), largeIt + splitIndex, largeEnd, smallIt, smallEnd);
        z1.resizeToFit();
        addViaIterators(resultIt + splitIndex, resultEnd, z1.rlcBegin(), z1.rlcEnd());
    }

    void BigUInt::multiplySortedViaIterators(rlIterator        resultIt,
                                             const rlIterator  resultEnd,
                                             rlcIterator       smallIt,
                                             const rlcIterator smallEnd,
                                             rlcIterator       largeIt,
                                             const rlcIterator largeEnd) {
        const auto smallSize = static_cast<size_t>(smallEnd - smallIt);
        const auto largeSize = static_cast<size_t>(largeEnd - largeIt);
        assert(largeSize >= smallSize);

        if (largeSize < s_karatsubaLowerLimit) {
            schoolMultiply(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        } else if (smallSize < s_karatsubaLowerLimit) {
            splitOneMultiplicationViaIterators(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        } else if (2ul * smallSize <= largeSize) {
            splitOneMultiplicationViaIterators(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        } else if (smallSize >= s_toomCook4LowerLimit) {
            toomCook_4(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        } else if (smallSize >= s_toomCook3LowerLimit) {
            toomCook_3(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        } else {
            karatsubaMultiplyViaIterators(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        }
    }

    void BigUInt::multiplyViaIterators(rlIterator  resultIt,
                                       rlIterator  resultEnd,
                                       rlcIterator smallIt,
                                       rlcIterator smallEnd,
                                       rlcIterator largeIt,
                                       rlcIterator largeEnd) {
        const auto copySize = static_cast<size_t>(largeEnd - largeIt);
        const auto rhsSize  = static_cast<size_t>(smallEnd - smallIt);
        if (copySize > rhsSize) {
            multiplySortedViaIterators(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        } else {
            multiplySortedViaIterators(resultIt, resultEnd, largeIt, largeEnd, smallIt, smallEnd);
        }
    }

    void BigUInt::toomCook_3(rlIterator  resultIt,
                             rlIterator  resultEnd,
                             rlcIterator smallIt,
                             rlcIterator smallEnd,
                             rlcIterator largeIt,
                             rlcIterator largeEnd) {
        assert((largeEnd - largeIt) >= (smallEnd - smallIt));
        const size_t i = (smallEnd - smallIt) / 3ul;

        const auto [m0, m1, m2] = static_cast<std::tuple<BigInt, BigInt, BigInt>>(splitThree(smallIt, smallEnd, i));
        const auto [n0, n1, n2] = static_cast<std::tuple<BigInt, BigInt, BigInt>>(splitThree(largeIt, largeEnd, i));

        const BigInt p_aux      = m0 + m2;
        const BigInt p_minusOne = p_aux - m1;
        const BigInt p_minusTwo = 2ul * (p_minusOne + m2) - m0;

        const BigInt q_aux      = n0 + n2;
        const BigInt q_minusOne = q_aux - n1;
        const BigInt q_minusTwo = 2ul * (q_minusOne + n2) - n0;

        const BigInt a0         = m0 * n0;
        const BigInt r_one      = (p_aux + m1) * (q_aux + n1);
        const BigInt r_minusOne = p_minusOne * q_minusOne;
        const BigInt r_minusTwo = q_minusTwo * p_minusTwo;
        const BigInt a4         = m2 * n2;

        BigInt a1;
        BigInt a2;
        BigInt a3;

        a3 = (r_minusTwo - r_one) / 3ul;
        a1 = (r_one - r_minusOne) / 2ul;
        a2 = r_minusOne - a0;
        a3 = (a2 - a3) / 2ul + 2ul * a4;
        a2 = a2 + a1 - a4;
        a1 -= a3;

        addViaIterators(resultIt, resultEnd, a0.magnitude().rlcBegin(), a0.magnitude().rlcEnd());
        addViaIterators(resultIt + i, resultEnd, a1.magnitude().rlcBegin(), a1.magnitude().rlcEnd());
        addViaIterators(resultIt + 2ul * i, resultEnd, a2.magnitude().rlcBegin(), a2.magnitude().rlcEnd());
        addViaIterators(resultIt + 3ul * i, resultEnd, a3.magnitude().rlcBegin(), a3.magnitude().rlcEnd());
        addViaIterators(resultIt + 4ul * i, resultEnd, a4.magnitude().rlcBegin(), a4.magnitude().rlcEnd());
    }

    void BigUInt::toomCook_4(rlIterator  resultIt,
                             rlIterator  resultEnd,
                             rlcIterator smallIt,
                             rlcIterator smallEnd,
                             rlcIterator largeIt,
                             rlcIterator largeEnd) {
        assert((largeEnd - largeIt) >= (smallEnd - smallIt));
        const size_t i = (smallEnd - smallIt) / 4ul;

        const auto [m0, m1, m2, m3] = static_cast<std::tuple<BigInt, BigInt, BigInt, BigInt>>(splitFour(smallIt, smallEnd, i));
        auto [n0, n1, n2, n3] =
            static_cast<std::tuple<BigInt, const BigInt, const BigInt, BigInt>>(splitFour(largeIt, largeEnd, i));

        BigInt n02 = n0 + n2;
        BigInt n13 = n1 + n3;
        BigInt m02 = m0 + m2;
        BigInt m13 = m1 + m3;

        const big::BigInt r_one      = (n02 + n13) * (m02 + m13);
        const big::BigInt r_minusOne = (n02 - n13) * (m02 - m13);
        n02 += 3ul * n2;
        n13 += n1 + 7ul * n3;
        m02 += 3ul * m2;
        m13 += m1 + 7ul * m3;
        const big::BigInt r_two      = (n02 + n13) * (m02 + m13);
        const big::BigInt r_minusTwo = (n02 - n13) * (m02 - m13);
        const big::BigInt r_three    = (n0 + 3ul * n1 + 9ul * n2 + 27ul * n3) * (m0 + 3ul * m1 + 9ul * m2 + 27ul * m3);

        n0 *= m0;
        n3 *= m3;

        BigInt a1 =
            (-20ll * n0 - 720ll * n3 + 60ll * r_one - 30ll * r_minusOne - 15ll * r_two + 3ll * r_minusTwo + 2ll * r_three) / 60ul;
        BigInt a2 = (-30ll * n0 + 96ll * n3 + 16ll * r_one + 16ll * r_minusOne - r_two - r_minusTwo) / 24ul;
        BigInt a3 = (10ll * n0 + 360ll * n3 - 14ll * r_one - r_minusOne + 7ll * r_two - r_minusTwo - r_three) / 24ul;
        BigInt a4 = (6ll * n0 - 120ll * n3 - 4ll * r_one - 4ll * r_minusOne + r_two + r_minusTwo) / 24ul;
        BigInt a5 = (-10ll * n0 - 360ll * n3 + 10ll * r_one + 5ll * r_minusOne - 5ll * r_two - r_minusTwo + r_three) / 120ul;

        addViaIterators(resultIt, resultEnd, n0.magnitude().rlcBegin(), n0.magnitude().rlcEnd());
        addViaIterators(resultIt + i, resultEnd, a1.magnitude().rlcBegin(), a1.magnitude().rlcEnd());
        addViaIterators(resultIt + 2ul * i, resultEnd, a2.magnitude().rlcBegin(), a2.magnitude().rlcEnd());
        addViaIterators(resultIt + 3ul * i, resultEnd, a3.magnitude().rlcBegin(), a3.magnitude().rlcEnd());
        addViaIterators(resultIt + 4ul * i, resultEnd, a4.magnitude().rlcBegin(), a4.magnitude().rlcEnd());
        addViaIterators(resultIt + 5ul * i, resultEnd, a5.magnitude().rlcBegin(), a5.magnitude().rlcEnd());
        addViaIterators(resultIt + 6ul * i, resultEnd, n3.magnitude().rlcBegin(), n3.magnitude().rlcEnd());
    }

    void BigUInt::schoolMultiply(rlIterator  resultIt,
                                 rlIterator  resultEnd,
                                 rlcIterator smallIt,
                                 rlcIterator smallEnd,
                                 rlcIterator largeIt,
                                 rlcIterator largeEnd) {
        assert((largeEnd - largeIt) >= (smallEnd - smallIt));
        for (size_t i = 0; smallIt < smallEnd; ++i) {
            addMultipleViaIterators(resultIt + i, resultEnd, largeIt, largeEnd, *smallIt);
            ++smallIt;
        }
    }

    /** Division **/
    size_t BigUInt::divisionSubRoutine(const lrcIterator leftToRightConstIt,
                                       const lrcIterator leftToRightConstEnd,
                                       const rlIterator  rightToLeftIt,
                                       const rlIterator  rightToLeftEnd,
                                       const BigUInt &   divisor) {
        if (lessThanViaIterators(leftToRightConstIt, leftToRightConstEnd, divisor.lrcBegin(), divisor.lrcEnd())) {
            return 0ul;
        }
        assert(divisor != 0ul);
        assert(divisor.mostSignificantDigit() * 2ul >= BigUInt::s_base);
        const size_t n = divisor.digitCount();
        const auto   m = static_cast<size_t>(leftToRightConstEnd - leftToRightConstIt);
        assert(m <= n + 1ul);
        assert(m >= n);

        size_t correction = 0ul;
        while (not lessThanShiftedRhsViaIterators(
            leftToRightConstIt, leftToRightConstEnd, divisor.lrcBegin(), divisor.lrcEnd(), 1)) {
            subtractViaIterators(rightToLeftIt + 1ul, rightToLeftEnd, divisor.rlcBegin(), divisor.rlcEnd());
            correction += BigUInt::s_base;
        }

        size_t quotientEstimate;
        if (m == n) {
            quotientEstimate = (*leftToRightConstIt) / divisor.mostSignificantDigit();
        } else {
            quotientEstimate =
                (*leftToRightConstIt * BigUInt::s_base + *(leftToRightConstIt + 1ul)) / divisor.mostSignificantDigit();
        }
        quotientEstimate    = std::min(quotientEstimate, s_maxDigit);
        const size_t offset = *leftToRightConstIt == 0 ? 1ul : 0ul;

        BigUInt closestMultipleEstimate = createWithRoom(divisor.digitCount() + 1ul);
        closestMultipleEstimate         = divisor;
        closestMultipleEstimate.resize(closestMultipleEstimate.digitCount() + 1ul);
        multiplyBySingleDigitViaIterators(closestMultipleEstimate.rlBegin(), closestMultipleEstimate.rlEnd(), quotientEstimate);
        closestMultipleEstimate.resizeToFit();
        while (greaterThanViaIterators(closestMultipleEstimate.lrcBegin(),
                                       closestMultipleEstimate.lrcEnd(),
                                       leftToRightConstIt + offset,
                                       leftToRightConstEnd)) {
            --quotientEstimate;
            closestMultipleEstimate -= divisor;
        }
        subtractViaIterators(rightToLeftIt, rightToLeftEnd, closestMultipleEstimate.rlcBegin(), closestMultipleEstimate.rlcEnd());
        return quotientEstimate + correction;
    }

    BigUInt BigUInt::longDivision(BigUInt &dividend, const BigUInt &divisor) {
        if (dividend < divisor) {
            return 0ul;
        }
        const size_t factor = ceilingIntegerDivision(BigUInt::s_base, 2ul * divisor.mostSignificantDigit());
        if (factor > 1) {
            dividend *= factor;
            return longDivisionAfterAdjustingDivisor(dividend, factor * divisor);
            dividend /= factor;
        } else {
            return longDivisionAfterAdjustingDivisor(dividend, divisor);
        }
    }

    BigUInt BigUInt::longDivisionAfterAdjustingDivisor(BigUInt &dividend, const BigUInt &divisor) {
        assert(divisor <= dividend);
        assert(divisor.mostSignificantDigit() * 2ul >= BigUInt::s_base);

        size_t       m = dividend.digitCount();
        const size_t n = divisor.digitCount();
        if (m <= n + 1) {
            return BigUInt::divisionSubRoutine(
                dividend.lrcBegin(), dividend.lrcEnd(), dividend.rlBegin(), dividend.rlEnd(), divisor);
        }

        std::vector<size_t> divisorDigits(m - n + 1ul, 0ul);
        while (m > n + 1) {
            const size_t splitIndex   = m - n - 1ul;
            divisorDigits[splitIndex] = BigUInt::divisionSubRoutine(
                dividend.lrcBegin(), dividend.lrcBegin() + n + 1, dividend.rlBegin() + splitIndex, dividend.rlEnd(), divisor);
            dividend.resizeToFit();
            m = dividend.digitCount();
        }
        if (m > n || dividend.mostSignificantDigit() >= divisor.mostSignificantDigit()) {
            divisorDigits[0ul] = BigUInt::divisionSubRoutine(
                dividend.lrcBegin(), dividend.lrcEnd(), dividend.rlBegin(), dividend.rlEnd(), divisor);
        }
        bubbleViaIterators(divisorDigits.begin(), divisorDigits.end());
        BigUInt result(std::move(divisorDigits), false);
        return result;
    }

    /** Comparison **/
    bool BigUInt::lessThanShiftedRhsViaIterators(
        lrcIterator thisIt, lrcIterator thisEnd, lrcIterator rhsIt, lrcIterator rhsEnd, size_t trailingZeroesOfRhs) {
        if (static_cast<size_t>(thisEnd - thisIt) != static_cast<size_t>(rhsEnd - rhsIt) + trailingZeroesOfRhs) {
            return static_cast<size_t>(thisEnd - thisIt) < static_cast<size_t>(rhsEnd - rhsIt) + trailingZeroesOfRhs;
        }

        for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
            if (*thisIt != *rhsIt) {
                return *thisIt < *rhsIt;
            }
        }
        return false;
    }

} // namespace big