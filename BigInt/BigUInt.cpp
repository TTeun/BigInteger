#include "BigUInt.h"

#include "BigInt.h"

#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

namespace big {
    const size_t DigitVector::s_maxDigit; // So that we can use it in std::min in BigUInt::divisionSubRoutine

    static size_t ceilingIntegerDivision(size_t a, size_t b) {
        return (a + b - 1) / b;
    }

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
        auto next = thisIt + 1;

        for (; next != thisEnd; ++thisIt, ++next) {
            if (*thisIt >= DigitVector::s_base) {
                *next += *thisIt / DigitVector::s_base;
                *thisIt %= DigitVector::s_base;
            }
        }
    }

    size_t BigUInt::divisionSubRoutine(const lrcIterator leftToRightConstIt,
                                       const lrcIterator leftToRightConstEnd,
                                       const rlIterator rightToLeftIt,
                                       const rlIterator rightToLeftEnd,
                                       const BigUInt &divisor) {
        if (lessThanViaIterators(
                leftToRightConstIt, leftToRightConstEnd, divisor.leftToRightConstBegin(), divisor.leftToRightConstEnd())) {
            return 0ul;
        }
        assert(divisor != 0ul);
        assert(divisor.mostSignificantDigit() * 2ul >= BigUInt::s_base);
        const size_t n = divisor.digitCount();
        const auto m = static_cast<size_t>(leftToRightConstEnd - leftToRightConstIt);
        assert(m <= n + 1ul);
        assert(m >= n);
        size_t correction = 0ul;

        while (not lessThanShiftedRhsViaIterators(
            leftToRightConstIt, leftToRightConstEnd, divisor.leftToRightConstBegin(), divisor.leftToRightConstEnd(), 1)) {
            subtractViaIterators(
                rightToLeftIt + 1ul, rightToLeftEnd, divisor.rightToLeftConstBegin(), divisor.rightToLeftConstEnd());
            correction += BigUInt::s_base;
        }

        size_t quotientEstimate;
        if (m == n) {
            quotientEstimate = (*leftToRightConstIt) / divisor.mostSignificantDigit();
        } else {
            quotientEstimate =
                (*leftToRightConstIt * BigUInt::s_base + *(leftToRightConstIt + 1ul)) / divisor.mostSignificantDigit();
        }
        quotientEstimate = std::min(quotientEstimate, DigitVector::s_maxDigit);
        const size_t offset = *leftToRightConstIt == 0 ? 1ul : 0ul;

        auto closestMultipleEstimate = divisor * quotientEstimate;
        while (greaterThanViaIterators(closestMultipleEstimate.leftToRightConstBegin(),
                                       closestMultipleEstimate.leftToRightConstEnd(),
                                       leftToRightConstIt + offset,
                                       leftToRightConstEnd)) {
            --quotientEstimate;
            closestMultipleEstimate -= divisor;
        }
        subtractViaIterators(rightToLeftIt,
                             rightToLeftEnd,
                             closestMultipleEstimate.rightToLeftConstBegin(),
                             closestMultipleEstimate.rightToLeftConstEnd());
        return quotientEstimate + correction;
    }

    BigUInt BigUInt::longDivision(BigUInt &dividend, const BigUInt &divisor) {
        if (dividend < divisor) { return 0ul; }
        const size_t factor = ceilingIntegerDivision(BigUInt::s_base, 2ul * divisor.mostSignificantDigit());
        if (factor > 1) {
            dividend *= factor;
            return longDivisionAfterAdjustingDivisor(dividend, factor * divisor);
        } else {
            return longDivisionAfterAdjustingDivisor(dividend, divisor);
        }
    }

    BigUInt BigUInt::longDivisionAfterAdjustingDivisor(BigUInt &dividend, const BigUInt &divisor) {
        assert(divisor <= dividend);
        assert(divisor.mostSignificantDigit() * 2ul >= BigUInt::s_base);

        size_t m = dividend.digitCount();
        const size_t n = divisor.digitCount();

        if (m <= n + 1) {
            return BigUInt::divisionSubRoutine(dividend.leftToRightConstBegin(),
                                               dividend.leftToRightConstEnd(),
                                               dividend.rightToLeftBegin(),
                                               dividend.rightToLeftEnd(),
                                               divisor);
        }

        std::vector<size_t> divisorDigits(m - n + 1ul, 0ul);
        while (m > n + 1) {
            const size_t splitIndex = m - n - 1ul;
            divisorDigits[splitIndex] = BigUInt::divisionSubRoutine(dividend.leftToRightConstBegin(),
                                                                    dividend.leftToRightConstBegin() + n + 1,
                                                                    dividend.rightToLeftBegin() + splitIndex,
                                                                    dividend.rightToLeftEnd(),
                                                                    divisor);
            // The first n + 1 digits of dividend now contain their remainder after division by divisor, as per
            // the divisionSubRoutine
            dividend.resizeToFit();
            m = dividend.digitCount();
        }
        if (m > n || dividend.mostSignificantDigit() >= divisor.mostSignificantDigit()) {
            divisorDigits[0ul] = BigUInt::divisionSubRoutine(dividend.leftToRightConstBegin(),
                                                             dividend.leftToRightConstEnd(),
                                                             dividend.rightToLeftBegin(),
                                                             dividend.rightToLeftEnd(),
                                                             divisor);
        }
        bubbleViaIterators(divisorDigits.begin(), divisorDigits.end());
        BigUInt result(std::move(divisorDigits), false);
        return result;
    }

    std::ostream &operator<<(std::ostream &os, const BigUInt &bigUnsignedInt) {
        os << "( ";
        for (auto it = bigUnsignedInt.leftToRightConstBegin(); it != bigUnsignedInt.leftToRightConstEnd(); ++it) {
            os << *it << " ";
        }
        os << ")_" << BigUInt::s_base;
        return os;
    }

    BigUInt power(const BigUInt &base, size_t exponent) {
        if (exponent == 0ul) {
            return BigUInt(1);
        } else if (exponent == 1ul) {
            return BigUInt(base);
        }
        BigUInt aux(1);
        BigUInt copy(base);

        while (exponent > 1ul) {
            assert(aux.isWellFormed());
            assert(copy.isWellFormed());
            if (exponent % 2ul == 1ul) { aux *= copy; }
            exponent /= 2ul;
            copy.square();
        }
        assert(aux.isWellFormed());
        assert(copy.isWellFormed());

        copy *= aux;
        assert(aux.isWellFormed());
        assert(copy.isWellFormed());

        return copy;
    }

    /* =================== Member functions =================== */

    BigUInt::BigUInt(std::vector<size_t> &&digits, bool isAlreadyCorrectlySized) : DigitVector(std::move(digits)) {
        if (m_digits.empty()) { init(0); }
        if (isAlreadyCorrectlySized) {
            assert(isWellFormed());
        } else {
            resizeToFit();
        }
    }

    BigUInt::BigUInt(const std::string &val) {
        static const auto maximumDecimalDigitsInSizeType = static_cast<size_t>(std::log10(std::numeric_limits<size_t>::max()));
        static const size_t largestMultipleOfTenInSizeType =
            modPower(10ul, maximumDecimalDigitsInSizeType, std::numeric_limits<size_t>::max());

        m_digits = {0ul};
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

    void BigUInt::bubble(size_t startIndex) {
        assert(not m_digits.empty());
        assert(startIndex < digitCount());

        auto it = rightToLeftBegin() + startIndex;
        auto next = it + 1;

        for (; next != rightToLeftEnd(); ++it, ++next) {
            if (*it >= s_base) {
                *next += *it / s_base;
                *it %= s_base;
            }
        }

        if (mostSignificantDigit() >= s_base) {
            const size_t continueIndex = digitCount() - 1ul;
            resize(1ul + static_cast<size_t>(digitCount() + std::log(mostSignificantDigit()) / std::log(s_base)));
            bubble(continueIndex);
        }
        resizeToFit();
    }

    void BigUInt::init(size_t val) {
        static_assert(s_base <= std::numeric_limits<size_t>::max() / s_base, "s_base^2 should not exceed the maximum size_t");
        static_assert(s_base % 2 == 0, "s_base should be even");
        m_digits = {val};
        bubble();
    }

    bool BigUInt::operator==(const BigUInt &rhs) const {
        assert(isWellFormed());
        assert(rhs.isWellFormed());

        return m_digits == rhs.m_digits;
    }

    BigUInt &BigUInt::operator=(const BigUInt &rhs) {
        if (this == &rhs) { return *this; }
        m_digits = rhs.m_digits;
        return *this;
    }

    BigUInt &BigUInt::operator+=(const BigUInt &rhs) {
        assert(isWellFormed() && rhs.isWellFormed());

        resize(std::max(digitCount(), rhs.digitCount()) + 1ul);
        addViaIterators(rightToLeftBegin(), rightToLeftEnd(), rhs.rightToLeftConstBegin(), rhs.rightToLeftConstEnd());
        if (mostSignificantDigit() == 0ul) { resize(digitCount() - 1ul); }

        return *this;
    }

    bool BigUInt::operator<(const BigUInt &rhs) const {
        assert(isWellFormed() && rhs.isWellFormed());

        return lessThanViaIterators(
            leftToRightConstBegin(), leftToRightConstEnd(), rhs.leftToRightConstBegin(), rhs.leftToRightConstEnd());
    }

    bool BigUInt::operator<=(const BigUInt &rhs) const {
        assert(isWellFormed() && rhs.isCorrectlySized());

        if (digitCount() != rhs.digitCount()) { return digitCount() < rhs.digitCount(); }

        auto thisIt = leftToRightConstBegin();
        auto rhsIt = rhs.leftToRightConstBegin();
        for (; thisIt != leftToRightConstEnd(); ++thisIt, ++rhsIt) {
            if (*thisIt != *rhsIt) { return *thisIt <= *rhsIt; }
        }
        return true;
    }

    BigUInt BigUInt::multiply(const BigUInt &smaller, const BigUInt &larger) {
        BigUInt result;
        result.resize(smaller.digitCount() + larger.digitCount() + 1);

        multiplyViaIterators(result.rightToLeftBegin(),
                             result.rightToLeftEnd(),
                             smaller.rightToLeftConstBegin(),
                             smaller.rightToLeftConstEnd(),
                             larger.rightToLeftConstBegin(),
                             larger.rightToLeftConstEnd());
        result.resizeToFit();
        assert(result.isWellFormed());
        return result;
    }

    BigUInt &BigUInt::operator*=(const BigUInt &rhs) {
        assert(isWellFormed());
        assert(rhs.isWellFormed());

        if (this == &rhs) {
            square();
            return *this;
        }

        *this = rhs.digitCount() > digitCount() ? multiply(*this, rhs) : multiply(rhs, *this);
        return *this;
    }

    BigUInt BigUInt::operator*(size_t rhs) const {
        auto copy = *this;
        return copy *= rhs;
    }

    void BigUInt::square() {
        const auto copy = *this;
        *this *= copy;
        bubble(0);
    }

    BigUInt &BigUInt::operator*=(const size_t rhs) {
        if (rhs == 0 || *this == BigUInt(0ul)) {
            init(0);
            return *this;
        }
        reserve(digitCount() + 3ul);
        if (rhs != 1ul) {
            if (rhs < s_base) {
                for (auto &it : m_digits) { it *= rhs; }
                bubble();
            } else {
                *this *= BigUInt(rhs);
            }
        }
        return *this;
    }

    BigUInt &BigUInt::operator+=(size_t rhs) {
        static const size_t additionRoom = std::numeric_limits<size_t>::max() - s_base;
        if (rhs <= additionRoom) {
            *rightToLeftBegin() += rhs;
            bubble();
        } else {
            *this += BigUInt(rhs);
        }
        return *this;
    }

    size_t BigUInt::operator%(size_t mod) const {
        if (mod == 2ul) { return leastSignificantDigit() % 2ul; }
        if (mod >= std::numeric_limits<size_t>::max() / mod) {
            auto copy = *this % BigUInt(mod);
            return copy.leastSignificantDigit() + s_base * *(copy.rightToLeftConstBegin() + 1) +
                   s_base * s_base * copy.mostSignificantDigit();
        }
        size_t result = 0ul;
        for (size_t i = 0; i != digitCount(); ++i) {
            result += (digitAt(i) % mod) * modPower(s_base, i, mod);
            result %= mod;
        }
        return result;
    }

    BigUInt &BigUInt::operator%=(size_t mod) {
        *this = *this % mod;
        return *this;
    }

    BigUInt BigUInt::operator/(const BigUInt &divisor) const {
        assert(divisor != 0ul);
        if (&divisor == this) { return BigUInt(1); }
        if (divisor == 1ul) { return BigUInt(*this); }
        if (divisor < s_base) {
            auto copy = *this;
            copy.divideByLessThanBase(divisor.mostSignificantDigit());
            return copy;
        }
        if (*this < divisor) { return 0; }
        const size_t m = this->digitCount();
        const size_t n = divisor.digitCount();

        if (m < n) { return 0ul; }
        if (m == 1) { return mostSignificantDigit() / divisor.mostSignificantDigit(); }
        BigUInt copy(*this);
        return longDivision(copy, divisor);
    }

    BigUInt BigUInt::operator*(const BigUInt &rhs) const {
        return multiply(*this, rhs);
    }

    BigUInt &BigUInt::operator-=(const BigUInt &rhs) {
        assert(rhs <= *this);

        subtractViaIterators(rightToLeftBegin(), rightToLeftEnd(), rhs.rightToLeftConstBegin(), rhs.rightToLeftConstEnd());
        resizeToFit();
        return *this;
    }

    BigUInt BigUInt::operator-(const BigUInt &rhs) const {
        auto copy = *this;
        return copy -= rhs;
    }

    BigUInt BigUInt::operator+(size_t rhs) const {
        auto copy = *this;
        return copy += rhs;
    }

    BigUInt BigUInt::operator+(const BigUInt &rhs) const {
        auto copy = *this;
        return copy += rhs;
    }

    BigUInt &BigUInt::operator%=(const BigUInt &mod) {
        assert(mod != 0ul);

        if (*this < mod) { return *this; }
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

    BigUInt BigUInt::createRandom(size_t numberOfDigits) {
        assert(numberOfDigits > 0ul);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> dis(0, DigitVector::s_maxDigit);

        BigUInt result(0ul);
        result.m_digits[0] = dis(gen);
        for (size_t i = 1; i != numberOfDigits; ++i) {
            result.shift(1);
            result.m_digits[0] = dis(gen);
        }
        result.resizeToFit();
        return result;
    }

    BigUInt BigUInt::createRandomFromDecimalDigits(size_t orderOfMagnitude) {
        assert(orderOfMagnitude > 0);
        const auto numberOfDigits = orderOfMagnitude * (std::log(10) / log(DigitVector::s_base)) + 1ul;
        auto result = createRandom((size_t)numberOfDigits);
        assert(result.isWellFormed());
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

    void BigUInt::divideByLessThanBase(size_t factor) {
        assert(factor < s_base);

        size_t carry = 0ul;
        for (auto thisIt = leftToRightBegin(); thisIt != leftToRightEnd(); ++thisIt) {
            *thisIt += carry * s_base;
            carry = *thisIt % factor;
            *thisIt /= factor;
        }
        resizeToFit();
    }

    void BigUInt::schoolMultiply(rlIterator resultIt,
                                 rlIterator resultEnd,
                                 rlcIterator smallIt,
                                 rlcIterator smallEnd,
                                 rlcIterator largeIt,
                                 rlcIterator largeEnd) {
        for (size_t i = 0; largeIt != largeEnd; ++i) {
            addMultipleViaIterators(resultIt + i, resultEnd, smallIt, smallEnd, *largeIt);
            ++largeIt;
        }
    }

    void BigUInt::toomCook_3(rlIterator resultIt,
                             rlIterator resultEnd,
                             rlcIterator rhsIt,
                             rlcIterator rhsEnd,
                             rlcIterator copyIt,
                             rlcIterator copyEnd) {
        assert((copyEnd - copyIt) >= (rhsEnd - rhsIt));
        const size_t i = (copyEnd - copyIt) / 3ul;

        const BigInt m0 = BigUInt(std::vector<size_t>{rhsIt, rhsIt + i}, false);
        const BigInt m1 = BigUInt(std::vector<size_t>{rhsIt + i, rhsIt + 2ul * i}, false);
        const BigInt m2 = BigUInt(std::vector<size_t>{rhsIt + 2ul * i, rhsEnd}, false);

        const BigInt n0 = BigUInt(std::vector<size_t>{copyIt, copyIt + i}, false);
        const BigInt n1 = BigUInt(std::vector<size_t>{copyIt + i, copyIt + 2ul * i}, false);
        const BigInt n2 = BigUInt(std::vector<size_t>{copyIt + 2ul * i, copyEnd}, false);

        // Towards Optimal Toom-Cook Multiplication for Univariate and Multivariate Polynomials in Characteristic 2 and 0
        // By Bodrato
        const BigInt p_aux = m0 + m2;
        const BigInt p_minusOne = p_aux - m1;
        const BigInt p_minusTwo = 2ul * (p_minusOne + m2) - m0;
        const BigInt q_aux = n0 + n2;
        const BigInt q_minusOne = q_aux - n1;
        const BigInt q_minusTwo = 2ul * (q_minusOne + n2) - n0;

        const BigInt a0 = m0 * n0;
        const BigInt r_one = (p_aux + m1) * (q_aux + n1);
        const BigInt r_minusOne = p_minusOne * q_minusOne;
        const BigInt r_minusTwo = q_minusTwo * p_minusTwo;
        const BigInt a4 = m2 * n2;

        BigInt a1;
        BigInt a2;
        BigInt a3;

        a3 = (r_minusTwo - r_one) / 3ul;
        a1 = (r_one - r_minusOne) / 2ul;
        a2 = r_minusOne - a0;
        a3 = (a2 - a3) / 2ul + 2ul * a4;
        a2 = a2 + a1 - a4;
        a1 -= a3;

        addViaIterators(resultIt, resultEnd, a0.magnitude().rightToLeftConstBegin(), a0.magnitude().rightToLeftConstEnd());
        addViaIterators(resultIt + i, resultEnd, a1.magnitude().rightToLeftConstBegin(), a1.magnitude().rightToLeftConstEnd());
        addViaIterators(
            resultIt + 2ul * i, resultEnd, a2.magnitude().rightToLeftConstBegin(), a2.magnitude().rightToLeftConstEnd());
        addViaIterators(
            resultIt + 3ul * i, resultEnd, a3.magnitude().rightToLeftConstBegin(), a3.magnitude().rightToLeftConstEnd());
        addViaIterators(
            resultIt + 4ul * i, resultEnd, a4.magnitude().rightToLeftConstBegin(), a4.magnitude().rightToLeftConstEnd());
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

    BigUInt &BigUInt::operator/=(const BigUInt &divisor) {
        *this = *this / divisor;
        return *this;
    }

    BigUInt &BigUInt::operator/=(size_t divisor) {
        if (divisor < s_base) {
            divideByLessThanBase(divisor);
        } else {
            *this = *this / divisor;
        }
        return *this;
    }

    void BigUInt::karatsubaMultiplyViaIterators(rlIterator resultIt,
                                                rlIterator resultEnd,
                                                rlcIterator smallIt,
                                                rlcIterator smallEnd,
                                                rlcIterator largeIt,
                                                rlcIterator largeEnd) {
        assert((largeEnd - largeIt) >= (smallEnd - smallIt));

        const auto m = static_cast<size_t>(smallEnd - smallIt);
        assert(m >= s_karatsubaLowerLimit);
        const size_t splitIndex = m / 2ul;

        BigUInt high1(largeIt + splitIndex, largeEnd);
        BigUInt high2(smallIt + splitIndex, smallEnd);

        BigUInt z0;
        z0.resize(2ul * splitIndex + 1ul);
        multiplyViaIterators(
            z0.rightToLeftBegin(), z0.rightToLeftEnd(), smallIt, smallIt + splitIndex, largeIt, largeIt + splitIndex);
        z0.resizeToFit();

        BigUInt z2;
        z2.resize((smallEnd - smallIt - splitIndex) + (largeEnd - largeIt - splitIndex) + 1ul);
        multiplyViaIterators(
            z2.rightToLeftBegin(), z2.rightToLeftEnd(), smallIt + splitIndex, smallEnd, largeIt + splitIndex, largeEnd);
        z2.resizeToFit();

        high1.resize(high1.digitCount() + 1ul);
        high2.resize(high2.digitCount() + 1ul);

        addViaIterators(high1.rightToLeftBegin(), high1.rightToLeftEnd(), largeIt, largeIt + splitIndex);
        addViaIterators(high2.rightToLeftBegin(), high2.rightToLeftEnd(), smallIt, smallIt + splitIndex);

        BigUInt z1;
        z1.resize((smallEnd - smallIt - splitIndex) + (largeEnd - largeIt - splitIndex) + 3ul);

        high1.resizeToFit();
        high2.resizeToFit();

        multiplyViaIterators(z1.rightToLeftBegin(),
                             z1.rightToLeftEnd(),
                             high1.rightToLeftConstBegin(),
                             high1.rightToLeftConstEnd(),
                             high2.rightToLeftConstBegin(),
                             high2.rightToLeftConstEnd()); // z1 = (high1 + low1) * (high2 + low2);

        subtractViaIterators(z1.rightToLeftBegin(), z1.rightToLeftEnd(), z2.rightToLeftConstBegin(), z2.rightToLeftConstEnd());
        subtractViaIterators(z1.rightToLeftBegin(), z1.rightToLeftEnd(), z0.rightToLeftConstBegin(), z0.rightToLeftConstEnd());

        addViaIterators(resultIt, resultEnd, z0.rightToLeftConstBegin(), z0.rightToLeftConstEnd());
        addViaIterators(resultIt + 2ul * splitIndex, resultEnd, z2.rightToLeftConstBegin(), z2.rightToLeftConstEnd());
        addViaIterators(resultIt + splitIndex, resultEnd, z1.rightToLeftConstBegin(), z1.rightToLeftConstEnd());
    }

    void BigUInt::splitOneMultiplicationViaIterators(rlIterator resultIt,
                                                     rlIterator resultEnd,
                                                     rlcIterator smallIt,
                                                     rlcIterator smallEnd,
                                                     rlcIterator largeIt,
                                                     rlcIterator largeEnd) {
        const auto m = static_cast<size_t>(largeEnd - largeIt);
        assert(m >= s_karatsubaLowerLimit);
        const size_t splitIndex = m / 2ul;

        BigUInt z0;
        z0.resize(splitIndex + (smallEnd - smallIt) + 1ul);
        multiplyViaIterators(z0.rightToLeftBegin(), z0.rightToLeftEnd(), largeIt, largeIt + splitIndex, smallIt, smallEnd);
        z0.resizeToFit();

        addViaIterators(resultIt, resultEnd, z0.rightToLeftConstBegin(), z0.rightToLeftConstEnd());

        BigUInt z1;
        z1.resize((largeEnd - (largeIt + splitIndex)) + (smallEnd - smallIt) + 1ul);
        multiplyViaIterators(z1.rightToLeftBegin(), z1.rightToLeftEnd(), largeIt + splitIndex, largeEnd, smallIt, smallEnd);
        z1.resizeToFit();

        addViaIterators(resultIt + splitIndex, resultEnd, z1.rightToLeftConstBegin(), z1.rightToLeftConstEnd());
    }

    void BigUInt::multiplySortedViaIterators(rlIterator resultIt,
                                             const rlIterator resultEnd,
                                             rlcIterator smallIt,
                                             const rlcIterator smallEnd,
                                             rlcIterator largeIt,
                                             const rlcIterator largeEnd) {
        const auto smallSize = static_cast<size_t>(smallEnd - smallIt);
        const auto largeSize = static_cast<size_t>(largeEnd - largeIt);
        assert(largeSize >= smallSize);

        if (largeSize < s_karatsubaLowerLimit) {
            schoolMultiply(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        } else if (smallSize >= s_toomCookLowerLimit) {
            toomCook_3(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        } else if (smallSize >= s_karatsubaLowerLimit) {
            karatsubaMultiplyViaIterators(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        } else {
            splitOneMultiplicationViaIterators(resultIt, resultEnd, smallIt, smallEnd, largeIt, largeEnd);
        }
    }

    void BigUInt::multiplyViaIterators(rlIterator resultIt,
                                       rlIterator resultEnd,
                                       rlcIterator rhsIt,
                                       rlcIterator rhsEnd,
                                       rlcIterator copyIt,
                                       rlcIterator copyEnd) {
        const auto copySize = static_cast<size_t>(copyEnd - copyIt);
        const auto rhsSize = static_cast<size_t>(rhsEnd - rhsIt);
        if (copySize > rhsSize) {
            multiplySortedViaIterators(resultIt, resultEnd, rhsIt, rhsEnd, copyIt, copyEnd);
        } else {
            multiplySortedViaIterators(resultIt, resultEnd, copyIt, copyEnd, rhsIt, rhsEnd);
        }
    }

    bool BigUInt::lessThanShiftedRhsViaIterators(
        lrcIterator thisIt, lrcIterator thisEnd, lrcIterator rhsIt, lrcIterator rhsEnd, size_t trailingZeroesOfRhs) {
        if (static_cast<size_t>(thisEnd - thisIt) != static_cast<size_t>(rhsEnd - rhsIt) + trailingZeroesOfRhs) {
            return static_cast<size_t>(thisEnd - thisIt) < static_cast<size_t>(rhsEnd - rhsIt) + trailingZeroesOfRhs;
        }

        for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
            if (*thisIt != *rhsIt) { return *thisIt < *rhsIt; }
        }
        return false;
    }
    void BigUInt::carryAdditionViaIterators(rlIterator thisIt, rlIterator thisEnd, size_t carry) {
        assert(carry != 0ul);
        assert(carry + DigitVector::s_base < std::numeric_limits<size_t>::max());
        for (; thisIt != thisEnd; ++thisIt) {
            *thisIt += carry;
            if (*thisIt < DigitVector::s_base) { return; }
            carry = (*thisIt) / DigitVector::s_base;
            *thisIt %= DigitVector::s_base;
        }
        assert(false);
    }

    void BigUInt::addViaIterators(rlIterator thisIt, rlIterator thisEnd, rlcIterator rhsIt, rlcIterator rhsEnd) {
        assert(std::distance(thisIt, thisEnd) >= std::distance(rhsIt, rhsEnd) + 1l);
        size_t carry = 0ul;
        for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
            *thisIt += *rhsIt + carry;
            if (*thisIt >= DigitVector::s_base) {
                carry = (*thisIt) / DigitVector::s_base;
                *thisIt %= DigitVector::s_base;
            } else {
                carry = 0ul;
            }
        }
        if (carry == 1ul) { carryAdditionViaIterators(thisIt, thisEnd, 1ul); }
    }

    void BigUInt::addMultipleViaIterators(
        rlIterator thisIt, rlIterator thisEnd, rlcIterator rhsIt, rlcIterator rhsEnd, size_t multiplier) {
        if (multiplier == 0ul) { return; }
        assert(multiplier < DigitVector::s_base);
        assert(std::distance(thisIt, thisEnd) >= std::distance(rhsIt, rhsEnd) + 1l);
        size_t carry = 0ul;
        for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
            *thisIt += (*rhsIt * multiplier) + carry;
            if (*thisIt >= DigitVector::s_base) {
                carry = (*thisIt) / DigitVector::s_base;
                *thisIt %= DigitVector::s_base;
            } else {
                carry = 0ul;
            }
        }
        if (carry > 0UL) { carryAdditionViaIterators(thisIt, thisEnd, carry); }
    }
    void BigUInt::subtractViaIterators(rlIterator thisIt, rlIterator thisEnd, rlcIterator rhsIt, rlcIterator rhsEnd) {
        assert(std::distance(thisIt, thisEnd) >= std::distance(rhsIt, rhsEnd));
        size_t carry = 0ul;
        for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
            if (*thisIt >= *rhsIt + carry) {
                *thisIt -= *rhsIt + carry;
                carry = 0;
            } else {
                *thisIt += DigitVector::s_base;
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
                *thisIt = DigitVector::s_maxDigit;
                ++thisIt;
            }
        }
    }

    std::string BigUInt::toString() const {
        std::stringstream ss;
        ss << *leftToRightConstBegin();
        for (auto it = leftToRightConstBegin() + 1ul; it != leftToRightConstEnd(); ++it) {
            ss << std::setfill('0') << std::setw(s_digitsPerLimb) << *it;
        }
        return ss.str();
    }
} // namespace big