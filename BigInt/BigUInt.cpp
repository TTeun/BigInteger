#include "BigUInt.h"

#include "BigUIntBase.h"

#include <cassert>
#include <iostream>
#include <random>
#include <sstream>

/* =================== Static functions =================== */

static void carryAdditionViaIterators(rightToLeftIterator thisIt, const rightToLeftIterator &thisEnd, size_t carry) {
    assert(carry != 0ul);
    assert(carry + DigitVector::s_base < std::numeric_limits<size_t>::max());
    for (; thisIt != thisEnd; ++thisIt) {
        *thisIt += carry;
        if (*thisIt < DigitVector::s_base) {
            return;
        }
        carry = (*thisIt) / DigitVector::s_base;
        *thisIt %= DigitVector::s_base;
    }
    assert(false);
}

static void addViaIterators(rightToLeftIterator                       thisIt,
                            const rightToLeftIterator                 thisEnd,
                            std::vector<size_t>::const_iterator       rhsIt,
                            const std::vector<size_t>::const_iterator rhsEnd) {
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
    if (carry == 1ul) {
        carryAdditionViaIterators(thisIt, thisEnd, 1ul);
    }
}

static void addMultipleViaIterators(rightToLeftIterator                       thisIt,
                                    const rightToLeftIterator                 thisEnd,
                                    std::vector<size_t>::const_iterator       rhsIt,
                                    const std::vector<size_t>::const_iterator rhsEnd,
                                    const size_t                              multiplier) {
    if (multiplier == 0ul) {
        return;
    }
    assert(multiplier < DigitVector::s_base);
    assert(std::distance(thisIt, thisEnd) >= std::distance(rhsIt, rhsEnd) + 1l); // Always enough space for adding two proper digit vectors
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
    if (carry > 0UL) {
        carryAdditionViaIterators(thisIt, thisEnd, carry);
    }
}

static void subtractViaIterators(rightToLeftIterator                       thisIt,
                                 const rightToLeftIterator                 thisEnd,
                                 std::vector<size_t>::const_iterator       rhsIt,
                                 const std::vector<size_t>::const_iterator rhsEnd) {
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

static bool lessThanShiftedRhsViaIterators(leftToRightConstIterator        thisIt,
                                           const leftToRightConstIterator &thisEnd,
                                           leftToRightConstIterator        rhsIt,
                                           const leftToRightConstIterator &rhsEnd,
                                           const size_t                    trailingZeroesOfRhs) {
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

static bool lessThanViaIterators(const leftToRightConstIterator &thisIt,
                                 const leftToRightConstIterator &thisEnd,
                                 const leftToRightConstIterator &rhsIt,
                                 const leftToRightConstIterator &rhsEnd) {
    return lessThanShiftedRhsViaIterators(thisIt, thisEnd, rhsIt, rhsEnd, 0ul);
}

static bool greaterThanViaIterators(const leftToRightConstIterator &thisIt,
                                    const leftToRightConstIterator &thisEnd,
                                    const leftToRightConstIterator &rhsIt,
                                    const leftToRightConstIterator &rhsEnd) {
    return lessThanViaIterators(rhsIt, rhsEnd, thisIt, thisEnd);
}

static size_t ceilingIntegerDivision(size_t a, size_t b) {
    return (a + b - 1) / b;
}

template <typename T>
static size_t modPower(T base, T exponent, size_t mod) {
    if (exponent == 0ul) {
        return 1;
    } else if (exponent == 1ul) {
        return base % mod;
    }
    T aux(1);
    T copy(base);
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

void bubbleViaIterators(rightToLeftIterator thisIt, const rightToLeftIterator &thisEnd) {
    auto next = thisIt + 1;

    for (; next != thisEnd; ++thisIt, ++next) {
        if (*thisIt >= DigitVector::s_base) {
            *next += *thisIt / DigitVector::s_base;
            *thisIt %= DigitVector::s_base;
        }
    }
}

/* =================== Friend functions =================== */

size_t divisionSubRoutine(const leftToRightConstIterator &leftToRightConstIt,
                          const leftToRightConstIterator &leftToRightConstEnd,
                          const rightToLeftIterator &     rightToLeftIt,
                          const rightToLeftIterator &     rightToLeftEnd,
                          const BigUInt &                 divisor) {
    if (lessThanViaIterators(leftToRightConstIt, leftToRightConstEnd, divisor.leftToRightConstBegin(), divisor.leftToRightConstEnd())) {
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
        leftToRightConstIt, leftToRightConstEnd, divisor.leftToRightConstBegin(), divisor.leftToRightConstEnd(), 1)) {
        // Divisor fits into dividend at least s_base times, so we subtract the divisor shifted by one digit
        subtractViaIterators(rightToLeftIt + 1ul, rightToLeftEnd, divisor.rightToLeftConstBegin(), divisor.rightToLeftConstEnd());
        correction += BigUInt::s_base;
    }

    size_t quotientEstimate;
    if (m == n) {
        quotientEstimate = (*leftToRightConstIt) / divisor.mostSignificantDigit();
    } else {
        quotientEstimate = (*leftToRightConstIt * BigUInt::s_base + *(leftToRightConstIt + 1ul)) / divisor.mostSignificantDigit();
    }
    quotientEstimate    = std::min(quotientEstimate, DigitVector::s_base - 1ul);
    const size_t offset = *leftToRightConstIt == 0 ? 1ul : 0ul;

    auto closestMultipleEstimate = divisor * quotientEstimate;
    while (greaterThanViaIterators(closestMultipleEstimate.leftToRightConstBegin(),
                                   closestMultipleEstimate.leftToRightConstEnd(),
                                   leftToRightConstIt + offset,
                                   leftToRightConstEnd)) {
        --quotientEstimate;
        closestMultipleEstimate -= divisor;
    }
    subtractViaIterators(
        rightToLeftIt, rightToLeftEnd, closestMultipleEstimate.rightToLeftConstBegin(), closestMultipleEstimate.rightToLeftConstEnd());
    return quotientEstimate + correction;
}

BigUInt longDivision(BigUInt &dividend, const BigUInt &divisor) {
    if (dividend < divisor) {
        return 0ul;
    }
    const size_t t = ceilingIntegerDivision(BigUInt::s_base, 2ul * divisor.mostSignificantDigit());
    if (t > 1) {
        return longDivisionAfterAdjustingDivisor(dividend *= t, divisor * t);
    } else {
        return longDivisionAfterAdjustingDivisor(dividend, divisor);
    }
}

static void resizeToFitVector(std::vector<size_t> &digits) {
    auto it = digits.rbegin();
    for (; it != digits.rend() && *it == 0ul; ++it)
        ;

    digits.resize(static_cast<unsigned long>(std::distance(it, digits.rend())));

    if (digits.empty()) {
        digits = {0ul};
    }
}

BigUInt longDivisionAfterAdjustingDivisor(BigUInt &dividend, const BigUInt &divisor) {
    assert(divisor <= dividend);
    assert(divisor.mostSignificantDigit() * 2ul >= BigUInt::s_base);

    size_t       m = dividend.digitCount();
    const size_t n = divisor.digitCount();

    if (m <= n + 1) {
        return divisionSubRoutine(dividend.leftToRightConstBegin(),
                                  dividend.leftToRightConstEnd(),
                                  dividend.rightToLeftBegin(),
                                  dividend.rightToLeftEnd(),
                                  divisor);
    }

    std::vector<size_t> divisorDigits(m - n + 1ul, 0ul);
    while (m > n + 1) {
        const size_t splitIndex   = m - n - 1ul;
        divisorDigits[splitIndex] = divisionSubRoutine(dividend.leftToRightConstBegin(),
                                                       dividend.leftToRightConstBegin() + n + 1,
                                                       dividend.rightToLeftBegin() + splitIndex,
                                                       dividend.rightToLeftEnd(),
                                                       divisor);
        // The first n + 1 digits of dividend now contain their remainder after division by divisor, as per the divisionSubRoutine
        dividend.resizeToFit();
        m = dividend.digitCount();
    }
    if (m > n || dividend.mostSignificantDigit() >= divisor.mostSignificantDigit()) {
        divisorDigits[0ul] = divisionSubRoutine(dividend.leftToRightConstBegin(),
                                                dividend.leftToRightConstEnd(),
                                                dividend.rightToLeftBegin(),
                                                dividend.rightToLeftEnd(),
                                                divisor);
    }
    bubbleViaIterators(divisorDigits.begin(), divisorDigits.end());
    resizeToFitVector(divisorDigits);
    return BigUInt(std::move(divisorDigits));
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
        if (exponent % 2ul == 1ul) {
            aux *= copy;
        }
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

void karatsubaMultiplyViaIterators(rightToLeftIterator             resultIt,
                                   const rightToLeftIterator &     resultEnd,
                                   rightToLeftConstIterator        rhsIt,
                                   const rightToLeftConstIterator &rhsEnd,
                                   rightToLeftConstIterator        copyIt,
                                   const rightToLeftConstIterator &copyEnd);

void splitOneMultiplicationViaIterators(rightToLeftIterator             resultIt,
                                        const rightToLeftIterator &     resultEnd,
                                        rightToLeftConstIterator        rhsIt,
                                        const rightToLeftConstIterator &rhsEnd,
                                        rightToLeftConstIterator        largeIt,
                                        const rightToLeftConstIterator &largeEnd);

void multiplyViaIterators(rightToLeftIterator             resultIt,
                          const rightToLeftIterator &     resultEnd,
                          rightToLeftConstIterator        rhsIt,
                          const rightToLeftConstIterator &rhsEnd,
                          rightToLeftConstIterator        copyIt,
                          const rightToLeftConstIterator &copyEnd) {
    const auto copySize = static_cast<size_t>(copyEnd - copyIt);
    const auto rhsSize  = static_cast<size_t>(rhsEnd - rhsIt);
    if (std::min(copySize, rhsSize) >= BigUInt::s_karatsubaLowerLimit) {
        karatsubaMultiplyViaIterators(resultIt, resultEnd, rhsIt, rhsEnd, copyIt, copyEnd);
    } else if (copySize >= BigUInt::s_karatsubaLowerLimit) {
        splitOneMultiplicationViaIterators(resultIt, resultEnd, rhsIt, rhsEnd, copyIt, copyEnd);
    } else {
        for (size_t i = 0; rhsIt != rhsEnd; ++i) {
            addMultipleViaIterators(resultIt + i, resultEnd, copyIt, copyEnd, *rhsIt);
            ++rhsIt;
        }
    }
}

void splitOneMultiplicationViaIterators(rightToLeftIterator             resultIt,
                                        const rightToLeftIterator &     resultEnd,
                                        rightToLeftConstIterator        rhsIt,
                                        const rightToLeftConstIterator &rhsEnd,
                                        rightToLeftConstIterator        largeIt,
                                        const rightToLeftConstIterator &largeEnd) {
    const auto m = static_cast<size_t>(largeEnd - largeIt);
    assert(m >= BigUInt::s_karatsubaLowerLimit);
    const size_t splitIndex = m / 2ul;

    BigUIntBase z0;
    z0.resize(splitIndex + (rhsEnd - rhsIt));
    multiplyViaIterators(z0.rightToLeftBegin(), z0.rightToLeftEnd(), largeIt, largeIt + splitIndex, rhsIt, rhsEnd);
    resizeToFitVector(z0.m_digits);

    addViaIterators(resultIt, resultEnd, z0.rightToLeftConstBegin(), z0.rightToLeftConstEnd());

    BigUIntBase z1;
    z1.resize((largeEnd - (largeIt + splitIndex)) + (rhsEnd - rhsIt));
    multiplyViaIterators(z1.rightToLeftBegin(), z1.rightToLeftEnd(), largeIt + splitIndex, largeEnd, rhsIt, rhsEnd);
    resizeToFitVector(z1.m_digits);

    addViaIterators(resultIt + splitIndex, resultEnd, z1.rightToLeftConstBegin(), z1.rightToLeftConstEnd());
}

void karatsubaMultiplyViaIterators(rightToLeftIterator             resultIt,
                                   const rightToLeftIterator &     resultEnd,
                                   rightToLeftConstIterator        rhsIt,
                                   const rightToLeftConstIterator &rhsEnd,
                                   rightToLeftConstIterator        copyIt,
                                   const rightToLeftConstIterator &copyEnd) {
    const auto m = static_cast<size_t>(std::min(copyEnd - copyIt, rhsEnd - rhsIt));
    assert(m >= BigUInt::s_karatsubaLowerLimit);
    const size_t splitIndex = m / 2ul;

    BigUIntBase high1(copyIt + splitIndex, copyEnd);
    BigUIntBase high2(rhsIt + splitIndex, rhsEnd);

    BigUIntBase z0;
    z0.resize(2ul * splitIndex + 1ul);
    multiplyViaIterators(
        z0.rightToLeftBegin(), z0.rightToLeftEnd(), rhsIt, rhsIt + splitIndex, copyIt, copyIt + splitIndex); // z0 = low1 * low2;
    resizeToFitVector(z0.m_digits);

    BigUIntBase z2;
    z2.resize((rhsEnd - rhsIt - splitIndex) + (copyEnd - copyIt - splitIndex) + 1ul);
    multiplyViaIterators(
        z2.rightToLeftBegin(), z2.rightToLeftEnd(), rhsIt + splitIndex, rhsEnd, copyIt + splitIndex, copyEnd); // z2 = high1 * high2;
    resizeToFitVector(z2.m_digits);

    high1.resize(high1.digitCount() + 1ul);
    high2.resize(high2.digitCount() + 1ul);

    addViaIterators(high1.rightToLeftBegin(), high1.rightToLeftEnd(), copyIt, copyIt + splitIndex);
    addViaIterators(high2.rightToLeftBegin(), high2.rightToLeftEnd(), rhsIt, rhsIt + splitIndex);

    BigUIntBase z1;
    z1.resize((rhsEnd - rhsIt - splitIndex) + (copyEnd - copyIt - splitIndex) + 3ul);

    resizeToFitVector(high1.m_digits);
    resizeToFitVector(high2.m_digits);

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

/* =================== Member functions =================== */

BigUInt::BigUInt() {
    init(0);
    assert(isWellFormed());
}

BigUInt::BigUInt(size_t val) {
    init(val);
    assert(isWellFormed());
}

BigUInt::BigUInt(std::vector<size_t> &&digits) : BigUIntBase(std::move(digits)) {
    if (m_digits.empty()) {
        init(0);
    }
    assert(isWellFormed());
}

BigUInt::BigUInt(const BigUInt &other) : BigUIntBase(std::vector<size_t>(other.m_digits.begin(), other.m_digits.end())) {
    assert(isWellFormed());
}

BigUInt::BigUInt(const std::string &val) {
    static const auto   maximumDecimalDigitsInSizeType = static_cast<size_t>(std::log10(std::numeric_limits<size_t>::max()));
    static const size_t largestMultipleOfTenInSizeType = modPower(10ul, maximumDecimalDigitsInSizeType, std::numeric_limits<size_t>::max());

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

BigUInt::BigUInt(rightToLeftConstIterator it, rightToLeftConstIterator endIt) : BigUIntBase(std::vector<size_t>{it, endIt}) {
    if (it == endIt) {
        m_digits = {0ul};
    } else {
        resizeToFit();
    }
}

void BigUInt::bubble(size_t startIndex) {
    assert(not m_digits.empty());
    assert(startIndex < digitCount());
    auto it   = rightToLeftBegin() + startIndex;
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

bool BigUInt::operator!=(const BigUInt &rhs) const {
    return !(rhs == *this);
}

void BigUInt::resizeToFit() {
    resizeToFitVector(m_digits);
}

BigUInt &BigUInt::operator=(const BigUInt &rhs) {
    if (this == &rhs) {
        return *this;
    }
    m_digits = rhs.m_digits;
    return *this;
}

BigUInt &BigUInt::operator+=(const BigUInt &rhs) {
    assert(isWellFormed() && rhs.isWellFormed());

    resize(std::max(digitCount(), rhs.digitCount()) + 1ul);
    addViaIterators(rightToLeftBegin(), rightToLeftEnd(), rhs.rightToLeftConstBegin(), rhs.rightToLeftConstEnd());
    if (mostSignificantDigit() == 0ul) {
        resize(digitCount() - 1ul);
    }

    return *this;
}

bool BigUInt::operator<(const BigUInt &rhs) const {
    assert(isWellFormed() && rhs.isWellFormed());

    return lessThanViaIterators(leftToRightConstBegin(), leftToRightConstEnd(), rhs.leftToRightConstBegin(), rhs.leftToRightConstEnd());
}

bool BigUInt::operator>(const BigUInt &rhs) const {
    return rhs < *this;
}

bool BigUInt::operator<=(const BigUInt &rhs) const {
    assert(isWellFormed() && rhs.isCorrectlySized());

    if (digitCount() != rhs.digitCount()) {
        return digitCount() < rhs.digitCount();
    }

    auto thisIt = leftToRightConstBegin();
    auto rhsIt  = rhs.leftToRightConstBegin();
    for (; thisIt != leftToRightConstEnd(); ++thisIt, ++rhsIt) {
        if (*thisIt != *rhsIt) {
            return *thisIt <= *rhsIt;
        }
    }
    return true;
}

bool BigUInt::operator>=(const BigUInt &rhs) const {
    return !(*this < rhs);
}

BigUInt &BigUInt::operator*=(const BigUInt &rhs) {
    assert(isWellFormed());
    assert(rhs.isWellFormed());
    if (this == &rhs) {
        square();
        return *this;
    }
    auto copy = *this;

    this->resize(digitCount() + rhs.digitCount() + 1);
    for (auto thisIt = leftToRightBegin(); thisIt != leftToRightEnd(); ++thisIt) {
        *thisIt = 0ul;
    }

    multiplyViaIterators(rightToLeftBegin(),
                         rightToLeftEnd(),
                         rhs.rightToLeftConstBegin(),
                         rhs.rightToLeftConstEnd(),
                         copy.rightToLeftBegin(),
                         copy.rightToLeftEnd());
    resizeToFit();
    assert(isWellFormed());
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
    if (rhs == 0) {
        init(0);
        return *this;
    }
    reserve(digitCount() + 3ul);
    if (rhs != 1ul) {
        if (rhs < s_base) {
            for (auto &it : m_digits) {
                it *= rhs;
            }
            bubble();
        } else {
            *this *= BigUInt(rhs);
        }
    }
    resizeToFit();
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

BigUInt &BigUInt::operator=(size_t rhs) {
    init(rhs);
    return *this;
}

size_t BigUInt::operator%(size_t mod) const {
    if (mod == 2ul) {
        return leastSignificantDigit() % 2ul;
    }
    if (mod >= std::numeric_limits<size_t>::max() / mod) {
        std::cout << "large modulo moet nog\n";
        return 0;
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
    if (&divisor == this) {
        return BigUInt(1);
    }
    if (divisor == 1ul) {
        return BigUInt(*this);
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

BigUInt BigUInt::operator*(const BigUInt &rhs) const {
    BigUInt copy(*this);
    return copy *= rhs;
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

BigUInt BigUInt::copyPrefix(size_t length) const {
    assert(length <= digitCount());
    return BigUInt(std::vector<size_t>(rightToLeftConstEnd() - length, rightToLeftConstEnd()));
}

BigUInt BigUInt::copySuffix(size_t length) const {
    assert(length <= digitCount());
    return BigUInt(std::vector<size_t>(rightToLeftConstBegin(), rightToLeftConstBegin() + length));
}

BigUInt BigUInt::createFromMostSignificantDigit(size_t digit, size_t position) {
    std::vector<size_t> digits(position, 0);
    digits[position] = digit;
    return BigUInt(std::move(digits));
}

BigUInt BigUInt::shiftedCopy(size_t shiftAmount) const {
    auto copy = *this;
    copy.shift(shiftAmount);
    return copy;
}

BigUInt &BigUInt::operator%=(const BigUInt &mod) {
    assert(mod != 0ul);
    BigUInt copy(*this);
    *this -= longDivision(copy, mod) * mod;
    return *this;
}

BigUInt BigUInt::createRandom(size_t numberOfDigits) {
    assert(numberOfDigits > 0ul);
    std::random_device                    rd;        // Will be used to obtain a seed for the random number engine
    std::mt19937                          gen(rd()); // Standard mersenne_twister_engine seeded with rd()
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

size_t BigUInt::approximatePowerOfTen() const {
    return static_cast<size_t>(std::log10(mostSignificantDigit() + 1ul) + (digitCount() - 1ul) * std::log10(DigitVector::s_base));
}

BigUInt BigUInt::createRandomFromDecimalDigits(size_t orderOfMagnitude) {
    assert(orderOfMagnitude > 0);
    const auto numberOfDigits = orderOfMagnitude * (std::log(10) / log(DigitVector::s_base)) + 1ul;
    auto       result         = createRandom((size_t)numberOfDigits);
    assert(result.isWellFormed());
    return result;
}

BigUInt BigUInt::operator%(const BigUInt &mod) const {
    auto copy = *this;
    if (mod.digitCount() == 1ul) {
        copy %= mod.leastSignificantDigit();
    } else if (mod.digitCount() == 2ul) {
        copy %= mod.mostSignificantDigit() * s_base + mod.leastSignificantDigit();
    }
    return copy;
}