#include "BigUnsignedInt.h"

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

static void multiplyViaIterators(rightToLeftIterator             thisIt,
                                 const rightToLeftIterator &     thisEnd,
                                 rightToLeftConstIterator        rhsIt,
                                 const rightToLeftConstIterator &rhsEnd,
                                 rightToLeftConstIterator        copyIt,
                                 const rightToLeftConstIterator &copyEnd) {
    for (size_t i = 0; rhsIt != rhsEnd; ++i) {
        addMultipleViaIterators(thisIt + i, thisEnd, copyIt, copyEnd, *rhsIt);
        ++rhsIt;
    }
}

static bool lessThanShiftedRhsViaIterators(leftToRightConstIterator        thisIt,
                                           const leftToRightConstIterator &thisEnd,
                                           leftToRightConstIterator        rhsIt,
                                           const leftToRightConstIterator &rhsEnd,
                                           const size_t                    trailingZeroesOfRhs) {
    if (std::distance(thisIt, thisEnd) != std::distance(rhsIt, rhsEnd) + trailingZeroesOfRhs) {
        return std::distance(thisIt, thisEnd) < std::distance(rhsIt, rhsEnd) + trailingZeroesOfRhs;
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
                          const BigUnsignedInt &          divisor) {
    if (lessThanViaIterators(leftToRightConstIt, leftToRightConstEnd, divisor.leftToRightConstBegin(), divisor.leftToRightConstEnd())) {
        return 0ul;
    }
    assert(divisor != 0);
    assert(divisor.mostSignificantDigit() * 2 >= BigUnsignedInt::s_base);
    const size_t n = divisor.digitCount();
    assert(std::distance(leftToRightConstIt, leftToRightConstEnd) <= n + 1);
    assert(std::distance(leftToRightConstIt, leftToRightConstEnd) >= n);
    size_t correction = 0ul;

    while (not lessThanShiftedRhsViaIterators(
        leftToRightConstIt, leftToRightConstEnd, divisor.leftToRightConstBegin(), divisor.leftToRightConstEnd(), 1)) {
        subtractViaIterators(rightToLeftIt + 1ul, rightToLeftEnd, divisor.rightToLeftConstBegin(), divisor.rightToLeftConstEnd());
        correction += BigUnsignedInt::s_base;
    }

    size_t quotientEstimate;
    if (std::distance(leftToRightConstIt, leftToRightConstEnd) == n) {
        quotientEstimate = (*leftToRightConstIt) / divisor.mostSignificantDigit();
    } else {
        quotientEstimate = (*leftToRightConstIt * BigUnsignedInt::s_base + *(leftToRightConstIt + 1ul)) / divisor.mostSignificantDigit();
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

BigUnsignedInt longDivision(BigUnsignedInt &dividend, const BigUnsignedInt &divisor) {
    if (dividend < divisor) {
        return 0ul;
    }
    const size_t t = ceilingIntegerDivision(BigUnsignedInt::s_base, 2ul * divisor.mostSignificantDigit());
    if (t > 1) {
        return longDivisionAfterAdjustingDivisor(dividend *= t, divisor * t);
    } else {
        return longDivisionAfterAdjustingDivisor(dividend, divisor);
    }
}

void resizeToFitVector(std::vector<size_t> &digits) {
    auto it = digits.rbegin();
    for (; it != digits.rend() && *it == 0ul; ++it)
        ;

    digits.resize(std::distance(it, digits.rend()));

    if (digits.empty()) {
        digits = {0ul};
    }
}

BigUnsignedInt longDivisionAfterAdjustingDivisor(BigUnsignedInt &dividend, const BigUnsignedInt &divisor) {
    assert(divisor <= dividend);
    assert(divisor.mostSignificantDigit() * 2ul >= BigUnsignedInt::s_base);

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
        const size_t splitIndex = m - n - 1ul;

        auto prefix = dividend.prefix(n + 1);

        auto quotientRemainder = divisionSubRoutine(
            prefix.leftToRightConstBegin(), prefix.leftToRightConstEnd(), prefix.rightToLeftBegin(), prefix.rightToLeftEnd(), divisor);
        divisorDigits[splitIndex] = quotientRemainder;
        prefix.resizeToFit();

        dividend.m_digits.resize(splitIndex);
        dividend.m_digits.insert(dividend.m_digits.end(), prefix.m_digits.begin(), prefix.m_digits.end());
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
    return BigUnsignedInt(std::move(divisorDigits));
    ;
}

void swap(BigUnsignedInt &a, BigUnsignedInt &b) {
    std::swap(a.m_digits, b.m_digits);
}

std::ostream &operator<<(std::ostream &os, const BigUnsignedInt &bigUnsignedInt) {
    os << "( ";
    for (auto it = bigUnsignedInt.leftToRightConstBegin(); it != bigUnsignedInt.leftToRightConstEnd(); ++it) {
        os << *it << " ";
    }
    os << ")_" << BigUnsignedInt::s_base;
    return os;
}

BigUnsignedInt power(const BigUnsignedInt &base, size_t exponent) {
    if (exponent == 0ul) {
        return BigUnsignedInt(1);
    } else if (exponent == 1ul) {
        return BigUnsignedInt(base);
    }
    BigUnsignedInt aux(1);
    BigUnsignedInt copy(base);

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

/* =================== Member functions =================== */

BigUnsignedInt::BigUnsignedInt() {
    init(0);
    assert(isWellFormed());
}

BigUnsignedInt::BigUnsignedInt(size_t val) {
    init(val);
    assert(isWellFormed());
}

BigUnsignedInt::BigUnsignedInt(std::vector<size_t> &&digits) : DigitVector(std::move(digits)) {
    if (m_digits.empty()) {
        init(0);
    }
    assert(isWellFormed());
}

BigUnsignedInt::BigUnsignedInt(const BigUnsignedInt &other)
    : DigitVector(std::vector<size_t>(other.m_digits.begin(), other.m_digits.end())) {
    assert(isWellFormed());
}

BigUnsignedInt::BigUnsignedInt(const std::string &val) {
    static const size_t maximumDecimalDigitsInSizeType = std::log10(std::numeric_limits<size_t>::max());
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

BigUnsignedInt::BigUnsignedInt(rightToLeftConstIterator it, rightToLeftConstIterator endIt) : DigitVector({it, endIt}) {
    if (it == endIt) {
        m_digits = {0ul};
    } else {
        resizeToFit();
    }
}

void BigUnsignedInt::bubble(size_t startIndex) {
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
        m_digits.resize(1ul + static_cast<size_t>(digitCount() + std::log(mostSignificantDigit()) / std::log(s_base)));
        bubble(continueIndex);
    }
    resizeToFit();
}

void BigUnsignedInt::init(size_t val) {
    static_assert(s_base <= std::numeric_limits<size_t>::max() / s_base, "s_base^2 should not exceed the maximum size_t");
    static_assert(s_base % 2 == 0, "s_base should be even");
    m_digits = {val};
    bubble();
}

bool BigUnsignedInt::operator==(const BigUnsignedInt &rhs) const {
    assert(isWellFormed());
    assert(rhs.isWellFormed());

    return m_digits == rhs.m_digits;
}

bool BigUnsignedInt::operator!=(const BigUnsignedInt &rhs) const {
    return !(rhs == *this);
}

void BigUnsignedInt::resizeToFit() {
    resizeToFitVector(m_digits);
}

BigUnsignedInt &BigUnsignedInt::operator=(const BigUnsignedInt &rhs) {
    if (this == &rhs) {
        return *this;
    }
    m_digits = rhs.m_digits;
    return *this;
}

BigUnsignedInt &BigUnsignedInt::operator+=(const BigUnsignedInt &rhs) {
    m_digits.resize(std::max(digitCount(), rhs.digitCount()) + 1ul);
    addViaIterators(rightToLeftBegin(), rightToLeftEnd(), rhs.rightToLeftConstBegin(), rhs.rightToLeftConstEnd());
    if (m_digits.back() == 0ul) {
        m_digits.resize(m_digits.size() - 1ul);
    }

    return *this;
}

bool BigUnsignedInt::operator<(const BigUnsignedInt &rhs) const {
    assert(isWellFormed() && rhs.isWellFormed());

    return lessThanViaIterators(leftToRightConstBegin(), leftToRightConstEnd(), rhs.leftToRightConstBegin(), rhs.leftToRightConstEnd());
}

bool BigUnsignedInt::operator>(const BigUnsignedInt &rhs) const {
    return rhs < *this;
}

bool greaterOrEqualViaIterators(const leftToRightConstIterator &thisIt,
                                const leftToRightConstIterator &thisEnd,
                                const leftToRightConstIterator &rhsIt,
                                const leftToRightConstIterator &rhsEnd) {
    return !(lessThanViaIterators(thisIt, thisEnd, rhsIt, rhsEnd));
}

bool BigUnsignedInt::operator<=(const BigUnsignedInt &rhs) const {
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

bool BigUnsignedInt::operator>=(const BigUnsignedInt &rhs) const {
    return !(*this < rhs);
}

void karatsubaMultiplyViaIterators(rightToLeftIterator             thisIt,
                                   const rightToLeftIterator &     thisEnd,
                                   rightToLeftConstIterator        rhsIt,
                                   const rightToLeftConstIterator &rhsEnd,
                                   rightToLeftIterator             copyIt,
                                   const rightToLeftIterator &     copyEnd) {
    const size_t m = std::min(copyEnd - copyIt, rhsEnd - rhsIt);
    if (m < 200ul) {
        multiplyViaIterators(thisIt, thisEnd, rhsIt, rhsEnd, copyIt, copyEnd);
        return;
    }
    const size_t splitIndex = m / 2ul;

    BigUnsignedInt high1(copyIt + splitIndex, copyEnd);
    BigUnsignedInt high2(rhsIt + splitIndex, rhsEnd);

    const auto     z2 = high1 * high2;
    BigUnsignedInt z0;
    z0.m_digits.resize(2ul * splitIndex + 1);
    karatsubaMultiplyViaIterators(z0.rightToLeftBegin(), z0.rightToLeftEnd(), rhsIt, rhsIt + splitIndex, copyIt, copyIt + splitIndex);
    z0.resizeToFit();
    addViaIterators(thisIt, thisEnd, z0.rightToLeftConstBegin(), z0.rightToLeftConstEnd());

    addViaIterators(thisIt + 2ul * splitIndex, thisEnd, z2.rightToLeftConstBegin(), z2.rightToLeftConstEnd());

    high1.m_digits.push_back(0ul);
    addViaIterators(high1.rightToLeftBegin(), high1.rightToLeftEnd(), copyIt, copyIt + splitIndex);
    high1.resizeToFit();

    high2.m_digits.push_back(0ul);
    addViaIterators(high2.rightToLeftBegin(), high2.rightToLeftEnd(), rhsIt, rhsIt + splitIndex);
    high2.resizeToFit();

    auto z1 = high1 * high2 - (z2 + z0);

    addViaIterators(thisIt + splitIndex, thisEnd, z1.rightToLeftConstBegin(), z1.rightToLeftConstEnd());
}

BigUnsignedInt &BigUnsignedInt::operator*=(const BigUnsignedInt &rhs) {
    assert(isWellFormed());
    assert(rhs.isWellFormed());
    if (this == &rhs) {
        square();
        return *this;
    }
    auto copy = *this;

    this->m_digits.resize(digitCount() + rhs.digitCount() + 1);
    for (auto thisIt = leftToRightBegin(); thisIt != leftToRightEnd(); ++thisIt) {
        *thisIt = 0ul;
    }

    karatsubaMultiplyViaIterators(rightToLeftBegin(),
                                  rightToLeftEnd(),
                                  rhs.rightToLeftConstBegin(),
                                  rhs.rightToLeftConstEnd(),
                                  copy.rightToLeftBegin(),
                                  copy.rightToLeftEnd());
    resizeToFit();
    assert(isWellFormed());
    return *this;
}

BigUnsignedInt BigUnsignedInt::operator*(size_t rhs) const {
    auto copy = *this;
    return copy *= rhs;
}

void BigUnsignedInt::square() {
    const auto copy = *this;
    *this *= copy;
    bubble(0);
}

BigUnsignedInt &BigUnsignedInt::operator*=(const size_t rhs) {
    if (rhs == 0) {
        init(0);
        return *this;
    }
    m_digits.reserve(m_digits.size() + 3ul);
    if (rhs != 1ul) {
        if (rhs < s_base) {
            for (auto &it : m_digits) {
                it *= rhs;
            }
            bubble();
        } else {
            *this *= BigUnsignedInt(rhs);
        }
    }
    resizeToFit();
    return *this;
}

BigUnsignedInt &BigUnsignedInt::operator+=(size_t rhs) {
    static const size_t additionRoom = std::numeric_limits<size_t>::max() - s_base;
    if (rhs <= additionRoom) {
        *rightToLeftBegin() += rhs;
        bubble();
    } else {
        *this += BigUnsignedInt(rhs);
    }
    return *this;
}

BigUnsignedInt &BigUnsignedInt::operator=(size_t rhs) {
    init(rhs);
    return *this;
}

size_t BigUnsignedInt::operator%(size_t mod) const {
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

BigUnsignedInt &BigUnsignedInt::operator%=(size_t mod) {
    *this = *this % mod;
    return *this;
}

BigUnsignedInt BigUnsignedInt::operator/(const BigUnsignedInt &divisor) const {
    assert(divisor != 0ul);
    if (&divisor == this) {
        return BigUnsignedInt(1);
    }
    if (divisor == 1ul) {
        return BigUnsignedInt(*this);
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
    BigUnsignedInt copy(*this);
    return longDivision(copy, divisor);
}

BigUnsignedInt BigUnsignedInt::operator*(const BigUnsignedInt &rhs) const {
    BigUnsignedInt copy(*this);
    return copy *= rhs;
}

BigUnsignedInt &BigUnsignedInt::operator-=(const BigUnsignedInt &rhs) {
    assert(rhs <= *this);

    subtractViaIterators(rightToLeftBegin(), rightToLeftEnd(), rhs.rightToLeftConstBegin(), rhs.rightToLeftConstEnd());
    resizeToFit();
    return *this;
}

BigUnsignedInt BigUnsignedInt::operator-(const BigUnsignedInt &rhs) const {
    auto copy = *this;
    return copy -= rhs;
}

BigUnsignedInt BigUnsignedInt::operator+(size_t rhs) const {
    auto copy = *this;
    return copy += rhs;
}

BigUnsignedInt BigUnsignedInt::operator+(const BigUnsignedInt &rhs) const {
    auto copy = *this;
    return copy += rhs;
}

BigUnsignedInt BigUnsignedInt::prefix(size_t length) const {
    assert(length <= digitCount());
    return BigUnsignedInt(std::vector<size_t>(rightToLeftConstEnd() - length, rightToLeftConstEnd()));
}

BigUnsignedInt BigUnsignedInt::suffix(size_t length) const {
    assert(length <= digitCount());
    return BigUnsignedInt(std::vector<size_t>(rightToLeftConstBegin(), rightToLeftConstBegin() + length));
}

BigUnsignedInt BigUnsignedInt::createFromMostSignificantDigit(size_t digit, size_t position) {
    std::vector<size_t> digits(position, 0);
    digits[position] = digit;
    return BigUnsignedInt(std::move(digits));
}

BigUnsignedInt BigUnsignedInt::shiftedCopy(size_t shiftAmount) const {
    auto copy = *this;
    copy.shift(shiftAmount);
    return copy;
}

BigUnsignedInt &BigUnsignedInt::operator%=(const BigUnsignedInt &mod) {
    assert(mod != 0ul);
    BigUnsignedInt copy(*this);
    *this -= longDivision(copy, mod) * mod;
    return *this;
}

BigUnsignedInt BigUnsignedInt::createRandom(size_t numberOfDigits) {
    assert(numberOfDigits > 0ul);
    std::random_device                    rd;        // Will be used to obtain a seed for the random number engine
    std::mt19937                          gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<size_t> dis(0, DigitVector::s_maxDigit);

    BigUnsignedInt result(0ul);
    result.m_digits[0] = dis(gen);
    for (size_t i = 1; i != numberOfDigits; ++i) {
        result.shift(1);
        result.m_digits[0] = dis(gen);
    }
    result.resizeToFit();
    return result;
}

size_t BigUnsignedInt::approximatePowerOfTen() const {
    return std::log10(mostSignificantDigit() + 1ul) + (digitCount() - 1ul) * std::log10(DigitVector::s_base);
}

BigUnsignedInt BigUnsignedInt::createRandomFromDecimalDigits(size_t orderOfMagnitude) {
    assert(orderOfMagnitude > 0);
    const auto numberOfDigits = orderOfMagnitude * (std::log(10) / log(DigitVector::s_base)) + 1ul;
    auto       result         = createRandom(numberOfDigits);
    assert(result.isWellFormed());
    return result;
}
