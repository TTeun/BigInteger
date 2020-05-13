#include "BigUnsignedInt.h"

#include <cassert>
#include <iostream>

/* =================== Static functions =================== */

static void subtractViaIterators(std::vector<size_t>::iterator             thisIt,
                                 const std::vector<size_t>::iterator       thisEnd,
                                 std::vector<size_t>::const_iterator       rhsIt,
                                 const std::vector<size_t>::const_iterator rhsEnd) {
    unsigned short carry = 0ul;
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
            *thisIt = DigitVector::s_base - 1ul;
            ++thisIt;
        }
    }
}

static void addViaIterators(std::vector<size_t>::iterator             thisIt,
                            const std::vector<size_t>::iterator       thisEnd,
                            std::vector<size_t>::const_iterator       rhsIt,
                            const std::vector<size_t>::const_iterator rhsEnd) {
    unsigned short carry = 0ul;
    for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
        *thisIt += *rhsIt + carry;
        if (*thisIt >= DigitVector::s_base) {
            carry = 1ul;
            *thisIt %= DigitVector::s_base;
        } else {
            carry = 0ul;
        }
    }
    if (carry == 1ul) {
        while (true) {
            assert(thisIt != thisEnd);
            ++*thisIt;
            if (*thisIt < DigitVector::s_base) {
                break;
            }
            *thisIt = 0ul;
            ++thisIt;
        }
    }
}

static bool lessThanShiftedRhsViaIterators(std::vector<size_t>::const_reverse_iterator        thisIt,
                                           const std::vector<size_t>::const_reverse_iterator &thisEnd,
                                           std::vector<size_t>::const_reverse_iterator        rhsIt,
                                           const std::vector<size_t>::const_reverse_iterator &rhsEnd,
                                           const size_t                                       trailingZeroesOfRhs) {
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

static bool lessThanViaIterators(const std::vector<size_t>::const_reverse_iterator &thisIt,
                                 const std::vector<size_t>::const_reverse_iterator &thisEnd,
                                 const std::vector<size_t>::const_reverse_iterator &rhsIt,
                                 const std::vector<size_t>::const_reverse_iterator &rhsEnd) {
    return lessThanShiftedRhsViaIterators(thisIt, thisEnd, rhsIt, rhsEnd, 0ul);
}

static bool greaterThanViaIterators(const std::vector<size_t>::const_reverse_iterator &thisIt,
                                    const std::vector<size_t>::const_reverse_iterator &thisEnd,
                                    const std::vector<size_t>::const_reverse_iterator &rhsIt,
                                    const std::vector<size_t>::const_reverse_iterator &rhsEnd) {
    return lessThanViaIterators(rhsIt, rhsEnd, thisIt, thisEnd);
}

/* =================== Member functions =================== */

BigUnsignedInt::BigUnsignedInt() {
    init(0);
}

BigUnsignedInt::BigUnsignedInt(size_t val) {
    init(val);
}

BigUnsignedInt::BigUnsignedInt(std::vector<size_t> &&digits) : DigitVector(std::move(digits)) {
    if (m_digits.empty()) {
        init(0);
    }
    bubble();
}

BigUnsignedInt::BigUnsignedInt(const std::string &val) {
    init(0);
    for (char charIt : val) {
        *this *= 10;
        *this += charIt - '0';
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

std::ostream &operator<<(std::ostream &os, const BigUnsignedInt &bigUnsignedInt) {
    os << "( ";
    for (auto it = bigUnsignedInt.leftToRightConstBegin(); it != bigUnsignedInt.leftToRightConstEnd(); ++it) {
        os << *it << " ";
    }
    os << ")_" << BigUnsignedInt::s_base;
    return os;
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
    auto it = leftToRightConstBegin();
    for (; it != leftToRightConstEnd() && *it == 0ul; ++it)
        ;

    m_digits.resize(std::distance(it, leftToRightConstEnd()));

    if (m_digits.empty()) {
        m_digits = {0ul};
    }
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

bool greaterOrEqualViaIterators(const std::vector<size_t>::const_reverse_iterator &thisIt,
                                const std::vector<size_t>::const_reverse_iterator &thisEnd,
                                const std::vector<size_t>::const_reverse_iterator &rhsIt,
                                const std::vector<size_t>::const_reverse_iterator &rhsEnd) {
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

BigUnsignedInt &BigUnsignedInt::operator*=(const BigUnsignedInt &rhs) {
    if (this == &rhs) {
        square();
        return *this;
    }
    BigUnsignedInt copy = *this * rhs.leastSignificantDigit();
    copy.m_digits.reserve(digitCount() + rhs.digitCount());
    auto rhsIt = rhs.rightToLeftConstBegin() + 1ul;
    for (size_t i = 1; i != rhs.digitCount(); ++i) {
        copy.addShiftedMultiplied(*this, i, *rhsIt);
        //        std::cout << copy << '\n';
        ++rhsIt;
    }
    swap(copy, *this);
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
    } else if (rhs != 1ul) {
        if (rhs < s_base) {
            for (auto &it : m_digits) {
                it *= rhs;
            }
            bubble();
        } else {
            *this *= BigUnsignedInt(rhs);
        }
    }
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

BigUnsignedInt power(const BigUnsignedInt &base, size_t exponent) {
    if (exponent == 0ul) {
        return BigUnsignedInt(1);
    } else if (exponent == 1ul) {
        return BigUnsignedInt(base);
    }
    BigUnsignedInt aux(1);
    BigUnsignedInt copy(base);
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

BigUnsignedInt &BigUnsignedInt::operator=(size_t rhs) {
    init(rhs);
    return *this;
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

size_t divisionSubRoutine(const std::vector<size_t>::const_reverse_iterator &lrcb,
                          const std::vector<size_t>::const_reverse_iterator &lrce,
                          const std::vector<size_t>::iterator &              rlb,
                          std::vector<size_t>::iterator                      rle,
                          const BigUnsignedInt &                             divisor) {
    if (lessThanViaIterators(lrcb, lrce, divisor.leftToRightConstBegin(), divisor.leftToRightConstEnd())) {
        return 0ul;
    }
    assert(divisor != 0);
    assert(divisor.mostSignificantDigit() * 2 >= BigUnsignedInt::s_base);
    const size_t n = divisor.digitCount();
    assert(std::distance(lrcb, lrce) <= n + 1);
    assert(std::distance(lrcb, lrce) >= n);
    size_t correction = 0ul;

    while (not lessThanShiftedRhsViaIterators(lrcb, lrce, divisor.leftToRightConstBegin(), divisor.leftToRightConstEnd(), 1)) {
        subtractViaIterators(rlb + 1ul, rle, divisor.rightToLeftConstBegin(), divisor.rightToLeftConstEnd());
        correction += BigUnsignedInt::s_base;
    }

    size_t q;
    if (std::distance(lrcb, lrce) == n) {
        q = (*lrcb) / divisor.mostSignificantDigit();
    } else {
        q = (*lrcb * BigUnsignedInt::s_base + *(lrcb + 1ul)) / divisor.mostSignificantDigit();
    }
    q                   = std::min(q, BigUnsignedInt::s_base - 1ul);
    const size_t offset = *lrcb == 0 ? 1ul : 0ul;

    BigUnsignedInt T = divisor * q;
    while (greaterThanViaIterators(T.leftToRightConstBegin(), T.leftToRightConstEnd(), lrcb + offset, lrce)) {
        --q;
        T -= divisor;
    }
    subtractViaIterators(rlb, rle, T.rightToLeftConstBegin(), T.rightToLeftConstEnd());
    return q + correction;
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
    return longDivision(*this, divisor).first;
}

BigUnsignedInt BigUnsignedInt::operator*(const BigUnsignedInt &rhs) const {
    auto copy = *this;
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

static size_t ceilingIntegerDivision(size_t a, size_t b) {
    return (a + b - 1) / b;
}

std::pair<BigUnsignedInt, BigUnsignedInt> longDivision(const BigUnsignedInt &dividend, const BigUnsignedInt &divisor) {
    if (dividend < divisor) {
        return {0, dividend};
    }
    const size_t t = ceilingIntegerDivision(BigUnsignedInt::s_base, 2ul * divisor.mostSignificantDigit());
    if (t > 1) {
        return longDivisionAfterAdjustingDivisor(dividend * t, divisor * t);
    } else {
        return longDivisionAfterAdjustingDivisor(dividend, divisor);
    }
}

std::pair<BigUnsignedInt, BigUnsignedInt> longDivisionAfterAdjustingDivisor(BigUnsignedInt dividend, const BigUnsignedInt &divisor) {
    assert(divisor <= dividend);
    assert(divisor.mostSignificantDigit() * 2ul >= BigUnsignedInt::s_base);

    size_t       m = dividend.digitCount();
    const size_t n = divisor.digitCount();

    if (m <= n + 1) {
        return {divisionSubRoutine(dividend.leftToRightConstBegin(),
                                   dividend.leftToRightConstEnd(),
                                   dividend.rightToLeftBegin(),
                                   dividend.rightToLeftEnd(),
                                   divisor),
                dividend};
    }

    std::vector<size_t> divisorDigits(m - n, 0ul);
    while (m > n + 1) {
        const size_t splitIndex = m - n - 1ul;

        auto prefix = dividend.prefix(n + 1);

        auto quotientRemainder = divisionSubRoutine(
            prefix.leftToRightConstBegin(), prefix.leftToRightConstEnd(), prefix.rightToLeftBegin(), prefix.rightToLeftEnd(), divisor);
        divisorDigits[splitIndex] = quotientRemainder;
        prefix.resizeToFit();
        //        prefix -= quotientRemainder.second;

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

    return {BigUnsignedInt(std::move(divisorDigits)), 0ul};
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

void swap(BigUnsignedInt &a, BigUnsignedInt &b) {
    std::swap(a.m_digits, b.m_digits);
}

BigUnsignedInt BigUnsignedInt::shiftedCopy(size_t shiftAmount) const {
    auto copy = *this;
    copy.shift(shiftAmount);
    return copy;
}

void BigUnsignedInt::addShiftedMultiplied(const BigUnsignedInt &rhs, size_t shiftAmount, size_t multiplier) {
    assert(multiplier < s_base);
    assert(multiplier * (s_base - 1ul) <= std::numeric_limits<size_t>::max() - s_base);

    m_digits.resize(std::max(rhs.digitCount() + shiftAmount, digitCount()) + 1);

    auto thisIt = rightToLeftBegin() + shiftAmount;
    auto rhsIt  = rhs.rightToLeftConstBegin();

    size_t carry = 0ul;
    for (; rhsIt != rhs.rightToLeftConstEnd(); ++thisIt, ++rhsIt) {
        *thisIt += (*rhsIt * multiplier) + carry;
        if (*thisIt >= s_base) {
            carry = *thisIt / s_base;
            *thisIt %= s_base;
        } else {
            carry = 0ul;
        }
    }
    if (carry > 0ul) {
        *thisIt += carry;
    }
    if (m_digits.back() == 0ul) {
        m_digits.resize(m_digits.size() - 1);
    }
    bubble(0);
}

size_t BigUnsignedInt::twoPrefix() const {
    return s_base * mostSignificantDigit() + secondMostSignificantDigit();
}
