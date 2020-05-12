#include "BigUnsignedInt.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

BigUnsignedInt::BigUnsignedInt()
{
    init(0);
}

BigUnsignedInt::BigUnsignedInt(size_t val)
{
    init(val);
}

BigUnsignedInt::BigUnsignedInt(std::vector<size_t> && digits) : m_digits(std::move(digits))
{
    if (m_digits.empty()) {
        init(0);
    }
    bubble();
}

void BigUnsignedInt::bubble(size_t startIndex)
{
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
        size_t continueIndex = digitCount() - 1ul;
        m_digits.resize(digitCount() + std::log(mostSignificantDigit()) / std::log(s_base) + 1ul);
        bubble(continueIndex);
    }
    resizeToFit();
}

std::ostream & operator<<(std::ostream & os, const BigUnsignedInt & bigUnsignedInt)
{
    os << "( ";
    for (auto it = bigUnsignedInt.leftToRightConstBegin(); it != bigUnsignedInt.leftToRightConstEnd(); ++it) {
        os << *it << " ";
    }
    os << ")_" << BigUnsignedInt::s_base;
    return os;
}

void BigUnsignedInt::init(size_t val)
{
    static_assert(s_base <= std::numeric_limits<size_t>::max() / s_base,
                  "s_base^2 should not exceed the maximum size_t");
    m_digits = {val};
    bubble(0);
}

bool BigUnsignedInt::operator==(const BigUnsignedInt & rhs) const
{
    assert(isWellFormed());
    assert(rhs.isWellFormed());

    return m_digits == rhs.m_digits;
}

bool BigUnsignedInt::operator!=(const BigUnsignedInt & rhs) const
{
    return !(rhs == *this);
}

void BigUnsignedInt::resizeToFit()
{
    auto it = leftToRightConstBegin();
    for (; it != leftToRightConstEnd(); ++it) {
        if (*it != 0ul) {
            break;
        }
    }
    m_digits.resize(digitCount() - std::distance(leftToRightConstBegin(), it));

    if (m_digits.empty()) {
        m_digits.resize(1);
        *rightToLeftBegin() = 0ul;
    }
}

BigUnsignedInt & BigUnsignedInt::operator=(const BigUnsignedInt & rhs)
{
    if (this == &rhs) {
        return *this;
    }
    m_digits = rhs.m_digits;
    return *this;
}

BigUnsignedInt & BigUnsignedInt::operator+=(const BigUnsignedInt & rhs)
{
    if (rhs.digitCount() > digitCount()) {
        m_digits.resize(rhs.digitCount());
    }

    auto rhsIt  = rhs.rightToLeftConstBegin();
    auto thisIt = rightToLeftBegin();
    for (; rhsIt != rhs.rightToLeftConstEnd(); ++thisIt, ++rhsIt) {
        *thisIt += *rhsIt;
    }
    bubble(0);
    return *this;
}

bool BigUnsignedInt::isCorrectlySized() const
{
    assert(digitCount() > 0);
    if (mostSignificantDigit() == 0ul) {
        return digitCount() == 1; // Is the value zero
    }

    return true;
}

bool BigUnsignedInt::operator<(const BigUnsignedInt & rhs) const
{
    assert(isWellFormed() && rhs.isCorrectlySized());

    if (digitCount() != rhs.digitCount()) {
        return digitCount() < rhs.digitCount();
    }
    if (digitCount() == 1) {
        return m_digits.front() < m_digits.front();
    }

    auto thisIt = leftToRightConstBegin();
    auto rhsIt  = rhs.leftToRightConstBegin();

    for (; thisIt != leftToRightConstEnd(); ++thisIt, ++rhsIt) {
        if (*thisIt != *rhsIt) {
            return *thisIt < *rhsIt;
        }
    }
    return false;
}

bool BigUnsignedInt::operator>(const BigUnsignedInt & rhs) const
{
    return rhs < *this;
}

bool BigUnsignedInt::operator<=(const BigUnsignedInt & rhs) const
{
    assert(isWellFormed() && rhs.isCorrectlySized());

    if (digitCount() != rhs.digitCount()) {
        return digitCount() < rhs.digitCount();
    }
    if (digitCount() == 1) {
        return m_digits.front() <= m_digits.front();
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

bool BigUnsignedInt::operator>=(const BigUnsignedInt & rhs) const
{
    return !(*this < rhs);
}

size_t BigUnsignedInt::digitCount() const
{
    return m_digits.size();
}

BigUnsignedInt & BigUnsignedInt::operator*=(const BigUnsignedInt & rhs)
{
    if (this == &rhs) {
        square();
        return *this;
    }

    BigUnsignedInt copy = 0ul;
    std::swap(copy.m_digits, m_digits);
    auto rhsIt = rhs.rightToLeftConstBegin();
    for (size_t i = 0; i != rhs.digitCount(); ++i) {
        shiftAdd(copy * (*rhsIt), i);
        ++rhsIt;
    }
    return *this;
}

BigUnsignedInt & BigUnsignedInt::shift(size_t shiftAmount)
{
    m_digits.insert(rightToLeftBegin(), shiftAmount, 0);
    return *this;
}

BigUnsignedInt BigUnsignedInt::operator*(size_t rhs) const
{
    if (rhs == 0) {
        return BigUnsignedInt(0);
    }
    auto copy = *this;
    copy *= rhs;
    return copy;
}

void BigUnsignedInt::square()
{
    const auto copy = *this;
    *this *= copy;
    bubble(0);
}

BigUnsignedInt::BigUnsignedInt(const std::string & val)
{
    init(0);
    for (auto charIt = val.cbegin(); charIt != val.cend(); ++charIt) {
        *this *= 10;
        *this += *charIt - '0';
    }
    bubble(0);
}

BigUnsignedInt & BigUnsignedInt::operator*=(const size_t rhs)
{
    if (rhs < s_base) {
        for (auto & it : m_digits) {
            it = it * rhs;
        }
        bubble();
        return *this;
    } else {
        return *this *= BigUnsignedInt(rhs);
    }
}

BigUnsignedInt & BigUnsignedInt::operator+=(size_t rhs)
{
    static const size_t additionRoom = std::numeric_limits<size_t>::max() - s_base;
    if (rhs <= additionRoom) {
        *rightToLeftBegin() += rhs;
        bubble();
        return *this;
    } else {
        return *this += BigUnsignedInt(rhs);
    }
}

BigUnsignedInt power(const BigUnsignedInt & base, size_t exponent)
{
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

void BigUnsignedInt::shiftAdd(const BigUnsignedInt & rhs, size_t shiftAmount)
{
    if (rhs.digitCount() + shiftAmount > digitCount()) {
        m_digits.resize(rhs.digitCount() + shiftAmount);
    }

    auto thisIt = rightToLeftBegin() + shiftAmount;
    auto rhsIt  = rhs.rightToLeftConstBegin();

    for (; rhsIt != rhs.m_digits.end(); ++thisIt, ++rhsIt) {
        *thisIt += *rhsIt;
    }
    bubble(0);
}

BigUnsignedInt & BigUnsignedInt::operator=(size_t rhs)
{
    init(rhs);
    return *this;
}

template <typename T>
static size_t modPower(T base, T exponent, size_t mod)
{
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

size_t BigUnsignedInt::operator%(size_t mod) const
{
    if (mod >= std::numeric_limits<size_t>::max() / mod) {
        std::cout << "large modulo moet nog\n";
        return 0;
    }
    size_t result = 0ul;
    for (size_t i = 0; i != digitCount(); ++i) {
        result += (m_digits.at(i) % mod) * (modPower(s_base, i, mod) % mod);
        result %= mod;
    }
    return result;
}

BigUnsignedInt & BigUnsignedInt::operator%=(size_t mod)
{
    *this = *this % mod;
    return *this;
}

std::pair<size_t, BigUnsignedInt> divisionSubRoutine(const BigUnsignedInt & dividend, const BigUnsignedInt & divisor)
{
    assert(divisor != 0);
    assert(divisor.mostSignificantDigit() * 2 >= BigUnsignedInt::s_base);
    const size_t n = divisor.digitCount();
    assert(dividend.digitCount() <= n + 1);
    if (dividend >= divisor * BigUnsignedInt::s_base) {
        auto result = divisionSubRoutine(dividend - divisor * BigUnsignedInt::s_base, divisor);
        result.first += BigUnsignedInt::s_base;
        return result;
    }
    size_t q;
    if (dividend.digitCount() == n + 1) {
        q = (BigUnsignedInt::s_base * dividend.mostSignificantDigit() + dividend.secondMostSignificantDigit()) /
            divisor.mostSignificantDigit();
    } else {
        q = divisor.mostSignificantDigit() / divisor.mostSignificantDigit();
    }
    q                = std::min(q, BigUnsignedInt::s_base - 1ul);
    BigUnsignedInt T = divisor * q;
    while (T > dividend) {
        --q;
        T -= divisor;
    }
    return {q, dividend - T};
}

BigUnsignedInt BigUnsignedInt::operator/(const BigUnsignedInt & divisor) const
{
    if (&divisor == this) {
        return BigUnsignedInt(0);
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

size_t BigUnsignedInt::mostSignificantDigit() const
{
    return m_digits.back();
}

BigUnsignedInt BigUnsignedInt::operator*(const BigUnsignedInt & rhs) const
{
    auto copy = *this;
    copy *= rhs;
    return copy;
}

BigUnsignedInt & BigUnsignedInt::operator-=(const BigUnsignedInt & rhs)
{
    assert(rhs <= *this);
    auto rhsIt           = rhs.rightToLeftConstBegin();
    auto thisIt          = rightToLeftBegin();
    unsigned short carry = 0ul;
    for (; rhsIt != rhs.rightToLeftConstEnd(); ++thisIt, ++rhsIt) {
        if (*thisIt >= *rhsIt + carry) {
            *thisIt -= *rhsIt + carry;
            carry = 0;
        } else {
            *thisIt += s_base;
            *thisIt -= *rhsIt + carry;
            carry = 1ul;
        }
    }
    while (carry != 0ul) {
        assert(thisIt != rightToLeftEnd());
        if (*thisIt > 0) {
            --*thisIt;
            break;
        } else {
            *thisIt = s_base - 1ul;
        }
        ++thisIt;
    }
    resizeToFit();
    return *this;
}

BigUnsignedInt BigUnsignedInt::operator-(const BigUnsignedInt & rhs) const
{
    auto copy = *this;
    copy -= rhs;
    return copy;
}

BigUnsignedInt BigUnsignedInt::operator+(size_t rhs) const
{
    auto copy = *this;
    copy += rhs;
    return copy;
}

BigUnsignedInt BigUnsignedInt::operator+(const BigUnsignedInt & rhs) const
{
    auto copy = *this;
    copy += rhs;
    return copy;
}

size_t BigUnsignedInt::secondMostSignificantDigit() const
{
    assert(digitCount() > 1);
    return *(m_digits.rbegin() + 1);
}

static size_t ceilingIntegerDivision(size_t a, size_t b)
{
    return (a + b - 1) / b;
}

std::pair<BigUnsignedInt, BigUnsignedInt> longDivision(const BigUnsignedInt & dividend, const BigUnsignedInt & divisor)
{
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

std::pair<BigUnsignedInt, BigUnsignedInt> longDivisionAfterAdjustingDivisor(const BigUnsignedInt & dividend,
                                                                            const BigUnsignedInt & divisor)
{
    assert(divisor <= dividend);
    assert(divisor.mostSignificantDigit() * 2 > BigUnsignedInt::s_base);

    const size_t m = dividend.digitCount();
    const size_t n = divisor.digitCount();

    if (m <= n + 1) {
        return divisionSubRoutine(dividend, divisor);
    }
    const auto prefix = dividend.prefix(n + 1);
    const auto suffix = dividend.suffix(m - n - 1ul);

    auto quotientRemainder = divisionSubRoutine(prefix, divisor);

    auto quotient = BigUnsignedInt(quotientRemainder.first).shift(m - n - 1);

    const auto rec = longDivisionAfterAdjustingDivisor(suffix + quotientRemainder.second.shift(m - n - 1), divisor);
    quotient += rec.first;
    return {quotient, rec.second};
}

BigUnsignedInt BigUnsignedInt::prefix(size_t length) const
{
    assert(length <= m_digits.size());
    return BigUnsignedInt(std::vector<size_t>(m_digits.end() - length, m_digits.end()));
}

BigUnsignedInt BigUnsignedInt::suffix(size_t length) const
{
    assert(length <= m_digits.size());
    return BigUnsignedInt(std::vector<size_t>(rightToLeftConstBegin(), rightToLeftConstBegin() + length));
}

BigUnsignedInt BigUnsignedInt::createFromMostSignificantDigit(size_t digit, size_t position)
{
    std::vector<size_t> digits(position, 0);
    digits[position] = digit;
    return BigUnsignedInt(std::move(digits));
}

std::vector<size_t>::reverse_iterator BigUnsignedInt::leftToRightBegin()
{
    return m_digits.rbegin();
}
std::vector<size_t>::reverse_iterator BigUnsignedInt::leftToRightEnd()
{
    return m_digits.rend();
}
std::vector<size_t>::iterator BigUnsignedInt::rightToLeftBegin()
{
    return m_digits.begin();
}
std::vector<size_t>::iterator BigUnsignedInt::rightToLeftEnd()
{
    return m_digits.end();
}
std::vector<size_t>::const_reverse_iterator BigUnsignedInt::leftToRightConstBegin() const
{
    return m_digits.crbegin();
}
std::vector<size_t>::const_reverse_iterator BigUnsignedInt::leftToRightConstEnd() const
{
    return m_digits.crend();
}
std::vector<size_t>::const_iterator BigUnsignedInt::rightToLeftConstBegin() const
{
    return m_digits.cbegin();
}
std::vector<size_t>::const_iterator BigUnsignedInt::rightToLeftConstEnd() const
{
    return m_digits.cend();
}
bool BigUnsignedInt::isWellFormed() const
{
    if (not isCorrectlySized()) {
        return false;
    }
    for (const auto it : m_digits) {
        if (it >= s_base) {
            return false;
        }
    }
    return true;
}
