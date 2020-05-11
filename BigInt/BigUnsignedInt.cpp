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

void BigUnsignedInt::bubble(size_t startIndex)
{
    assert(not m_digits.empty());
    assert(startIndex < digitCount());
    auto it   = m_digits.begin() + startIndex;
    auto next = it + 1;

    for (; next != m_digits.end(); ++it, ++next) {
        if (*it >= s_base) {
            *next += *it / s_base;
            *it %= s_base;
        }
    }

    if (m_digits.back() >= s_base) {
        size_t continueIndex = digitCount() - 1ul;
        m_digits.resize(digitCount() + std::log(m_digits.back()) / std::log(s_base) + 1ul);
        bubble(continueIndex);
    }
    resizeToFit();
}

std::ostream & operator<<(std::ostream & os, const BigUnsignedInt & bigUnsignedInt)
{
    os << "( ";
    for (auto it = bigUnsignedInt.m_digits.rbegin(); it != bigUnsignedInt.m_digits.rend(); ++it) {
        os << *it << " ";
    }
    os << ")_" << BigUnsignedInt::s_base;
    return os;
}

void BigUnsignedInt::init(size_t val)
{
    static_assert(s_base <= std::numeric_limits<size_t>::max() / s_base,
                  "s_base^2 should not exceed the maximum size_t");
    if (val < s_base) {
        m_digits = {val};
    } else if (val < s_base * s_base) {
        m_digits = {val % s_base, val / s_base};
    } else {
        m_digits = {val % s_base, (val / s_base) % s_base, val / (s_base * s_base)};
    }
}

bool BigUnsignedInt::operator==(const BigUnsignedInt & rhs) const
{
    assert(isCorrectlySized());
    assert(rhs.isCorrectlySized());

    return m_digits == rhs.m_digits;
}

bool BigUnsignedInt::operator!=(const BigUnsignedInt & rhs) const
{
    return !(rhs == *this);
}

void BigUnsignedInt::resizeToFit()
{
    for (auto it = m_digits.crbegin(); it != m_digits.crend(); ++it) {
        if (*it != 0ul) {
            m_digits.resize(digitCount() - std::distance(m_digits.crbegin(), it));
            break;
        }
    }

    if (m_digits.empty()) {
        m_digits.resize(1);
        *m_digits.begin() = 0ul;
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

    auto rhsIt = rhs.m_digits.cbegin();
    for (auto thisIt = m_digits.begin(); rhsIt != rhs.m_digits.end(); ++thisIt, ++rhsIt) {
        *thisIt += *rhsIt;
    }
    bubble(0);
    return *this;
}

bool BigUnsignedInt::isCorrectlySized() const
{
    assert(digitCount() > 0);
    if (m_digits.back() == 0ul) {
        return digitCount() == 1; // Is the value zero
    }

    return true;
}

bool BigUnsignedInt::operator<(const BigUnsignedInt & rhs) const
{
    assert(isCorrectlySized() && rhs.isCorrectlySized());

    if (digitCount() != rhs.digitCount()) {
        return digitCount() < rhs.digitCount();
    }
    if (digitCount() == 1) {
        return m_digits.front() < m_digits.front();
    }

    auto thisIt = m_digits.crbegin();
    auto rhsIt  = rhs.m_digits.crbegin();

    for (; thisIt != m_digits.crend(); ++thisIt, ++rhsIt) {
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
    assert(isCorrectlySized() && rhs.isCorrectlySized());

    if (digitCount() != rhs.digitCount()) {
        return digitCount() < rhs.digitCount();
    }
    if (digitCount() == 1) {
        return m_digits.front() <= m_digits.front();
    }

    auto thisIt = m_digits.crbegin();
    auto rhsIt  = rhs.m_digits.crbegin();

    for (; thisIt != m_digits.crend(); ++thisIt, ++rhsIt) {
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
    auto rhsIt = rhs.m_digits.cbegin();
    for (size_t i = 0; i != rhs.digitCount(); ++i) {
        shiftAdd(copy * *rhsIt, i);
        ++rhsIt;
    }
    return *this;
}

BigUnsignedInt & BigUnsignedInt::shift(size_t shiftAmount)
{
    m_digits.insert(m_digits.begin(), shiftAmount, 0);
    return *this;
}

BigUnsignedInt BigUnsignedInt::operator*(const size_t & rhs)
{
    assert(rhs < s_base);
    auto result = *this;
    for (auto & it : result.m_digits) {
        it *= rhs;
    }
    result.bubble(0);
    return result;
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
    if (rhs <= std::numeric_limits<size_t>::max() - m_digits.front()) {
        *m_digits.begin() += rhs;
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

    auto thisIt = m_digits.begin() + shiftAmount;
    auto rhsIt  = rhs.m_digits.cbegin();

    for (; rhsIt != rhs.m_digits.end(); ++thisIt, ++rhsIt) {
        *thisIt += *rhsIt;
    }
    bubble(0);
}

BigUnsignedInt & BigUnsignedInt::operator=(size_t rhs)
{
    init(rhs);
}

static size_t modPower(size_t base, size_t exponent, size_t mod)
{
    if (exponent == 0ul) {
        return 1;
    } else if (exponent == 1ul) {
        return base % mod;
    }
    size_t aux(1);
    size_t copy(base);
    while (exponent > 1ul) {
        if (exponent % 2ul == 0ul) {
            exponent /= 2ul;
        } else {
            aux *= copy;
            aux %= mod;
            exponent = (exponent - 1) / 2ul;
        }
        copy *= copy;
        copy %= mod;
    }
    copy *= aux;
    return copy % mod;
}

size_t BigUnsignedInt::operator%(size_t mod) const
{
    if (mod >= s_base) {
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

BigUnsignedInt BigUnsignedInt::operator/(size_t mod) const
{
    return 0;
}
