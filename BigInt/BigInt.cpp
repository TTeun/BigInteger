#include "BigInt.h"

#include <cassert>

BigInt::BigInt() : m_magnitude() {
}

BigInt::BigInt(BigUInt &&magnitude) : m_magnitude(std::move(magnitude)) {
}

BigInt &BigInt::operator=(size_t rhs) {
    m_magnitude  = rhs;
    m_isNegative = false;
    return *this;
}

BigInt::BigInt(const BigUInt &other) : m_magnitude(other) {
}

BigInt::BigInt(BigInt &&other) noexcept : m_magnitude(other.m_magnitude), m_isNegative(other.m_isNegative) {
}

BigInt::BigInt(const BigInt &other) : m_magnitude(other.m_magnitude), m_isNegative(other.m_isNegative) {
}

BigInt::BigInt(const std::string &val) {
    if (val.front() == '-') {
        m_isNegative = true;
        m_magnitude  = BigUInt(val.substr(1));
    } else {
        m_magnitude = BigUInt(val);
    }
}

BigInt &BigInt::operator+=(const BigInt &rhs) {
    if ((!m_isNegative && !rhs.m_isNegative) || (m_isNegative && rhs.m_isNegative)) {
        m_magnitude += rhs.m_magnitude;
    } else if (rhs.m_isNegative) {
        if (rhs.m_magnitude <= m_magnitude) {
            m_magnitude -= rhs.m_magnitude;
        } else {
            m_isNegative = true;
            m_magnitude  = rhs.m_magnitude - m_magnitude;
        }
    } else {
        assert(m_isNegative);
        if (m_magnitude > rhs.m_magnitude) {
            m_magnitude -= rhs.m_magnitude;
        } else {
            m_isNegative = false;
            m_magnitude  = rhs.m_magnitude - m_magnitude;
        }
    }

    return *this;
}

BigInt &BigInt::operator+=(const BigUInt &rhs) {
    if (!m_isNegative) {
        m_magnitude += rhs;
    } else {
        if (m_magnitude > rhs) {
            m_magnitude -= rhs;
        } else {
            m_isNegative = false;
            m_magnitude  = rhs - m_magnitude;
        }
    }
    return *this;
}

BigInt &BigInt::operator-=(const BigInt &rhs) {
    BigInt copy       = rhs;
    copy.m_isNegative = !rhs.m_isNegative;
    return *this += copy;
}

BigInt &BigInt::operator-=(const BigUInt &rhs) {
    BigInt copy;
    copy.m_magnitude  = rhs;
    copy.m_isNegative = true;
    return *this += copy;
}

BigInt &BigInt::operator+=(size_t rhs) {
    // ToDo optimize
    BigInt copy;
    copy.m_magnitude = rhs;
    return *this += copy;
}

BigInt &BigInt::operator=(const BigUInt &rhs) {
    m_magnitude = rhs;
    return *this;
}

BigInt &BigInt::operator=(BigInt &&rhs) noexcept {
    std::swap(m_magnitude.m_digits, rhs.m_magnitude.m_digits);
    m_isNegative = rhs.m_isNegative;
    return *this;
}

BigInt &BigInt::operator=(BigUInt &&rhs) {
    std::swap(m_magnitude.m_digits, rhs.m_digits);
    return *this;
}

std::ostream &operator<<(std::ostream &os, const BigInt &anInt) {
    os << (anInt.m_isNegative ? "-" : "") << "" << anInt.m_magnitude;
    return os;
}

BigInt &BigInt::operator=(const BigInt &rhs) {
    m_magnitude  = rhs.m_magnitude;
    m_isNegative = rhs.m_isNegative;
    return *this;
}

BigInt BigInt::operator+(size_t rhs) const {
    auto result = *this;
    return result += rhs;
}

BigInt BigInt::operator-(const BigInt &rhs) const {
    auto result = *this;
    return result -= rhs;
}

BigInt BigInt::operator-(const BigUInt &rhs) const {
    auto result = *this;
    return result -= rhs;
}

BigInt &BigInt::operator*=(const BigInt &rhs) {
    m_magnitude *= rhs.m_magnitude;
    m_isNegative = m_isNegative ^ rhs.m_isNegative;
    return *this;
}

BigInt &BigInt::operator*=(size_t rhs) {
    m_magnitude *= rhs;
    return *this;
}

BigInt BigInt::operator+(const BigInt &rhs) const {
    auto result = *this;
    return result += rhs;
}

BigInt BigInt::operator*(size_t rhs) const {
    auto result = *this;
    return result *= rhs;
}

BigInt BigInt::operator*(const BigInt &rhs) const {
    auto result = *this;
    return result *= rhs;
}

BigInt BigInt::operator/(const BigInt &divisor) const {
    auto result = *this;
    return result /= divisor;
}

BigInt BigInt::operator/(const BigUInt &divisor) const {
    auto result = *this;
    return result /= divisor;
}

BigInt BigInt::operator/(size_t divisor) const {
    auto result = *this;
    return result /= divisor;
}

BigInt &BigInt::operator/=(const BigInt &divisor) {
    m_magnitude  = m_magnitude / divisor.m_magnitude;
    m_isNegative = m_isNegative ^ divisor.m_isNegative;
    return *this;
}

BigInt &BigInt::operator/=(const BigUInt &divisor) {
    m_magnitude = m_magnitude / divisor;
    return *this;
}

BigInt &BigInt::operator/=(size_t divisor) {
    m_magnitude = m_magnitude / divisor;
    return *this;
}

void BigInt::resize(size_t size) {
    m_magnitude.resize(size);
}
void BigInt::reserve(size_t size) {
    m_magnitude.reserve(size);
}
BigInt operator+(size_t lhs, const BigInt &rhs) {
    return rhs + lhs;
}
BigInt operator-(size_t lhs, const BigInt &rhs) {
    auto result = rhs - lhs;
    result.flipSign();
    return result;
}
BigInt operator*(size_t lhs, const BigInt &rhs) {
    return rhs * lhs;
}
