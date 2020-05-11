//
// Created by pc on 5/11/20.
//

#include "BigInt.h"

BigInt::BigInt()
{
}

BigInt::BigInt(size_t val) : m_magnitude(val)
{
}

BigInt::BigInt(long long int val) : m_magnitude(std::abs(val))
{
    if (val < 0) {
        m_isPositive = false;
    }
}

BigInt::BigInt(const std::string & val)
{
    if (val.front() == '-') {
        m_isPositive = false;
        m_magnitude  = BigUnsignedInt(val.substr(1));
    } else {
        m_magnitude = BigUnsignedInt(val);
    }
}

bool BigInt::operator==(const BigInt & rhs) const
{
    return m_magnitude == rhs.m_magnitude && m_isPositive == rhs.m_isPositive;
}
bool BigInt::operator!=(const BigInt & rhs) const
{
    return !(rhs == *this);
}

BigInt::BigInt(const BigUnsignedInt & unsignedInt) : m_magnitude(unsignedInt)
{
}

BigInt::BigInt(BigUnsignedInt && unsignedInt) : m_magnitude(std::move(unsignedInt))
{
}
