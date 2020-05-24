#ifndef __BIG_INT__H__
#define __BIG_INT__H__

#include "BigUInt.h"

#include <ostream>

class BigInt {

public:
    BigInt();

    BigInt(BigInt &&other) noexcept;

    BigInt(const BigInt &other);

    BigInt(BigUInt &&magnitude);

    explicit BigInt(const BigUInt &other);

    explicit BigInt(const std::string &val);

    BigInt &operator=(const BigUInt &rhs);

    BigInt &operator=(const BigInt &rhs);

    BigInt &operator=(BigInt &&rhs) noexcept;

    BigInt &operator=(BigUInt &&rhs);

    friend std::ostream &operator<<(std::ostream &os, const BigInt &anInt);

    BigInt &operator=(size_t rhs);

    BigInt &operator+=(const BigInt &rhs);

    BigInt &operator+=(const BigUInt &rhs);

    BigInt &operator-=(const BigInt &rhs);

    BigInt &operator-=(const BigUInt &rhs);

    BigInt operator-(const BigInt &rhs) const;

    BigInt operator-(const BigUInt &rhs) const;

    BigInt &operator+=(size_t rhs);

    BigInt &operator*=(const BigInt &rhs);

    BigInt &operator/=(const BigInt &divisor);

    BigInt &operator/=(const BigUInt &divisor);

    BigInt &operator/=(size_t divisor);

    BigInt &operator*=(size_t rhs);

    BigInt operator+(size_t rhs) const;

    BigInt operator+(const BigInt &rhs) const;

    BigInt operator*(size_t rhs) const;

    BigInt operator*(const BigInt &rhs) const;

    BigInt operator/(const BigInt &divisor) const;

    BigInt operator/(const BigUInt &divisor) const;

    BigInt operator/(size_t divisor) const;

    friend BigInt operator+(size_t lhs, const BigInt &rhs);

    friend BigInt operator-(size_t lhs, const BigInt &rhs);

    friend BigInt operator*(size_t lhs, const BigInt &rhs);

    void resize(size_t size);

    void reserve(size_t size);

    friend BigUInt;

    const BigUInt &magnitude() const {
        return m_magnitude;
    }

private:
    void negate() {
        m_isNegative = !m_isNegative;
    }

    BigUInt &magnitude() {
        return m_magnitude;
    }

    BigUInt m_magnitude;

    bool m_isNegative = false;
};

#endif // __BIG_INT__H__
