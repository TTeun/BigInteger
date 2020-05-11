#ifndef TEUN_GAME_BIGUNSIGNEDINT_H
#define TEUN_GAME_BIGUNSIGNEDINT_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <vector>

class BigUnsignedInt {

private:
    static const size_t s_base = std::sqrt(std::numeric_limits<size_t>::max()) - 1ul;

    static constexpr bool baseIsEven()
    {
        return s_base % 2 == 0;
    }

public:
    BigUnsignedInt();

    BigUnsignedInt(size_t val);

    BigUnsignedInt(const std::string & val);

    BigUnsignedInt & operator=(const BigUnsignedInt & rhs);

    BigUnsignedInt & operator=(size_t rhs);

    BigUnsignedInt & operator+=(const BigUnsignedInt & rhs);

    BigUnsignedInt & operator+=(size_t rhs);

    BigUnsignedInt & operator*=(const BigUnsignedInt & rhs);

    BigUnsignedInt & operator*=(size_t rhs);

    BigUnsignedInt operator%(const BigUnsignedInt & mod) const;

    size_t operator%(size_t mod) const;

    BigUnsignedInt & operator%=(size_t mod);

    BigUnsignedInt operator/(size_t dividend) const;

    bool operator==(const BigUnsignedInt & rhs) const;

    bool operator!=(const BigUnsignedInt & rhs) const;

    bool operator<(const BigUnsignedInt & rhs) const;

    bool operator>(const BigUnsignedInt & rhs) const;

    bool operator<=(const BigUnsignedInt & rhs) const;

    bool operator>=(const BigUnsignedInt & rhs) const;

    friend std::ostream & operator<<(std::ostream & os, const BigUnsignedInt & anInt);

    friend BigUnsignedInt power(const BigUnsignedInt & base, size_t exponent);

protected:
    void shiftAdd(const BigUnsignedInt & rhs, size_t shiftAmount);

    void square();

    BigUnsignedInt operator*(const size_t & rhs);

    bool isCorrectlySized() const;

    size_t digitCount() const;

    void resizeToFit();

    void init(size_t val);

    void bubble(size_t startIndex = 0ul);

    BigUnsignedInt & shift(size_t shiftAmount);

    std::vector<size_t> m_digits;
};

#endif // TEUN_GAME_BIGUNSIGNEDINT_H
