#ifndef TEUN_GAME_BIGUNSIGNEDINT_H
#define TEUN_GAME_BIGUNSIGNEDINT_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <vector>

class BigUnsignedInt {

    // ToDo static const additionRoom = size_t::max() - s_base;
    // subtraction using carry

public:
    //    static const size_t s_base = std::sqrt(std::numeric_limits<size_t>::max()) - 1ul;
    static const size_t s_base = 10ul;

public:
    BigUnsignedInt();

    explicit BigUnsignedInt(std::vector<size_t> && digits);

    BigUnsignedInt(size_t val);

    BigUnsignedInt(const std::string & val);

    BigUnsignedInt & operator=(const BigUnsignedInt & rhs);

    BigUnsignedInt & operator=(size_t rhs);

    BigUnsignedInt & operator+=(const BigUnsignedInt & rhs);

    BigUnsignedInt & operator-=(const BigUnsignedInt & rhs);

    BigUnsignedInt operator-(const BigUnsignedInt & rhs) const;

    BigUnsignedInt & operator+=(size_t rhs);

    BigUnsignedInt & operator*=(const BigUnsignedInt & rhs);

    BigUnsignedInt & operator*=(size_t rhs);

    BigUnsignedInt operator+(size_t rhs) const;

    BigUnsignedInt operator+(const BigUnsignedInt & rhs) const;

    BigUnsignedInt operator*(size_t rhs) const;

    BigUnsignedInt operator*(const BigUnsignedInt & rhs) const;

    size_t operator%(size_t mod) const;

    BigUnsignedInt & operator%=(size_t mod);

    BigUnsignedInt operator/(const BigUnsignedInt & divisor) const;

    bool operator==(const BigUnsignedInt & rhs) const;

    bool operator!=(const BigUnsignedInt & rhs) const;

    bool operator<(const BigUnsignedInt & rhs) const;

    bool operator>(const BigUnsignedInt & rhs) const;

    bool operator<=(const BigUnsignedInt & rhs) const;

    bool operator>=(const BigUnsignedInt & rhs) const;

    size_t mostSignificantDigit() const;

    size_t secondMostSignificantDigit() const;

    friend std::ostream & operator<<(std::ostream & os, const BigUnsignedInt & anInt);

    friend BigUnsignedInt power(const BigUnsignedInt & base, size_t exponent);

    size_t digitCount() const;

    BigUnsignedInt prefix(size_t length) const;

    BigUnsignedInt suffix(size_t length) const;

protected:
    friend std::pair<size_t, BigUnsignedInt> divisionSubRoutine(const BigUnsignedInt & dividend, const BigUnsignedInt & divisor);

    friend std::pair<BigUnsignedInt, BigUnsignedInt> longDivision(const BigUnsignedInt & dividend, const BigUnsignedInt & divisor);

    friend std::pair<BigUnsignedInt, BigUnsignedInt> longDivisionAfterAdjustingDivisor(const BigUnsignedInt & dividend,
                                                            const BigUnsignedInt & divisor);

    bool isWellFormed() const;

    void shiftAdd(const BigUnsignedInt & rhs, size_t shiftAmount);

    void square();

    bool isCorrectlySized() const;

    void resizeToFit();

    void init(size_t val);

    static BigUnsignedInt createFromMostSignificantDigit(size_t digit, size_t position);

    void bubble(size_t startIndex = 0ul);

    BigUnsignedInt & shift(size_t shiftAmount);

    std::vector<size_t>::reverse_iterator leftToRightBegin();

    std::vector<size_t>::reverse_iterator leftToRightEnd();

    std::vector<size_t>::iterator rightToLeftBegin();

    std::vector<size_t>::iterator rightToLeftEnd();

    std::vector<size_t>::const_reverse_iterator leftToRightConstBegin() const;

    std::vector<size_t>::const_reverse_iterator leftToRightConstEnd() const;

    std::vector<size_t>::const_iterator rightToLeftConstBegin() const;

    std::vector<size_t>::const_iterator rightToLeftConstEnd() const;

    std::vector<size_t> m_digits;
};

#endif // TEUN_GAME_BIGUNSIGNEDINT_H
