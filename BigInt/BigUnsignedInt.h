#ifndef TEUN_GAME_BIGUNSIGNEDINT_H
#define TEUN_GAME_BIGUNSIGNEDINT_H

#include "DigitVector.h"

#include <cmath>
#include <ostream>

class BigUnsignedInt : private DigitVector {

public:
    BigUnsignedInt();

    explicit BigUnsignedInt(std::vector<size_t> &&digits);

    BigUnsignedInt(size_t val);

    BigUnsignedInt(const std::string &val);

    BigUnsignedInt &operator=(const BigUnsignedInt &rhs);

    BigUnsignedInt &operator=(size_t rhs);

    BigUnsignedInt &operator+=(const BigUnsignedInt &rhs);

    BigUnsignedInt &operator-=(const BigUnsignedInt &rhs);

    BigUnsignedInt operator-(const BigUnsignedInt &rhs) const;

    BigUnsignedInt &operator+=(size_t rhs);

    BigUnsignedInt &operator*=(const BigUnsignedInt &rhs);

    BigUnsignedInt &operator*=(size_t rhs);

    BigUnsignedInt operator+(size_t rhs) const;

    BigUnsignedInt operator+(const BigUnsignedInt &rhs) const;

    BigUnsignedInt operator*(size_t rhs) const;

    BigUnsignedInt operator*(const BigUnsignedInt &rhs) const;

    size_t operator%(size_t mod) const;

    BigUnsignedInt &operator%=(size_t mod);

    BigUnsignedInt operator/(const BigUnsignedInt &divisor) const;

    bool operator==(const BigUnsignedInt &rhs) const;

    bool operator!=(const BigUnsignedInt &rhs) const;

    bool operator<(const BigUnsignedInt &rhs) const;

    bool operator>(const BigUnsignedInt &rhs) const;

    bool operator<=(const BigUnsignedInt &rhs) const;

    bool operator>=(const BigUnsignedInt &rhs) const;

    friend std::ostream &operator<<(std::ostream &        os,
                                    const BigUnsignedInt &anInt);

    friend BigUnsignedInt power(const BigUnsignedInt &base, size_t exponent);

    BigUnsignedInt prefix(size_t length) const;

    size_t twoPrefix() const;

    BigUnsignedInt suffix(size_t length) const;

protected:
    friend void swap(BigUnsignedInt &a, BigUnsignedInt &b);

    friend std::pair<size_t, BigUnsignedInt>
    divisionSubRoutine(BigUnsignedInt dividend,
                       const BigUnsignedInt &divisor);

    friend std::pair<BigUnsignedInt, BigUnsignedInt>
    longDivision(const BigUnsignedInt &dividend, const BigUnsignedInt &divisor);

    friend std::pair<BigUnsignedInt, BigUnsignedInt>
    longDivisionAfterAdjustingDivisor(BigUnsignedInt        dividend,
                                      const BigUnsignedInt &divisor);

    void addShifted(const BigUnsignedInt &rhs, size_t shiftAmount);

    void addShiftedMultiplied(const BigUnsignedInt &rhs, size_t shiftAmount,
                              size_t multiplier);

    void square();

    void resizeToFit();

    void init(size_t val);

    static BigUnsignedInt createFromMostSignificantDigit(size_t digit,
                                                         size_t position);

    void bubble(size_t startIndex = 0ul);

    //    void shiftedCopy(size_t shiftAmount);

    BigUnsignedInt shiftedCopy(size_t shiftAmount) const;
};

#endif // TEUN_GAME_BIGUNSIGNEDINT_H
