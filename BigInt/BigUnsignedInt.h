#ifndef TEUN_GAME_BIGUNSIGNEDINT_H
#define TEUN_GAME_BIGUNSIGNEDINT_H

#include "DigitVector.h"

#include <ostream>

class BigUnsignedInt : private DigitVector {

public:
    BigUnsignedInt();

    BigUnsignedInt(const BigUnsignedInt &other);

    explicit BigUnsignedInt(std::vector<size_t> &&digits);

    BigUnsignedInt(size_t val);

    BigUnsignedInt(const std::string &val);

    static BigUnsignedInt createRandom(size_t numberOfDigits);

    static BigUnsignedInt createRandomFromDecimalDigits(size_t orderOfMagnitude);

    size_t approximatePowerOfTen() const;

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

    BigUnsignedInt &operator%=(const BigUnsignedInt &mod);

    BigUnsignedInt operator/(const BigUnsignedInt &divisor) const;

    bool operator==(const BigUnsignedInt &rhs) const;

    bool operator!=(const BigUnsignedInt &rhs) const;

    bool operator<(const BigUnsignedInt &rhs) const;

    bool operator>(const BigUnsignedInt &rhs) const;

    bool operator<=(const BigUnsignedInt &rhs) const;

    bool operator>=(const BigUnsignedInt &rhs) const;

    friend std::ostream &operator<<(std::ostream &os, const BigUnsignedInt &anInt);

    friend BigUnsignedInt power(const BigUnsignedInt &base, size_t exponent);

    BigUnsignedInt prefix(size_t length) const;

    BigUnsignedInt suffix(size_t length) const;

private:
    friend void swap(BigUnsignedInt &a, BigUnsignedInt &b);

    friend size_t divisionSubRoutine(const std::vector<size_t>::const_reverse_iterator &leftToRightConstIt,
                                     const std::vector<size_t>::const_reverse_iterator &leftToRightConstEnd,
                                     const std::vector<size_t>::iterator &              rightToLeftIt,
                                     const std::vector<size_t>::iterator &              rightToLeftEnd,
                                     const BigUnsignedInt &                             divisor);

    friend BigUnsignedInt longDivision(BigUnsignedInt &dividend, const BigUnsignedInt &divisor);

    friend BigUnsignedInt longDivisionAfterAdjustingDivisor(BigUnsignedInt &dividend, const BigUnsignedInt &divisor);

    void square();

    void resizeToFit();

    void init(size_t val);

    static BigUnsignedInt createFromMostSignificantDigit(size_t digit, size_t position);

    void bubble(size_t startIndex = 0ul);

    //    void shiftedCopy(size_t shiftAmount);

    BigUnsignedInt shiftedCopy(size_t shiftAmount) const;
};

#endif // TEUN_GAME_BIGUNSIGNEDINT_H