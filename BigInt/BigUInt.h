#ifndef __BIG_U_INT__H__
#define __BIG_U_INT__H__

#include "BigUIntBase.h"
#include "DigitVector.h"

#include <ostream>

class BigUInt : public BigUIntBase {

public:
    BigUInt();

    BigUInt(const BigUInt &other);

    explicit BigUInt(std::vector<size_t> &&digits);

    BigUInt(size_t val);

    BigUInt(const std::string &val);

    static BigUInt createRandom(size_t numberOfDigits);

    static BigUInt createRandomFromDecimalDigits(size_t orderOfMagnitude);

    size_t approximatePowerOfTen() const;

    BigUInt &operator=(const BigUInt &rhs);

    BigUInt &operator=(size_t rhs);

    BigUInt &operator+=(const BigUInt &rhs);

    BigUInt &operator-=(const BigUInt &rhs);

    BigUInt operator-(const BigUInt &rhs) const;

    BigUInt &operator+=(size_t rhs);

    BigUInt &operator*=(const BigUInt &rhs);

    BigUInt &operator*=(size_t rhs);

    BigUInt operator+(size_t rhs) const;

    BigUInt operator+(const BigUInt &rhs) const;

    BigUInt operator*(size_t rhs) const;

    BigUInt operator*(const BigUInt &rhs) const;

    size_t operator%(size_t mod) const;

    BigUInt operator%(const BigUInt &mod) const;

    BigUInt &operator%=(size_t mod);

    BigUInt &operator%=(const BigUInt &mod);

    BigUInt operator/(const BigUInt &divisor) const;

    bool operator==(const BigUInt &rhs) const;

    bool operator!=(const BigUInt &rhs) const;

    bool operator<(const BigUInt &rhs) const;

    bool operator>(const BigUInt &rhs) const;

    bool operator<=(const BigUInt &rhs) const;

    bool operator>=(const BigUInt &rhs) const;

    friend std::ostream &operator<<(std::ostream &os, const BigUInt &anInt);

    friend BigUInt power(const BigUInt &base, size_t exponent);

    BigUInt copyPrefix(size_t length) const;

    BigUInt copySuffix(size_t length) const;

private:
    static size_t divisionSubRoutine(const std::vector<size_t>::const_reverse_iterator &leftToRightConstIt,
                                     const std::vector<size_t>::const_reverse_iterator &leftToRightConstEnd,
                                     const std::vector<size_t>::iterator &              rightToLeftIt,
                                     const std::vector<size_t>::iterator &              rightToLeftEnd,
                                     const BigUInt &                                    divisor);

    friend BigUInt longDivision(BigUInt &dividend, const BigUInt &divisor);

    friend BigUInt longDivisionAfterAdjustingDivisor(BigUInt &dividend, const BigUInt &divisor);

    void square();


    void init(size_t val);

    void bubble(size_t startIndex = 0ul);

    //    void shiftedCopy(size_t shiftAmount);

    BigUInt shiftedCopy(size_t shiftAmount) const;
};

#endif // __BIG_U_INT__H__