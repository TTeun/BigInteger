//
// Created by pc on 5/11/20.
//

#ifndef TEUN_GAME_BIGINT_H
#define TEUN_GAME_BIGINT_H

#include "BigUnsignedInt.h"

class BigInt {

public:
    BigInt();

    BigInt(size_t val);

    BigInt(long long val);

    BigInt(const std::string & val);

    BigInt(const BigUnsignedInt & unsignedInt);

    BigInt(BigUnsignedInt && unsignedInt);

    bool operator==(const BigInt & rhs) const;
    bool operator!=(const BigInt & rhs) const;

private:
    BigUnsignedInt m_magnitude;

    bool m_isPositive = true;
};

#endif // TEUN_GAME_BIGINT_H
