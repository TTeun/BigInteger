//
// Created by pc on 5/15/20.
//

#ifndef TEUN_GAME_BIGUINTBASE_H
#define TEUN_GAME_BIGUINTBASE_H

#include "DigitVector.h"

#include <ostream>

class BigUIntBase : public DigitVector {

private:
    static const size_t s_karatsubaLowerLimit = 200ul;

public:
    BigUIntBase();

    BigUIntBase(std::vector<size_t> &&digits);
};

#endif // TEUN_GAME_BIGUINTBASE_H
