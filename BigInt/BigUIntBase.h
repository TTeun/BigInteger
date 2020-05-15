#ifndef TEUN_GAME_BIGUINTBASE_H
#define TEUN_GAME_BIGUINTBASE_H

#include "DigitVector.h"

#include <ostream>

class BigUIntBase : public DigitVector {

protected:
    static const size_t s_karatsubaLowerLimit = 200ul;

    friend void karatsubaMultiplyViaIterators(rightToLeftIterator             resultIt,
                                              const rightToLeftIterator &     resultEnd,
                                              rightToLeftConstIterator        rhsIt,
                                              const rightToLeftConstIterator &rhsEnd,
                                              rightToLeftConstIterator        copyIt,
                                              const rightToLeftConstIterator &copyEnd);

    friend void splitOneMultiplicationViaIterators(rightToLeftIterator             resultIt,
                                                   const rightToLeftIterator &     resultEnd,
                                                   rightToLeftConstIterator        rhsIt,
                                                   const rightToLeftConstIterator &rhsEnd,
                                                   rightToLeftConstIterator        largeIt,
                                                   const rightToLeftConstIterator &largeEnd);

public:
    BigUIntBase();

    BigUIntBase(std::vector<size_t> &&digits);

    BigUIntBase(rightToLeftConstIterator it, rightToLeftConstIterator endIt);
};

#endif // TEUN_GAME_BIGUINTBASE_H
