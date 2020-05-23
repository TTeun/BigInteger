#ifndef __BIG_U_INT_BASE__H__
#define __BIG_U_INT_BASE__H__

#include "DigitVector.h"

class BigUIntBase : public DigitVector {

public:
    BigUIntBase() {
    }

    BigUIntBase(std::vector<size_t> &&digits) : DigitVector(std::move(digits)) {
    }

    BigUIntBase(rightToLeftConstIterator it, rightToLeftConstIterator endIt) : DigitVector({it, endIt}) {
    }

protected:
    static const size_t s_karatsubaLowerLimit = 200ul;

    static void karatsubaMultiplyViaIterators(rightToLeftIterator            resultIt,
                                              const rightToLeftIterator      resultEnd,
                                              rightToLeftConstIterator       rhsIt,
                                              const rightToLeftConstIterator rhsEnd,
                                              rightToLeftConstIterator       copyIt,
                                              const rightToLeftConstIterator copyEnd);

    static void splitOneMultiplicationViaIterators(rightToLeftIterator            resultIt,
                                                   const rightToLeftIterator      resultEnd,
                                                   rightToLeftConstIterator       rhsIt,
                                                   const rightToLeftConstIterator rhsEnd,
                                                   rightToLeftConstIterator       largeIt,
                                                   const rightToLeftConstIterator largeEnd);

    static void multiplyViaIterators(rightToLeftIterator            resultIt,
                                     const rightToLeftIterator      resultEnd,
                                     rightToLeftConstIterator       rhsIt,
                                     const rightToLeftConstIterator rhsEnd,
                                     rightToLeftConstIterator       copyIt,
                                     const rightToLeftConstIterator copyEnd);

    static void addMultipleViaIterators(rightToLeftIterator      thisIt,
                                        rightToLeftIterator      thisEnd,
                                        rightToLeftConstIterator rhsIt,
                                        rightToLeftConstIterator rhsEnd,
                                        size_t                   multiplier);

    static void addViaIterators(rightToLeftIterator      thisIt,
                                rightToLeftIterator      thisEnd,
                                rightToLeftConstIterator rhsIt,
                                rightToLeftConstIterator rhsEnd);

    static void subtractViaIterators(rightToLeftIterator            thisIt,
                                     const rightToLeftIterator      thisEnd,
                                     rightToLeftConstIterator       rhsIt,
                                     const rightToLeftConstIterator rhsEnd);

    static bool lessThanViaIterators(const leftToRightConstIterator thisIt,
                                     const leftToRightConstIterator thisEnd,
                                     const leftToRightConstIterator rhsIt,
                                     const leftToRightConstIterator rhsEnd) {
        return lessThanShiftedRhsViaIterators(thisIt, thisEnd, rhsIt, rhsEnd, 0ul);
    }

    static bool greaterThanViaIterators(const leftToRightConstIterator thisIt,
                                        const leftToRightConstIterator thisEnd,
                                        const leftToRightConstIterator rhsIt,
                                        const leftToRightConstIterator rhsEnd) {
        return lessThanViaIterators(rhsIt, rhsEnd, thisIt, thisEnd);
    }

    static bool lessThanShiftedRhsViaIterators(leftToRightConstIterator       thisIt,
                                               const leftToRightConstIterator thisEnd,
                                               leftToRightConstIterator       rhsIt,
                                               const leftToRightConstIterator rhsEnd,
                                               size_t                         trailingZeroesOfRhs);

    static void
    carryAdditionViaIterators(rightToLeftIterator thisIt, const rightToLeftIterator thisEnd, size_t carry);
};

#endif // __BIG_U_INT_BASE__H__
