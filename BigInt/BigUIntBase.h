#ifndef __BIG_U_INT_BASE__H__
#define __BIG_U_INT_BASE__H__

#include "DigitVector.h"

class BigUIntBase : public DigitVector {

public:
    BigUIntBase();

    BigUIntBase(rightToLeftConstIterator it, rightToLeftConstIterator endIt);

    explicit BigUIntBase(std::vector<size_t> &&digits);

protected:
    static const size_t s_karatsubaLowerLimit = 200ul;

    static void karatsubaMultiplyViaIterators(rightToLeftIterator             resultIt,
                                              const rightToLeftIterator &     resultEnd,
                                              rightToLeftConstIterator        rhsIt,
                                              const rightToLeftConstIterator &rhsEnd,
                                              rightToLeftConstIterator        copyIt,
                                              const rightToLeftConstIterator &copyEnd);

    static void splitOneMultiplicationViaIterators(rightToLeftIterator             resultIt,
                                                   const rightToLeftIterator &     resultEnd,
                                                   rightToLeftConstIterator        rhsIt,
                                                   const rightToLeftConstIterator &rhsEnd,
                                                   rightToLeftConstIterator        largeIt,
                                                   const rightToLeftConstIterator &largeEnd);

    static void multiplyViaIterators(rightToLeftIterator             resultIt,
                                     const rightToLeftIterator &     resultEnd,
                                     rightToLeftConstIterator        rhsIt,
                                     const rightToLeftConstIterator &rhsEnd,
                                     rightToLeftConstIterator        copyIt,
                                     const rightToLeftConstIterator &copyEnd);

public:
    static void multiplyBySmallNumberViaIterators(leftToRightIterator resultIt, leftToRightIterator resultEnd, const size_t rhs);

protected:
    static void addMultipleViaIterators(rightToLeftIterator                 thisIt,
                                        rightToLeftIterator                 thisEnd,
                                        std::vector<size_t>::const_iterator rhsIt,
                                        std::vector<size_t>::const_iterator rhsEnd,
                                        size_t                              multiplier);

    static void addViaIterators(rightToLeftIterator                 thisIt,
                                rightToLeftIterator                 thisEnd,
                                std::vector<size_t>::const_iterator rhsIt,
                                std::vector<size_t>::const_iterator rhsEnd);

    static void subtractViaIterators(rightToLeftIterator                 thisIt,
                                     rightToLeftIterator                 thisEnd,
                                     std::vector<size_t>::const_iterator rhsIt,
                                     std::vector<size_t>::const_iterator rhsEnd);

    static void subtractViaIterators(leftToRightIterator            thisIt,
                                     const leftToRightIterator      thisEnd,
                                     leftToRightConstIterator       rhsIt,
                                     const leftToRightConstIterator rhsEnd);

    static void subtractShiftedViaIterators(leftToRightIterator            thisIt,
                                            const leftToRightIterator      thisEnd,
                                            leftToRightConstIterator       rhsIt,
                                            const leftToRightConstIterator rhsEnd,
                                            const size_t                   shiftAmount);

    static bool lessThanViaIterators(const leftToRightConstIterator &thisIt,
                                     const leftToRightConstIterator &thisEnd,
                                     const leftToRightConstIterator &rhsIt,
                                     const leftToRightConstIterator &rhsEnd);

    static bool greaterThanViaIterators(const leftToRightConstIterator &thisIt,
                                        const leftToRightConstIterator &thisEnd,
                                        const leftToRightConstIterator &rhsIt,
                                        const leftToRightConstIterator &rhsEnd);

    static bool lessThanShiftedRhsViaIterators(leftToRightConstIterator        thisIt,
                                               const leftToRightConstIterator &thisEnd,
                                               leftToRightConstIterator        rhsIt,
                                               const leftToRightConstIterator &rhsEnd,
                                               size_t                          trailingZeroesOfRhs);

    static void carryAdditionViaIterators(rightToLeftIterator thisIt, const rightToLeftIterator &thisEnd, size_t carry);
};

#endif // __BIG_U_INT_BASE__H__
