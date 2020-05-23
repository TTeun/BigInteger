//
// Created by pc on 5/15/20.
//

#include "BigUIntBase.h"

#include "BigUInt.h"

#include <cassert>
#include <iostream>
#include <random>
#include <sstream>
BigUIntBase::BigUIntBase() {
}

BigUIntBase::BigUIntBase(std::vector<size_t> &&digits) : DigitVector(std::move(digits)) {
}

BigUIntBase::BigUIntBase(rightToLeftConstIterator it, rightToLeftConstIterator endIt) : DigitVector({it, endIt}) {
}

bool BigUIntBase::lessThanShiftedRhsViaIterators(leftToRightConstIterator        thisIt,
                                                 const leftToRightConstIterator &thisEnd,
                                                 leftToRightConstIterator        rhsIt,
                                                 const leftToRightConstIterator &rhsEnd,
                                                 const size_t                    trailingZeroesOfRhs) {
    if (static_cast<size_t>(thisEnd - thisIt) != static_cast<size_t>(rhsEnd - rhsIt) + trailingZeroesOfRhs) {
        return static_cast<size_t>(thisEnd - thisIt) < static_cast<size_t>(rhsEnd - rhsIt) + trailingZeroesOfRhs;
    }

    for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
        if (*thisIt != *rhsIt) {
            return *thisIt < *rhsIt;
        }
    }
    return false;
}

bool BigUIntBase::lessThanViaIterators(const leftToRightConstIterator &thisIt,
                                       const leftToRightConstIterator &thisEnd,
                                       const leftToRightConstIterator &rhsIt,
                                       const leftToRightConstIterator &rhsEnd) {
    return lessThanShiftedRhsViaIterators(thisIt, thisEnd, rhsIt, rhsEnd, 0ul);
}

bool BigUIntBase::greaterThanViaIterators(const leftToRightConstIterator &thisIt,
                                          const leftToRightConstIterator &thisEnd,
                                          const leftToRightConstIterator &rhsIt,
                                          const leftToRightConstIterator &rhsEnd) {
    return lessThanViaIterators(rhsIt, rhsEnd, thisIt, thisEnd);
}

void BigUIntBase::carryAdditionViaIterators(rightToLeftIterator thisIt, const rightToLeftIterator &thisEnd, size_t carry) {
    assert(carry != 0ul);
    assert(carry + DigitVector::s_base < std::numeric_limits<size_t>::max());
    for (; thisIt != thisEnd; ++thisIt) {
        *thisIt += carry;
        if (*thisIt < DigitVector::s_base) {
            return;
        }
        carry = (*thisIt) / DigitVector::s_base;
        *thisIt %= DigitVector::s_base;
    }
    assert(false);
}

void BigUIntBase::addViaIterators(rightToLeftIterator                       thisIt,
                                  const rightToLeftIterator                 thisEnd,
                                  std::vector<size_t>::const_iterator       rhsIt,
                                  const std::vector<size_t>::const_iterator rhsEnd) {
    assert(std::distance(thisIt, thisEnd) >= std::distance(rhsIt, rhsEnd) + 1l);
    size_t carry = 0ul;
    for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
        *thisIt += *rhsIt + carry;
        if (*thisIt >= DigitVector::s_base) {
            carry = (*thisIt) / DigitVector::s_base;
            *thisIt %= DigitVector::s_base;
        } else {
            carry = 0ul;
        }
    }
    if (carry == 1ul) {
        carryAdditionViaIterators(thisIt, thisEnd, 1ul);
    }
}
void BigUIntBase::addMultipleViaIterators(rightToLeftIterator                       thisIt,
                                          const rightToLeftIterator                 thisEnd,
                                          std::vector<size_t>::const_iterator       rhsIt,
                                          const std::vector<size_t>::const_iterator rhsEnd,
                                          const size_t                              multiplier) {
    if (multiplier == 0ul) {
        return;
    }
    assert(multiplier < DigitVector::s_base);
    assert(std::distance(thisIt, thisEnd) >= std::distance(rhsIt, rhsEnd) + 1l); // Always enough space for adding two proper digit vectors
    size_t carry = 0ul;
    for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
        *thisIt += (*rhsIt * multiplier) + carry;
        if (*thisIt >= DigitVector::s_base) {
            carry = (*thisIt) / DigitVector::s_base;
            *thisIt %= DigitVector::s_base;
        } else {
            carry = 0ul;
        }
    }
    if (carry > 0UL) {
        carryAdditionViaIterators(thisIt, thisEnd, carry);
    }
}

void BigUIntBase::multiplyViaIterators(rightToLeftIterator             resultIt,
                                       const rightToLeftIterator &     resultEnd,
                                       rightToLeftConstIterator        rhsIt,
                                       const rightToLeftConstIterator &rhsEnd,
                                       rightToLeftConstIterator        copyIt,
                                       const rightToLeftConstIterator &copyEnd) {
    const auto copySize = static_cast<size_t>(copyEnd - copyIt);
    const auto rhsSize  = static_cast<size_t>(rhsEnd - rhsIt);
    if (std::max(copySize, rhsSize) < BigUIntBase::s_karatsubaLowerLimit) {
        for (size_t i = 0; rhsIt != rhsEnd; ++i) {
            addMultipleViaIterators(resultIt + i, resultEnd, copyIt, copyEnd, *rhsIt);
            ++rhsIt;
        }
    } else if (std::min(copySize, rhsSize) >= BigUIntBase::s_karatsubaLowerLimit) {
        karatsubaMultiplyViaIterators(resultIt, resultEnd, rhsIt, rhsEnd, copyIt, copyEnd);
    } else if (copySize >= BigUIntBase::s_karatsubaLowerLimit) {
        splitOneMultiplicationViaIterators(resultIt, resultEnd, rhsIt, rhsEnd, copyIt, copyEnd);
    } else {
        splitOneMultiplicationViaIterators(resultIt, resultEnd, copyIt, copyEnd, rhsIt, rhsEnd);
    }
}

void BigUIntBase::splitOneMultiplicationViaIterators(rightToLeftIterator             resultIt,
                                                     const rightToLeftIterator &     resultEnd,
                                                     rightToLeftConstIterator        rhsIt,
                                                     const rightToLeftConstIterator &rhsEnd,
                                                     rightToLeftConstIterator        largeIt,
                                                     const rightToLeftConstIterator &largeEnd) {
    const auto m = static_cast<size_t>(largeEnd - largeIt);
    assert(m >= BigUIntBase::s_karatsubaLowerLimit);
    const size_t splitIndex = m / 2ul;

    BigUIntBase z0;
    z0.resize(splitIndex + (rhsEnd - rhsIt) + 1ul);
    multiplyViaIterators(z0.rightToLeftBegin(), z0.rightToLeftEnd(), largeIt, largeIt + splitIndex, rhsIt, rhsEnd);
    z0.resizeToFit();

    addViaIterators(resultIt, resultEnd, z0.rightToLeftConstBegin(), z0.rightToLeftConstEnd());

    BigUIntBase z1;
    z1.resize((largeEnd - (largeIt + splitIndex)) + (rhsEnd - rhsIt) + 1ul);
    multiplyViaIterators(z1.rightToLeftBegin(), z1.rightToLeftEnd(), largeIt + splitIndex, largeEnd, rhsIt, rhsEnd);
    z1.resizeToFit();

    addViaIterators(resultIt + splitIndex, resultEnd, z1.rightToLeftConstBegin(), z1.rightToLeftConstEnd());
}

void BigUIntBase::subtractViaIterators(rightToLeftIterator                       thisIt,
                                       const rightToLeftIterator                 thisEnd,
                                       std::vector<size_t>::const_iterator       rhsIt,
                                       const std::vector<size_t>::const_iterator rhsEnd) {
    assert(std::distance(thisIt, thisEnd) >= std::distance(rhsIt, rhsEnd));
    size_t carry = 0ul;
    for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
        if (*thisIt >= *rhsIt + carry) {
            *thisIt -= *rhsIt + carry;
            carry = 0;
        } else {
            *thisIt += DigitVector::s_base;
            *thisIt -= *rhsIt + carry;
            carry = 1ul;
        }
    }
    if (carry != 0ul) {
        while (true) {
            assert(thisIt != thisEnd);
            if (*thisIt > 0) {
                --*thisIt;
                break;
            }
            *thisIt = DigitVector::s_maxDigit;
            ++thisIt;
        }
    }
}

void BigUIntBase::karatsubaMultiplyViaIterators(rightToLeftIterator             resultIt,
                                                const rightToLeftIterator &     resultEnd,
                                                rightToLeftConstIterator        rhsIt,
                                                const rightToLeftConstIterator &rhsEnd,
                                                rightToLeftConstIterator        copyIt,
                                                const rightToLeftConstIterator &copyEnd) {
    const auto m = static_cast<size_t>(std::min(copyEnd - copyIt, rhsEnd - rhsIt));
    assert(m >= BigUIntBase::s_karatsubaLowerLimit);
    const size_t splitIndex = m / 2ul;

    BigUIntBase high1(copyIt + splitIndex, copyEnd);
    BigUIntBase high2(rhsIt + splitIndex, rhsEnd);

    BigUIntBase z0;
    z0.resize(2ul * splitIndex + 1ul);
    multiplyViaIterators(
        z0.rightToLeftBegin(), z0.rightToLeftEnd(), rhsIt, rhsIt + splitIndex, copyIt, copyIt + splitIndex); // z0 = low1 * low2;
    z0.resizeToFit();

    BigUIntBase z2;
    z2.resize((rhsEnd - rhsIt - splitIndex) + (copyEnd - copyIt - splitIndex) + 1ul);
    multiplyViaIterators(
        z2.rightToLeftBegin(), z2.rightToLeftEnd(), rhsIt + splitIndex, rhsEnd, copyIt + splitIndex, copyEnd); // z2 = high1 * high2;
    z2.resizeToFit();

    high1.resize(high1.digitCount() + 1ul);
    high2.resize(high2.digitCount() + 1ul);

    addViaIterators(high1.rightToLeftBegin(), high1.rightToLeftEnd(), copyIt, copyIt + splitIndex);
    addViaIterators(high2.rightToLeftBegin(), high2.rightToLeftEnd(), rhsIt, rhsIt + splitIndex);

    BigUIntBase z1;
    z1.resize((rhsEnd - rhsIt - splitIndex) + (copyEnd - copyIt - splitIndex) + 3ul);

    high1.resizeToFit();
    high2.resizeToFit();

    multiplyViaIterators(z1.rightToLeftBegin(),
                         z1.rightToLeftEnd(),
                         high1.rightToLeftConstBegin(),
                         high1.rightToLeftConstEnd(),
                         high2.rightToLeftConstBegin(),
                         high2.rightToLeftConstEnd()); // z1 = (high1 + low1) * (high2 + low2);

    subtractViaIterators(z1.rightToLeftBegin(), z1.rightToLeftEnd(), z2.rightToLeftConstBegin(), z2.rightToLeftConstEnd());
    subtractViaIterators(z1.rightToLeftBegin(), z1.rightToLeftEnd(), z0.rightToLeftConstBegin(), z0.rightToLeftConstEnd());

    addViaIterators(resultIt, resultEnd, z0.rightToLeftConstBegin(), z0.rightToLeftConstEnd());
    addViaIterators(resultIt + 2ul * splitIndex, resultEnd, z2.rightToLeftConstBegin(), z2.rightToLeftConstEnd());
    addViaIterators(resultIt + splitIndex, resultEnd, z1.rightToLeftConstBegin(), z1.rightToLeftConstEnd());
}
