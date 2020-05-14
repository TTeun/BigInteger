#ifndef TEUN_GAME_DIGITVECTOR_H
#define TEUN_GAME_DIGITVECTOR_H

#include <cmath>
#include <limits>
#include <vector>

typedef std::vector<size_t>::reverse_iterator       leftToRightIterator;
typedef std::vector<size_t>::const_reverse_iterator leftToRightConstIterator;
typedef std::vector<size_t>::iterator               rightToLeftIterator;
typedef std::vector<size_t>::const_iterator         rightToLeftConstIterator;

class DigitVector {

public:
    static const size_t s_base =
        static_cast<size_t>(std::sqrt(std::numeric_limits<size_t>::max()) - 1ul) & (std::numeric_limits<size_t>::max() - 1ul);
    //    static const size_t s_base = 20ul;

    static const size_t s_maxDigit = s_base - 1ul;

    size_t digitCount() const;

protected:
    size_t digitAt(size_t index) const;

    DigitVector();

    DigitVector(std::vector<size_t> &&digits);

    void shift(size_t shiftAmount);

    size_t mostSignificantDigit() const;

    size_t leastSignificantDigit() const;

    bool isWellFormed() const;

    leftToRightIterator leftToRightBegin();

    leftToRightIterator leftToRightEnd();

    rightToLeftIterator rightToLeftBegin();

    rightToLeftIterator rightToLeftEnd();

    leftToRightConstIterator leftToRightConstBegin() const;

    leftToRightConstIterator leftToRightConstEnd() const;

    rightToLeftConstIterator rightToLeftConstBegin() const;

    rightToLeftConstIterator rightToLeftConstEnd() const;

    bool isCorrectlySized() const;

    std::vector<size_t> m_digits;
};

#endif // TEUN_GAME_DIGITVECTOR_H