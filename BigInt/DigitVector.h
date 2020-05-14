#ifndef TEUN_GAME_DIGITVECTOR_H
#define TEUN_GAME_DIGITVECTOR_H

#include <cmath>
#include <limits>
#include <vector>

class DigitVector {

public:
    static const size_t s_base =
        static_cast<size_t>(std::sqrt(std::numeric_limits<size_t>::max()) - 1ul) & (std::numeric_limits<size_t>::max() - 1ul);
    //    static const size_t s_base = 100000ul;

    static const size_t s_maxDigit = s_base - 1ul;

protected:
    size_t digitAt(size_t index) const;

    DigitVector();

    DigitVector(std::vector<size_t> &&digits);

    size_t digitCount() const;

    void shift(size_t shiftAmount);

    size_t mostSignificantDigit() const;

    size_t leastSignificantDigit() const;

    bool isWellFormed() const;

    std::vector<size_t> m_digits;

    std::vector<size_t>::reverse_iterator leftToRightBegin();

    std::vector<size_t>::reverse_iterator leftToRightEnd();

    std::vector<size_t>::iterator rightToLeftBegin();

    std::vector<size_t>::iterator rightToLeftEnd();

    std::vector<size_t>::const_reverse_iterator leftToRightConstBegin() const;

    std::vector<size_t>::const_reverse_iterator leftToRightConstEnd() const;

    std::vector<size_t>::const_iterator rightToLeftConstBegin() const;

    std::vector<size_t>::const_iterator rightToLeftConstEnd() const;

    bool isCorrectlySized() const;
};

#endif // TEUN_GAME_DIGITVECTOR_H