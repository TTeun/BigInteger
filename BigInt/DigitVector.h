#ifndef TEUN_GAME_DIGITVECTOR_H
#define TEUN_GAME_DIGITVECTOR_H

#include <cassert>
#include <cmath>
#include <limits>
#include <ostream>
#include <vector>

typedef std::vector<size_t>::reverse_iterator       leftToRightIterator;
typedef std::vector<size_t>::const_reverse_iterator leftToRightConstIterator;
typedef std::vector<size_t>::iterator               rightToLeftIterator;
typedef std::vector<size_t>::const_iterator         rightToLeftConstIterator;

class DigitVector {

public:
    static const size_t s_base = 1000000000ul;

    static const size_t s_maxDigit = s_base - 1ul;

    size_t digitCount() const {
        return m_digits.size();
    }

    friend class BigInt;

public:
    DigitVector() = default;

    explicit DigitVector(std::vector<size_t> &&digits) : m_digits(std::move(digits)) {
    }

    size_t digitAt(size_t index) const {
        assert(index < m_digits.size());
        return m_digits.at(index);
    }

    void shift(size_t shiftAmount) {
        m_digits.insert(rightToLeftBegin(), shiftAmount, 0);
    }

    size_t leastSignificantDigit() const {
        return m_digits.front();
    }

    bool isWellFormed() const;

    size_t mostSignificantDigit() const {
        return m_digits.back();
    }

    leftToRightIterator leftToRightBegin() {
        return m_digits.rbegin();
    }

    leftToRightIterator leftToRightEnd() {
        return m_digits.rend();
    }

    rightToLeftIterator rightToLeftEnd() {
        return m_digits.end();
    }

    leftToRightConstIterator leftToRightConstBegin() const {
        return m_digits.crbegin();
    }

    leftToRightConstIterator leftToRightConstEnd() const {
        return m_digits.crend();
    }

    rightToLeftConstIterator rightToLeftConstBegin() const {
        return m_digits.cbegin();
    }

    rightToLeftConstIterator rightToLeftConstEnd() const {
        return m_digits.cend();
    }

    bool isCorrectlySized() const;

    void reserve(size_t size) {
        m_digits.reserve(size);
    }

    void resizeToFit() {
        resizeToFitVector(m_digits);
    }

    static void resizeToFitVector(std::vector<size_t> &digits);

    std::vector<size_t> m_digits;
    rightToLeftIterator rightToLeftBegin() {
        return m_digits.begin();
    }
    void resize(size_t size) {
        m_digits.resize(size);
    }
};

#endif // TEUN_GAME_DIGITVECTOR_H