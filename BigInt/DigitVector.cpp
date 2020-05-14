#include "DigitVector.h"

#include <cassert>

DigitVector::DigitVector() {
}

DigitVector::DigitVector(std::vector<size_t> &&digits) : m_digits(std::move(digits)) {
}

bool DigitVector::isCorrectlySized() const {
    assert(digitCount() > 0);
    if (mostSignificantDigit() == 0ul) {
        return digitCount() == 1; // Is the value zero
    }
    return true;
}

size_t DigitVector::mostSignificantDigit() const {
    return m_digits.back();
}

leftToRightIterator DigitVector::leftToRightBegin() {
    return m_digits.rbegin();
}

leftToRightIterator DigitVector::leftToRightEnd() {
    return m_digits.rend();
}

rightToLeftIterator DigitVector::rightToLeftBegin() {
    return m_digits.begin();
}

rightToLeftIterator DigitVector::rightToLeftEnd() {
    return m_digits.end();
}

leftToRightConstIterator DigitVector::leftToRightConstBegin() const {
    return m_digits.crbegin();
}

leftToRightConstIterator DigitVector::leftToRightConstEnd() const {
    return m_digits.crend();
}

rightToLeftConstIterator DigitVector::rightToLeftConstBegin() const {
    return m_digits.cbegin();
}

rightToLeftConstIterator DigitVector::rightToLeftConstEnd() const {
    return m_digits.cend();
}

bool DigitVector::isWellFormed() const {
    if (not isCorrectlySized()) {
        return false;
    }
    for (const auto &it : m_digits) {
        if (it >= s_base) {
            return false;
        }
    }
    return true;
}

size_t DigitVector::digitCount() const {
    return m_digits.size();
}

size_t DigitVector::digitAt(size_t index) const {
    assert(index < m_digits.size());
    return m_digits.at(index);
}

void DigitVector::shift(size_t shiftAmount) {
    m_digits.insert(rightToLeftBegin(), shiftAmount, 0);
}

size_t DigitVector::leastSignificantDigit() const {
    return m_digits.front();
}
