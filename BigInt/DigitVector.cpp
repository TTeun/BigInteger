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

size_t DigitVector::secondMostSignificantDigit() const {
    assert(digitCount() > 1);
    return *(m_digits.rbegin() + 1);
}

std::vector<size_t>::reverse_iterator DigitVector::leftToRightBegin() {
    return m_digits.rbegin();
}

std::vector<size_t>::reverse_iterator DigitVector::leftToRightEnd() {
    return m_digits.rend();
}

std::vector<size_t>::iterator DigitVector::rightToLeftBegin() {
    return m_digits.begin();
}

std::vector<size_t>::iterator DigitVector::rightToLeftEnd() {
    return m_digits.end();
}

std::vector<size_t>::const_reverse_iterator DigitVector::leftToRightConstBegin() const {
    return m_digits.crbegin();
}

std::vector<size_t>::const_reverse_iterator DigitVector::leftToRightConstEnd() const {
    return m_digits.crend();
}

std::vector<size_t>::const_iterator DigitVector::rightToLeftConstBegin() const {
    return m_digits.cbegin();
}

std::vector<size_t>::const_iterator DigitVector::rightToLeftConstEnd() const {
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
