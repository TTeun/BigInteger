#include "DigitVector.h"

#include <cassert>

bool DigitVector::isCorrectlySized() const {
    assert(digitCount() > 0);
    if (mostSignificantDigit() == 0ul) {
        return digitCount() == 1; // Is the value zero
    }
    return true;
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

void DigitVector::resizeToFitVector(std::vector<size_t> &digits) {
    auto it = digits.rbegin();
    for (; it != digits.rend() && *it == 0ul; ++it)
        ;

    digits.resize(static_cast<unsigned long>(std::distance(it, digits.rend())));

    if (digits.empty()) {
        digits = {0ul};
    }
}
