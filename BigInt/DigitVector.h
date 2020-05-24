#ifndef __DIGIT_VECTOR__H__
#define __DIGIT_VECTOR__H__

#include <cassert>
#include <ostream>
#include <vector>

typedef std::vector<size_t>::reverse_iterator lrIterator;
typedef std::vector<size_t>::const_reverse_iterator lrcIterator;
typedef std::vector<size_t>::iterator rlIterator;
typedef std::vector<size_t>::const_iterator rlcIterator;

namespace big {

    class DigitVector {

    public:
        static const size_t s_base = 1000000000ul;
        static const size_t s_maxDigit = 999999999ul;
        static const size_t s_digitsPerLimb = 9ul;

        friend class BigInt;

    public:
        /***************** Constructors *****************/
        DigitVector() = default;
        explicit DigitVector(std::vector<size_t> &&digits) : m_digits(std::move(digits)) {
        }

        /***************** Const functions *****************/
        size_t digitAt(size_t index) const {
            assert(index < m_digits.size());
            return m_digits.at(index);
        }

        size_t leastSignificantDigit() const {
            return m_digits.front();
        }

        bool isWellFormed() const;

        size_t digitCount() const {
            return m_digits.size();
        }

        size_t mostSignificantDigit() const {
            return m_digits.back();
        }
        bool isCorrectlySized() const;

        /***************** Iterators *****************/
        lrIterator leftToRightBegin() {
            return m_digits.rbegin();
        }
        lrIterator leftToRightEnd() {
            return m_digits.rend();
        }
        rlIterator rightToLeftEnd() {
            return m_digits.end();
        }
        lrcIterator leftToRightConstBegin() const {
            return m_digits.crbegin();
        }
        lrcIterator leftToRightConstEnd() const {
            return m_digits.crend();
        }
        rlcIterator rightToLeftConstBegin() const {
            return m_digits.cbegin();
        }
        rlcIterator rightToLeftConstEnd() const {
            return m_digits.cend();
        }
        rlIterator rightToLeftBegin() {
            return m_digits.begin();
        }

        /***************** Vector functions *****************/
        void reserve(size_t size) {
            m_digits.reserve(size);
        }
        void resizeToFit() {
            resizeToFitVector(m_digits);
        }
        void resize(size_t size) {
            m_digits.resize(size);
        }
        void shift(size_t shiftAmount) {
            m_digits.insert(rightToLeftBegin(), shiftAmount, 0);
        }
        static void resizeToFitVector(std::vector<size_t> &digits);

        /***************** Data members *****************/
        std::vector<size_t> m_digits;
    };
} // namespace big

#endif // __DIGIT_VECTOR__H__