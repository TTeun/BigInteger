#ifndef __BIG_INT__H__
#define __BIG_INT__H__

#include "BigUInt.h"

#include <ostream>

namespace big {

    class BigInt {

    public:
        /***************** Constructors *****************/
        BigInt();

        BigInt(BigInt &&other) noexcept;

        BigInt(const BigInt &other);

        explicit BigInt(BigUInt &&magnitude);

        explicit BigInt(const BigUInt &other);

        explicit BigInt(const std::string &val);

        /***************** Operators *****************/
        BigInt &operator=(size_t rhs);
        BigInt &operator=(const BigUInt &rhs);
        BigInt &operator=(const BigInt &rhs);
        BigInt &operator=(BigInt &&rhs);
        BigInt &operator=(BigUInt &&rhs);

        /** Addition **/
        BigInt &operator+=(size_t rhs);
        BigInt &operator+=(const BigInt &rhs);
        BigInt &operator+=(const BigUInt &rhs);
        BigInt operator+(size_t rhs) const;
        BigInt operator+(const BigInt &rhs) const;

        /** Subtraction **/
        BigInt &operator-=(const BigInt &rhs);
        BigInt &operator-=(const BigUInt &rhs);
        BigInt operator-(const BigInt &rhs) const;
        BigInt operator-(const BigUInt &rhs) const;

        /** Multiplication **/
        BigInt &operator*=(size_t rhs);
        BigInt &operator*=(const BigInt &rhs);
        BigInt operator*(size_t rhs) const;
        BigInt operator*(const BigInt &rhs) const;

        /** Division **/
        BigInt &operator/=(size_t divisor);
        BigInt &operator/=(const BigInt &divisor);
        BigInt &operator/=(const BigUInt &divisor);
        BigInt operator/(const BigInt &divisor) const;
        BigInt operator/(const BigUInt &divisor) const;
        BigInt operator/(size_t divisor) const;

        /** Comparison **/

        /** Friends **/
        friend BigInt operator+(size_t lhs, const BigInt &rhs);
        friend BigInt operator-(size_t lhs, const BigInt &rhs);
        friend BigInt operator*(size_t lhs, const BigInt &rhs);

        /***************** Output *****************/
        std::string toString() const {
            return m_isNegative ? "=" : "" + m_magnitude.toString();
        }
        friend std::ostream &operator<<(std::ostream &os, const BigInt &anInt);

        /***************** Vector functions *****************/
        void resize(size_t size);
        void reserve(size_t size);

        friend BigUInt;

        /***************** Accessors *****************/
        const BigUInt &magnitude() const {
            return m_magnitude;
        }

    private:
        void negate() {
            m_isNegative = !m_isNegative;
        }
        BigUInt &magnitude() {
            return m_magnitude;
        }

        /***************** Data members *****************/
        BigUInt m_magnitude;
        bool m_isNegative = false;
    };
} // namespace big
#endif // __BIG_INT__H__
