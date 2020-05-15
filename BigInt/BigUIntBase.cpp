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
