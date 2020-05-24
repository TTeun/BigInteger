#ifndef __BENCHMARK__H__
#define __BENCHMARK__H__

#include <cstddef>

namespace BenchMark {

    void run();

    double multiply(size_t powerOfTen = 40ul, size_t numberOfRepetitions = 50ul);
    double add(size_t powerOfTen = 40ul, size_t numberOfRepetitions = 50ul);
    double divide(size_t powerOfTen = 40ul, size_t numberOfRepetitions = 50ul);
    double modulo(size_t powerOfTen = 40ul, size_t numberOfRepetitions = 50ul);

} // namespace BenchMark

#endif // __BENCHMARK__H__
