#define NDEBUG

#include "Time/BenchMark.h"

int main() {
    BenchMark::run({{10, 100}, {100, 100}, {1000, 100}, {10000, 50}, {100000, 10}, {1000000, 5}});
}
