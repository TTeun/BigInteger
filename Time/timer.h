#ifndef __BIG_INT__H__
#define __BIG_INT__H__

#include <chrono>
#include <iostream>
#include <utility>

class Timer {

public:
    Timer();
    ~Timer();
    double elapsed() const;

private:
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
};

#endif // __BIG_INT__H__
