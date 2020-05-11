#include "timer.h"

Timer::Timer(std::string string) : m_string(std::move(string))
{
}

Timer::~Timer()
{
    if (m_string.length() != 0) {
        std::cout << m_string << " took ";
    }
    std::cout << std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() -
                                                                           t1)
                     .count()
              << " seconds\n";
}
