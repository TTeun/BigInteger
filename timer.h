#ifndef TEUN_GAME_TIMER_H
#define TEUN_GAME_TIMER_H

#include <chrono>
#include <iostream>
#include <utility>
class Timer {

public:
    explicit Timer(std::string string = "");

    ~Timer();

private:
    std::string m_string;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
};

#endif // TEUN_GAME_TIMER_H
