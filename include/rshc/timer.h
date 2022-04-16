#ifndef RSHC_TIMER_H_
#define RSHC_TIMER_H_

#include <chrono>

namespace rshc {
class Timer {
 public:
  Timer() { time_ = std::chrono::high_resolution_clock::now(); }

  ~Timer() = default;

  void start() { time_ = std::chrono::high_resolution_clock::now(); }

  std::chrono::duration<double> elapsed() {
    return std::chrono::high_resolution_clock::now() - time_;
  }

 protected:
  std::chrono::time_point<std::chrono::high_resolution_clock> time_;
};
}  // namespace rshc

#endif  // RSHC_TIMER_H_