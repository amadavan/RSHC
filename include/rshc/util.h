#ifndef RSHC_UTIL_H_
#define RSHC_UTIL_H_

#include <iostream>
#include <iomanip>

#include <Eigen/Core>
#include <vector>

#include "timer.h"

namespace rshc::util {
template <typename T, typename U>
std::vector<T> getEigenSTL(const U &op) {
  Eigen::VectorXd vec = op;
  std::vector<T> data = std::vector<T>(vec.size(), 0);
  std::transform(vec.data(), vec.data() + vec.size(), data.begin(),
                 [](const double &arg) { return (T)arg; });
  return data;
}

class ProgressBar {
 public:
  ProgressBar(size_t max) : max_(max) { timer_.start(); }
  ~ProgressBar() {}

  void setValue(size_t value) {
    double cp = (100. * value) / max_;
    if (cp > cval_ + 0.1) {
      cval_ = cp;
      double elapsed = timer_.elapsed().count();

      std::cout << std::right << std::setw(5) << std::setprecision(4) << cp
                << "%"
                << "\t"
                << "Elapsed: " << std::left << std::setw(10)
                << std::setprecision(10) << elapsed << "s"
                << "\t"
                << "Estimated: " << std::left << std::setw(10)
                << std::setprecision(10) << (100. - cp) * (elapsed / cp) << "s"
                << std::endl;
    }
  }

  void finalize() {
    double elapsed = timer_.elapsed().count();
    std::cout << std::right << std::setw(5) << std::setprecision(4) << 100
              << "%"
              << "\t"
              << "Total time: " << std::left << std::setw(10)
              << std::setprecision(10) << elapsed << "s" << std::endl;
  }

 private:
  Timer timer_;

  size_t max_ = 0;
  double cval_ = 0;
};

}  // namespace rshc::util

#endif  // RSHC_UTIL_H_