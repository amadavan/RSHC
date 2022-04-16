#ifndef RSHC_UTIL_H_
#define RSHC_UTIL_H_

#include <Eigen/Core>
#include <vector>

namespace rshc::util {
template <typename T, typename U>
std::vector<T> getEigenSTL(const U &op) {
  Eigen::VectorXd vec = op;
  std::vector<T> data = std::vector<T>(vec.size(), 0);
  std::transform(vec.data(), vec.data() + vec.size(), data.begin(),
                 [](const double &arg) { return (T)arg; });
  return data;
}
}  // namespace rshc::util

#endif  // RSHC_UTIL_H_