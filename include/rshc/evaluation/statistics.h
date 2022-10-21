#ifndef RSHC_EVALUATION_STATISTICS_H_
#define RSHC_EVALUATION_STATISTICS_H_

#include <vector>

namespace rshc::eval {
struct Statistics {
  double P_joint;            ///< Joint probability of violation
  
  std::vector<double> P_wl;  ///< Probability of violating voltage lower limit
  std::vector<double> P_wu;  ///< Probability of violating voltage upper limit
  std::vector<double> P_f;   ///< Probability of violating line flow limit

  std::vector<double> CVaR_wl; ///< CVaR of voltage lower limit
  std::vector<double> CVaR_wu; ///< CVaR of voltage upper limit
  std::vector<double> CVaR_f;  ///< CVaR of line flow limit

  std::vector<double> E_wl;   ///< Expected voltage lower limit violation
  std::vector<double> E_wu;   ///< Expected voltage upper limit violation
  std::vector<double> E_f;    ///< Expected line flow limit violation

  std::vector<double> max_wl;   ///< Worst case voltage (lower)
  std::vector<double> max_wu;   ///< Worst case voltage (upper)
  std::vector<double> max_f;    ///< Worst case line flow
};
}  // namespace rshc::eval

#endif  // RSHC_EVALUATION_STATISTICS_H_