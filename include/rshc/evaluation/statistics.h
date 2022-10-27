#ifndef RSHC_EVALUATION_STATISTICS_H_
#define RSHC_EVALUATION_STATISTICS_H_

#include <iostream>
#include <vector>

namespace rshc::eval {
struct Statistics {
  double P_joint;  ///< Joint probability of violation

  std::vector<double> P_wl;  ///< Probability of violating voltage lower limit
  std::vector<double> P_wu;  ///< Probability of violating voltage upper limit
  std::vector<double> P_f;   ///< Probability of violating line flow limit

  std::vector<double> CVaR_wl;  ///< CVaR of voltage lower limit
  std::vector<double> CVaR_wu;  ///< CVaR of voltage upper limit
  std::vector<double> CVaR_f;   ///< CVaR of line flow limit

  std::vector<double> E_wl;  ///< Expected voltage lower limit violation
  std::vector<double> E_wu;  ///< Expected voltage upper limit violation
  std::vector<double> E_f;   ///< Expected line flow limit violation

  std::vector<double> max_wl;  ///< Worst case voltage (lower)
  std::vector<double> max_wu;  ///< Worst case voltage (upper)
  std::vector<double> max_f;   ///< Worst case line flow

  friend std::ostream &operator<<(std::ostream &output,
                                  const Statistics &stats) {
    output << "P[fail]:\t" << stats.P_joint << std::endl;
    output << std::endl;

    for (size_t b = 0; b < stats.P_wl.size(); ++b)
      output << "P[w(" << b << ")"
             << "<wl(" << b << ")]:\t" << stats.P_wl[b] << std::endl;
    for (size_t b = 0; b < stats.P_wu.size(); ++b)
      output << "P[w(" << b << ")"
             << ">wu(" << b << ")]:\t" << stats.P_wu[b] << std::endl;
    for (size_t l = 0; l < stats.P_f.size(); ++l)
      output << "P[S2(" << l << ")"
             << "<f(" << l << ")]:\t" << stats.P_f[l] << std::endl;
    output << std::endl;

    for (size_t b = 0; b < stats.CVaR_wl.size(); ++b)
      output << "CVaR[w(" << b << ")]:\t" << stats.CVaR_wl[b] << std::endl;
    for (size_t b = 0; b < stats.CVaR_wu.size(); ++b)
      output << "CVaR[w(" << b << ")]:\t" << stats.CVaR_wu[b] << std::endl;
    for (size_t l = 0; l < stats.CVaR_f.size(); ++l)
      output << "CVaR[S2(" << l << ")]:\t" << stats.CVaR_f[l] << std::endl;
    output << std::endl;

    for (size_t b = 0; b < stats.E_wl.size(); ++b)
      output << "E[w(" << b << ") | w(" << b << ") < wl(" << b << ")]:\t"
             << stats.E_wl[b] << std::endl;
    for (size_t b = 0; b < stats.E_wu.size(); ++b)
      output << "E[w(" << b << ") | w(" << b << ") > wu(" << b << ")]:\t"
             << stats.E_wu[b] << std::endl;
    for (size_t l = 0; l < stats.E_f.size(); ++l)
      output << "E[S2(" << l << ") | S2(" << l << ") < f(" << l << ")]:\t"
             << stats.E_f[l] << std::endl;
    output << std::endl;

    for (size_t b = 0; b < stats.max_wl.size(); ++b)
      output << "max[w(" << b << ")]:\t" << stats.max_wl[b] << std::endl;
    for (size_t b = 0; b < stats.max_wu.size(); ++b)
      output << "max[w(" << b << ")]:\t" << stats.max_wu[b] << std::endl;
    for (size_t l = 0; l < stats.max_f.size(); ++l)
      output << "max[S2(" << l << ")]:\t" << stats.max_f[l] << std::endl;
    output << std::endl;
  }
};
}  // namespace rshc::eval

#endif  // RSHC_EVALUATION_STATISTICS_H_