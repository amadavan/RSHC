#ifndef RSHC_ACCEPTABILITY_ASSESSOR_H_
#define RSHC_ACCEPTABILITY_ASSESSOR_H_

#include <fusion.h>

#include <numeric>
#include <vector>

#include "../model.h"
#include "unacceptability_check.h"

namespace rshc::acc {
class Assessor {
 public:
  /**
   * @brief Assesses whether a hosting capacity is feasible.
   *
   * Solves the mosek model to evaluate whether a hosting capacity is feasible.
   * Constructs the infeasibility certificate if infeasible and uses it to
   * construct a dual cut that can be used to verify infeasibility.
   *
   * @param model Network model.
   */
  Assessor(const Model &model);
  ~Assessor();

  /**
   * @brief Verify whether a hosting capacity is acceptable.
   *
   * Solves the mosek model and evaluates its feasibility.
   *
   * @param psi Hosting capacity to verify.
   * @return true Hosting capacity is acceptable.
   * @return false Hosting capacity is unacceptable.
   */
  bool verify(const std::vector<double> &psi);

  /**
   * @brief Get a infeasibility certificate.
   *
   * Returns a set of constraints, whose violation guarantees that the hosting
   * capacity is unacceptable.
   *
   * @return InfeasibilityCut Dual infeasibility certificate.
   */
  InfeasibilityCut getInfeasibilityCut();

 private:
  const Model &network_;  ///< Network model

  mosek::fusion::Model::t model_;  ///< Mosek model

  mosek::fusion::Parameter::t psi_;  ///< Hosting capacity

  mosek::fusion::Variable::t w_l_;  ///< VaR of voltage lower limit
  mosek::fusion::Variable::t w_u_;  ///< VaR of voltage upper limit
  mosek::fusion::Variable::t s_;    ///< VaR of line capacity limit

  mosek::fusion::Variable::t t_l_;  ///< Lower voltage limit epigraph
  mosek::fusion::Variable::t t_u_;  ///< Upper voltage limit epigraph
  mosek::fusion::Variable::t t_s_;  ///< Line capacity limit epigraph

  struct ScenarioVars {
    mosek::fusion::Variable::t P_;  ///< Real power flows
    mosek::fusion::Variable::t Q_;  ///< Reactive power flows
    mosek::fusion::Variable::t w_;  ///< Squared bus voltages
    mosek::fusion::Variable::t l_;  ///< Line
  };

  std::vector<ScenarioVars> k_vars_;  ///< List of scenario variables

  mosek::fusion::Constraint::t c_v_l_;  ///< CVaR voltage lower limit
  mosek::fusion::Constraint::t c_v_u_;  ///< CVaR voltage upper limit
  mosek::fusion::Constraint::t c_f_;    ///< CVaR flow capacity limit

  struct ScenarioCons {
    mosek::fusion::Constraint::t v_ref_;  ///< Reference (feeder) bus voltage
    mosek::fusion::Constraint::t pb_p_;   ///< Real power balance constraint
    mosek::fusion::Constraint::t pb_q_;   ///< Reactive power balance constraint
    mosek::fusion::Constraint::t v_b_;    ///< Bus voltage balance
    mosek::fusion::Constraint::t flow_;   ///< SOC line flow constraint
    mosek::fusion::Constraint::t v_l_;    ///< Lower voltage limit
    mosek::fusion::Constraint::t v_u_;    ///< Upper voltage limit
    mosek::fusion::Constraint::t f_;      ///< Flow capacity limit
    mosek::fusion::Constraint::t f_c_;    ///< Flow capacity limit (constants)
  };

  std::vector<ScenarioCons> k_cons_;  ///< List of scenario constraints

  std::shared_ptr<monty::ndarray<double, 1>> eta_;      ///< DER pf
  std::shared_ptr<monty::ndarray<double, 1>> n_w_min_;  ///< Negative min w
  std::shared_ptr<monty::ndarray<double, 1>> wl_min_;   ///< Min wl
  std::shared_ptr<monty::ndarray<double, 1>> w_max_;    ///< Max w
  std::shared_ptr<monty::ndarray<double, 1>> s_lim_;    ///< Flow limit

  size_t n_b_;  ///< Number of buses
  size_t n_l_;  ///< Number of lines
  size_t n_d_;  ///< Number of DERs
  size_t n_k_;  ///< Number of scenarios
};
}  // namespace rshc::acc

#endif  // RSHC_ACCEPTABILITY_ASSESSOR_H_