#ifndef RSHC_DETERMINISTIC_MOSEK_DETERMINISTIC_H_
#define RSHC_DETERMINISTIC_MOSEK_DETERMINISTIC_H_

#include <fusion.h>

#include <algorithm>
#include <memory>
#include <vector>

#include "../model.h"
#include "../solver_options.h"

namespace rshc::det {
class MosekDeterministic {
 public:
  /**
   * @brief Return value for solver.
   *
   * Specifies the problem status and the optimal hosting capacity.
   */
  struct RetVal {
    mosek::fusion::ProblemStatus status;   ///< Solver status code.
    std::vector<double> hosting_capacity;  ///< Optimal hosting capacity.
  };

  MosekDeterministic(const Model &model,
                     const SolverOptions &opts = SolverOptions());
  ~MosekDeterministic();

  RetVal solve();

 private:
  const Model &network_;  ///< Network parameters

  mosek::fusion::Model::t model_;  ///< Mosek model

  mosek::fusion::Variable::t psi_;  ///< Hosting capacity

  mosek::fusion::Variable::t w_l_;  ///< VaR of voltage lower limit
  mosek::fusion::Variable::t w_u_;  ///< VaR of voltage upper limit
  mosek::fusion::Variable::t s_;    ///< VaR of line capacity limit

  struct ScenarioVars {
    mosek::fusion::Variable::t P_;  ///< Real power flows
    mosek::fusion::Variable::t Q_;  ///< Reactive power flows
    mosek::fusion::Variable::t w_;  ///< Squared bus voltages
    mosek::fusion::Variable::t l_;  ///< Line
  };

  std::vector<ScenarioVars> k_vars_;  ///< List of scenario variables

  struct ScenarioCons {
    mosek::fusion::Constraint::t v_ref_;  ///< Reference (feeder) bus voltage
    mosek::fusion::Constraint::t pb_p_;   ///< Real power balance constraint
    mosek::fusion::Constraint::t pb_q_;   ///< Reactive power balance constraint
    mosek::fusion::Constraint::t v_b_;    ///< Bus voltage balance
    mosek::fusion::Constraint::t flow_;   ///< SOC line flow constraint
    mosek::fusion::Constraint::t v_l_;    ///< Lower voltage limit
    mosek::fusion::Constraint::t v_u_;    ///< Upper voltage limit
    mosek::fusion::Constraint::t f_;      ///< Flow capacity limit
  };

  std::vector<ScenarioCons> k_cons_;  ///< List of scenario constraints

  size_t n_b_;  ///< Number of buses
  size_t n_l_;  ///< Number of lines
  size_t n_d_;  ///< Number of DERs
  size_t n_k_;  ///< Number of scenarios

  const SolverOptions opts_;  ///< Solver options
};
}  // namespace rshc::det

#endif  // RSHC_DETERMINISTIC_MOSEK_DETERMINISTIC_H_