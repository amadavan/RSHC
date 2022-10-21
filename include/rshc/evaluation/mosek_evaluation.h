#ifndef RSHC_EVALUATION_MOSEK_EVALUATION_H_
#define RSHC_EVALUATION_MOSEK_EVALUATION_H_

#include <fusion.h>

#include <algorithm>
#include <numeric>
#include <vector>

#include <Eigen/Core>

#include "../model.h"
#include "../solver_options.h"
#include "../util.h"
#include "statistics.h"

namespace rshc::eval {
class MosekEvaluation {
 public:
  MosekEvaluation(const Model &model,
                  const SolverOptions &opts = SolverOptions());
  ~MosekEvaluation();

  Statistics evaluate(const std::vector<double> &psi);

 private:
  const Model &network_;  ///< Network parameters

  mosek::fusion::Model::t model_;  ///< Mosek model

  mosek::fusion::Parameter::t psi_;  ///< Hosting capacity

  mosek::fusion::Parameter::t alpha_;  ///< Solar irradiance profile
  mosek::fusion::Parameter::t d_;      ///< Real power demand
  mosek::fusion::Parameter::t e_;      ///< Reactive power demand

  mosek::fusion::Variable::t P_;  ///< Real power flows
  mosek::fusion::Variable::t Q_;  ///< Reactive power flows
  mosek::fusion::Variable::t w_;  ///< Squared voltage magnitudes
  mosek::fusion::Variable::t l_;  ///< Line losses

  mosek::fusion::Variable::t wl_slack_;  ///< Slack for lower limit
  mosek::fusion::Variable::t wu_slack_;  ///< Slack for upper limit
  mosek::fusion::Variable::t f_slack_;   ///< Slack for line flow

  mosek::fusion::Constraint::t v_ref_;  ///< Reference (feeder) bus voltage
  mosek::fusion::Constraint::t pb_p_;   ///< Real power balance constraint
  mosek::fusion::Constraint::t pb_q_;   ///< Reactive power balance constraint
  mosek::fusion::Constraint::t v_b_;    ///< Bus voltage balance
  mosek::fusion::Constraint::t flow_;   ///< SOC line flow constraint
  mosek::fusion::Constraint::t v_l_;    ///< Lower voltage limit
  mosek::fusion::Constraint::t v_u_;    ///< Upper voltage limit
  mosek::fusion::Constraint::t f_;      ///< Flow capacity limit

  size_t n_b_;  ///< Number of buses
  size_t n_l_;  ///< Number of lines
  size_t n_d_;  ///< Number of DERs
  size_t n_k_;  ///< Number of scenarios

  const SolverOptions opts_;  ///< Solver options
};
}  // namespace rshc::eval

#endif  // RSHC_EVALUATION_MOSEK_EVALUATION_H_