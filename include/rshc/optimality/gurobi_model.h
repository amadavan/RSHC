#ifndef RSHC_OPTIMALITY_GUROBI_MODEL_H_
#define RSHC_OPTIMALITY_GUROBI_MODEL_H_

#include <gurobi_c++.h>

#include <vector>

#include "../model.h"
#include "../solver_options.h"

namespace rshc::opt {
class GurobiModel {
 public:
  /**
   * @brief Return value for solver.
   *
   * Specifies the problem status and the optimal hosting capacity.
   */
  struct RetVal {
    int status;                            ///< Solver status code
    std::vector<double> hosting_capacity;  ///< Optimal hosting capacity
  };

  /**
   * @brief Construct a new Gurobi-based solver for the Optimality problem.
   *
   * Defines the system properties that will be used to construct the problem in
   * the formulate routine. This model seeks to determine the maximum allowable
   * hosting capacity for the prescribed network.
   *
   * @param model Network parameters
   * @param opts Solver options
   */
  GurobiModel(const Model &model, const SolverOptions &opts = SolverOptions());
  ~GurobiModel();

  /**
   * @brief Solve the Optimality problem.
   *
   * Solves the generated Gurobi model and return the optimal hosting capacity.
   */
  RetVal solve();

 private:
  const Model &network_;  ///< Network model

  GRBEnv env_;      ///< Gurobi environment
  GRBModel model_;  ///< Gurobi model of main problem

  GRBVar *psi_;                ///< Hosting capacity
  GRBVar *w_l_;                ///< VaR of voltage lower limit
  GRBVar *w_u_;                ///< VaR of voltage upper limit
  GRBVar *s_;                  ///< VaR of line capacity limit
  std::vector<GRBVar *> t_l_;  ///< Lower voltage limit epigraph
  std::vector<GRBVar *> t_u_;  ///< Upper voltage limit epigraph
  std::vector<GRBVar *> t_s_;  ///< Line capacity limit epigraph
  std::vector<GRBVar *> P_;    ///< Real power flows
  std::vector<GRBVar *> Q_;    ///< Reactive power flows
  std::vector<GRBVar *> w_;    ///< Squared bus voltages
  std::vector<GRBVar *> l_;    ///< Line losses

  GRBConstr *c_v_l_;  ///< CVaR voltage lower limit
  GRBConstr *c_v_u_;  ///< CVaR voltage upper limit
  GRBConstr *c_f_;    ///< CVaR flow capacity limit

  // std::vector<GRBConstr> v_ref_;  ///< Reference (feeder) bus voltage
  std::vector<GRBConstr *> pb_p_;  ///< Real power balance constraint
  std::vector<GRBConstr *> pb_q_;  ///< Reactive power balance constraint
  std::vector<GRBConstr *> v_b_;   ///< Bus voltage balance
  // std::vector<GRBConstr*> flow_;   ///< SOC line flow constraint
  std::vector<GRBConstr *> v_l_;  ///< Lower voltage limit
  std::vector<GRBConstr *> v_u_;  ///< Upper voltage limit
  // std::vector<GRBConstr*> f_;      ///< Flow capacity limit

  size_t n_b_;  ///< Number of buses
  size_t n_l_;  ///< Number of lines
  size_t n_d_;  ///< Number of DERs
  size_t n_k_;  ///< Number of scenarios

  const SolverOptions opts_;  ///< Solver options
};
}  // namespace rshc::opt

#endif  // RSHC_OPTIMALITY_GUROBI_MODEL_H_