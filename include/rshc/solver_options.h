#ifndef RSHC_SOLVER_OPTIONS_H_
#define RSHC_SOLVER_OPTIONS_H_

namespace rshc {

/**
 * @brief Options for solver.
 *
 * Solver parameters that can be modified. Note: not all parameters will apply
 * to all solvers.
 */
struct SolverOptions {
  double tolerance = 1e-8;  ///< Solver tolerance
  double obj_scale = 1;     ///< Scaling parameter for objective heuristic
  bool verbose = false;     ///< Specify whether to output log
  bool save = false;        ///< Save the resulting variables on solution
};

}  // namespace rshc

#endif  // RSHC_SOLVER_OPTIONS_H_