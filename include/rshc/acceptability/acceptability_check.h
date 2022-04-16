#ifndef RSHC_ACCEPTABILITY_ACCEPTABILITY_CHECK_H_
#define RSHC_ACCEPTABILITY_ACCEPTABILITY_CHECK_H_

#include <fusion.h>
#include <hdf5.h>

#include <string>

#ifdef ENABLE_FS
#include <filesystem>
#endif

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <highfive/H5Easy.hpp>

#include "../model.h"

namespace rshc::acc {
class AcceptabilityCheck {
 public:
  /**
   * @brief Construct a new Acceptability Check object.
   *
   * Creates a class that verifies whether hosting capacities are acceptable
   * using prior acceptable hosting capacities for a specified model.
   *
   * @param model Network model on which capacities are evaluated
   * @param data_file_name HDF5 file containing feasible psis
   */
  AcceptabilityCheck(const Model &model, const std::string &data_file_name);
  ~AcceptabilityCheck();

  /**
   * @brief Verify whether capacity is acceptable.
   *
   * Uses the prior acceptable hosting capacities to evaluate whether the
   * solution lies within the convex hull.
   *
   * @param psi Hosting capacity to verify.
   * @return true Capacity is acceptable.
   * @return false Capacity is unacceptable.
   */
  bool verify(const std::vector<double> &psi);

  /**
   * @brief Add acceptable capacity.
   *
   * Expand the state information to include the following capacity as
   * acceptable.
   *
   * @param psi Acceptable hosting capacity.
   */
  void addCapacity(const std::vector<double> &psi);

  /**
   * @brief Save data to file.
   *
   * Save set of feasible hosting capacities to a HDF5 file.
   *
   * @param filename Save file name.
   */
  void save(const std::string &filename = "");

 private:
  const Model &network_;               ///< Network model to initialize sizes
  const std::string &data_file_name_;  ///< Data file to load data

  mosek::fusion::Model::t model_;  ///< Mosek model

  std::vector<mosek::fusion::Variable::t> alpha_;  ///< Combination weights
  mosek::fusion::Parameter::t psi_;                ///< Test capacity

  mosek::fusion::Constraint::t con_;  ///< Determine that psi can be recreated
  mosek::fusion::Constraint::t balance_;  ///< Enforce combination is convex

  Eigen::MatrixXd psis_;  ///< Matrix of feasible psis

  const size_t n_d_;  ///< Number of DER isntallations
};
}  // namespace rshc::acc

#endif  // RSHC_ACCEPTABILITY_ACCEPTABILITY_CHECK_H_