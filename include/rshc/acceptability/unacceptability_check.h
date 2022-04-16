#ifndef RSHC_ACCEPTABILITY_UNACCEPTABILITY_CHECK_H_
#define RSHC_ACCEPTABILITY_UNACCEPTABILITY_CHECK_H_

#include <hdf5.h>

#include <Eigen/Core>
#include <Eigen/StdVector>

#ifdef ENABLE_FS
#include <filesystem>
#endif

#include <highfive/H5Easy.hpp>
#include <string>

#include "../model.h"

namespace rshc::acc {

/**
 * @brief Infeasibility certificate.
 *
 * Constraint describing the infeasibility certificate for  a specified hosting
 * capacity. Violation of the constraint implies that the hosting capacity is
 * infeasible, where violation refers to A * psi < b.
 */
struct InfeasibilityCut {
  Eigen::VectorXd A;  ///< Constraint coefficients
  double b;           ///< RHS of constraint
};

class UnacceptabilityCheck {
 public:
  /**
   * @brief  Construct a new Unacceptability Check object.
   *
   * Creates a class that verifies whether hosting capacities are unacceptable
   * using infeasibility certificates of prior hosting capacities for a
   * specified model.
   *
   * @param model Network model on which capacities are evaluated.
   * @param data_file_name HDF5 file containing infeasibility certificates.
   */
  UnacceptabilityCheck(const Model &model, const std::string &data_file_name);
  ~UnacceptabilityCheck();

  /**
   * @brief Verify whether capacity is unacceptable.
   *
   * Uses the prior infeasibility certificates of unacceptable hosting
   * capacities to evaluate whether the capacity violates any of the specified
   * constraints.
   *
   * @param psi Hosting capacity to verify.
   * @return true Capacity is unacceptable.
   * @return false Capacity is acceptable.
   */
  bool verify(const std::vector<double> &psi);

  /**
   * @brief Add infeasibility certificate.
   *
   * Adds an infeasibility certificate for an unacceptable hosting capacity.
   *
   * @param cut Infeasibility certificate.
   */
  void addCapacity(const InfeasibilityCut &cut);

  /**
   * @brief Save data to file.
   *
   * Save set of unacceptable certificates to a HDF5 file.
   *
   * @param filename Save file name.
   */
  void save(const std::string &filename);

 private:
  const Model &network_;               ///< Network model to initialize sizes
  const std::string &data_file_name_;  ///< Data file from which data is read

  Eigen::MatrixXd A_;  ///< Matrix of cut coefficients
  Eigen::VectorXd b_;  ///< Vector of cut RHS
};
}  // namespace rshc::acc

#endif  // RSHC_ACCEPTABILITY_UNACCEPTABILITY_CHECK_H_