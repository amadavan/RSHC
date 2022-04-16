#ifndef RSHC_MODEL_H_
#define RSHC_MODEL_H_

#include <phasor/phasor.h>

#include <fstream>
#include <string>
#include <vector>
#include <random>

#include "util.h"

namespace rshc {

struct Scenario {
  std::vector<double> alpha;  ///< DER scaling
  std::vector<double> d;      ///< Real demand (absent 0-bus)
  std::vector<double> e;      ///< Reactive demand (abset 0-bus)

  Scenario(std::vector<double> a, std::vector<double> d, std::vector<double> e)
      : alpha(a), d(d), e(e) {}
};

struct Model {
  std::vector<double> wl;  ///< Lower voltage squared limit
  std::vector<double> wu;  ///< Upper voltage squared limit

  std::vector<double> r;      ///< Line resistance (p.u.)
  std::vector<double> x;      ///< Line reactances (p.u.)
  std::vector<double> S2;     ///< Squared line flow capacity limit
  std::vector<int> from_bus;  ///< Line sending end index
  std::vector<int> to_bus;    ///< Line receiving end index

  std::vector<double> der_cap;  ///< DER installation capacity
  std::vector<double> eta;      ///< DER reactive power scale
  std::vector<int> der_ind;     ///< Bus index of DER installation

  std::vector<Scenario> scenarios;  ///< Scenarios to consider

  double nu;     ///< Voltage limit CVaR parameter
  double gamma;  ///< Line capacity CVaR parameter

  size_t n_b;  ///< Number of buses
  size_t n_l;  ///< Number of lines
  size_t n_d;  ///< Number of DER installations
  size_t n_k;  ///< Number of scenarios

  // Properties required for data loading
  std::vector<int> der_file_ind;  ///< Column index of DER in scenario file
  std::vector<double> d;          ///< Nodal real demand (p.u.)
  std::vector<double> e;          ///< Nodal reactive demand (p.u.)

  /**
   * @brief Load network case file.
   *
   * Populate the network parameters based on the specified matpower data file.
   * Also converts line parmaeters to per unit if necessary.
   *
   * @param casefile Matpower network case filename.
   * @param perunit Whether to convert line parameters to per unit (default
   * false).
   */
  void load_case(const std::string &casefile, bool perunit = false);

  /**
   * @brief Load DER data file.
   *
   * Reads the specified DER allocation file. Contains information regarding
   * location, maximum allowable value, and the power factor.
   *
   * @param derfile CSV file containing DER data.
   *
   * @return int 0 for success, -1 for failure
   */
  int load_der(const std::string &derfile);

  /**
   * @brief Load scenario (load/pv) data.
   *
   * Reads the two files and populates the scenario information.
   *
   * NOTE: must be called after specifying the DER allocation variable der_ind.
   *
   * @param loadfile CSV file containing load information.
   * @param pvfile CSV file containing PV data.
   * @param count Number of scenarios to include (0 = all)
   * @param randomize Whether to randomize the order of the scenarios
   *
   * @return int 0 for success, -1 for failure
   */
  int load_scenarios(const std::string &loadfile, const std::string &pvfile, int count = 0, bool randomize = false);
};

}  // namespace rshc

#endif  // RSHC_MODEL_H_
