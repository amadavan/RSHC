#include <rshc/model.h>

void rshc::Model::load_case(const std::string &casefile, bool perunit) {
  phasor::matpower::Parser network_parser = phasor::matpower::Parser(casefile);
  phasor::matpower::Network network = network_parser.parseMatpower();
  network.toInternalOrder();

  n_b = network.bus.getCol(0).size();
  n_l = network.branch.getCol(0).size();

  wl = rshc::util::getEigenSTL<double>(
      network.bus[phasor::matpower::BusIndex::VMIN].cwiseProduct(
          network.bus[phasor::matpower::BusIndex::VMIN]));
  wu = rshc::util::getEigenSTL<double>(
      network.bus[phasor::matpower::BusIndex::VMAX].cwiseProduct(
          network.bus[phasor::matpower::BusIndex::VMAX]));

  if (perunit) {
    double baseS = network.baseMVA * 1e6;
    double baseV = network.bus[phasor::matpower::BusIndex::BASE_KV][0] * 1e3;
    r = rshc::util::getEigenSTL<double>(
        network.branch[phasor::matpower::BranchIndex::BR_R] * baseS /
        (baseV * baseV));
    x = rshc::util::getEigenSTL<double>(
        network.branch[phasor::matpower::BranchIndex::BR_X] * baseS /
        (baseV * baseV));
  } else {
    r = rshc::util::getEigenSTL<double>(
        network.branch[phasor::matpower::BranchIndex::BR_R]);
    x = rshc::util::getEigenSTL<double>(
        network.branch[phasor::matpower::BranchIndex::BR_X]);
  }
  S2 = rshc::util::getEigenSTL<double>(
      network.branch[phasor::matpower::BranchIndex::RATE_A].cwiseProduct(
          network.branch[phasor::matpower::BranchIndex::RATE_A]));
  from_bus = rshc::util::getEigenSTL<int>(
      network.branch[phasor::matpower::BranchIndex::F_BUS]);
  to_bus = rshc::util::getEigenSTL<int>(
      network.branch[phasor::matpower::BranchIndex::T_BUS]);

  Eigen::VectorXd PD = network.bus[phasor::matpower::BusIndex::PD];
  Eigen::VectorXd QD = network.bus[phasor::matpower::BusIndex::QD];

  d = std::vector<double>(PD.data(), PD.data() + n_l);
  e = std::vector<double>(QD.data(), QD.data() + n_l);
}

int rshc::Model::load_der(const std::string &derfile) {
  std::ifstream der_stream(derfile);
  std::string line;

  der_file_ind.clear();

  der_ind.clear();
  der_cap.clear();
  eta.clear();

  while (std::getline(der_stream, line)) {
    std::istringstream ss(line);
    int bus_index;
    int pv_col_index;
    double capacity;
    double power_factor;
    char _;

    ss >> bus_index >> _ >> pv_col_index >> _ >> capacity >> _ >>
        power_factor >> _;

    der_ind.emplace_back(bus_index);
    der_cap.emplace_back(capacity);
    eta.emplace_back(sqrt(1. / pow(power_factor, 2) - 1));

    der_file_ind.emplace_back(pv_col_index);
  }

  n_d = der_ind.size();

  return 0;
}

int rshc::Model::load_scenarios(const std::string &loadfile,
                                const std::string &pvfile, int count,
                                bool randomize) {
  std::ifstream load_stream(loadfile);
  std::ifstream pv_stream(pvfile);

  if (!load_stream.good() || !pv_stream.good()) {
    std::cout << "Unable to load data files." << std::endl;
    return -1;
  }

  std::string load_line, pv_line;

  std::vector<double> load_data;
  std::vector<double> pv_data;

  load_data.reserve(56);
  pv_data.reserve(60);
  char _;

  scenarios.clear();

  while (std::getline(load_stream, load_line) &&
         std::getline(pv_stream, pv_line)) {
    std::istringstream load_line_stream(load_line);
    std::istringstream pv_line_stream(pv_line);

    for (int i = 0; i < 56; ++i) load_line_stream >> load_data[i] >> _;
    for (int i = 0; i < 60; ++i) pv_line_stream >> pv_data[i] >> _;

    std::vector<double> d =
        std::vector<double>(load_data.begin(), load_data.begin() + n_l);
    std::vector<double> e = d;

    for (size_t i = 0; i < n_l; ++i) {
      d[i] *= this->d[i];
      e[i] *= this->e[i];
    }

    std::vector<double> alpha = std::vector<double>();
    for (int index : der_file_ind) alpha.push_back(pv_data[index]);

    scenarios.emplace_back(alpha, d, e);
  }

  // Randomize order of scenarios
  if (randomize) {
    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(scenarios.begin(), scenarios.end(), g);
  }

  // Get a subset of scenarios
  if (count > 0 && count < scenarios.size())
    scenarios.erase(scenarios.begin() + count, scenarios.end());

  n_k = scenarios.size();

  return 0;
}