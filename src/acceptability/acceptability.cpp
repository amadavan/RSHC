#include <phasor/phasor.h>
#include <rshc/acceptability/acceptability_check.h>
#include <rshc/acceptability/assessor.h>
#include <rshc/acceptability/unacceptability_check.h>
#include <rshc/model.h>
#include <rshc/timer.h>
#include <rshc/util.h>

#include <cxxopts.hpp>
#include <regex>

#ifdef ENABLE_FS
#include <filesystem>
#endif

#ifdef ENABLE_GUROBI
#endif

#include "tools/cpp/runfiles/runfiles.h"
using bazel::tools::cpp::runfiles::Runfiles;

int main(int argc, char **argv) {
  std::string error;
  std::unique_ptr<Runfiles> runfiles(Runfiles::Create(argv[0], &error));
  if (runfiles == nullptr) {
    std::cout << "Runfiles not found." << std::endl;
  }

  cxxopts::Options options("rshca",
                           "Risk-sensitive hosting capacity analysis.");

  options.add_options()
      ("n,network", "Network name",cxxopts::value<std::string>()->default_value("case56_sce.m"))
      ("d,der", "DER allocation", cxxopts::value<std::string>()->default_value("case56_der.csv"))
      ("l,load", "Load data",cxxopts::value<std::string>()->default_value("load_data.csv"))
      ("p,pv", "Solar PV data", cxxopts::value<std::string>()->default_value("pv_data.csv"))
      ("v,nu", "Voltage limit CVaR parameter",cxxopts::value<double>()->default_value("0"))
      ("g,gamma", "Line flow capacity CVaR parameter",cxxopts::value<double>()->default_value("0"))
      ("f,file", "File of DERs to test",cxxopts::value<std::string>()->default_value(runfiles->Rlocation("rshc/data/psis_case56.csv")))
      ("k,scenarios", "Scenarios", cxxopts::value<int>()->default_value("0"))
      ("u,unit", "Per-unit rescaling", cxxopts::value<bool>()->default_value("false"))
  ;

  auto result = options.parse(argc, argv);

  std::string filestring = result["network"].as<std::string>();
  std::string derstring = result["der"].as<std::string>();
  std::string loadstring = result["load"].as<std::string>();
  std::string pvstring = result["pv"].as<std::string>();
  std::string psistring = result["file"].as<std::string>();
  double nu = result["nu"].as<double>();
  double gamma = result["gamma"].as<double>();
  int scenarios = result["scenarios"].as<int>();
  bool perunit = result["unit"].as<bool>();

  std::string case_file = runfiles->Rlocation("rshc/data/" + filestring);
  std::string der_file = runfiles->Rlocation("rshc/data/" + derstring);
  std::string load_file = runfiles->Rlocation("rshc/data/" + loadstring);
  std::string pv_file = runfiles->Rlocation("rshc/data/" + pvstring);
  std::string psi_file = runfiles->Rlocation("rshc/data/" + psistring);
#ifdef ENABLE_FS
  if (!std::filesystem::exists(case_file)) case_file = filestring;
  if (!std::filesystem::exists(der_file)) der_file = derstring;
  if (!std::filesystem::exists(load_file)) load_file = loadstring;
  if (!std::filesystem::exists(pv_file)) pv_file = pvstring;
  if (!std::filesystem::exists(psi_file)) psi_file = psistring;
#endif

  rshc::Model model;
  model.load_case(case_file, perunit);
  model.load_der(der_file);
  model.load_scenarios(load_file, pv_file, scenarios);

  model.nu = nu;
  model.gamma = gamma;

  // Define output data files based on the casename
  std::string casename = "case";
  std::smatch m;
  if (std::regex_search(case_file, m, std::regex("(case[0-9]+.*).m")))
    casename = m[1];
  std::string data_filename = casename + std::to_string((int)(nu * 100)) + "_" +
                              std::to_string((int)(gamma * 100)) + ".h5";
  std::string run_filename = casename + std::to_string((int)(nu * 100)) + "_" +
                             std::to_string((int)(gamma * 100)) + ".csv";

  // Set up the test objects
  // #ifdef ENABLE_GUROBI
  //   std::cout << "Using gurobi feasibility test" << std::endl;
  //   rshc::GRBFeasibilityTest feasibility_test(model, "");
  // #else
  std::cout << "Using mosek feasibility test" << std::endl;
  rshc::acc::AcceptabilityCheck acceptability_check(model, "");
  // #endif
  rshc::acc::UnacceptabilityCheck unacceptability_check(model, "");
  rshc::acc::Assessor assessor(model);

  // Load the data file and iterate over it
  rshc::Timer timer;

  std::ifstream psi_stream(psi_file);

  std::ofstream output_stream(run_filename);
  char _;
  int count = 0;

  if (psi_stream.good()) {
    std::string psi_line;
    std::vector<double> psi(model.n_d);
    while (std::getline(psi_stream, psi_line)) {
      std::istringstream psi_line_stream(psi_line);
      for (int i = 0; i < model.n_d; ++i) psi_line_stream >> psi[i] >> _;

      // Check for feasibility
      bool isAcceptable = true;
      bool assessed = false;
      timer.start();
      if (unacceptability_check.verify(psi)) {
        isAcceptable = false;
      } else if (acceptability_check.verify(psi)) {
        isAcceptable = true;
      } else {
        assessed = true;
        isAcceptable = assessor.verify(psi);
        if (isAcceptable)
          acceptability_check.addCapacity(psi);
        else
          unacceptability_check.addCapacity(assessor.getInfeasibilityCut());
      }
      output_stream << timer.elapsed().count() << "," << isAcceptable << ","
                    << assessed << std::endl;
      std::cout << count++ << "\t" << isAcceptable << "\t" << assessed << std::endl;
    }

    std::cout << "Saving feasibility test" << std::endl;
    acceptability_check.save(data_filename);
    std::cout << "Saving infeasibility test" << std::endl;
    unacceptability_check.save(data_filename);
    std::cout << "[DONE]" << std::endl;
  } else {
    std::cout << "Unable to read file." << std::endl;
  }
}