#include <phasor/phasor.h>
#include <rshc/model.h>
#include <rshc/optimality/mosek_model.h>
#include <rshc/timer.h>
#include <rshc/util.h>

#include <algorithm>
#include <cxxopts.hpp>

#ifdef ENABLE_FS
#include <filesystem>
#endif

#ifdef ENABLE_GUROBI
#include <rshc/optimality/gurobi_model.h>
#endif

#include "tools/cpp/runfiles/runfiles.h"
using bazel::tools::cpp::runfiles::Runfiles;

int main(int argc, char **argv) {
  std::string error;
  std::unique_ptr<Runfiles> runfiles(Runfiles::Create(argv[0], &error));
  if (runfiles == nullptr) {
    std::cout << "Runfiles not found." << std::endl;
  }

  cxxopts::Options options(
      "rshc", "Risk-sensitive maximum hosting capacity problem solver.");

  options.add_options()
      ("n,network", "Network name",cxxopts::value<std::string>()->default_value("case56_sce.m"))
      ("d,der", "DER allocation", cxxopts::value<std::string>()->default_value("case56_der.csv"))
      ("l,load", "Load data",cxxopts::value<std::string>()->default_value("load_data.csv"))
      ("p,pv", "Solar PV data", cxxopts::value<std::string>()->default_value("pv_data.csv"))
      ("v,nu", "Voltage limit CVaR parameter", cxxopts::value<double>()->default_value("0"))
      ("g,gamma", "Line flow capacity CVaR parameter", cxxopts::value<double>()->default_value("0"))
      ("k,scenarios", "Scenarios", cxxopts::value<int>()->default_value("0"))
      ("r,randomize", "Randomize scenarios", cxxopts::value<bool>()->default_value("false"))
      ("u,unit", "Per-unit rescaling of line parameters", cxxopts::value<bool>()->default_value("false"))
      ("t,tolerance", "Tolerance", cxxopts::value<double>()->default_value("1e-8"))
      ("z,scale", "Scaling factor", cxxopts::value<double>()->default_value("1"))
      ("s,solver", "Solver", cxxopts::value<std::string>()->default_value("mosek"))
      ("b,verbose", "Verbose solver output/save", cxxopts::value<bool>()->default_value("false"))
  ;

  auto result = options.parse(argc, argv);

  std::string filestring = result["network"].as<std::string>();
  std::string derstring = result["der"].as<std::string>();
  std::string loadstring = result["load"].as<std::string>();
  std::string pvstring = result["pv"].as<std::string>();
  double nu = result["nu"].as<double>();
  double gamma = result["gamma"].as<double>();
  int scenarios = result["scenarios"].as<int>();
  bool randomize = result["randomize"].as<bool>();
  bool perunit = result["unit"].as<bool>();
  double tolerance = result["tolerance"].as<double>();
  double scale = result["scale"].as<double>();
  std::string solver = result["solver"].as<std::string>();
  bool verbose = result["verbose"].as<bool>();

  std::string case_file = runfiles->Rlocation("rshc/data/" + filestring);
  std::string der_file = runfiles->Rlocation("rshc/data/" + derstring);
  std::string load_file = runfiles->Rlocation("rshc/data/" + loadstring);
  std::string pv_file = runfiles->Rlocation("rshc/data/" + pvstring);

#ifdef ENABLE_FS
  if (!std::filesystem::exists(case_file)) case_file = filestring;
  if (!std::filesystem::exists(der_file)) der_file = derstring;
  if (!std::filesystem::exists(load_file)) load_file = loadstring;
  if (!std::filesystem::exists(pv_file)) pv_file = pvstring;
#endif

  rshc::Model model;
  model.load_case(case_file, perunit);
  model.load_der(der_file);
  model.load_scenarios(load_file, pv_file, scenarios, randomize);

  model.nu = nu;
  model.gamma = gamma;

  rshc::Timer timer;
  double formulate_time;
  double solve_time;
  std::vector<double> hosting_capacity;

  rshc::SolverOptions opts;
  opts.tolerance = tolerance;
  opts.obj_scale = scale;
  opts.verbose = verbose;
  opts.save = verbose;

#ifdef ENABLE_GUROBI
  if (solver[0] == 'g' || solver[0] == 'G') {
    timer.start();
    rshc::opt::GurobiModel model_solver(model);
    formulate_time = timer.elapsed().count();

    timer.start();
    rshc::opt::GurobiModel::RetVal ret = model_solver.solve();
    solve_time = timer.elapsed().count();

    hosting_capacity = ret.hosting_capacity;

    std::cout << ret.status << std::endl;
  } else
#endif
      if (solver[0] == 'm' || solver[0] == 'M') {
    timer.start();
    rshc::opt::MosekModel model_solver(model, opts);
    formulate_time = timer.elapsed().count();

    timer.start();
    rshc::opt::MosekModel::RetVal ret = model_solver.solve();
    solve_time = timer.elapsed().count();

    hosting_capacity = ret.hosting_capacity;

    std::cout << ret.status << std::endl;
  }

  std::cout << "Hosting capacity: [";
  for (size_t i = 0; i < hosting_capacity.size() - 1; ++i)
    std::cout << hosting_capacity.at(i) << ", ";
  std::cout << hosting_capacity.at(hosting_capacity.size() - 1) << "]"
            << std::endl;

  std::cout << "Formulate time: " << formulate_time << std::endl;
  std::cout << "Solve time: " << solve_time << std::endl;
  std::cout << "Total time: " << solve_time + formulate_time << std::endl;
}