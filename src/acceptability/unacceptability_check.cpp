#include <rshc/acceptability/unacceptability_check.h>

rshc::acc::UnacceptabilityCheck::UnacceptabilityCheck(const Model &model, const std::string &data_file_name)
    : network_(model), data_file_name_(data_file_name) {
#ifdef ENABLE_FS
  if (std::filesystem::exists(data_file_name_)) {
    H5Easy::File file(data_file_name_, H5Easy::File::ReadOnly);

    A_ = H5Easy::load<Eigen::MatrixXd>(file, "/infeasible/A");
    b_ = H5Easy::load<Eigen::VectorXd>(file, "/infeasible/b");
  } else {
#endif
    A_ = Eigen::MatrixXd(0, network_.eta.size());
    b_ = Eigen::VectorXd(0);
#ifdef ENABLE_FS
  }
#endif
}

rshc::acc::UnacceptabilityCheck::~UnacceptabilityCheck() {}

bool rshc::acc::UnacceptabilityCheck::verify(const std::vector<double> &psi) {
  if (b_.size() == 0) return false;
  
  Eigen::VectorXd _psi = Eigen::VectorXd::Map(psi.data(), psi.size());

  auto diff = b_ - A_ * _psi;

  return (diff.array() < 0).any();
}

void rshc::acc::UnacceptabilityCheck::addCapacity(const InfeasibilityCut &cut) {
  A_.conservativeResize(A_.rows() + 1, A_.cols());
  b_.conservativeResize(b_.size() + 1);

  A_.row(A_.rows() - 1) = cut.A;
  b_.coeffRef(b_.size() - 1) = cut.b;
}

void rshc::acc::UnacceptabilityCheck::save(const std::string &filename) {
  std::string file_name = filename;
  if (filename == "") file_name = data_file_name_;

  H5Easy::File file(file_name, H5Easy::File::ReadWrite | H5Easy::File::Create);

  H5Easy::dump(file, "/infeasible/A", A_, H5Easy::DumpMode::Overwrite);
  H5Easy::dump(file, "/infeasible/b", b_, H5Easy::DumpMode::Overwrite);
}