#include <rshc/acceptability/acceptability_check.h>

using namespace mosek::fusion;

rshc::acc::AcceptabilityCheck::AcceptabilityCheck(
    const Model &model, const std::string &data_file_name)
    : network_(model), data_file_name_(data_file_name), n_d_(model.eta.size()) {
  model_ = new mosek::fusion::Model();

  model_->objective(ObjectiveSense::Maximize, 0);

  psi_ = model_->parameter(n_d_);

  con_ = model_->constraint(Expr::constTerm(n_d_, 0), Domain::equalsTo(0.));
  balance_ = model_->constraint(Expr::constTerm(1, 0), Domain::equalsTo(1.));

#ifdef ENABLE_FS
  if (std::filesystem::exists(data_file_name_)) {
    H5Easy::File file(data_file_name_, H5Easy::File::ReadOnly);

    psis_ = H5Easy::load<Eigen::MatrixXd>(file, "/feasible/psis");
    alpha_.reserve(psis_.rows());
    for (size_t k = 0; k < psis_.rows(); k++) {
      auto psi = psis_.row(k);
      Variable::t alpha = model_->variable(Domain::inRange(0, 1));
      alpha_.emplace_back(alpha);

      con_->update(Expr::mulElm(std::make_shared<monty::ndarray<double, 1>>(
                                    psi.data(), psi.size()),
                                Var::vrepeat(alpha, n_d_)),
                   alpha);
      balance_->update(alpha, alpha);
    }
  } else {
#endif
    psis_ = Eigen::MatrixXd(0, n_d_);
#ifdef ENABLE_FS
  }
#endif
}

rshc::acc::AcceptabilityCheck::~AcceptabilityCheck() { model_->dispose(); }

bool rshc::acc::AcceptabilityCheck::verify(const std::vector<double> &psi) {
  std::vector<double> n_psi(psi);
  std::for_each(n_psi.begin(), n_psi.end(), [](double &d) { d *= -1; });
  con_->update(monty::new_array_ptr(n_psi));

  model_->solve();

  return (model_->getProblemStatus() == ProblemStatus::PrimalAndDualFeasible);
}

void rshc::acc::AcceptabilityCheck::addCapacity(
    const std::vector<double> &psi) {
  psis_.conservativeResize(psis_.rows() + 1, psis_.cols());
  psis_.row(psis_.rows() - 1) = Eigen::VectorXd::Map(psi.data(), psi.size());

  Variable::t alpha = model_->variable(Domain::inRange(0, 1));
  alpha_.emplace_back(alpha);

  con_->update(
      Expr::mulElm(monty::new_array_ptr(psi), Var::vrepeat(alpha, n_d_)), alpha,
      false);
  balance_->update(alpha, alpha);
}

void rshc::acc::AcceptabilityCheck::save(const std::string &filename) {
  std::string file_name = filename;
  if (filename == "") file_name = data_file_name_;

  H5Easy::File file(file_name, H5Easy::File::ReadWrite | H5Easy::File::Create);

  H5Easy::dump(file, "/feasible/psis", psis_, H5Easy::DumpMode::Overwrite);
}