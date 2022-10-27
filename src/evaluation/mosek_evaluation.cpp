#include <rshc/evaluation/mosek_evaluation.h>

using namespace mosek::fusion;
typedef std::shared_ptr<monty::ndarray<double, 1>> ndarray_ptr_t;

rshc::eval::MosekEvaluation::MosekEvaluation(const Model &model,
                                             const SolverOptions &opts)
    : network_(model), opts_(opts) {
  model_ = new mosek::fusion::Model();

  n_b_ = network_.wl.size();
  n_l_ = network_.r.size();
  n_d_ = network_.der_ind.size();
  n_k_ = network_.scenarios.size();

  assert(n_b_ == n_l_ + 1);

  /*******************************************************************
   * Solver parameters
   ******************************************************************/
  model_->acceptedSolutionStatus(AccSolutionStatus::Anything);
  model_->setSolverParam("intpntCoTolDfeas", opts_.tolerance);
  model_->setSolverParam("intpntCoTolInfeas", 1e-12);
  model_->setSolverParam("intpntCoTolMuRed", opts_.tolerance);
  model_->setSolverParam("intpntCoTolNearRel", 1000);
  model_->setSolverParam("intpntCoTolPfeas", opts_.tolerance);
  model_->setSolverParam("intpntCoTolRelGap", opts_.tolerance);
  if (opts_.verbose)
    model_->setLogHandler(
        [=](const std::string &msg) { std::cout << msg << std::flush; });

  /*******************************************************************
   * Useful arrays and matrices
   ******************************************************************/
  ndarray_ptr_t psi_max = monty::new_array_ptr(network_.der_cap);

  // Note: we ignore voltage of feeder bus in CVaR constraints
  ndarray_ptr_t w_min = monty::new_array_ptr(
      std::vector<double>(network_.wl.begin(), network_.wl.end()));
  ndarray_ptr_t w_max = monty::new_array_ptr(
      std::vector<double>(network_.wu.begin(), network_.wu.end()));
  ndarray_ptr_t s_lim = monty::new_array_ptr(network_.S2);

  ndarray_ptr_t r = monty::new_array_ptr(network_.r);
  ndarray_ptr_t x = monty::new_array_ptr(network_.x);
  ndarray_ptr_t r2_x2 =
      std::make_shared<monty::ndarray<double, 1>>(n_l_, [this](ptrdiff_t i) {
        return network_.r[i] * network_.r[i] + network_.x[i] * network_.x[i];
      });

  ndarray_ptr_t eta = monty::new_array_ptr(network_.eta);

  // Adjacency matrices
  std::vector<int> der_adj_col(n_d_);
  std::iota(der_adj_col.begin(), der_adj_col.end(), 0);
  Matrix::t der_loc =
      Matrix::sparse(n_l_, n_d_, monty::new_array_ptr(network_.der_ind),
                     monty::new_array_ptr(der_adj_col),
                     monty::new_array_ptr(std::vector<double>(n_d_, 1.)));

  std::vector<int> adj_pos_row;
  std::vector<int> adj_pos_col;
  std::vector<double> adj_pos_val;

  std::vector<int> adj_row;
  std::vector<int> adj_col;
  std::vector<double> adj_val;
  for (size_t l = 0; l < n_l_; l++) {
    adj_row.emplace_back(l);
    adj_col.emplace_back(network_.from_bus[l] - 1);
    adj_val.emplace_back(1);

    adj_pos_row.emplace_back(l);
    adj_pos_col.emplace_back(network_.from_bus[l] - 1);
    adj_pos_val.emplace_back(1);

    adj_row.emplace_back(l);
    adj_col.emplace_back(network_.to_bus[l] - 1);
    adj_val.emplace_back(-1);
  }
  Matrix::t adj_pos = Matrix::sparse(
      n_l_, n_b_, monty::new_array_ptr(adj_pos_row),
      monty::new_array_ptr(adj_pos_col), monty::new_array_ptr(adj_pos_val));
  Matrix::t adj = Matrix::sparse(n_l_, n_b_, monty::new_array_ptr(adj_row),
                                 monty::new_array_ptr(adj_col),
                                 monty::new_array_ptr(adj_val));

  // Assumes that the first index corresponds to the line from the feeder
  std::vector<int> pi_adj_row(adj_row.begin() + 1, adj_row.end());
  std::vector<int> pi_adj_col(adj_col.begin() + 1, adj_col.end());
  std::vector<double> pi_adj_val(adj_val.begin() + 1, adj_val.end());
  Matrix::t pi_adj = Matrix::sparse(
      n_l_, n_l_, monty::new_array_ptr(pi_adj_row),
      monty::new_array_ptr(pi_adj_col), monty::new_array_ptr(pi_adj_val));

  /*******************************************************************
   * Build the model
   ******************************************************************/
  // Variables
  psi_ = model_->parameter(n_d_);

  alpha_ = model_->parameter(n_d_);
  d_ = model_->parameter(n_l_);
  e_ = model_->parameter(n_l_);

  P_ = model_->variable(n_l_, Domain::unbounded());
  Q_ = model_->variable(n_l_, Domain::unbounded());
  w_ = model_->variable(n_b_, Domain::greaterThan(0.));
  l_ = model_->variable(n_l_, Domain::greaterThan(0.));

  wl_slack_ = model_->variable(n_b_, Domain::greaterThan(0.));
  wu_slack_ = model_->variable(n_b_, Domain::greaterThan(0.));
  f_slack_ = model_->variable(n_l_, Domain::greaterThan(0.));

  // Set the objective
  model_->objective(
      ObjectiveSense::Minimize,
      Expr::add(Expr::sum(f_slack_),
                Expr::add(Expr::sum(wl_slack_), Expr::sum(wu_slack_))));

  // Constraints
  v_ref_ = model_->constraint(w_->index(n_b_ - 1), Domain::equalsTo(1.));
  pb_p_ = model_->constraint(
      Expr::sub(Expr::mul(der_loc, Expr::mulElm((Expression::t)alpha_, psi_)),
                Expr::add(Expr::mul(pi_adj->transpose(), P_),
                          Expr::add(Expr::mulElm(r, l_), d_))),
      Domain::equalsTo(0.));
  pb_q_ = model_->constraint(
      Expr::sub(Expr::mul(der_loc,
                          Expr::mulElm(
                              eta, Expr::mulElm((Expression::t)alpha_, psi_))),
                Expr::add(Expr::mul(pi_adj->transpose(), Q_),
                          Expr::add(Expr::mulElm(x, l_), e_))),
      Domain::equalsTo(0.));
  v_b_ = model_->constraint(
      Expr::sub(
          Expr::add(Expr::mul(adj, w_), Expr::mulElm(r2_x2, l_)),
          Expr::mul(2, Expr::add(Expr::mulElm(r, P_), Expr::mulElm(x, Q_)))),
      Domain::equalsTo(0.));

  std::shared_ptr<monty::ndarray<Expression::t, 1>> flow_p =
      monty::new_array_ptr<Expression::t, 1>(
          {Expr::mul(0.5, Expr::mul(adj_pos, w_)), l_, P_, Q_});
  flow_ = model_->constraint(Expr::hstack(flow_p), Domain::inRotatedQCone());

  v_l_ =
      model_->constraint(Expr::add(w_, wl_slack_), Domain::greaterThan(w_min));
  v_u_ = model_->constraint(Expr::sub(w_, wu_slack_), Domain::lessThan(w_max));
  std::shared_ptr<monty::ndarray<Expression::t, 1>> f_p =
      monty::new_array_ptr<Expression::t, 1>(
          {Expr::constTerm(n_l_, 0.5), Expr::add(s_lim, f_slack_), P_, Q_});
  f_ = model_->constraint(Expr::hstack(f_p), Domain::inRotatedQCone());
}

rshc::eval::MosekEvaluation::~MosekEvaluation() { model_->dispose(); }

rshc::eval::Statistics rshc::eval::MosekEvaluation::evaluate(
    const std::vector<double> &psi) {
  psi_->setValue(monty::new_array_ptr(psi));

  Statistics stats;
  stats.P_joint = 0;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> w =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                    Eigen::RowMajor>::Zero(n_b_, n_k_);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> S2 =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                    Eigen::RowMajor>::Zero(n_l_, n_k_);

  util::ProgressBar progressBar(n_k_);
  for (size_t k = 0; k < n_k_; ++k) {
    progressBar.setValue(k);
    
    const Scenario &scenario = network_.scenarios[k];
    alpha_->setValue(monty::new_array_ptr(scenario.alpha));
    d_->setValue(monty::new_array_ptr(scenario.d));
    e_->setValue(monty::new_array_ptr(scenario.e));

    model_->solve();

    // Verify feasibility
    if (model_->getProblemStatus() != ProblemStatus::PrimalAndDualFeasible) {
      std::cout << "Houston, we have a problem." << std::endl;
    }

    if (model_->primalObjValue() != 0) {
      // Problem was infeasible
      stats.P_joint += 1. / k;
    }

    ndarray_ptr_t w_val = w_->level();
    ndarray_ptr_t P_val = P_->level();
    ndarray_ptr_t Q_val = Q_->level();

    Eigen::ArrayXd P = Eigen::Map<Eigen::ArrayXd>(P_val->begin(), n_l_);
    Eigen::ArrayXd Q = Eigen::Map<Eigen::ArrayXd>(Q_val->begin(), n_l_);

    w.col(k) = Eigen::Map<Eigen::VectorXd>(w_val->begin(), n_b_);
    S2.col(k) = P.square() + Q.square();
  }
  progressBar.finalize();

  for (auto row : w.rowwise()) std::sort(row.begin(), row.end());
  for (auto row : S2.rowwise()) std::sort(row.begin(), row.end());

  // Identify violation indices
  std::vector<size_t> wl_ind(n_b_, 0);
  std::vector<size_t> wu_ind(n_b_, 0);
  std::vector<size_t> f_ind(n_l_, 0);

  for (size_t b = 0; b < n_b_; ++b) {
    for (size_t k = 0; k < n_k_ && w.coeff(b, k) < network_.wl[b]; ++k)
      wl_ind[b] = k;
    for (size_t k = n_k_ - 1; k >= 0 && w.coeff(b, k) > network_.wu[b]; ++k)
      wu_ind[b] = k;
  }

  for (size_t l = 0; l < n_l_; ++l) {
    for (size_t k = n_k_ - 1; k >= 0 && S2.coeff(l, k) > network_.S2[l]; ++k)
      f_ind[l] = k;
  }

  // Probability of violation (given by indices)
  // Can be computed directly above, but kept separate for clarity
  stats.P_wl = std::vector<double>(wl_ind.begin(), wl_ind.end());
  stats.P_wu = std::vector<double>(wu_ind.begin(), wu_ind.end());
  stats.P_f = std::vector<double>(f_ind.begin(), f_ind.end());

  std::for_each(stats.P_wl.begin(), stats.P_wl.end(),
                [this](double &a) { a /= this->n_k_; });
  std::for_each(stats.P_wu.begin(), stats.P_wu.end(),
                [this](double &a) { a /= this->n_k_; });
  std::for_each(stats.P_f.begin(), stats.P_f.end(),
                [this](double &a) { a /= this->n_k_; });

  // CVaR of limits
  // Ignore the component of the expectation between extremal values (i.e.
  // immediately less than/greater than; results in <2% error)
  size_t nu_idx = static_cast<size_t>((1 - network_.nu) * n_k_);
  size_t gamma_idx = static_cast<size_t>((1 - network_.gamma) * n_k_);

  stats.CVaR_wl =
      util::getEigenSTL<double>(w.leftCols(nu_idx).rowwise().mean());
  stats.CVaR_wu =
      util::getEigenSTL<double>(w.rightCols(nu_idx).rowwise().mean());
  stats.CVaR_f =
      util::getEigenSTL<double>(S2.rightCols(gamma_idx).rowwise().mean());

  // Expected violation
  stats.E_wl = std::vector<double>(n_b_, 0);
  stats.E_wu = std::vector<double>(n_b_, 0);
  stats.E_f = std::vector<double>(n_l_, 0);

  for (size_t b = 0; b < n_b_; ++b) {
    if (wl_ind[b] < n_k_ - 1) stats.E_wl[b] = w.row(b).head(wl_ind[b]).mean();
    if (wu_ind[b] > 0) stats.E_wu[b] = w.row(b).tail(wu_ind[b]).mean();
  }

  for (size_t l = 0; l < n_l_; ++l) {
    stats.E_f[l] = S2.row(l).tail(f_ind[l]).mean();
  }

  // Worst case violations
  stats.max_wl = util::getEigenSTL<double>(w.col(0));
  stats.max_wu = util::getEigenSTL<double>(w.col(n_k_ - 1));
  stats.max_f = util::getEigenSTL<double>(S2.col(n_k_ - 1));

  return stats;
}
