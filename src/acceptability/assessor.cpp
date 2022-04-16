#include <rshc/acceptability/assessor.h>

using namespace mosek::fusion;
typedef std::shared_ptr<monty::ndarray<double, 1>> ndarray_ptr_t;

rshc::acc::Assessor::Assessor(const Model &model) : network_(model) {
  model_ = new mosek::fusion::Model();

  n_b_ = network_.wl.size();
  n_l_ = network_.r.size();
  n_d_ = network_.eta.size();
  n_k_ = network_.scenarios.size();

  /*******************************************************************
   * Solver parameters
   ******************************************************************/
  // Ensure that infeasible rays are generated
  model_->acceptedSolutionStatus(AccSolutionStatus::Anything);

  /*******************************************************************
   * Useful arrays and matrices
   ******************************************************************/
  ndarray_ptr_t psi_max = monty::new_array_ptr(network_.der_cap);

  // Note: we ignore voltage of feeder bus in CVaR constraints
  std::vector<double> neg_w_min(network_.wl.begin(), network_.wl.end() - 1);
  std::for_each(neg_w_min.begin(), neg_w_min.end(), [](double &d) { d *= -1; });
  n_w_min_ = monty::new_array_ptr(neg_w_min);
  w_max_ = monty::new_array_ptr(
      std::vector<double>(network_.wu.begin(), network_.wu.end() - 1));
  s_lim_ = monty::new_array_ptr(network_.S2);

  ndarray_ptr_t r = monty::new_array_ptr(network_.r);
  ndarray_ptr_t x = monty::new_array_ptr(network_.x);
  ndarray_ptr_t r2_x2 =
      std::make_shared<monty::ndarray<double, 1>>(n_l_, [this](ptrdiff_t i) {
        return network_.r[i] * network_.r[i] + network_.x[i] * network_.x[i];
      });

  eta_ = monty::new_array_ptr(network_.eta);

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

  w_l_ = model_->variable(n_l_, Domain::lessThan(0));
  w_u_ = model_->variable(n_l_, Domain::greaterThan(0));
  s_ = model_->variable(n_l_, Domain::greaterThan(0));

  t_l_ = model_->variable(monty::new_array_ptr<int, 1>({(int)n_l_, (int)n_k_}),
                          Domain::greaterThan(0.));
  t_u_ = model_->variable(monty::new_array_ptr<int, 1>({(int)n_l_, (int)n_k_}),
                          Domain::greaterThan(0.));
  t_s_ = model_->variable(monty::new_array_ptr<int, 1>({(int)n_l_, (int)n_k_}),
                          Domain::greaterThan(0.));

  // CVaR constraints
  c_v_l_ = model_->constraint(
      Expr::add(w_l_, Expr::mul(1. / (1. - network_.nu) * 1. / n_k_,
                                Expr::sum(t_l_, 1))),
      Domain::lessThan(n_w_min_));
  c_v_u_ = model_->constraint(
      Expr::add(w_u_, Expr::mul(1. / (1. - network_.nu) * 1. / n_k_,
                                Expr::sum(t_u_, 1))),
      Domain::lessThan(w_max_));
  c_f_ = model_->constraint(
      Expr::add(s_, Expr::mul(1. / (1. - network_.gamma) * 1. / n_k_,
                              Expr::sum(t_s_, 1))),
      Domain::lessThan(s_lim_));

  // Set objective
  model_->objective(ObjectiveSense::Maximize, 0);

  for (size_t k = 0; k < n_k_; ++k) {
    const Scenario &scenario = network_.scenarios[k];

    // Create variable for each scenario
    ScenarioVars vars_k;

    vars_k.P_ = model_->variable(n_l_, Domain::unbounded());
    vars_k.Q_ = model_->variable(n_l_, Domain::unbounded());
    vars_k.w_ = model_->variable(n_b_, Domain::greaterThan(0.));
    vars_k.l_ = model_->variable(n_l_, Domain::greaterThan(0.));

    k_vars_.emplace_back(vars_k);

    // Create constraints for each scenario
    ScenarioCons cons_k;

    ndarray_ptr_t d_k = monty::new_array_ptr(scenario.d);
    ndarray_ptr_t e_k = monty::new_array_ptr(scenario.e);

    ndarray_ptr_t alpha = monty::new_array_ptr(scenario.alpha);

    // Assumes that the lines are ordered such that the k-th line is received by
    // the k-th bus
    cons_k.v_ref_ =
        model_->constraint(vars_k.w_->index(n_b_ - 1), Domain::equalsTo(1.));
    cons_k.pb_p_ = model_->constraint(
        Expr::sub(Expr::mul(der_loc, Expr::mulElm(alpha, psi_)),
                  Expr::add(Expr::mul(pi_adj->transpose(), vars_k.P_),
                            Expr::add(Expr::mulElm(r, vars_k.l_), d_k))),
        Domain::equalsTo(0.));
    cons_k.pb_q_ = model_->constraint(
        Expr::sub(
            Expr::mul(der_loc, Expr::mulElm(eta_, Expr::mulElm(alpha, psi_))),
            Expr::add(Expr::mul(pi_adj->transpose(), vars_k.Q_),
                      Expr::add(Expr::mulElm(x, vars_k.l_), e_k))),
        Domain::equalsTo(0.));
    cons_k.v_b_ = model_->constraint(
        Expr::sub(Expr::add(Expr::mul(adj, vars_k.w_),
                            Expr::mulElm(r2_x2, vars_k.l_)),
                  Expr::mul(2, Expr::add(Expr::mulElm(r, vars_k.P_),
                                         Expr::mulElm(x, vars_k.Q_)))),
        Domain::equalsTo(0.));

    std::shared_ptr<monty::ndarray<Expression::t, 1>> flow_p =
        monty::new_array_ptr<Expression::t, 1>(
            {Expr::mul(0.5, Expr::mul(adj_pos, vars_k.w_)), vars_k.l_,
             vars_k.P_, vars_k.Q_});
    cons_k.flow_ =
        model_->constraint(Expr::hstack(flow_p), Domain::inRotatedQCone());

    cons_k.v_l_ = model_->constraint(
        Expr::add(vars_k.w_->slice(0, n_b_ - 1),
                  Expr::add(w_l_, t_l_->slice(
                                      monty::new_array_ptr<int, 1>({0, (int)k}),
                                      monty::new_array_ptr<int, 1>(
                                          {(int)n_l_, (int)k + 1})))),
        Domain::greaterThan(0.));
    cons_k.v_u_ = model_->constraint(
        Expr::sub(vars_k.w_->slice(0, n_b_ - 1),
                  Expr::add(w_u_, t_u_->slice(
                                      monty::new_array_ptr<int, 1>({0, (int)k}),
                                      monty::new_array_ptr<int, 1>(
                                          {(int)n_l_, (int)k + 1})))),
        Domain::lessThan(0.));

    Variable::t t_s_k =
        t_s_->slice(monty::new_array_ptr<int>({0, (int)k}),
                    monty::new_array_ptr<int>({(int)n_l_, (int)k + 1}));
    std::shared_ptr<monty::ndarray<Expression::t, 1>> f_p =
        monty::new_array_ptr<Expression::t, 1>({Expr::constTerm(n_l_, 0.5),
                                                Expr::add(t_s_k, s_), vars_k.P_,
                                                vars_k.Q_});
    cons_k.f_ = model_->constraint(Expr::hstack(f_p), Domain::inRotatedQCone());

    std::shared_ptr<monty::ndarray<int, 1>> idx_start2 =
        monty::new_array_ptr<int, 1>({0, 0});
    std::shared_ptr<monty::ndarray<int, 1>> idx_end2 =
        monty::new_array_ptr<int, 1>({(int)n_l_, 1});

    cons_k.f_c_ = cons_k.f_->slice(idx_start2, idx_end2);

    k_cons_.emplace_back(cons_k);
  }
}

rshc::acc::Assessor::~Assessor() { model_->dispose(); }

bool rshc::acc::Assessor::verify(const std::vector<double> &psi) {
  psi_->setValue(monty::new_array_ptr(psi));

  model_->solve();

  return (model_->getProblemStatus() == ProblemStatus::PrimalAndDualFeasible);
}

rshc::acc::InfeasibilityCut rshc::acc::Assessor::getInfeasibilityCut() {
  double tmp, constant = 0;

  InfeasibilityCut cut{Eigen::VectorXd::Zero(n_d_), 0};

  for (size_t k = 0; k < network_.scenarios.size(); ++k) {
    const Scenario &scenario = network_.scenarios[k];
    ndarray_ptr_t alpha = monty::new_array_ptr(scenario.alpha);
    ndarray_ptr_t d = monty::new_array_ptr(scenario.d);
    ndarray_ptr_t e = monty::new_array_ptr(scenario.e);

    ndarray_ptr_t pb_p_dual = k_cons_[k].pb_p_->dual();
    ndarray_ptr_t pb_q_dual = k_cons_[k].pb_q_->dual();

    mosek::LinAlg::dot(n_l_, pb_p_dual, d, tmp);
    cut.b += tmp;

    mosek::LinAlg::dot(n_l_, pb_q_dual, e, tmp);
    cut.b += tmp;

    for (size_t l = 0; l < n_l_; ++l)
      cut.b -= 0.5 * (*k_cons_[k].f_c_->dual())[l];

    cut.b += (*k_cons_[k].v_ref_->dual())[0];

    for (size_t d = 0; d < n_d_; ++d) {
      cut.A.coeffRef(d) += (*pb_p_dual)[network_.der_ind[d]] * (*alpha)[d];
      cut.A.coeffRef(d) += (*pb_q_dual)[network_.der_ind[d]] * (*alpha)[d] * (*eta_)[d];
    }
  }

  mosek::LinAlg::dot(n_l_, w_max_, c_v_u_->dual(), tmp);
  cut.b += tmp;

  mosek::LinAlg::dot(n_l_, n_w_min_, c_v_l_->dual(), tmp);
  cut.b += tmp;

  mosek::LinAlg::dot(n_l_, s_lim_, c_f_->dual(), tmp);
  cut.b += tmp;

  return cut;
}