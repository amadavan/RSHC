#include <rshc/optimality/mosek_model.h>

using namespace mosek::fusion;
typedef std::shared_ptr<monty::ndarray<double, 1>> ndarray_ptr_t;

rshc::opt::MosekModel::MosekModel(const Model &model, const SolverOptions &opts)
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
  std::vector<double> neg_w_min(network_.wl.begin(), network_.wl.end() - 1);
  std::for_each(neg_w_min.begin(), neg_w_min.end(), [](double &d) { d *= -1; });
  ndarray_ptr_t n_w_min = monty::new_array_ptr(neg_w_min);
  ndarray_ptr_t w_max = monty::new_array_ptr(
      std::vector<double>(network_.wu.begin(), network_.wu.end() - 1));
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
  psi_ = model_->variable(n_d_, Domain::inRange(0, psi_max));

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
      Domain::lessThan(n_w_min));
  c_v_u_ = model_->constraint(
      Expr::add(w_u_, Expr::mul(1. / (1. - network_.nu) * 1. / n_k_,
                                Expr::sum(t_u_, 1))),
      Domain::lessThan(w_max));
  c_f_ = model_->constraint(
      Expr::add(s_, Expr::mul(1. / (1. - network_.gamma) * 1. / n_k_,
                              Expr::sum(t_s_, 1))),
      Domain::lessThan(s_lim));

  // Set the objective
  model_->objective(
      ObjectiveSense::Maximize,
      Expr::mul(opts_.obj_scale * 1e12 / (sqrt((double)n_k_) * n_d_),
                Expr::sum(psi_)));

  // Add each scenario
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
            Expr::mul(der_loc, Expr::mulElm(eta, Expr::mulElm(alpha, psi_))),
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
        monty::new_array_ptr<Expression::t, 1>(
            {Expr::mul(0.5, Expr::add(Expr::add(t_s_k, s_), 1)),
             Expr::mul(0.5, Expr::sub(Expr::add(t_s_k, s_), 1)), vars_k.P_,
             vars_k.Q_});
    cons_k.f_ = model_->constraint(Expr::hstack(f_p), Domain::inQCone());

    k_cons_.emplace_back(cons_k);
  }
}

rshc::opt::MosekModel::~MosekModel() {
  model_->dispose();
}

rshc::opt::MosekModel::RetVal rshc::opt::MosekModel::solve() {
  model_->solve();

  ndarray_ptr_t psi_val = psi_->level();

  RetVal ret = {
    model_->getProblemStatus(),
    monty::new_vector_from_array_ptr<double>(psi_val)
  };

  if (opts_.save) {
    // Write each of the variables to file
    std::ofstream output("variables.txt");
    output << "psi: " << *(psi_->level()) << std::endl;
    output << "wl: " << *(w_l_->level()) << std::endl;
    output << "wu: " << *(w_u_->level()) << std::endl;
    output << "s: " << *(s_->level()) << std::endl;
    output << "tl: " << *(t_l_->level()) << std::endl;
    output << "tu: " << *(t_u_->level()) << std::endl;
    output << "ts: " << *(t_s_->level()) << std::endl;
    for (size_t k = 0; k < n_k_; ++k) {
      const ScenarioVars &v = k_vars_[k];
      output << "Scenario " << k << std::endl;
      output << "P: " << *(v.P_->level()) << std::endl;
      output << "Q: " << *(v.Q_->level()) << std::endl;
      output << "w: " << *(v.w_->level()) << std::endl;
      output << "l: " << *(v.l_->level()) << std::endl;
    }
  }

  return ret;
}