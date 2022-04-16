#include <rshc/optimality/gurobi_model.h>

rshc::opt::GurobiModel::GurobiModel(const Model &model,
                                    const SolverOptions &opts)
    : network_(model), env_(), model_(env_), opts_(opts) {
  n_b_ = network_.wl.size();
  n_l_ = network_.r.size();
  n_d_ = network_.der_ind.size();
  n_k_ = network_.scenarios.size();
  // Set model to maximize
  model_.setObjective(GRBLinExpr(), GRB_MAXIMIZE);

  // Create variables
  {
    std::vector<double> psi_min(n_d_, 0);
    std::vector<double> psi_obj(n_d_, opts.obj_scale);
    psi_ = model_.addVars(psi_min.data(), network_.der_cap.data(),
                          psi_obj.data(), nullptr, nullptr, n_d_);

    std::vector<double> w_zero(n_b_, 0.);
    std::vector<double> neg_w_min(network_.wl);
    std::for_each(neg_w_min.begin(), neg_w_min.end(),
                  [](double &d) { d *= -1; });
    std::vector<double> n_inf(n_b_, -GRB_INFINITY);

    w_l_ = model_.addVars(n_inf.data(), neg_w_min.data(), nullptr, nullptr,
                          nullptr, n_b_);
    w_u_ = model_.addVars(w_zero.data(), network_.wu.data(), nullptr, nullptr,
                          nullptr, n_b_);

    std::vector<double> s_zero(n_l_, 0);
    s_ = model_.addVars(s_zero.data(), network_.S2.data(), nullptr, nullptr,
                        nullptr, n_l_);

    std::vector<double> l_inv(n_l_, -GRB_INFINITY);
    t_l_.reserve(n_k_);
    t_u_.reserve(n_k_);
    t_s_.reserve(n_k_);
    P_.reserve(n_k_);
    Q_.reserve(n_k_);
    w_.reserve(n_k_);
    l_.reserve(n_k_);
    for (size_t k = 0; k < n_k_; ++k) {
      // Note that Gurobi has default bound of 0
      t_l_.emplace_back(
          model_.addVars(nullptr, nullptr, nullptr, nullptr, nullptr, n_b_));
      t_u_.emplace_back(
          model_.addVars(nullptr, nullptr, nullptr, nullptr, nullptr, n_b_));
      t_s_.emplace_back(
          model_.addVars(nullptr, nullptr, nullptr, nullptr, nullptr, n_l_));

      P_.emplace_back(model_.addVars(l_inv.data(), nullptr, nullptr, nullptr,
                                     nullptr, n_l_));
      Q_.emplace_back(model_.addVars(l_inv.data(), nullptr, nullptr, nullptr,
                                     nullptr, n_l_));
      w_.emplace_back(
          model_.addVars(nullptr, nullptr, nullptr, nullptr, nullptr, n_b_));
      l_.emplace_back(
          model_.addVars(nullptr, nullptr, nullptr, nullptr, nullptr, n_l_));
    }

    // Create CVaR constraints
    c_v_l_ = model_.addConstrs(n_b_);
    c_v_u_ = model_.addConstrs(n_b_);
    c_f_ = model_.addConstrs(n_l_);

    for (size_t b = 0; b < n_b_; ++b) {
      (c_v_l_ + b)->set(GRB_CharAttr_Sense, '<');
      (c_v_l_ + b)->set(GRB_DoubleAttr_RHS, neg_w_min[b]);

      (c_v_u_ + b)->set(GRB_CharAttr_Sense, '<');
      (c_v_u_ + b)->set(GRB_DoubleAttr_RHS, network_.wu[b]);

      if (b < n_l_) {
        (c_f_ + b)->set(GRB_CharAttr_Sense, '<');
        (c_f_ + b)->set(GRB_DoubleAttr_RHS, network_.S2[b]);
      }
    }
  }

  // Set CVaR constraint coefficients
  std::vector<double> ones(n_b_, 1);
  model_.chgCoeffs(c_v_l_, w_l_, ones.data(), n_b_);
  model_.chgCoeffs(c_v_u_, w_u_, ones.data(), n_b_);
  model_.chgCoeffs(c_f_, s_, ones.data(), n_l_);

  for (size_t k = 0; k < n_k_; ++k) {
    std::vector<double> nu_param(n_b_, 1. / (n_k_ * (1. - network_.nu)));
    std::vector<double> gamma_param(n_l_, 1. / (n_k_ * (1. - network_.gamma)));
    model_.chgCoeffs(c_v_l_, t_l_[k], nu_param.data(), n_b_);
    model_.chgCoeffs(c_v_u_, t_u_[k], nu_param.data(), n_b_);
    model_.chgCoeffs(c_f_, t_s_[k], gamma_param.data(), n_l_);

    // Create linear constraints
    model_.addConstr(w_[k][n_b_ - 1], '=', 1);

    pb_p_.emplace_back(model_.addConstrs(n_l_));
    pb_q_.emplace_back(model_.addConstrs(n_l_));
    v_b_.emplace_back(model_.addConstrs(n_l_));
    v_l_.emplace_back(model_.addConstrs(n_b_));
    v_u_.emplace_back(model_.addConstrs(n_b_));

    for (size_t b = 0; b < n_b_; ++b) {
      if (b < n_l_) {
        (pb_p_[k] + b)->set(GRB_CharAttr_Sense, '=');
        (pb_q_[k] + b)->set(GRB_CharAttr_Sense, '=');
        (v_b_[k] + b)->set(GRB_CharAttr_Sense, '=');
      }
      (v_l_[k] + b)->set(GRB_CharAttr_Sense, '<');
      (v_u_[k] + b)->set(GRB_CharAttr_Sense, '<');
    }

    std::vector<double> n_ones(n_b_, -1);
    model_.chgCoeffs(v_l_[k], w_l_, n_ones.data(), n_b_);
    model_.chgCoeffs(v_l_[k], t_l_[k], n_ones.data(), n_b_);
    model_.chgCoeffs(v_l_[k], w_[k], n_ones.data(), n_b_);
    model_.chgCoeffs(v_u_[k], w_u_, n_ones.data(), n_b_);
    model_.chgCoeffs(v_u_[k], t_u_[k], n_ones.data(), n_b_);
    model_.chgCoeffs(v_u_[k], w_[k], ones.data(), n_b_);

    std::vector<double> n_r(n_l_, 0);
    std::vector<double> n_x(n_l_, 0);
    std::vector<double> n_two_r(n_l_, 0);
    std::vector<double> n_two_x(n_l_, 0);
    std::vector<double> r2_x2(n_l_, 0);
    for (size_t l = 0; l < n_l_; ++l) {
      n_r[l] = -network_.r[l];
      n_x[l] = -network_.x[l];
      n_two_r[l] = -2 * network_.r[l];
      n_two_x[l] = -2 * network_.x[l];
      r2_x2[l] = network_.r[l] * network_.r[l] + network_.x[l] * network_.x[l];

      (pb_p_[k] + l)->set(GRB_DoubleAttr_RHS, network_.scenarios[k].d[l]);
      (pb_q_[k] + l)->set(GRB_DoubleAttr_RHS, network_.scenarios[k].e[l]);
    }

    model_.chgCoeffs(v_b_[k], P_[k], n_two_r.data(), n_l_);
    model_.chgCoeffs(v_b_[k], Q_[k], n_two_x.data(), n_l_);
    model_.chgCoeffs(v_b_[k], l_[k], r2_x2.data(), n_l_);

    model_.chgCoeffs(pb_p_[k], l_[k], n_r.data(), n_l_);
    model_.chgCoeffs(pb_q_[k], l_[k], n_x.data(), n_l_);

    for (size_t l = 0; l < n_l_; ++l) {
      const int fbus = network_.from_bus[l] - 1;
      const int tbus = network_.to_bus[l] - 1;

      if (fbus < n_l_) {
        model_.chgCoeff(pb_p_[k][fbus], P_[k][l], -1);
        model_.chgCoeff(pb_q_[k][fbus], Q_[k][l], -1);
      }
      model_.chgCoeff(pb_p_[k][tbus], P_[k][l], 1);
      model_.chgCoeff(pb_q_[k][tbus], Q_[k][l], 1);

      model_.chgCoeff(v_b_[k][l], w_[k][fbus], 1);
      model_.chgCoeff(v_b_[k][l], w_[k][tbus], -1);

      model_.addQConstr(
          P_[k][l] * P_[k][l] + Q_[k][l] * Q_[k][l] - w_[k][fbus] * l_[k][l],
          '<', 0);
      model_.addQConstr(
          P_[k][l] * P_[k][l] + Q_[k][l] * Q_[k][l] - s_[l] - t_s_[k][l], '<',
          0);
    }

    for (size_t d = 0; d < n_d_; ++d) {
      const int bus = network_.der_ind[d];
      model_.chgCoeff(pb_p_[k][bus], psi_[d], network_.scenarios[k].alpha[d]);
      model_.chgCoeff(pb_q_[k][bus], psi_[d],
                      network_.eta[d] * network_.scenarios[k].alpha[d]);
    }
  }

  // Set parameters
  model_.set(GRB_IntParam_OutputFlag, opts_.verbose);
  model_.set(GRB_IntParam_Method, GRB_METHOD_DUAL);
  model_.update();
}

rshc::opt::GurobiModel::~GurobiModel() {
  delete[] psi_;
  delete[] w_l_;
  delete[] w_u_;
  delete[] s_;

  for (size_t k = 0; k < n_k_; ++k) {
    delete[] t_l_[k];
    delete[] t_u_[k];
    delete[] t_s_[k];
    delete[] P_[k];
    delete[] Q_[k];
    delete[] w_[k];
    delete[] l_[k];
  }

  delete[] c_v_l_;
  delete[] c_v_u_;
  delete[] c_f_;

  for (size_t k = 0; k < n_k_; ++k) {
    delete[] pb_p_[k];
    delete[] pb_q_[k];
    delete[] v_b_[k];
    delete[] v_l_[k];
    delete[] v_u_[k];
  }
}

rshc::opt::GurobiModel::RetVal rshc::opt::GurobiModel::solve() {
  model_.optimize();

  double *_psi = model_.get(GRB_DoubleAttr_X, psi_, n_d_);

  RetVal ret = {model_.get(GRB_IntAttr_Status),
                std::vector<double>(_psi, _psi + n_d_)};

  delete[] _psi;

  return ret;
}