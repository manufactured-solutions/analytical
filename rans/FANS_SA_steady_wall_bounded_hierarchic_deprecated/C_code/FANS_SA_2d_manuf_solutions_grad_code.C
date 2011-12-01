// Manufactured solutions, gradients and aux relations, FANS-SA
void Manufactured_solutions (
  double x,
  double y,
  double kappa,
  double _C,
  double eta1,
  double b,
  double r_T,
  double Gamma,
  double M_inf,
  double C_cf,
  double eta_v,
  double T_inf,
  double p_0,
  double R,
  double mu,
  double alpha)
{
  double Re_x;
  double c_f;
  double u_tau;
  double y_plus;
  double u_eq_plus;
  double u_eq;
  double d_ueqplus_yplus;
  double C1;
  double u_inf;
  double rho_inf;
  double T_aw;
  double rho_w;
  double F_c;
  double nu_w;
  double A;
  double rho_an;
  double T_an;
  double nu_sa_an;
  double u_an;
  double v_an;
  double grad_rho_an[2];
  double grad_nu_sa_an[2];
  double grad_u_an[2];
  double grad_v_an[2];
  double grad_T_an[2];
    // "constants computed from parameters -------------------------"
;
  C1 = -0.1e1 / kappa * log(kappa) + _C;
  u_inf = M_inf * sqrt(Gamma * R * T_inf);
  rho_inf = p_0 / R / T_inf;
  T_aw = T_inf * (0.1e1 + r_T * (Gamma - 0.1e1) * M_inf * M_inf / 0.2e1);
  rho_w = p_0 / R / T_aw;
  A = sqrt(0.1e1 - T_inf / T_aw);
  F_c = (T_aw / T_inf - 0.1e1) * pow(asin(A), -0.2e1);
  nu_w = mu / rho_w;
    // " ---------------------------------------------------------"
;
  Re_x = rho_inf * u_inf * x / mu;
  c_f = C_cf / F_c * pow(0.1e1 / F_c * Re_x, -0.1e1 / 0.7e1);
  u_tau = u_inf * sqrt(c_f / 0.2e1);
  y_plus = y * u_tau / nu_w;
  u_eq_plus = 0.1e1 / kappa * log(0.1e1 + kappa * y_plus) + C1 * (0.1e1 - exp(-y_plus / eta1) - y_plus / eta1 * exp(-y_plus * b));
  u_eq = u_tau * u_eq_plus;
  d_ueqplus_yplus = 0.1e1 / (0.1e1 + kappa * y_plus) + C1 * (exp(-y_plus / eta1) / eta1 - exp(-y_plus * b) / eta1 + y_plus * b * exp(-y_plus * b) / eta1);
    // "Manufactured solutions -------------------------------------"
;
  u_an = u_inf / A * sin(A / u_inf * u_eq);
  v_an = eta_v * u_tau * y / x / 0.14e2;
  T_an = T_inf * (0.1e1 + r_T * (Gamma - 0.1e1) * M_inf * M_inf * (0.1e1 - u_an * u_an * pow(u_inf, -0.2e1)) / 0.2e1);
  rho_an = p_0 / R / T;
  nu_sa_an = kappa * u_tau * y - alpha * y * y;
    // "Gradients of Manufactured solutions -------------------------"
;
  grad_u_an[0] = -cos(A / u_inf * u_eq) * u_tau * (u_eq_plus + y_plus * d_ueqplus_yplus) / x / 0.14e2;
  grad_u_an[1] = cos(A / u_inf * u_eq) * y_plus * d_ueqplus_yplus * u_tau / y;
  grad_v_an[0] = -0.15e2 / 0.196e3 * eta_v * u_tau * y * pow(x, -0.2e1);
  grad_v_an[1] = eta_v * u_tau / x / 0.14e2;
  grad_T_an[0] = T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * u_an * cos(A / u_inf * u_eq) * u_tau * (u_eq_plus + y_plus * d_ueqplus_yplus) * pow(u_inf, -0.2e1) / x / 0.14e2;
  grad_T_an[1] = -T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * u_an * cos(A / u_inf * u_eq) * y_plus * d_ueqplus_yplus * u_tau * pow(u_inf, -0.2e1) / y;
  grad_rho_an[0] = -p_0 * T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * u_an * cos(A / u_inf * u_eq) * u_tau * (u_eq_plus + y_plus * d_ueqplus_yplus) / R * pow(T_an, -0.2e1) * pow(u_inf, -0.2e1) / x / 0.14e2;
  grad_rho_an[1] = p_0 * T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * u_an * cos(A / u_inf * u_eq) * y_plus * d_ueqplus_yplus * u_tau / R * pow(T_an, -0.2e1) * pow(u_inf, -0.2e1) / y;
  grad_nu_sa_an[0] = -kappa * u_tau * y / x / 0.14e2;
  grad_nu_sa_an[1] = kappa * u_tau - 0.2e1 * alpha * y;
}
