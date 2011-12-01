// Source term for the continuity equation, FANS-SA
double SourceQ_rho (
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
  double mu)
{
  double RHO;
  double U;
  double V;
  double T;
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
  double Q_rho;
  double Q_rho_convection;
    // "constants computed from parameters -------------"
;
  C1 = -0.1e1 / kappa * log(kappa) + _C;
  u_inf = M_inf * sqrt(Gamma * R * T_inf);
  rho_inf = p_0 / R / T_inf;
  T_aw = T_inf * (0.1e1 + r_T * (Gamma - 0.1e1) * M_inf * M_inf / 0.2e1);
  rho_w = p_0 / R / T_aw;
  A = sqrt(0.1e1 - T_inf / T_aw);
  F_c = (T_aw / T_inf - 0.1e1) * pow(asin(A), -0.2e1);
  nu_w = mu / rho_w;
    // "---------------------------------------------"
;
  Re_x = rho_inf * u_inf * x / mu;
  c_f = C_cf / F_c * pow(0.1e1 / F_c * Re_x, -0.1e1 / 0.7e1);
  u_tau = u_inf * sqrt(c_f / 0.2e1);
  y_plus = y * u_tau / nu_w;
  u_eq_plus = 0.1e1 / kappa * log(0.1e1 + kappa * y_plus) + C1 * (0.1e1 - exp(-y_plus / eta1) - y_plus / eta1 * exp(-y_plus * b));
  u_eq = u_tau * u_eq_plus;
  d_ueqplus_yplus = 0.1e1 / (0.1e1 + kappa * y_plus) + C1 * (exp(-y_plus / eta1) / eta1 - exp(-y_plus * b) / eta1 + y_plus * b * exp(-y_plus * b) / eta1);
  U = u_inf / A * sin(A / u_inf * u_eq);
  V = eta_v * u_tau * y / x / 0.14e2;
  T = T_inf * (0.1e1 + r_T * (Gamma - 0.1e1) * M_inf * M_inf * (0.1e1 - U * U * pow(u_inf, -0.2e1)) / 0.2e1);
  RHO = p_0 / R / T;
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_rho_convection = r_T * (Gamma - 0.1e1) * d_ueqplus_yplus * M_inf * M_inf * T_inf * y_plus * u_tau * RHO * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T - r_T * (Gamma - 0.1e1) * (d_ueqplus_yplus * y_plus + u_eq_plus) * M_inf * M_inf * T_inf * u_tau * RHO * U * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 - (d_ueqplus_yplus * y_plus + u_eq_plus) * u_tau * RHO * cos(A / u_inf * u_eq) / x / 0.14e2 + RHO * V / y;
    // "Total source term ----------------------------------------------------------------------------------"
;
  Q_rho = Q_rho_convection;
  return(Q_rho);
}
