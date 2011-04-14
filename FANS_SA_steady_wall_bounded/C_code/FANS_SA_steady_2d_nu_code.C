#include <math.h>

double SourceQ_nu (double x, double y)
{
  double Q_nu;
  double RHO;
  double U;
  double V;
  double T;
  double NU_SA;
  double Re_x;
  double c_f;
  double u_tau;
  double y_plus;
  double u_eq_plus;
  double u_eq;
  double d_eqplus_yplus;
  double mu_t;
  double chi;
  double f_v1;
  double f_v2;
  double Sm_orig;
  double Sm1;
  double Sm2;
  double Sm;
  double S_sa;
  double Omega;
  double f_w;
  double g;
  double r;
  Re_x = rho_inf * u_inf * x / mu;
  c_f = C_cf / F_c * pow(0.1e1 / F_c * Re_x, -0.1e1 / 0.7e1);
  u_tau = u_inf * sqrt(c_f / 0.2e1);
  y_plus = y * u_tau / nu_w;
  u_eq_plus = 0.1e1 / kappa * log(0.1e1 + kappa * y_plus) + C1 * (0.1e1 - exp(-y_plus / eta1) - y_plus / eta1 * exp(-y_plus * b));
  u_eq = u_tau * u_eq_plus;
  d_eqplus_yplus = 0.1e1 / (0.1e1 + kappa * y_plus) + C1 * (exp(-y_plus / eta1) / eta1 - exp(-y_plus * b) / eta1 + y_plus * b * exp(-y_plus * b) / eta1);
  U = u_inf / A * sin(A / u_inf * u_eq);
  V = eta_v * u_tau * y / x / 0.14e2;
  T = T_inf * (0.1e1 - r_T * (double) (Gamma - 1) * M_inf * M_inf * (0.1e1 - U * U * pow(u_inf, -0.2e1)) / 0.2e1);
  RHO = p_0 / R / T;
  NU_SA = kappa * u_tau * y - alpha * y * y;
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  f_v2 = 0.1e1 - chi / (0.1e1 + chi * f_v1);
  Omega = sqrt(pow(0.196e3 * x * x * y_plus * d_eqplus_yplus * cos(A / u_inf * u_eq) + 0.15e2 * eta_v * y * y, 0.2e1) * u_tau * u_tau * pow(x, -0.4e1) * pow(y, -0.2e1)) / 0.196e3;
  Sm_orig = NU_SA * pow(kappa, -0.2e1) * pow(d, -0.2e1) * f_v2;
  Sm1 = Sm_orig;
  Sm2 = Omega * (c_v2 * c_v2 * Omega + c_v3 * Sm_orig) / ((c_v3 + (-0.1e1) * 0.20e1 * c_v2) * Omega - Sm_orig);
  if (-c_v2 * Omega <= Sm_orig)
    Sm = Sm1;
  else
    Sm = Sm2;
  S_sa = Sm + Omega;
  r = NU_SA / S_sa * pow(kappa, -0.2e1) * pow(d, -0.2e1);
  g = r + c_w2 * (pow(r, 0.6e1) - r);
  f_w = g * pow((0.1e1 + pow(c_w3, 0.6e1)) / (pow(g, 0.6e1) + pow(c_w3, 0.6e1)), 0.1e1 / 0.6e1);
  Q_nu = -r_T * (double) (Gamma - 1) * M_inf * M_inf * T_inf * y_plus * u_tau * d_eqplus_yplus * NU_SA * RHO * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T + (double) (Gamma - 1) * (y_plus * d_eqplus_yplus + u_eq_plus) * y * kappa * r_T * M_inf * M_inf * T_inf * u_tau * u_tau * NU_SA * RHO * U * cos(A / u_inf * u_eq) / sigma * pow(u_inf, -0.2e1) * pow(x, -0.2e1) / T / 0.196e3 - (double) (Gamma - 1) * (0.2e1 * alpha * y - kappa * u_tau) * r_T * M_inf * M_inf * T_inf * y_plus * u_tau * d_eqplus_yplus * NU_SA * RHO * U * cos(A / u_inf * u_eq) / sigma * pow(u_inf, -0.2e1) / y / T + (double) (Gamma - 1) * (y_plus * d_eqplus_yplus + u_eq_plus) * r_T * M_inf * M_inf * T_inf * u_tau * NU_SA * RHO * U * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 - (double) (1 + c_b2) * y * y * kappa * kappa * u_tau * u_tau * RHO / sigma * pow(x, -0.2e1) / 0.196e3 - kappa * u_tau * y * RHO * U / x / 0.14e2 - (y_plus * d_eqplus_yplus + u_eq_plus) * u_tau * NU_SA * RHO * cos(A / u_inf * u_eq) / x / 0.14e2 + c_w1 * f_w * NU_SA * NU_SA * RHO * pow(d, -0.2e1) + RHO * NU_SA * V / y - S_sa * c_b1 * NU_SA * RHO + (-0.2e1 * alpha * y + kappa * u_tau) * RHO * V - pow(0.2e1 * alpha * y - kappa * u_tau, 0.2e1) * (double) (1 + c_b2) * RHO / sigma + (RHO * NU_SA + mu) * (0.392e3 * alpha * x * x - 0.15e2 * kappa * u_tau * y) / sigma * pow(x, -0.2e1) / 0.196e3;
  return(Q_nu);
}
