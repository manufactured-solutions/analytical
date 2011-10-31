#include <math.h>

double SourceQ_E (double x, double y)
{
  double Qe;
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
  double D2TDx2;
  double D2TDy2;
  double D2uDx2;
  double D2uDy2;
  double D2uDxy;
  double D2vDx2;
  double D2vDy2;
  double D2vDxy;
  double d_ueqx2;
  double d_ueqy2;
  double mu_t;
  double chi;
  double f_v1;
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
  mu_t = RHO * NU_SA * f_v1;
  d_ueqx2 = -y_plus * y_plus * u_tau * d_eqplus_yplus / eta1 * pow(x, -0.2e1) / 0.196e3 - (b * b * eta1 * y_plus - y_plus * b + 0.1e1 - 0.2e1 * b * eta1) * exp(-y_plus * b) * C1 * y_plus * y_plus * u_tau * pow(eta1, -0.2e1) * pow(x, -0.2e1) / 0.196e3 - (-kappa * y_plus + kappa * eta1 - 0.1e1) * y_plus * y_plus * u_tau * pow(0.1e1 + kappa * y_plus, -0.2e1) / eta1 * pow(x, -0.2e1) / 0.196e3 + 0.17e2 / 0.196e3 * y_plus * u_tau * d_eqplus_yplus * pow(x, -0.2e1) + 0.15e2 / 0.196e3 * u_tau * u_eq_plus * pow(x, -0.2e1);
  d_ueqy2 = -y_plus * y_plus * u_tau * d_eqplus_yplus * pow(y, -0.2e1) / eta1 + y_plus * y_plus * u_tau * (-(b * b * eta1 * y_plus - y_plus * b + 0.1e1 - 0.2e1 * b * eta1) * exp(-y_plus * b) * C1 * pow(y, -0.2e1) * pow(eta1, -0.2e1) - (-kappa * y_plus + kappa * eta1 - 0.1e1) * pow(0.1e1 + kappa * y_plus, -0.2e1) / eta1 * pow(y, -0.2e1));
  D2uDx2 = -U * A * A * u_tau * u_tau * pow(u_eq_plus + y_plus * d_eqplus_yplus, 0.2e1) * pow(u_inf, -0.2e1) * pow(x, -0.2e1) / 0.196e3 + cos(A / u_inf * u_eq) * d_ueqx2;
  D2uDy2 = -U * A * A * y_plus * y_plus * d_eqplus_yplus * d_eqplus_yplus * u_tau * u_tau * pow(u_inf, -0.2e1) * pow(y, -0.2e1) + cos(A / u_inf * u_eq) * d_ueqy2;
  D2uDxy = U * A * A * u_tau * u_tau * (u_eq_plus + y_plus * d_eqplus_yplus) * y_plus * d_eqplus_yplus * pow(u_inf, -0.2e1) / x / y / 0.14e2;
  D2vDx2 = 0.435e3 / 0.196e3 * V * pow(x, -0.2e1);
  D2vDy2 = 0.0e0;
  D2vDxy = -0.15e2 / 0.14e2 * V / y / x;
  D2TDx2 = T_inf * r_T * (double) (Gamma - 1) * M_inf * M_inf * u_tau * u_tau * pow(u_eq_plus + y_plus * d_eqplus_yplus, 0.2e1) * pow(u_inf, -0.2e1) * pow(x, -0.2e1) * (pow(cos(A / u_inf * u_eq), 0.2e1) - pow(sin(A / u_inf * u_eq), 0.2e1)) / 0.196e3 + T_inf * r_T * (double) (Gamma - 1) * M_inf * M_inf * sin(A / u_inf * u_eq) * cos(A / u_inf * u_eq) * d_ueqx2 / u_inf / A;
  D2TDy2 = T_inf * r_T * (double) (Gamma - 1) * M_inf * M_inf * y_plus * y_plus * d_eqplus_yplus * d_eqplus_yplus * u_tau * u_tau * pow(u_inf, -0.2e1) * pow(y, -0.2e1) * (pow(cos(A / u_inf * u_eq), 0.2e1) - pow(sin(A / u_inf * u_eq), 0.2e1)) + T_inf * r_T * (double) (Gamma - 1) * M_inf * M_inf * sin(A / u_inf * u_eq) * cos(A / u_inf * u_eq) * d_ueqy2 / u_inf / A;
  Qe = 0.4e1 / 0.3e1 * (0.2e1 * alpha * y - kappa * u_tau) * f_v1 * RHO * V * V / y + (-0.2e1 * mu - 0.2e1 * mu_t) * (y_plus * y_plus * u_tau * u_tau * d_eqplus_yplus * d_eqplus_yplus * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(y, -0.2e1) / 0.2e1 + pow(u_eq_plus + y_plus * d_eqplus_yplus, 0.2e1) * u_tau * u_tau * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(x, -0.2e1) / 0.294e3 - (0.43e2 * y_plus * d_eqplus_yplus - 0.2e1 * u_eq_plus) * u_tau * V * cos(A / u_inf * u_eq) / y / x / 0.42e2 + (0.2e1 / 0.3e1 * D2uDx2 + D2vDxy / 0.6e1 + D2uDy2 / 0.2e1) * U + (D2uDxy / 0.6e1 + D2vDx2 / 0.2e1 + 0.2e1 / 0.3e1 * D2vDy2) * V + 0.225e3 / 0.392e3 * V * V * pow(x, -0.2e1) + 0.2e1 / 0.3e1 * V * V * pow(y, -0.2e1)) + (double) (Gamma - 1) * mu_t * r_T * M_inf * M_inf * T_inf * y_plus * y_plus * u_tau * u_tau * d_eqplus_yplus * d_eqplus_yplus * U * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) / T - (double) (Gamma - 1) * (u_eq_plus + y_plus * d_eqplus_yplus) * mu_t * r_T * M_inf * M_inf * T_inf * y_plus * u_tau * u_tau * d_eqplus_yplus * U * V * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) / x / y / T / 0.42e2 - (double) (Gamma - 1) * (u_eq_plus + y_plus * d_eqplus_yplus) * y * f_v1 * kappa * cp * r_T * M_inf * M_inf * T_inf * u_tau * u_tau * RHO * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) * pow(x, -0.2e1) / Pr_t / 0.196e3 + (0.2e1 * alpha * y - kappa * u_tau) * (double) (Gamma - 1) * f_v1 * cp * r_T * M_inf * M_inf * T_inf * y_plus * u_tau * d_eqplus_yplus * RHO * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / Pr_t - (double) (Gamma - 1) * (U * U + V * V) * r_T * M_inf * M_inf * T_inf * y_plus * u_tau * d_eqplus_yplus * RHO * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T / 0.2e1 + 0.4e1 / 0.3e1 * (double) (Gamma - 1) * mu_t * r_T * M_inf * M_inf * T_inf * y_plus * u_tau * d_eqplus_yplus * U * V * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) / T + (double) (Gamma - 1) * pow(u_eq_plus + y_plus * d_eqplus_yplus, 0.2e1) * mu_t * r_T * M_inf * M_inf * T_inf * u_tau * u_tau * U * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) * pow(x, -0.2e1) / T / 0.147e3 - (double) (Gamma - 1) * (0.43e2 * y_plus * d_eqplus_yplus - 0.2e1 * u_eq_plus) * mu_t * r_T * M_inf * M_inf * T_inf * u_tau * U * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / y / T / 0.42e2 + (double) (Gamma - 1) * (u_eq_plus + y_plus * d_eqplus_yplus) * (U * U + V * V) * r_T * M_inf * M_inf * T_inf * u_tau * RHO * U * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.28e2 + 0.15e2 / 0.196e3 * (double) (Gamma - 1) * (u_eq_plus + y_plus * d_eqplus_yplus) * mu_t * r_T * M_inf * M_inf * T_inf * u_tau * U * V * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) * pow(x, -0.2e1) / T + (double) (int) pow((double) (Gamma - 1), (double) 2) * mu_t * cp * r_T * r_T * pow(M_inf, 0.4e1) * T_inf * T_inf * y_plus * y_plus * u_tau * u_tau * d_eqplus_yplus * d_eqplus_yplus * U * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.4e1) * pow(y, -0.2e1) / Pr_t / T + (double) (int) pow((double) (Gamma - 1), (double) 2) * pow(u_eq_plus + y_plus * d_eqplus_yplus, 0.2e1) * mu_t * cp * r_T * r_T * pow(M_inf, 0.4e1) * T_inf * T_inf * u_tau * u_tau * U * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.4e1) * pow(x, -0.2e1) / Pr_t / T / 0.196e3 + 0.2e1 / 0.21e2 * (u_eq_plus + y_plus * d_eqplus_yplus) * alpha * y * f_v1 * u_tau * RHO * V * cos(A / u_inf * u_eq) / x - (u_eq_plus + y_plus * d_eqplus_yplus) * y * f_v1 * kappa * u_tau * u_tau * RHO * U * cos(A / u_inf * u_eq) * pow(x, -0.2e1) / 0.147e3 + (0.2e1 * alpha * y - kappa * u_tau) * f_v1 * y_plus * u_tau * d_eqplus_yplus * RHO * U * cos(A / u_inf * u_eq) / y - 0.15e2 / 0.196e3 * y * f_v1 * kappa * u_tau * RHO * V * V * pow(x, -0.2e1) - (0.90e2 * alpha * y - 0.43e2 * kappa * u_tau) * f_v1 * RHO * U * V / x / 0.42e2 - (u_eq_plus + y_plus * d_eqplus_yplus) * (0.2e1 * cp * T + 0.3e1 * U * U + V * V) * u_tau * RHO * cos(A / u_inf * u_eq) / x / 0.28e2 + (y_plus * d_eqplus_yplus - 0.2e1 * u_eq_plus) * f_v1 * kappa * u_tau * u_tau * RHO * V * cos(A / u_inf * u_eq) / x / 0.42e2 + y_plus * u_tau * d_eqplus_yplus * RHO * U * V * cos(A / u_inf * u_eq) / y + (-mu / Pr - mu_t / Pr_t) * (D2TDx2 + D2TDy2) * cp - 0.15e2 / 0.14e2 * RHO * U * V * V / x + (0.2e1 * cp * T + U * U + 0.3e1 * V * V) * RHO * V / y / 0.2e1;
  return(Qe);
}
