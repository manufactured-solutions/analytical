#include <math.h>

double SourceQ_u (
  double x,
  double y,
  double t,
  double mu,
  double c_v1)
{
  double Q_u;
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
  double D2uDx2;
  double D2uDy2;
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
  D2vDxy = -0.15e2 / 0.14e2 * V / y / x;
  Q_u = r_T * (double) (Gamma - 1) * mu_t * M_inf * M_inf * T_inf * y_plus * y_plus * u_tau * u_tau * d_eqplus_yplus * d_eqplus_yplus * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) / T - r_T * (double) (Gamma - 1) * M_inf * M_inf * T_inf * y_plus * u_tau * d_eqplus_yplus * RHO * U * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T + r_T * (double) (Gamma - 1) * (u_eq_plus + y_plus * d_eqplus_yplus) * M_inf * M_inf * T_inf * u_tau * RHO * pow(U, 0.3e1) * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 + r_T * (double) (Gamma - 1) * pow(u_eq_plus + y_plus * d_eqplus_yplus, 0.2e1) * mu_t * M_inf * M_inf * T_inf * u_tau * u_tau * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) * pow(x, -0.2e1) / T / 0.147e3 + y_plus * u_tau * d_eqplus_yplus * RHO * V * cos(A / u_inf * u_eq) / y - (u_eq_plus + y_plus * d_eqplus_yplus) * u_tau * RHO * U * cos(A / u_inf * u_eq) / x / 0.7e1 + RHO * U * V / y + (-0.2e1 * mu_t - 0.2e1 * mu) * (0.2e1 / 0.3e1 * D2uDx2 + D2vDxy / 0.6e1 + D2uDy2 / 0.2e1) - r_T * (double) (Gamma - 1) * mu_t * M_inf * M_inf * T_inf * u_tau * U * V * cos(A / u_inf * u_eq) * (0.43e2 * y_plus * d_eqplus_yplus - 0.2e1 * u_eq_plus) * pow(u_inf, -0.2e1) / x / y / T / 0.42e2;
  return(Q_u);
}
