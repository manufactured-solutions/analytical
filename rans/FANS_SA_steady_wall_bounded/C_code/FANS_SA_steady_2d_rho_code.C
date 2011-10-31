#include <math.h>

double SourceQ_rho (double x, double y)
{
  double Q_rho;
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
  double d_eqplus_yplus;
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
  Q_rho = -r_T * (double) (Gamma - 1) * M_inf * M_inf * T_inf * y_plus * u_tau * d_eqplus_yplus * RHO * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T + r_T * (double) (Gamma - 1) * (y_plus * d_eqplus_yplus + u_eq_plus) * M_inf * M_inf * T_inf * u_tau * RHO * U * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 - (y_plus * d_eqplus_yplus + u_eq_plus) * u_tau * RHO * cos(A / u_inf * u_eq) / x / 0.14e2 + RHO * V / y;
  return(Q_rho);
}
