#include <math.h>

double SourceQ_v (
  double x,
  double y,
  double t,
  double mu,
  double c_v1)
{
  double Q_v;
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
  double D2vDx2;
  double D2vDy2;
  double D2uDxy;
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
  D2uDxy = U * A * A * u_tau * u_tau * (u_eq_plus + y_plus * d_eqplus_yplus) * y_plus * d_eqplus_yplus * pow(u_inf, -0.2e1) / x / y / 0.14e2;
  D2vDx2 = 0.435e3 / 0.196e3 * V * pow(x, -0.2e1);
  D2vDy2 = 0.0e0;
  Q_v = -r_T * (double) (Gamma - 1) * M_inf * M_inf * T_inf * y_plus * u_tau * d_eqplus_yplus * RHO * U * V * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T - (double) (Gamma - 1) * (u_eq_plus + y_plus * d_eqplus_yplus) * mu_t * r_T * M_inf * M_inf * T_inf * y_plus * u_tau * u_tau * d_eqplus_yplus * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) / x / y / T / 0.42e2 + 0.4e1 / 0.3e1 * (double) (Gamma - 1) * mu_t * r_T * M_inf * M_inf * T_inf * y_plus * u_tau * d_eqplus_yplus * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) / T + (double) (Gamma - 1) * (u_eq_plus + y_plus * d_eqplus_yplus) * r_T * M_inf * M_inf * T_inf * u_tau * RHO * U * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 + 0.15e2 / 0.196e3 * (double) (Gamma - 1) * (u_eq_plus + y_plus * d_eqplus_yplus) * mu_t * r_T * M_inf * M_inf * T_inf * u_tau * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) * pow(x, -0.2e1) / T - (u_eq_plus + y_plus * d_eqplus_yplus) * u_tau * RHO * V * cos(A / u_inf * u_eq) / x / 0.14e2 - 0.15e2 / 0.14e2 * RHO * U * V / x + 0.2e1 * RHO * V * V / y + (-0.2e1 * mu - 0.2e1 * mu_t) * (D2uDxy / 0.6e1 + D2vDx2 / 0.2e1 + 0.2e1 / 0.3e1 * D2vDy2);
  return(Q_v);
}
