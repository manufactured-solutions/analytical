#include <math.h>

double Manufac_solution (void)
{
  double Q_nu;
  double Re_x;
  double c_f;
  double u_tau;
  double y_plus;
  double u_eq_plus;
  double u_eq;
  double d_utaux;
  double d_eqplus_yplus;
  double u_an;
  double v_an;
  double p_an;
  double rho_an;
  double T_an;
  double nu_sa_an;
  Re_x = rho_inf * u_inf * x / mu;
  c_f = C_cf / F_c * pow(0.1e1 / F_c * Re_x, -0.1e1 / 0.7e1);
  u_tau = u_inf * sqrt(c_f / 0.2e1);
  y_plus = y * u_tau / nu_w;
  u_eq_plus = 0.1e1 / kappa * log(0.1e1 + kappa * y_plus) + C1 * (0.1e1 - exp(-y_plus / eta1) - y_plus / eta1 * exp(-y_plus * b));
  u_eq = u_tau * u_eq_plus;
  d_eqplus_yplus = 0.1e1 / (0.1e1 + kappa * y_plus) + C1 * (exp(-y_plus / eta1) / eta1 - exp(-y_plus * b) / eta1 + y_plus * b * exp(-y_plus * b) / eta1);
  u_an = u_inf / A * sin(A / u_inf * u_eq);
  v_an = eta_v * u_tau * y / x / 0.14e2;
  p_an = p_0;
  T_an = T_inf * (0.1e1 - r_T * (double) (Gamma - 1) * M_inf * M_inf * (0.1e1 - u_an * u_an * pow(u_inf, -0.2e1)) / 0.2e1);
  nu_sa_an = kappa * u_tau * y - alpha * y * y;
  rho_an = p_0 / R / T_an;
  return(rho_an);
}
grad_rho_an[0] = p_0 * T_inf * r_T * (double) (Gamma - 1) * M_inf * M_inf * sin(A * u_eq / u_inf) * cos(A * u_eq / u_inf) * u_tau * (u_eq_plus + y_plus * d_eqplus_yplus) / R * pow(T, -0.2e1) / A / u_inf / x / 0.14e2;
grad_rho_an[1] = -p_0 * T_inf * r_T * (double) (Gamma - 1) * M_inf * M_inf * sin(A * u_eq / u_inf) * cos(A * u_eq / u_inf) * y_plus * d_eqplus_yplus * u_tau / R * pow(T, -0.2e1) / A / u_inf / y;
grad_p_an[0] = 0;
grad_p_an[1] = 0;
grad_u_an[0] = -cos(A * u_eq / u_inf) * u_tau * (u_eq_plus + y_plus * d_eqplus_yplus) / x / 0.14e2;
grad_u_an[1] = cos(A * u_eq / u_inf) * y_plus * d_eqplus_yplus * u_tau / y;
grad_v_an[0] = -0.15e2 / 0.196e3 * eta_v * u_tau * y * pow(x, -0.2e1);
grad_v_an[1] = eta_v * u_tau / x / 0.14e2;
grad_T_an[0] = -T_inf * r_T * (double) (Gamma - 1) * M_inf * M_inf * sin(A * u_eq / u_inf) * cos(A * u_eq / u_inf) * u_tau * (u_eq_plus + y_plus * d_eqplus_yplus) / A / u_inf / x / 0.14e2;
grad_T_an[1] = T_inf * r_T * (double) (Gamma - 1) * M_inf * M_inf * sin(A * u_eq / u_inf) * cos(A * u_eq / u_inf) * y_plus * d_eqplus_yplus * u_tau / A / u_inf / y;
grad_nu_sa_an[0] = -kappa * u_tau * y / x / 0.14e2;
grad_nu_sa_an[1] = kappa * u_tau - 0.2e1 * alpha * y;
