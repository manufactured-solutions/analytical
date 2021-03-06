// Source term for the y-momentum equation (wall normal velocity v), FANS-SA
double SourceQ_v (
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
  double c_v1,
  double alpha)
{
  double NU_SA;
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
  double chi;
  double mu_t;
  double f_v1;
  double C1;
  double u_inf;
  double rho_inf;
  double T_aw;
  double rho_w;
  double F_c;
  double nu_w;
  double A;
  double D2uDxy;
  double D2vDx2;
  double D2vDy2;
  double D2ueqDx2;
  double D2ueqDy2;
  double Q_v;
  double Q_v_convection;
  double Q_v_gradp;
  double Q_v_viscous;
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
  U = u_inf / A * sin(A / u_inf * u_eq);
  V = eta_v * u_tau * y / x / 0.14e2;
  T = T_inf * (0.1e1 + r_T * (Gamma - 0.1e1) * M_inf * M_inf * (0.1e1 - U * U * pow(u_inf, -0.2e1)) / 0.2e1);
  RHO = p_0 / R / T;
  NU_SA = kappa * u_tau * y - alpha * y * y;
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  mu_t = RHO * NU_SA * f_v1;
    // "---------------------------------------------"
;
  d_ueqplus_yplus = 0.1e1 / (0.1e1 + kappa * y_plus) + C1 * (exp(-y_plus / eta1) / eta1 - exp(-y_plus * b) / eta1 + y_plus * b * exp(-y_plus * b) / eta1);
  D2ueqDx2 = -u_tau * y_plus * y_plus * d_ueqplus_yplus / eta1 * pow(x, -0.2e1) / 0.196e3 + (0.17e2 * y_plus * d_ueqplus_yplus + 0.15e2 * u_eq_plus) * u_tau * pow(x, -0.2e1) / 0.196e3 - pow(x, -0.2e1) * ((b * b * eta1 * y_plus - 0.2e1 * b * eta1 - y_plus * b + 0.1e1) * C1 * u_tau * y_plus * y_plus * exp(-y_plus * b) * pow(eta1, -0.2e1) + (eta1 * kappa - kappa * y_plus - 0.1e1) * u_tau * y_plus * y_plus * pow(0.1e1 + kappa * y_plus, -0.2e1) / eta1) / 0.196e3;
  D2ueqDy2 = -u_tau * y_plus * y_plus * d_ueqplus_yplus / eta1 * pow(y, -0.2e1) - pow(y, -0.2e1) * ((b * b * eta1 * y_plus - 0.2e1 * b * eta1 - y_plus * b + 0.1e1) * C1 * u_tau * y_plus * y_plus * exp(-y_plus * b) * pow(eta1, -0.2e1) + (eta1 * kappa - kappa * y_plus - 0.1e1) * u_tau * y_plus * y_plus * pow(0.1e1 + kappa * y_plus, -0.2e1) / eta1);
  D2uDxy = A * A * (y_plus * d_ueqplus_yplus + u_eq_plus) * u_tau * u_tau * y_plus * d_ueqplus_yplus * U * pow(u_inf, -0.2e1) / x / y / 0.14e2 - u_tau * y_plus * d_ueqplus_yplus * cos(A / u_inf * u_eq) / x / y / 0.7e1 - y * cos(A / u_inf * u_eq) * D2ueqDy2 / x / 0.14e2;
  D2vDx2 = 0.435e3 / 0.196e3 * V * pow(x, -0.2e1);
  D2vDy2 = 0.0e0;
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_v_convection = r_T * (Gamma - 0.1e1) * M_inf * M_inf * T_inf * y_plus * u_tau * d_ueqplus_yplus * RHO * U * V * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T - r_T * (Gamma - 0.1e1) * (y_plus * d_ueqplus_yplus + u_eq_plus) * M_inf * M_inf * T_inf * u_tau * RHO * U * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 - (y_plus * d_ueqplus_yplus + u_eq_plus) * u_tau * RHO * V * cos(A / u_inf * u_eq) / x / 0.14e2 - 0.15e2 / 0.14e2 * RHO * U * V / x + 0.2e1 * RHO * V * V / y;
    // "Contribution from the gradient of pressure (body forces) to the total source term --------------------------------------------------"
;
  Q_v_gradp = 0.0e0;
    // "Contribution  from the viscous terms to the total source term ----------------------------------------------"
;
  Q_v_viscous = (-0.2e1 * mu - 0.2e1 * mu_t) * (D2uDxy / 0.6e1 + D2vDx2 / 0.2e1 + 0.2e1 / 0.3e1 * D2vDy2) + r_T * (Gamma - 0.1e1) * (y_plus * d_ueqplus_yplus + u_eq_plus) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * M_inf * M_inf * T_inf * y_plus * u_tau * u_tau * d_ueqplus_yplus * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) / x / y / NU_SA / RHO / T / 0.42e2 - 0.15e2 / 0.196e3 * r_T * (Gamma - 0.1e1) * (y_plus * d_ueqplus_yplus + u_eq_plus) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * M_inf * M_inf * T_inf * u_tau * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) * pow(x, -0.2e1) / NU_SA / RHO / T - 0.4e1 / 0.3e1 * r_T * (Gamma - 0.1e1) * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * M_inf * M_inf * T_inf * y_plus * u_tau * d_ueqplus_yplus * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) / NU_SA / RHO / T;
    // "Total source term ----------------------------------------------------------------------------------"
;
  Q_v = Q_v_convection + Q_v_gradp + Q_v_viscous;
  return(Q_v);
}
