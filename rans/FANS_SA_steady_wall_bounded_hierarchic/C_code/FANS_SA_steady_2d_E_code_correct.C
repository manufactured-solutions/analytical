// Source term for the total energy equation, FANS-SA
double SourceQ_E (
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
  double alpha,
  double Pr,
  double Pr_t)
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
  double cp;
  double D2uDx2;
  double D2uDy2;
  double D2uDxy;
  double D2vDx2;
  double D2vDy2;
  double D2vDxy;
  double D2TDx2;
  double D2TDy2;
  double D2ueqDx2;
  double D2ueqDy2;
  double Q_E;
  double Q_E_work;
  double Q_E_heat_flux;
  double Q_E_convection;
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
  cp = Gamma * R / (Gamma - 0.1e1);
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
  D2uDx2 = -u_tau * u_tau * A * A * pow(u_eq_plus + y_plus * d_ueqplus_yplus, 0.2e1) * U * pow(x, -0.2e1) * pow(u_inf, -0.2e1) / 0.196e3 + D2ueqDx2 * cos(A / u_inf * u_eq);
  D2uDxy = A * A * (u_eq_plus + y_plus * d_ueqplus_yplus) * u_tau * u_tau * y_plus * d_ueqplus_yplus * U * pow(u_inf, -0.2e1) / x / y / 0.14e2 - u_tau * y_plus * d_ueqplus_yplus * cos(A / u_inf * u_eq) / x / y / 0.7e1 - y * cos(A / u_inf * u_eq) * D2ueqDy2 / x / 0.14e2;
  D2uDy2 = -U * A * A * y_plus * y_plus * d_ueqplus_yplus * d_ueqplus_yplus * u_tau * u_tau * pow(u_inf, -0.2e1) * pow(y, -0.2e1) + cos(A / u_inf * u_eq) * D2ueqDy2;
  D2vDx2 = 0.435e3 / 0.196e3 * V * pow(x, -0.2e1);
  D2vDy2 = 0.0e0;
  D2vDxy = -0.15e2 / 0.14e2 * V / x / y;
  D2TDx2 = -T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * pow(u_eq_plus + y_plus * d_ueqplus_yplus, 0.2e1) * u_tau * u_tau * (pow(cos(A / u_inf * u_eq), 0.2e1) - pow(sin(A / u_inf * u_eq), 0.2e1)) * pow(x, -0.2e1) * pow(u_inf, -0.2e1) / 0.196e3 - T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * D2ueqDx2 * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1);
  D2TDy2 = -T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * y_plus * y_plus * d_ueqplus_yplus * d_ueqplus_yplus * u_tau * u_tau * (pow(cos(A / u_inf * u_eq), 0.2e1) - pow(sin(A / u_inf * u_eq), 0.2e1)) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) - T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * D2ueqDy2 * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1);
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_E_convection = T_inf * r_T * (Gamma - 0.1e1) * M_inf * M_inf * y_plus * u_tau * d_ueqplus_yplus * RHO * pow(U, 0.3e1) * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T / 0.2e1 - r_T * (Gamma - 0.1e1) * (u_eq_plus + y_plus * d_ueqplus_yplus) * (U * U + V * V) * M_inf * M_inf * T_inf * u_tau * RHO * U * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.28e2 + (r_T * Gamma * M_inf * M_inf * T_inf * V * V - r_T * M_inf * M_inf * T_inf * V * V + 0.2e1 * u_inf * u_inf * T) * y_plus * u_tau * d_ueqplus_yplus * RHO * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T / 0.2e1 - 0.15e2 / 0.14e2 * RHO * U * V * V / x - (u_eq_plus + y_plus * d_ueqplus_yplus) * (0.3e1 * U * U + V * V + 0.2e1 * cp * T) * u_tau * RHO * cos(A / u_inf * u_eq) / x / 0.28e2 + (U * U + 0.3e1 * V * V + 0.2e1 * cp * T) * RHO * V / y / 0.2e1;
    // "Contribution from the heat flux to the total source term --------------------------------------------------"
;
  Q_E_heat_flux = (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * pow(Gamma - 0.1e1, 0.2e1) * cp * r_T * r_T * mu_t * pow(M_inf, 0.4e1) * T_inf * T_inf * y_plus * y_plus * u_tau * u_tau * d_ueqplus_yplus * d_ueqplus_yplus * U * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.4e1) / Pr_t * pow(y, -0.2e1) / NU_SA / RHO / T + (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * pow(Gamma - 0.1e1, 0.2e1) * pow(u_eq_plus + y_plus * d_ueqplus_yplus, 0.2e1) * cp * r_T * r_T * mu_t * pow(M_inf, 0.4e1) * T_inf * T_inf * u_tau * u_tau * U * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.4e1) / Pr_t * pow(x, -0.2e1) / NU_SA / RHO / T / 0.196e3 + (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * (u_eq_plus + y_plus * d_ueqplus_yplus) * cp * r_T * mu_t * kappa * M_inf * M_inf * T_inf * u_tau * u_tau * y * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / Pr_t * pow(x, -0.2e1) * pow(NU_SA, -0.2e1) / RHO / 0.196e3 - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (-kappa * u_tau + 0.2e1 * alpha * y) * (Gamma - 0.1e1) * cp * r_T * mu_t * M_inf * M_inf * T_inf * y_plus * u_tau * d_ueqplus_yplus * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / Pr_t / y * pow(NU_SA, -0.2e1) / RHO + (-mu_t / Pr_t - mu / Pr) * (D2TDx2 + D2TDy2) * cp;
    // "Contribution  from the viscous/turbulence work to the total source term ----------------------------------------------"
;
  Q_E_work = -(0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * r_T * mu_t * M_inf * M_inf * T_inf * y_plus * y_plus * u_tau * u_tau * d_ueqplus_yplus * d_ueqplus_yplus * U * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) / NU_SA / RHO / T + (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * (u_eq_plus + y_plus * d_ueqplus_yplus) * r_T * mu_t * M_inf * M_inf * T_inf * y_plus * u_tau * u_tau * d_ueqplus_yplus * U * V * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(u_inf, -0.2e1) / x / y / NU_SA / RHO / T / 0.42e2 - 0.4e1 / 0.3e1 * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * r_T * mu_t * M_inf * M_inf * T_inf * y_plus * u_tau * d_ueqplus_yplus * U * V * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) * pow(y, -0.2e1) / NU_SA / RHO / T - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * pow(u_eq_plus + y_plus * d_ueqplus_yplus, 0.2e1) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * u_tau * U * U * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(x, -0.2e1) * pow(u_inf, -0.2e1) / NU_SA / RHO / T / 0.147e3 + (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * (0.43e2 * y_plus * d_ueqplus_yplus - 0.2e1 * u_eq_plus) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * U * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / y / NU_SA / RHO / T / 0.42e2 - 0.15e2 / 0.196e3 * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (Gamma - 0.1e1) * (u_eq_plus + y_plus * d_ueqplus_yplus) * r_T * mu_t * M_inf * M_inf * T_inf * u_tau * U * V * V * cos(A / u_inf * u_eq) * pow(x, -0.2e1) * pow(u_inf, -0.2e1) / NU_SA / RHO / T - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * mu_t * kappa * y_plus * u_tau * u_tau * d_ueqplus_yplus * U * cos(A / u_inf * u_eq) / y * pow(NU_SA, -0.2e1) / RHO + (0.8e1 * RHO * NU_SA - 0.6e1 * mu_t) * alpha * mu_t * y_plus * u_tau * d_ueqplus_yplus * U * cos(A / u_inf * u_eq) * pow(NU_SA, -0.2e1) / RHO + (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (y_plus * d_ueqplus_yplus - 0.2e1 * u_eq_plus) * mu_t * kappa * u_tau * u_tau * V * cos(A / u_inf * u_eq) / x * pow(NU_SA, -0.2e1) / RHO / 0.42e2 - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (u_eq_plus + y_plus * d_ueqplus_yplus) * mu_t * kappa * u_tau * u_tau * y * U * cos(A / u_inf * u_eq) * pow(x, -0.2e1) * pow(NU_SA, -0.2e1) / RHO / 0.147e3 + 0.2e1 / 0.21e2 * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * alpha * (u_eq_plus + y_plus * d_ueqplus_yplus) * mu_t * u_tau * y * V * cos(A / u_inf * u_eq) / x * pow(NU_SA, -0.2e1) / RHO + 0.8e1 / 0.3e1 * (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * alpha * mu_t * V * V * pow(NU_SA, -0.2e1) / RHO - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (0.784e3 * x * x + 0.45e2 * y * y) * mu_t * kappa * u_tau * V * V * pow(x, -0.2e1) / y * pow(NU_SA, -0.2e1) / RHO / 0.588e3 - (0.4e1 * RHO * NU_SA - 0.3e1 * mu_t) * (-0.43e2 * kappa * u_tau + 0.90e2 * alpha * y) * mu_t * U * V / x * pow(NU_SA, -0.2e1) / RHO / 0.42e2 + (-0.2e1 * mu_t - 0.2e1 * mu) * (y_plus * y_plus * u_tau * u_tau * d_ueqplus_yplus * d_ueqplus_yplus * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(y, -0.2e1) / 0.2e1 + pow(u_eq_plus + y_plus * d_ueqplus_yplus, 0.2e1) * u_tau * u_tau * pow(cos(A / u_inf * u_eq), 0.2e1) * pow(x, -0.2e1) / 0.294e3 - (0.43e2 * y_plus * d_ueqplus_yplus - 0.2e1 * u_eq_plus) * u_tau * V * cos(A / u_inf * u_eq) / x / y / 0.42e2 + (0.2e1 / 0.3e1 * D2uDx2 + D2vDxy / 0.6e1 + D2uDy2 / 0.2e1) * U + (D2uDxy / 0.6e1 + D2vDx2 / 0.2e1 + 0.2e1 / 0.3e1 * D2vDy2) * V + 0.225e3 / 0.392e3 * V * V * pow(x, -0.2e1) + 0.2e1 / 0.3e1 * V * V * pow(y, -0.2e1));
    // "Total source term ----------------------------------------------------------------------------------"
;
  Q_E = Q_E_work + Q_E_heat_flux + Q_E_convection;
  return(Q_E);
}
