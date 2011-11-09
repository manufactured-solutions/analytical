// Source term for the working variable nu_sa, FANS-SA
double SourceQ_nu_sa (
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
  double alpha,
  double sigma_sa,
  double c_w2,
  double c_w3,
  double c_b1,
  double c_b2,
  double c_v1,
  double c_v2,
  double c_v3)
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
  double f_v2;
  double C1;
  double u_inf;
  double rho_inf;
  double T_aw;
  double rho_w;
  double F_c;
  double nu_w;
  double A;
  double c_w1;
  double Sm_orig;
  double Sm1;
  double Sm2;
  double Sm;
  double S_sa;
  double Omega;
  double f_w;
  double g;
  double r;
  double d;
  double Q_nu_sa;
  double Q_nusa_convection;
  double Q_nusa_production;
  double Q_nusa_diffusion;
  double Q_nusa_gradsquare;
  double Q_nusa_dissipation;
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
  c_w1 = c_b1 * pow(kappa, -0.2e1) + (0.1e1 + c_b2) / sigma_sa;
    // "---------------------------------------------"
;
  Re_x = rho_inf * u_inf * x / mu;
  c_f = C_cf / F_c * pow(0.1e1 / F_c * Re_x, -0.1e1 / 0.7e1);
  u_tau = u_inf * sqrt(c_f / 0.2e1);
  y_plus = y * u_tau / nu_w;
  u_eq_plus = 0.1e1 / kappa * log(0.1e1 + kappa * y_plus) + C1 * (0.1e1 - exp(-y_plus / eta1) - y_plus / eta1 * exp(-y_plus * b));
  u_eq = u_tau * u_eq_plus;
  d_ueqplus_yplus = 0.1e1 / (0.1e1 + kappa * y_plus) + C1 * (exp(-y_plus / eta1) / eta1 - exp(-y_plus * b) / eta1 + y_plus * b * exp(-y_plus * b) / eta1);
    // "a suitable approximation for d------ -------------"
;
  d = y;
    // "---------------------------------------------"
;
  U = u_inf / A * sin(A / u_inf * u_eq);
  V = eta_v * u_tau * y / x / 0.14e2;
  T = T_inf * (0.1e1 + r_T * (Gamma - 0.1e1) * M_inf * M_inf * (0.1e1 - U * U * pow(u_inf, -0.2e1)) / 0.2e1);
  RHO = p_0 / R / T;
  NU_SA = kappa * u_tau * y - alpha * y * y;
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  mu_t = RHO * NU_SA * f_v1;
  f_v2 = 0.1e1 - chi / (0.1e1 + chi * f_v1);
  Omega = sqrt(pow(0.196e3 * x * x * y_plus * d_ueqplus_yplus * cos(A / u_inf * u_eq) + 0.15e2 * eta_v * y * y, 0.2e1) * u_tau * u_tau * pow(x, -0.4e1) * pow(y, -0.2e1)) / 0.196e3;
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
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_nusa_convection = r_T * (Gamma - 0.1e1) * d_ueqplus_yplus * M_inf * M_inf * T_inf * y_plus * u_tau * NU_SA * RHO * U * V * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / y / T - r_T * (Gamma - 0.1e1) * (d_ueqplus_yplus * y_plus + u_eq_plus) * M_inf * M_inf * T_inf * u_tau * NU_SA * RHO * U * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / x / T / 0.14e2 - kappa * u_tau * y * RHO * U / x / 0.14e2 - (d_ueqplus_yplus * y_plus + u_eq_plus) * u_tau * NU_SA * RHO * cos(A / u_inf * u_eq) / x / 0.14e2 + RHO * NU_SA * V / y + (-0.2e1 * alpha * y + kappa * u_tau) * RHO * V;
    // "Contribution from theproduction  to the total source term --------------------------------------------------"
;
  Q_nusa_production = -c_b1 * (Sm + Omega) * RHO * NU_SA;
    // "Contribution  from the diffusion to the total source term ----------------------------------------------"
;
  Q_nusa_diffusion = (0.2e1 * alpha * y - kappa * u_tau) * r_T * (Gamma - 0.1e1) * M_inf * M_inf * T_inf * y_plus * u_tau * d_ueqplus_yplus * NU_SA * RHO * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / sigma_sa / y / T - r_T * (Gamma - 0.1e1) * (d_ueqplus_yplus * y_plus + u_eq_plus) * kappa * y * M_inf * M_inf * T_inf * u_tau * u_tau * NU_SA * RHO * U * cos(A / u_inf * u_eq) * pow(u_inf, -0.2e1) / sigma_sa * pow(x, -0.2e1) / T / 0.196e3 - kappa * kappa * y * y * u_tau * u_tau * RHO / sigma_sa * pow(x, -0.2e1) / 0.196e3 - pow(0.2e1 * alpha * y - kappa * u_tau, 0.2e1) * RHO / sigma_sa - 0.15e2 / 0.196e3 * kappa * (RHO * NU_SA + mu) * y * u_tau / sigma_sa * pow(x, -0.2e1) + 0.2e1 * alpha * (RHO * NU_SA + mu) / sigma_sa;
    // "Contribution  from the grad squared  to the total source term ----------------------------------------------"
;
  Q_nusa_gradsquare = -c_b2 * kappa * kappa * y * y * u_tau * u_tau * RHO / sigma_sa * pow(x, -0.2e1) / 0.196e3 - c_b2 * pow(0.2e1 * alpha * y - kappa * u_tau, 0.2e1) * RHO / sigma_sa;
    // "Contribution  from the dissipation to the total source term ----------------------------------------------"
;
  Q_nusa_dissipation = c_w1 * f_w * RHO * NU_SA * NU_SA * pow(d, -0.2e1);
    // "Total source term ----------------------------------------------------------------------------------"
;
  Q_nu_sa = Q_nusa_convection + Q_nusa_production + Q_nusa_diffusion + Q_nusa_gradsquare + Q_nusa_dissipation;
  return(Q_nu_sa);
}
