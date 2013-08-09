#include <math.h>

double SourceQ_e (double x, double t)
{
  double RHO;
  double P;
  double U;
  double T;
  double MU;
  double DMu_Dx;
  double kappa;
  double Q_e;
  double Q_e_convection;
  double Q_e_work_pressure;
  double Q_e_work_viscous;
  double Q_e_conduction;
  double Q_e_time;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_t * sin(a_rhot * PI * t / Lt);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / Lt);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_t * cos(a_pt * PI * t / Lt);
  T = P / RHO / R;
  MU = A_mu * pow(T, 0.3e1 / 0.2e1) / (T + B_mu);
  DMu_Dx = a_rhox * PI * rho_x * MU * MU * cos(a_rhox * PI * x / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhox * PI * rho_x * MU * cos(a_rhox * PI * x / L) / L / RHO + a_px * PI * p_x * MU * MU * sin(a_px * PI * x / L) / A_mu / L / R / RHO * pow(T, -0.3e1 / 0.2e1) - 0.3e1 / 0.2e1 * a_px * PI * p_x * MU * sin(a_px * PI * x / L) / L / R / RHO / T;
  Q_e_convection = a_rhox * PI * rho_x * pow(U, 0.3e1) * cos(a_rhox * PI * x / L) / L / 0.2e1 + 0.3e1 / 0.2e1 * a_ux * PI * u_x * RHO * U * U * cos(a_ux * PI * x / L) / L - a_px * PI * p_x * U * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + a_ux * PI * u_x * P * cos(a_ux * PI * x / L) / (Gamma - 0.1e1) / L;
  Q_e_work_pressure = -a_px * PI * p_x * U * sin(a_px * PI * x / L) / L + a_ux * PI * u_x * P * cos(a_ux * PI * x / L) / L;
  Q_e_conduction = -0.2e1 * kappa * a_rhox * a_rhox * PI * PI * rho_x * rho_x * P * pow(cos(a_rhox * PI * x / L), 0.2e1) * pow(L, -0.2e1) / R * pow(RHO, -0.3e1) - 0.2e1 * kappa * a_rhox * a_px * PI * PI * rho_x * p_x * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) - kappa * a_rhox * a_rhox * PI * PI * rho_x * P * sin(a_rhox * PI * x / L) * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) + kappa * a_px * a_px * PI * PI * p_x * cos(a_px * PI * x / L) * pow(L, -0.2e1) / R / RHO;
  Q_e_work_viscous = -0.4e1 / 0.3e1 * a_ux * a_ux * MU * PI * PI * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * a_ux * a_ux * MU * PI * PI * u_x * U * sin(a_ux * PI * x / L) * pow(L, -0.2e1);
  Q_e_time = -a_ut * u_t * PI * RHO * U * sin(a_ut * PI * t / Lt) / Lt + a_rhot * rho_t * PI * U * U * cos(a_rhot * PI * t / Lt) / Lt / 0.2e1 - a_pt * p_t * PI * sin(a_pt * PI * t / Lt) / (Gamma - 0.1e1) / Lt;
  Q_e = Q_e_convection + Q_e_work_pressure + Q_e_work_viscous + Q_e_conduction + Q_e_time;
  return(Q_e);
}
