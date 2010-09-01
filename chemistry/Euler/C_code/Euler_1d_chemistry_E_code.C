#include <math.h>

double SourceQ_E (
  double x,
  double p_0,
  double p_x,
  double rho_N_0,
  double rho_N_x,
  double rho_N2_0,
  double rho_N2_x,
  double u_0,
  double u_x,
  double a_px,
  double a_rho_N_x,
  double a_rho_N2_x,
  double a_ux,
  double L,
  double M_N,
  double k_b1N,
  double k_f1N,
  double k_b1N2,
  double k_f1N2,
  double h0_N,
  double h0_N2)
{
  double Q_E;
  Q_E = -(0.2e1 * p_0 + 0.2e1 * p_x * cos(a_px * PI * x / L)) * (rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_rho_N_x * PI * rho_N_x * cos(a_rho_N_x * PI * x / L) / L * pow(0.2e1 * rho_N_0 + 0.2e1 * rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L), -0.2e1) + (pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) + 0.2e1 * h0_N) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_rho_N_x * PI * rho_N_x * cos(a_rho_N_x * PI * x / L) / L / 0.2e1 - (0.10e2 * rho_N_0 + 0.10e2 * rho_N_x * sin(a_rho_N_x * PI * x / L) + 0.7e1 * rho_N2_0 + 0.7e1 * rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_px * PI * p_x * sin(a_px * PI * x / L) / L / (0.2e1 * rho_N_0 + 0.2e1 * rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) / 0.2e1 + (0.2e1 * rho_N_0 + 0.2e1 * rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * h0_N * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L / 0.2e1 + (p_0 + p_x * cos(a_px * PI * x / L)) * (rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L / (0.2e1 * rho_N_0 + 0.2e1 * rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) - (-(0.3e1 * rho_N_0 + 0.3e1 * rho_N2_0 + 0.3e1 * rho_N_x * sin(a_rho_N_x * PI * x / L) + 0.3e1 * rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) + h0_N * (rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) - 0.2e1 * h0_N2 * (rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) - 0.5e1 * p_0 - 0.5e1 * p_x * cos(a_px * PI * x / L)) * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L / 0.2e1 - (p_0 + p_x * cos(a_px * PI * x / L)) * (u_0 + u_x * sin(a_ux * PI * x / L)) * sin(a_rho_N2_x * PI * x / L) * a_rho_N2_x * PI * rho_N2_x / L / (0.2e1 * rho_N_0 + 0.2e1 * rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) + (p_0 + p_x * cos(a_px * PI * x / L)) * (rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * (u_0 + u_x * sin(a_ux * PI * x / L)) * sin(a_rho_N2_x * PI * x / L) * a_rho_N2_x * PI * rho_N2_x / L * pow(0.2e1 * rho_N_0 + 0.2e1 * rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L), -0.2e1) - (pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) + 0.2e1 * h0_N2) * (u_0 + u_x * sin(a_ux * PI * x / L)) * sin(a_rho_N2_x * PI * x / L) * a_rho_N2_x * PI * rho_N2_x / L / 0.2e1;
  return(Q_E);
}
