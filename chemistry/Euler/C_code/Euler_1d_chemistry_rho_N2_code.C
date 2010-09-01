#include <math.h>

double SourceQ_rho_N2 (
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
  double k_f1N2)
{
  double Q_rho_N2;
  Q_rho_N2 = -rho_N2_x * sin(a_rho_N2_x * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_rho_N2_x * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * a_ux * PI / L - 0.2e1 * k_b1N * pow(rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L), 0.3e1) * pow(M_N, -0.2e1) - k_b1N2 * (rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * pow(rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L), 0.2e1) * pow(M_N, -0.2e1) + k_f1N * (rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * (rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L)) / M_N + k_f1N2 * pow(rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L), 0.2e1) / M_N / 0.2e1;
  return(Q_rho_N2);
}
