#include <math.h>

double SourceQ_u (
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
  double L)
{
  double Q_u;
  Q_u = -p_x * sin(a_px * PI * x / L) * a_px * PI / L - rho_N2_x * sin(a_rho_N2_x * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) * a_rho_N2_x * PI / L + rho_N_x * cos(a_rho_N_x * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) * a_rho_N_x * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_ux * PI / L;
  return(Q_u);
}
