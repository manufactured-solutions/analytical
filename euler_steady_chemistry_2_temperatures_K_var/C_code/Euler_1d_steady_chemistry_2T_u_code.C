#include <math.h>

double SourceQ_u (double x)
{
  double Q_u;
  double RHO;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  Q_u = a_rho_N_x * PI * rho_N_x * U * U * cos(a_rho_N_x * PI * x / L) / L - a_rho_N2_x * PI * rho_N2_x * U * U * sin(a_rho_N2_x * PI * x / L) / L + 0.2e1 * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L - (0.2e1 * RHO_N + RHO_N2) * a_Tx * PI * T_x * R * sin(a_Tx * PI * x / L) / L / M_N / 0.2e1 - (-0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * T / L / M_N / 0.2e1;
  return(Q_u);
}
