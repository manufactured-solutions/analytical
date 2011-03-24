#include <math.h>

double SourceQ_u (
  double x,
  double W_C,
  double W_C3,
  double mu)
{
  double Q_u;
  double RHO;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double P;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  Q_u = -a_rho_C3_x * PI * rho_C3_x * R * T * sin(a_rho_C3_x * PI * x / L) / L / W_C3 + a_rho_C_x * PI * rho_C_x * R * T * cos(a_rho_C_x * PI * x / L) / L / W_C + 0.4e1 / 0.3e1 * mu * a_ux * a_ux * PI * PI * u_x * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + 0.2e1 * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L;
  return(Q_u);
}
