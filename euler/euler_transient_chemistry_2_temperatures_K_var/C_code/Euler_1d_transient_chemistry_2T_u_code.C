#include <math.h>

double SourceQ_u (double x, double t)
{
  double Q_u;
  double RHO;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double TV;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_t * cos(a_rho_N_t * PI * t / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_t * sin(a_rho_N2_t * PI * t / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L) + T_t * cos(a_Tt * PI * t / L);
  TV = Tv_0 + Tv_x * cos(a_Tvx * PI * x / L) + Tv_t * sin(a_Tvt * PI * t / L);
  Q_u = 0.2e1 * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L - (0.2e1 * RHO_N + RHO_N2) * a_Tx * PI * T_x * R * sin(a_Tx * PI * x / L) / L / M_N / 0.2e1 - a_ut * PI * u_t * RHO * sin(a_ut * PI * t / L) / L - (-0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * T / L / M_N / 0.2e1 - (-a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * U * U / L + (-a_rho_N_t * rho_N_t * sin(a_rho_N_t * PI * t / L) + a_rho_N2_t * rho_N2_t * cos(a_rho_N2_t * PI * t / L)) * PI * U / L;
  return(Q_u);
}
