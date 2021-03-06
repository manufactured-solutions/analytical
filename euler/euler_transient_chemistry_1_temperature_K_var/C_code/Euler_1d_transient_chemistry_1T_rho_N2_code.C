#include <math.h>

double SourceQ_rho_N2 (
  double x,
  double t,
  double M_N,
  double h0_N,
  double h0_N2,
  double K,
  double Cf1_N,
  double Cf1_N2,
  double etaf1_N,
  double etaf1_N2,
  double Ea_N,
  double Ea_N2,
  double Function_to_Calculate_K)
{
  double Q_rho_N2;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double kf1_N;
  double kf1_N2;
  double K_eq;
  K_eq = Function_to_Calculate_K;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_t * cos(a_rho_N_t * PI * t / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_t * sin(a_rho_N2_t * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L) + T_t * cos(a_Tt * PI * t / L);
  kf1_N = Cf1_N * pow(T, etaf1_N) * exp(-Ea_N / R / T);
  kf1_N2 = Cf1_N2 * pow(T, etaf1_N2) * exp(-Ea_N2 / R / T);
  Q_rho_N2 = -a_rho_N2_x * PI * rho_N2_x * U * sin(a_rho_N2_x * PI * x / L) / L + a_ux * PI * u_x * RHO_N2 * cos(a_ux * PI * x / L) / L + a_rho_N2_t * PI * rho_N2_t * cos(a_rho_N2_t * PI * t / L) / L - (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N * RHO_N * pow(M_N, -0.2e1) / K_eq + (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N2 / M_N / 0.2e1;
  return(Q_rho_N2);
}
