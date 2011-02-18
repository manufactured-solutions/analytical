#include <math.h>

double SourceQ_e (
  double x,
  double t,
  double R_N,
  double R_N2,
  double h0_N,
  double h0_N2,
  double theta_v_N2)
{
  double Q_e;
  double RHO;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double alpha;
  double E_vib_N2;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_t * cos(a_rho_N_t * PI * t / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_t * sin(a_rho_N2_t * PI * t / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L) + T_t * cos(a_Tt * PI * t / L);
  alpha = exp(theta_v_N2 / T);
  E_vib_N2 = R_N * theta_v_N2 / (alpha - 0.1e1) / 0.2e1;
  Q_e = 0.3e1 / 0.2e1 * a_ux * PI * u_x * RHO * U * U * cos(a_ux * PI * x / L) / L - (0.2e1 * RHO_N + RHO_N2) * a_Tx * PI * T_x * R_N * U * sin(a_Tx * PI * x / L) / L / 0.2e1 + (0.10e2 * RHO_N + 0.7e1 * RHO_N2) * a_ux * PI * u_x * R_N * T * cos(a_ux * PI * x / L) / L / 0.4e1 + a_ux * PI * u_x * E_vib_N2 * RHO_N2 * cos(a_ux * PI * x / L) / L - a_ut * PI * u_t * RHO * U * sin(a_ut * PI * t / L) / L + (h0_N * RHO_N + h0_N2 * RHO_N2) * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L - (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R_N * RHO_N2 * U * T / L / RHO / 0.4e1 + (0.5e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - 0.4e1 * a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R_N * U * T / L / 0.2e1 + (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * pow(U, 0.3e1) / L / 0.2e1 - 0.3e1 / 0.2e1 * (a_rho_N_t * rho_N_t * sin(a_rho_N_t * PI * t / L) - a_rho_N2_t * rho_N2_t * cos(a_rho_N2_t * PI * t / L)) * PI * R_N * RHO_N * T / L / RHO - 0.5e1 / 0.4e1 * (a_rho_N_t * rho_N_t * sin(a_rho_N_t * PI * t / L) - a_rho_N2_t * rho_N2_t * cos(a_rho_N2_t * PI * t / L)) * PI * R_N * RHO_N2 * T / L / RHO + (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * E_vib_N2 * RHO_N2 * U / L / RHO - (a_rho_N_t * rho_N_t * sin(a_rho_N_t * PI * t / L) - a_rho_N2_t * rho_N2_t * cos(a_rho_N2_t * PI * t / L)) * PI * U * U / L / 0.2e1 - (a_rho_N_t * rho_N_t * sin(a_rho_N_t * PI * t / L) - a_rho_N2_t * rho_N2_t * cos(a_rho_N2_t * PI * t / L)) * PI * E_vib_N2 * RHO_N2 / L / RHO + (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * (h0_N * RHO_N + h0_N2 * RHO_N2) * PI * U / L / RHO - (a_rho_N_t * rho_N_t * sin(a_rho_N_t * PI * t / L) - a_rho_N2_t * rho_N2_t * cos(a_rho_N2_t * PI * t / L)) * (h0_N * RHO_N + h0_N2 * RHO_N2) * PI / L / RHO;
  return(Q_e);
}
