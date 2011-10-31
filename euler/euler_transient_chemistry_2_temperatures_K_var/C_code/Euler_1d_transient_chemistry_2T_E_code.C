#include <math.h>

double SourceQ_E (
  double x,
  double t,
  int energy_level_N,
  int energy_level_N2,
  double *g_N,
  double *g_N2,
  double *theta_e_N,
  double *theta_e_N2)
{
  double Q_E;
  double e_elec_N;
  double e_elec_N_num;
  double e_elec_N_den;
  double e_elec_N2;
  double e_elec_N2_num;
  double e_elec_N2_den;
  double e_vib_N2;
  double RHO;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double TV;
  double P;
  double De_elecN2_Dx;
  double De_elecN_Dx;
  double AUX1;
  double AUX2;
  int i;
  double De_elecN2_Dt;
  double De_elecN_Dt;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_t * cos(a_rho_N_t * PI * t / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_t * sin(a_rho_N2_t * PI * t / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L) + T_t * cos(a_Tt * PI * t / L);
  TV = Tv_0 + Tv_x * cos(a_Tvx * PI * x / L) + Tv_t * sin(a_Tvt * PI * t / L);
  P = RHO_N * R_N * T + RHO_N2 * R_N2 * T;
  e_elec_N_num = 0.0e0;
  e_elec_N2_num = 0.0e0;
  e_elec_N_den = g_N_0;
  e_elec_N2_den = g_N2_0;
  for (i = 1; i <= energy_level_N; i++)
  {
    e_elec_N_den = e_elec_N_den + g_N[i - 1] * exp(-theta_e_N[i - 1] / TV);
    e_elec_N_num = e_elec_N_num + theta_e_N[i - 1] * g_N[i - 1] * exp(-theta_e_N[i - 1] / TV);
    AUX1 = AUX1 + pow(theta_e_N[i - 1], 0.2e1) * g_N[i - 1] * exp(-theta_e_N[i - 1] / TV);
  }
  for (i = 1; i <= energy_level_N2; i++)
  {
    e_elec_N2_den = e_elec_N2_den + g_N2[i - 1] * exp(-theta_e_N2[i - 1] / TV);
    e_elec_N2_num = e_elec_N2_num + theta_e_N2[i - 1] * g_N2[i - 1] * exp(-theta_e_N2[i - 1] / TV);
    AUX2 = AUX2 + pow(theta_e_N2[i - 1], 0.2e1) * g_N2[i - 1] * exp(-theta_e_N2[i - 1] / TV);
  }
  e_elec_N = R * e_elec_N_num / M_N / e_elec_N_den;
  e_elec_N2 = R * e_elec_N2_num / M_N / e_elec_N2_den / 0.2e1;
  e_vib_N2 = theta_v_N2 * R / M_N / (exp(theta_v_N2 / TV) - 0.1e1) / 0.2e1;
  De_elecN_Dx = e_elec_N * e_elec_N * M_N * a_Tvx * PI * Tv_x * sin(a_Tvx * PI * x / L) / R * pow(TV, -0.2e1) / L - R * AUX1 * a_Tvx * PI * Tv_x * sin(a_Tvx * PI * x / L) / M_N / e_elec_N_den * pow(TV, -0.2e1) / L;
  De_elecN2_Dx = 0.2e1 * e_elec_N2 * e_elec_N2 * M_N * a_Tvx * PI * Tv_x * sin(a_Tvx * PI * x / L) / R * pow(TV, -0.2e1) / L - R * AUX2 * a_Tvx * PI * Tv_x * sin(a_Tvx * PI * x / L) / M_N / e_elec_N2_den * pow(TV, -0.2e1) / L / 0.2e1;
  De_elecN_Dt = -e_elec_N * e_elec_N * a_Tvt * PI * Tv_t * M_N * cos(a_Tvt * PI * t / L) / R * pow(TV, -0.2e1) / L + a_Tvt * PI * Tv_t * R * AUX1 * cos(a_Tvt * PI * t / L) / M_N / e_elec_N_den * pow(TV, -0.2e1) / L;
  De_elecN2_Dt = -0.2e1 * e_elec_N2 * e_elec_N2 * a_Tvt * PI * Tv_t * M_N * cos(a_Tvt * PI * t / L) / R * pow(TV, -0.2e1) / L + a_Tvt * PI * Tv_t * R * AUX2 * cos(a_Tvt * PI * t / L) / M_N / e_elec_N2_den * pow(TV, -0.2e1) / L / 0.2e1;
  Q_E = 0.3e1 / 0.2e1 * a_ux * PI * u_x * RHO * U * U * cos(a_ux * PI * x / L) / L - (0.10e2 * RHO_N + 0.7e1 * RHO_N2) * a_Tx * PI * T_x * R * U * sin(a_Tx * PI * x / L) / M_N / L / 0.4e1 + (0.10e2 * RHO_N + 0.7e1 * RHO_N2) * a_ux * PI * u_x * R * T * cos(a_ux * PI * x / L) / M_N / L / 0.4e1 - a_ut * PI * u_t * RHO * U * sin(a_ut * PI * t / L) / L + (h0_N + e_elec_N) * a_rho_N_x * PI * rho_N_x * U * cos(a_rho_N_x * PI * x / L) / L - (h0_N2 + e_vib_N2 + e_elec_N2) * a_rho_N2_x * PI * rho_N2_x * U * sin(a_rho_N2_x * PI * x / L) / L - (0.6e1 * RHO_N + 0.5e1 * RHO_N2) * a_Tt * PI * T_t * R * sin(a_Tt * PI * t / L) / M_N / L / 0.4e1 - e_vib_N2 * a_Tvx * PI * Tv_x * theta_v_N2 * RHO_N2 * U * sin(a_Tvx * PI * x / L) / L * pow(TV, -0.2e1) + (h0_N + e_elec_N) * a_ux * PI * u_x * RHO_N * cos(a_ux * PI * x / L) / L + (h0_N2 + e_vib_N2 + e_elec_N2) * a_ux * PI * u_x * RHO_N2 * cos(a_ux * PI * x / L) / L - (h0_N + e_elec_N) * a_rho_N_t * PI * rho_N_t * sin(a_rho_N_t * PI * t / L) / L + (h0_N2 + e_vib_N2 + e_elec_N2) * a_rho_N2_t * PI * rho_N2_t * cos(a_rho_N2_t * PI * t / L) / L - (-0.10e2 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + 0.7e1 * a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * U * T / M_N / L / 0.4e1 - (-a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * pow(U, 0.3e1) / L / 0.2e1 - 0.2e1 * M_N * e_vib_N2 * e_vib_N2 * a_Tvx * PI * Tv_x * RHO_N2 * U * sin(a_Tvx * PI * x / L) / R * pow(TV, -0.2e1) / L - (0.6e1 * a_rho_N_t * rho_N_t * sin(a_rho_N_t * PI * t / L) - 0.5e1 * a_rho_N2_t * rho_N2_t * cos(a_rho_N2_t * PI * t / L)) * PI * R * T / M_N / L / 0.4e1 + (-a_rho_N_t * rho_N_t * sin(a_rho_N_t * PI * t / L) + a_rho_N2_t * rho_N2_t * cos(a_rho_N2_t * PI * t / L)) * PI * U * U / L / 0.2e1 + e_vib_N2 * (theta_v_N2 * R + 0.2e1 * M_N * e_vib_N2) * a_Tvt * PI * Tv_t * RHO_N2 * cos(a_Tvt * PI * t / L) / R * pow(TV, -0.2e1) / L + De_elecN_Dx * RHO_N * U + De_elecN2_Dx * RHO_N2 * U + De_elecN_Dt * RHO_N + De_elecN2_Dt * RHO_N2;
  return(Q_E);
}
