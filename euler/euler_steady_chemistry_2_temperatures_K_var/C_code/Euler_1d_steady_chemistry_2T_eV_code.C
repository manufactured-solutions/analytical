#include <math.h>

double SourceQ_eV (
  double x,
  int energy_level_N,
  int energy_level_N2,
  double *g_N,
  double *g_N2,
  double *theta_e_N,
  double *theta_e_N2)
{
  double Q_eV;
  double e_elec_N;
  double e_elec_N_num;
  double e_elec_N_den;
  double e_elec_N2;
  double e_elec_N2_num;
  double e_elec_N2_den;
  double e_vib_N2;
  double e_vib_eq_N2;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double TV;
  double P;
  double tau_vib_N2_N;
  double tau_vib_N2_N2;
  double De_elecN2_Dx;
  double De_elecN_Dx;
  double AUX1;
  double AUX2;
  double w_dot_N;
  double w_dot_N2;
  double kf1_N;
  double kf1_N2;
  double R1;
  double K;
  double T_bar;
  int i;
  double w_dot_V;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  TV = Tv_0 + Tv_x * cos(a_Tvx * PI * x / L);
  P = RHO_N * R * T / M_N + RHO_N2 * R * T / M_N / 0.2e1;
  T_bar = pow(T, q) * pow(TV, 0.1e1 - q);
  K = calculate_equilibrium_constant_K(T);
  AUX1 = 0.0e0;
  AUX2 = 0.0e0;
  e_elec_N_num = 0.0e0;
  e_elec_N2_num = 0.0e0;
  e_elec_N_den = 0.0e0;
  e_elec_N2_den = 0.0e0;
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
  e_vib_N2 = R * theta_v_N2 / M_N / (exp(theta_v_N2 / TV) - 0.1e1) / 0.2e1;
  e_vib_eq_N2 = R * theta_v_N2 / M_N / (exp(theta_v_N2 / T) - 0.1e1) / 0.2e1;
  De_elecN_Dx = e_elec_N * e_elec_N * M_N * a_Tvx * PI * Tv_x * sin(a_Tvx * PI * x / L) / R * pow(TV, -0.2e1) / L - R * AUX1 * a_Tvx * PI * Tv_x * sin(a_Tvx * PI * x / L) / M_N / e_elec_N_den * pow(TV, -0.2e1) / L;
  De_elecN2_Dx = 0.2e1 * e_elec_N2 * e_elec_N2 * M_N * a_Tvx * PI * Tv_x * sin(a_Tvx * PI * x / L) / R * pow(TV, -0.2e1) / L - R * AUX2 * a_Tvx * PI * Tv_x * sin(a_Tvx * PI * x / L) / M_N / e_elec_N2_den * pow(TV, -0.2e1) / L / 0.2e1;
  tau_vib_N2_N = exp(0.29e2 / 0.75000e5 * sqrt(0.6e1) * sqrt(M_N) * pow(theta_v_N2, 0.4e1 / 0.3e1) * (pow(T, -0.1e1 / 0.3e1) - pow(0.54e2, 0.1e1 / 0.4e1) * pow(M_N, 0.1e1 / 0.4e1) / 0.200e3) - 0.921e3 / 0.50e2) / P;
  tau_vib_N2_N2 = exp(0.29e2 / 0.25000e5 * sqrt(M_N) * pow(theta_v_N2, 0.4e1 / 0.3e1) * (pow(T, -0.1e1 / 0.3e1) - 0.3e1 / 0.200e3 * pow(M_N, 0.1e1 / 0.4e1)) - 0.921e3 / 0.50e2) / P;
  kf1_N2 = Cf1_N2 * pow(T_bar, etaf1_N2) * exp(-Ea_N2 / R / T_bar);
  kf1_N = Cf1_N * pow(T_bar, etaf1_N) * exp(-Ea_N / R / T_bar);
  R1 = -RHO_N * RHO_N2 * kf1_N * pow(M_N, -0.2e1) / 0.2e1 - RHO_N2 * RHO_N2 * kf1_N2 * pow(M_N, -0.2e1) / 0.4e1 + RHO_N * RHO_N * RHO_N2 * kf1_N2 / K * pow(M_N, -0.3e1) / 0.2e1 + pow(RHO_N, 0.3e1) * kf1_N / K * pow(M_N, -0.3e1);
  w_dot_N = -0.2e1 * M_N * R1;
  w_dot_N2 = -w_dot_N;
  w_dot_V = e_elec_N * w_dot_N + (e_vib_N2 + e_elec_N2) * w_dot_N2 + (e_vib_eq_N2 - e_vib_N2) * RHO_N2 * (0.2e1 * RHO_N / (0.2e1 * RHO_N + RHO_N2) / tau_vib_N2_N + RHO_N2 / (0.2e1 * RHO_N + RHO_N2) / tau_vib_N2_N2);
  Q_eV = e_elec_N * a_rho_N_x * PI * rho_N_x * U * cos(a_rho_N_x * PI * x / L) / L - (e_vib_N2 + e_elec_N2) * a_rho_N2_x * PI * rho_N2_x * U * sin(a_rho_N2_x * PI * x / L) / L + e_elec_N * a_ux * PI * u_x * RHO_N * cos(a_ux * PI * x / L) / L + (e_vib_N2 + e_elec_N2) * a_ux * PI * u_x * RHO_N2 * cos(a_ux * PI * x / L) / L - e_vib_N2 * a_Tvx * PI * Tv_x * theta_v_N2 * RHO_N2 * U * sin(a_Tvx * PI * x / L) / L * pow(TV, -0.2e1) - 0.2e1 * e_vib_N2 * e_vib_N2 * M_N * a_Tvx * PI * Tv_x * RHO_N2 * U * sin(a_Tvx * PI * x / L) / R * pow(TV, -0.2e1) / L + De_elecN_Dx * RHO_N * U + De_elecN2_Dx * RHO_N2 * U - w_dot_V;
  return(Q_eV);
}
