#include <math.h>

double SourceQ_e (
  double x,
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
  alpha = exp(theta_v_N2 / T);
  E_vib_N2 = R_N2 * theta_v_N2 / (alpha - 0.1e1);
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  Q_e = -a_Tx * PI * T_x * alpha * theta_v_N2 * E_vib_N2 * RHO_N2 * U * sin(a_Tx * PI * x / L) / L / (alpha - 0.1e1) * pow(T, -0.2e1) - 0.5e1 / 0.2e1 * a_Tx * PI * T_x * R_N * RHO_N * U * sin(a_Tx * PI * x / L) / L - 0.7e1 / 0.4e1 * a_Tx * PI * T_x * R_N * RHO_N2 * U * sin(a_Tx * PI * x / L) / L + 0.5e1 / 0.2e1 * a_ux * PI * u_x * R_N * RHO_N * T * cos(a_ux * PI * x / L) / L + 0.7e1 / 0.4e1 * a_ux * PI * u_x * R_N * RHO_N2 * T * cos(a_ux * PI * x / L) / L + 0.3e1 / 0.2e1 * a_ux * PI * u_x * RHO * U * U * cos(a_ux * PI * x / L) / L - a_rho_N2_x * PI * rho_N2_x * E_vib_N2 * U * sin(a_rho_N2_x * PI * x / L) / L + a_ux * PI * u_x * E_vib_N2 * RHO_N2 * cos(a_ux * PI * x / L) / L + (h0_N * RHO_N + h0_N2 * RHO_N2) * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L - (-0.10e2 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + 0.7e1 * a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R_N * U * T / L / 0.4e1 - (-a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * pow(U, 0.3e1) / L / 0.2e1 - (-a_rho_N_x * rho_N_x * h0_N * cos(a_rho_N_x * PI * x / L) + a_rho_N2_x * rho_N2_x * h0_N2 * sin(a_rho_N2_x * PI * x / L)) * PI * U / L;
  return(Q_e);
}
