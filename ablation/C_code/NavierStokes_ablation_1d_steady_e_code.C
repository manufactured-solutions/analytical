#include <math.h>

double SourceQ_e (
  double x,
  double W_C,
  double W_C3,
  double mu,
  double Gamma,
  double k)
{
  double Q_e;
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
  P = R * T * (RHO_C / W_C + RHO_C3 / W_C3);
  Q_e = -0.4e1 / 0.3e1 * mu * a_ux * a_ux * PI * PI * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) * pow(L, -0.2e1) - a_rho_C3_x * PI * rho_C3_x * R * T * U * sin(a_rho_C3_x * PI * x / L) / W_C3 / L + a_rho_C_x * PI * rho_C_x * R * T * U * cos(a_rho_C_x * PI * x / L) / W_C / L - a_Tx * PI * T_x * R * RHO_C3 * U * sin(a_Tx * PI * x / L) / W_C3 / L - a_Tx * PI * T_x * R * RHO_C * U * sin(a_Tx * PI * x / L) / W_C / L - a_Tx * PI * T_x * R * RHO * U * sin(a_Tx * PI * x / L) / (Gamma - 0.1e1) / L + 0.4e1 / 0.3e1 * mu * a_ux * a_ux * PI * PI * u_x * U * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_ux * PI * u_x * R * RHO * T * cos(a_ux * PI * x / L) / (Gamma - 0.1e1) / L + 0.3e1 / 0.2e1 * a_ux * PI * u_x * RHO * U * U * cos(a_ux * PI * x / L) / L + k * a_Tx * a_Tx * PI * PI * T_x * cos(a_Tx * PI * x / L) * pow(L, -0.2e1) + a_ux * PI * u_x * P * cos(a_ux * PI * x / L) / L - (a_rho_C3_x * rho_C3_x * sin(a_rho_C3_x * PI * x / L) - a_rho_C_x * rho_C_x * cos(a_rho_C_x * PI * x / L)) * PI * R * T * U / (Gamma - 0.1e1) / L - (a_rho_C3_x * rho_C3_x * sin(a_rho_C3_x * PI * x / L) - a_rho_C_x * rho_C_x * cos(a_rho_C_x * PI * x / L)) * PI * pow(U, 0.3e1) / L / 0.2e1;
  return(Q_e);
}
