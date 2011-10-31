#include <math.h>

double SourceQ_e (double r, double z)
{
  double Q_e;
  double RHO;
  double P;
  double U;
  double W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L);
  P = p_0 + p_r * sin(a_pr * PI * r / L) + p_z * cos(a_pz * PI * z / L);
  U = u_r * u_z * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L);
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L);
  Q_e = -Gamma * a_ur * PI * u_r * u_z * P * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) / (Gamma - 0.1e1) / L + (cos(a_ur * PI * r / L) - 0.1e1) * a_uz * PI * u_r * u_z * RHO * U * W * cos(a_uz * PI * z / L) / L - (0.3e1 * U * U + W * W) * a_ur * PI * u_r * u_z * RHO * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) / L / 0.2e1 - a_wr * PI * w_r * RHO * U * W * sin(a_wr * PI * r / L) / L + Gamma * a_pr * PI * p_r * U * cos(a_pr * PI * r / L) / (Gamma - 0.1e1) / L - Gamma * a_pz * PI * p_z * W * sin(a_pz * PI * z / L) / (Gamma - 0.1e1) / L + Gamma * a_wz * PI * w_z * P * cos(a_wz * PI * z / L) / (Gamma - 0.1e1) / L - (U * U + W * W) * a_rhor * PI * rho_r * U * sin(a_rhor * PI * r / L) / L / 0.2e1 + (U * U + W * W) * a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L / 0.2e1 + (U * U + 0.3e1 * W * W) * a_wz * PI * w_z * RHO * cos(a_wz * PI * z / L) / L / 0.2e1 + Gamma * P * U / (Gamma - 0.1e1) / r + (U * U + W * W) * RHO * U / r / 0.2e1;
  return(Q_e);
}
