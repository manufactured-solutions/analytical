#include <math.h>

double SourceQ_e (double r, double z, double t)
{
  double Q_e;
  double RHO;
  double P;
  double U;
  double W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / L);
  P = p_0 + p_r * sin(a_pr * PI * r / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / L);
  U = u_r * (cos(a_ur * PI * r / L) - 0.1e1) * (u_z * sin(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / L));
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / L);
  Q_e = -(cos(a_ur * PI * r / L) - 0.1e1) * a_ut * PI * u_r * u_t * RHO * U * sin(a_ut * PI * t / L) / L + Gamma * a_pr * PI * p_r * U * cos(a_pr * PI * r / L) / (Gamma - 0.1e1) / L - Gamma * a_pz * PI * p_z * W * sin(a_pz * PI * z / L) / (Gamma - 0.1e1) / L - a_wt * PI * w_t * RHO * W * sin(a_wt * PI * t / L) / L - (U * U + W * W) * a_rhor * PI * rho_r * U * sin(a_rhor * PI * r / L) / L / 0.2e1 + (U * U + W * W) * a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L / 0.2e1 + (U * U + W * W) * a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L / 0.2e1 - a_pt * PI * p_t * sin(a_pt * PI * t / L) / (Gamma - 0.1e1) / L - (0.3e1 * a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + 0.3e1 * a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO * U * U / L / 0.2e1 + (a_uz * u_r * u_z * cos(a_ur * PI * r / L) * cos(a_uz * PI * z / L) - a_uz * u_r * u_z * cos(a_uz * PI * z / L) - a_wr * w_r * sin(a_wr * PI * r / L)) * PI * RHO * U * W / L - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - 0.3e1 * a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO * W * W / L / 0.2e1 - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * Gamma * PI * P / (Gamma - 0.1e1) / L + Gamma * P * U / (Gamma - 0.1e1) / r + (U * U + W * W) * RHO * U / r / 0.2e1;
  return(Q_e);
}
