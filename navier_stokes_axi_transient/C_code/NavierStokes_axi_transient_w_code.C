#include <math.h>

double SourceQ_w (double r, double z, double t)
{
  double Q_w;
  double RHO;
  double U;
  double W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_r * (cos(a_ur * PI * r / L) - 0.1e1) * (u_z * sin(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / L));
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / L);
  Q_w = -a_rhor * PI * rho_r * U * W * sin(a_rhor * PI * r / L) / L + a_rhoz * PI * rho_z * W * W * cos(a_rhoz * PI * z / L) / L - a_wr * PI * w_r * RHO * U * sin(a_wr * PI * r / L) / L + a_rhot * PI * rho_t * W * cos(a_rhot * PI * t / L) / L - a_wt * PI * w_t * RHO * sin(a_wt * PI * t / L) / L - a_pz * PI * p_z * sin(a_pz * PI * z / L) / L - (cos(a_ur * PI * r / L) - 0.1e1) * mu * a_uz * PI * u_r * u_z * cos(a_uz * PI * z / L) / L / r / 0.3e1 - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - 0.2e1 * a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO * W / L + (a_ur * a_uz * u_r * u_z * sin(a_ur * PI * r / L) * cos(a_uz * PI * z / L) + 0.4e1 * a_wz * a_wz * w_z * sin(a_wz * PI * z / L)) * mu * PI * PI * pow(L, -0.2e1) / 0.3e1 + RHO * U * W / r;
  return(Q_w);
}
