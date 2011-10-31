#include <math.h>

double SourceQ_rho (double r, double z, double t)
{
  double Q_rho;
  double RHO;
  double U;
  double W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_r * (cos(a_ur * PI * r / L) - 0.1e1) * (u_z * sin(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / L));
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / L);
  Q_rho = -a_rhor * PI * rho_r * U * sin(a_rhor * PI * r / L) / L + a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L + a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) + a_ur * u_r * u_t * sin(a_ur * PI * r / L) * cos(a_ut * PI * t / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO / L + RHO * U / r;
  return(Q_rho);
}
