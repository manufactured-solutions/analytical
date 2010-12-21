#include <math.h>

double SourceQ_u (double r, double z)
{
  double Q_u;
  double RHO;
  double U;
  double W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_r * u_z * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L);
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L);
  Q_u = (cos(a_ur * PI * r / L) - 0.1e1) * a_uz * PI * u_r * u_z * RHO * W * cos(a_uz * PI * z / L) / L - a_rhor * PI * rho_r * U * U * sin(a_rhor * PI * r / L) / L + a_rhoz * PI * rho_z * U * W * cos(a_rhoz * PI * z / L) / L + a_pr * PI * p_r * cos(a_pr * PI * r / L) / L - (0.2e1 * a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO * U / L + RHO * U * U / r;
  return(Q_u);
}
