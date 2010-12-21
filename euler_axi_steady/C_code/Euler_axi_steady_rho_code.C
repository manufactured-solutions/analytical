#include <math.h>

double SourceQ_rho (double r, double z)
{
  double Q_rho;
  double RHO;
  double U;
  double W;
  RHO = rho_0 + rho_r * cos(a_rhor * PI * r / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_r * u_z * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L);
  W = w_0 + w_r * cos(a_wr * PI * r / L) + w_z * sin(a_wz * PI * z / L);
  Q_rho = -a_rhor * PI * rho_r * U * sin(a_rhor * PI * r / L) / L + a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L - (a_ur * u_r * u_z * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) - a_wz * w_z * cos(a_wz * PI * z / L)) * PI * RHO / L + RHO * U / r;
  return(Q_rho);
}
