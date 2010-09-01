#include <math.h>

double SourceQ_rho (
  double r,
  double z,
  double p_0,
  double p_1,
  double rho_0,
  double rho_1,
  double u_1,
  double w_0,
  double w_1,
  double a_pr,
  double a_pz,
  double a_rhor,
  double a_rhoz,
  double a_ur,
  double a_uz,
  double a_wr,
  double a_wz,
  double PI,
  double L)
{
  double Q_rho;
  Q_rho = -u_1 * sin(a_uz * PI * z / L) * rho_1 * sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * (cos(a_ur * PI * r / L) - 0.1e1) * a_rhor * PI / L + rho_1 * cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_rhoz * PI / L - u_1 * sin(a_uz * PI * z / L) * sin(a_ur * PI * r / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * a_ur * PI / L + w_1 * cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * a_wz * PI / L + u_1 * sin(a_uz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (cos(a_ur * PI * r / L) - 0.1e1) / r;
  return(Q_rho);
}
