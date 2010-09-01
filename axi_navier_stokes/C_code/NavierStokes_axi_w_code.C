#include <math.h>

double SourceQ_w (
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
  double Q_w;
  Q_w = sin(a_ur * PI * r / L) * cos(a_uz * PI * z / L) * a_ur * a_uz * PI * PI * mu * u_1 * pow(L, -0.2e1) / 0.3e1 + 0.4e1 / 0.3e1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L) * a_wz * a_wz * PI * PI * mu * w_1 * pow(L, -0.2e1) - sin(a_pr * PI * r / L) * sin(a_pz * PI * z / L) * a_pz * PI * p_1 / L - sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_rhor * PI * rho_1 / L + cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) * pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) * a_rhoz * PI * rho_1 / L - sin(a_uz * PI * z / L) * sin(a_ur * PI * r / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_ur * PI * u_1 / L - sin(a_wr * PI * r / L) * sin(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L) * a_wr * PI * w_1 / L + 0.2e1 * cos(a_wz * PI * z / L) * cos(a_wr * PI * r / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_wz * PI * w_1 / L - (cos(a_ur * PI * r / L) - 0.1e1) * cos(a_uz * PI * z / L) * a_uz * PI * mu * u_1 / L / r / 0.3e1 + sin(a_uz * PI * z / L) * u_1 * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * (cos(a_ur * PI * r / L) - 0.1e1) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) / r;
  return(Q_w);
}
