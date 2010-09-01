#include <math.h>

double SourceQ_u (
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
  double Q_u;
  Q_u = 0.4e1 / 0.3e1 * mu * u_1 * cos(a_ur * PI * r / L) * sin(a_uz * PI * z / L) * a_ur * a_ur * PI * PI * pow(L, -0.2e1) + mu * u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L) * a_uz * a_uz * PI * PI * pow(L, -0.2e1) - 0.2e1 / 0.3e1 * mu * w_1 * sin(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * a_wr * a_wz * PI * PI * pow(L, -0.2e1) + p_1 * cos(a_pr * PI * r / L) * cos(a_pz * PI * z / L) * a_pr * PI / L - u_1 * u_1 * pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * rho_1 * sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * a_rhor * PI / L + u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * rho_1 * cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) * sin(a_uz * PI * z / L) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_rhoz * PI / L + u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * cos(a_uz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_uz * PI / L + 0.2e1 / 0.3e1 * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) * a_ur * PI * u_1 * mu / L / r - 0.2e1 * sin(a_ur * PI * r / L) * pow(sin(a_uz * PI * z / L), 0.2e1) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * u_1 * u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * a_ur * PI / L + 0.2e1 / 0.3e1 * cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * a_wz * PI * w_1 * mu / L / r + cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L) * a_wz * PI * w_1 / L + u_1 * u_1 * pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) / r;
  return(Q_u);
}
