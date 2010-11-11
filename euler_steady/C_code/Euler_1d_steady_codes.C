#include <math.h>

double SourceQ_rho (double x)
{
  double Q_rho;
  double RHO;
  double U;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  Q_rho = cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x / L + cos(a_rhox * PI * x / L) * U * a_rhox * PI * rho_x / L;
  return(Q_rho);
}
#include <math.h>

double SourceQ_u (double x)
{
  double Q_u;
  double RHO;
  double U;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  Q_u = 0.2e1 * cos(a_ux * PI * x / L) * RHO * U * a_ux * PI * u_x / L + cos(a_rhox * PI * x / L) * U * U * a_rhox * PI * rho_x / L - sin(a_px * PI * x / L) * a_px * PI * p_x / L;
  return(Q_u);
}
#include <math.h>

double SourceQ_e (double x)
{
  double Q_e;
  double RHO;
  double U;
  double P;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  P = p_0 + p_x * cos(a_px * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  Q_e = cos(a_rhox * PI * x / L) * pow(U, 0.3e1) * a_rhox * PI * rho_x / L / 0.2e1 + cos(a_ux * PI * x / L) * P * a_ux * PI * u_x * Gamma / L / (Gamma - 0.1e1) + 0.3e1 / 0.2e1 * cos(a_ux * PI * x / L) * RHO * U * U * a_ux * PI * u_x / L - sin(a_px * PI * x / L) * U * a_px * PI * p_x * Gamma / L / (Gamma - 0.1e1);
  return(Q_e);
}
