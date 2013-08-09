#include <math.h>

double SourceQ_rho (double x, double t)
{
  double RHO;
  double U;
  double Q_rho;
  double Q_rho_convection;
  double Q_rho_time;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / L);
  Q_rho_convection = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L + a_ux * PI * u_x * RHO * cos(a_ux * PI * x / L) / L;
  Q_rho_time = a_rhot * rho_t * cos(a_rhot * PI * t / Lt) * PI / Lt;
  Q_rho = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L + a_ux * PI * u_x * RHO * cos(a_ux * PI * x / L) / L + a_rhot * rho_t * cos(a_rhot * PI * t / Lt) * PI / Lt;
  return(Q_rho);
}
