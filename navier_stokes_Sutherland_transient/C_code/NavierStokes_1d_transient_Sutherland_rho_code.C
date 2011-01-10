#include <math.h>

double SourceQ_rho (double x, double t)
{
  double Q_rho;
  double RHO;
  double U;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / L);
  Q_rho = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L + a_ux * PI * u_x * RHO * cos(a_ux * PI * x / L) / L + a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L;
  return(Q_rho);
}
