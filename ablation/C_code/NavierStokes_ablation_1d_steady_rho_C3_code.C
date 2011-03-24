#include <math.h>

double SourceQ_rho_C3 (double x)
{
  double Q_rho_C3;
  double U;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  Q_rho_C3 = -a_rho_C3_x * PI * rho_C3_x * U * sin(a_rho_C3_x * PI * x / L) / L;
  return(Q_rho_C3);
}
