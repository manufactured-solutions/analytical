#include <math.h>

double SourceQ_rho_C (double x)
{
  double Q_rho_C;
  double U;
  U = u_0 + u_x * sin(a_ux * PI * x / L);
  Q_rho_C = a_rho_C_x * PI * rho_C_x * U * cos(a_rho_C_x * PI * x / L) / L;
  return(Q_rho_C);
}
