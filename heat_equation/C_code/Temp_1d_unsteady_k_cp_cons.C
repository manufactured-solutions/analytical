#include <math.h>

void SourceQ (
  double x,
  double t,
  double A_x,
  double A_t,
  double D_t,
  double k_0,
  double cp_0,
  double rho)
{
  double Q_T;
  double *gradT_an;
  double T_an;
  Q_T = cos(A_x * x + A_t * t) * cos(D_t * t) * k_0 * A_x * A_x - (sin(A_x * x + A_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(D_t * t) * D_t) * rho * cp_0;
  T_an = cos(A_x * x + A_t * t) * cos(D_t * t);
  gradT_an[0] = -A_x * cos(D_t * t) * sin(A_x * x + A_t * t);
}
