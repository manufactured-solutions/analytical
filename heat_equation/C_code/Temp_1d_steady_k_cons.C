#include <math.h>

void SourceQ (double x, double A_x, double k_0)
{
  double Q_T;
  double T_an;
  double *gradT_an;
  Q_T = A_x * A_x * k_0 * cos(A_x * x);
  T_an = cos(A_x * x);
  gradT_an[0] = -A_x * sin(A_x * x);
}
