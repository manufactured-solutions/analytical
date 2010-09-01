#include <math.h>

void SourceQ (
  double x,
  double A_x,
  double k_0,
  double k_1,
  double k_2)
{
  double Q_T;
  double T_an;
  double *gradT_an;
  Q_T = 0.3e1 * A_x * A_x * k_2 * pow(cos(A_x * x), 0.3e1) + 0.2e1 * A_x * A_x * k_1 * pow(cos(A_x * x), 0.2e1) - A_x * A_x * k_1 + (k_0 - 0.2e1 * k_2) * A_x * A_x * cos(A_x * x);
  T_an = cos(A_x * x);
  gradT_an[0] = -A_x * sin(A_x * x);
}
