#include <math.h>

void SourceQ (
  double x,
  double y,
  double A_x,
  double B_y,
  double k_0,
  double k_1,
  double k_2)
{
  double Q_T;
  double T_an;
  double *gradT_an;
  Q_T = (0.3e1 * A_x * A_x + 0.3e1 * B_y * B_y) * k_2 * pow(cos(A_x * x), 0.3e1) * pow(cos(B_y * y), 0.3e1) + (0.2e1 * A_x * A_x + 0.2e1 * B_y * B_y) * k_1 * pow(cos(A_x * x), 0.2e1) * pow(cos(B_y * y), 0.2e1) - (pow(cos(B_y * y), 0.2e1) * A_x * A_x + pow(cos(A_x * x), 0.2e1) * B_y * B_y) * k_1 + (k_0 * A_x * A_x + k_0 * B_y * B_y - 0.2e1 * pow(cos(B_y * y), 0.2e1) * k_2 * A_x * A_x - 0.2e1 * pow(cos(A_x * x), 0.2e1) * k_2 * B_y * B_y) * cos(A_x * x) * cos(B_y * y);
  T_an = cos(A_x * x) * cos(B_y * y);
  gradT_an[0] = -A_x * cos(B_y * y) * sin(A_x * x);
  gradT_an[1] = -B_y * cos(A_x * x) * sin(B_y * y);
}
