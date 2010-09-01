#include <math.h>

void SourceQ (
  double x,
  double y,
  double z,
  double A_x,
  double B_y,
  double C_z,
  double k_0)
{
  double Q_T;
  double T_an;
  double *gradT_an;
  Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * cos(C_z * z) * (A_x * A_x + B_y * B_y + C_z * C_z);
  T_an = cos(A_x * x) * cos(B_y * y) * cos(C_z * z);
  gradT_an[0] = -A_x * cos(B_y * y) * cos(C_z * z) * sin(A_x * x);
  gradT_an[1] = -B_y * cos(A_x * x) * cos(C_z * z) * sin(B_y * y);
  gradT_an[2] = -C_z * cos(A_x * x) * cos(B_y * y) * sin(C_z * z);
}
