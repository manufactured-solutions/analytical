#include <math.h>

void SourceQ (
  double x,
  double y,
  double t,
  double A_x,
  double A_t,
  double B_y,
  double B_t,
  double D_t,
  double rho,
  double k_0,
  double cp_0)
{
  double Q_T;
  double T_an;
  double *gradT_an;
  Q_T = -(sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t) * rho * cp_0 + (A_x * A_x + B_y * B_y) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * k_0;
  T_an = cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t);
  gradT_an[0] = -A_x * cos(D_t * t) * cos(B_y * y + B_t * t) * sin(A_x * x + A_t * t);
  gradT_an[1] = -B_y * cos(D_t * t) * cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t);
}
