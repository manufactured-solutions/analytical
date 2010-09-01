#include <math.h>

void SourceQ (
  double x,
  double t,
  double A_x,
  double A_t,
  double D_t,
  double rho,
  double k_0,
  double k_1,
  double k_2,
  double cp_0,
  double cp_1,
  double cp_2)
{
  double Q_T;
  double T_an;
  double *gradT_an;
  Q_T = -pow(cos(D_t * t), 0.2e1) * A_x * A_x * k_1 + 0.3e1 * pow(cos(A_x * x + A_t * t), 0.3e1) * pow(cos(D_t * t), 0.3e1) * A_x * A_x * k_2 - sin(A_x * x + A_t * t) * cos(D_t * t) * A_t * rho * cp_0 - cos(A_x * x + A_t * t) * sin(D_t * t) * D_t * rho * cp_0 + (A_x * A_x * k_0 - 0.2e1 * pow(cos(D_t * t), 0.2e1) * A_x * A_x * k_2 - sin(A_x * x + A_t * t) * cos(D_t * t) * A_t * rho * cp_1 - cos(A_x * x + A_t * t) * sin(D_t * t) * D_t * rho * cp_1) * cos(A_x * x + A_t * t) * cos(D_t * t) + (0.2e1 * A_x * A_x * k_1 - sin(A_x * x + A_t * t) * cos(D_t * t) * A_t * rho * cp_2 - cos(A_x * x + A_t * t) * sin(D_t * t) * D_t * rho * cp_2) * pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1);
  gradT_an[0] = -A_x * cos(D_t * t) * sin(A_x * x + A_t * t);
}
