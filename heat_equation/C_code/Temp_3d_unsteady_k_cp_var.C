#include <math.h>

void SourceQ (
  double x,
  double y,
  double z,
  double t,
  double A_x,
  double A_t,
  double B_y,
  double B_t,
  double C_z,
  double C_t,
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
  Q_T = -sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_0 * A_t - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_0 * B_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_0 * C_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * rho * cp_0 * D_t - pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(C_z * z + C_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_1 * A_x * A_x - pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(C_z * z + C_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_1 * B_y * B_y - pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_1 * C_z * C_z + 0.3e1 * pow(cos(A_x * x + A_t * t), 0.3e1) * pow(cos(B_y * y + B_t * t), 0.3e1) * pow(cos(C_z * z + C_t * t), 0.3e1) * pow(cos(D_t * t), 0.3e1) * (A_x * A_x + B_y * B_y + C_z * C_z) * k_2 + (-sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_1 * A_t - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_1 * B_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_1 * C_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * rho * cp_1 * D_t + k_0 * A_x * A_x + k_0 * B_y * B_y + k_0 * C_z * C_z - 0.2e1 * pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(C_z * z + C_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_2 * A_x * A_x - 0.2e1 * pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(C_z * z + C_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_2 * B_y * B_y - 0.2e1 * pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1) * k_2 * C_z * C_z) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) + (-sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_2 * A_t - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_2 * B_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_2 * C_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * rho * cp_2 * D_t + 0.2e1 * k_1 * A_x * A_x + 0.2e1 * k_1 * B_y * B_y + 0.2e1 * k_1 * C_z * C_z) * pow(cos(A_x * x + A_t * t), 0.2e1) * pow(cos(B_y * y + B_t * t), 0.2e1) * pow(cos(C_z * z + C_t * t), 0.2e1) * pow(cos(D_t * t), 0.2e1);
  T_an = cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t);
  gradT_an[0] = -A_x * cos(D_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(A_x * x + A_t * t);
  gradT_an[1] = -B_y * cos(D_t * t) * cos(A_x * x + A_t * t) * cos(C_z * z + C_t * t) * sin(B_y * y + B_t * t);
  gradT_an[2] = -C_z * cos(D_t * t) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t);
}
