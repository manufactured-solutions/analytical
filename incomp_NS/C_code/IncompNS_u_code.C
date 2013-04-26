#include <math.h>

double SourceQ_u (
  double x,
  double y,
  double z,
  double t,
  double nu,
  double Ly,
  double k_l,
  double k_n,
  double a_12,
  double b_12)
{
  double Qu;
  double Q_time;
  double Q_gradp;
  double Q_convection;
  double Q_dissipation;
  double g;
  double dgdy;
  double d2gdy2;
  double d3gdy3;
  double alpha;
  g = Ly * Ly * cos(PI * y / Ly) / PI + y * y - Ly * Ly / 0.4e1 + pow(y - Ly / 0.2e1, 0.2e1) * pow(y, 0.4e1) * pow(y + Ly / 0.2e1, 0.2e1);
  dgdy = pow(Ly, 0.4e1) * pow(y, 0.3e1) / 0.4e1 - 0.3e1 * Ly * Ly * pow(y, 0.5e1) + 0.8e1 * pow(y, 0.7e1) - Ly * sin(PI * y / Ly) + 0.2e1 * y;
  d2gdy2 = 0.3e1 / 0.4e1 * pow(Ly, 0.4e1) * y * y - 0.15e2 * Ly * Ly * pow(y, 0.4e1) + 0.56e2 * pow(y, 0.6e1) - PI * cos(PI * y / Ly) + 0.2e1;
  d3gdy3 = 0.3e1 / 0.2e1 * pow(Ly, 0.4e1) * y - 0.60e2 * Ly * Ly * pow(y, 0.3e1) + 0.336e3 * pow(y, 0.5e1) + PI * PI * sin(PI * y / Ly) / Ly;
  alpha = -k_l * Lx / (k_l * k_l + k_n * k_n);
  Q_time = (a_12 * cos(k_l * x / Lx + a_12 * t / Lt) * cos(k_n * z / Lz + b_12 * t / Lt) - b_12 * sin(k_l * x / Lx + a_12 * t / Lt) * sin(k_n * z / Lz + b_12 * t / Lt)) * alpha * dgdy / Lt;
  Q_gradp = (-0.2e1 * y * y + 0.2e1 * Ly) * k_l * cos(k_l * x / Lx + a_12 * t / Lt) * sin(k_l * x / Lx + a_12 * t / Lt) * pow(sin(k_n * z / Lz + b_12 * t / Lt), 0.2e1) / Lx;
  Q_convection = g * d2gdy2 * k_l * Lx * cos(k_l * x / Lx + a_12 * t / Lt) * sin(k_l * x / Lx + a_12 * t / Lt) * pow(cos(k_n * z / Lz + b_12 * t / Lt), 0.2e1) / (k_l * k_l + k_n * k_n) - (k_l * k_l * pow(cos(k_n * z / Lz + b_12 * t / Lt), 0.2e1) / (k_l * k_l + k_n * k_n) - k_n * k_n * pow(sin(k_n * z / Lz + b_12 * t / Lt), 0.2e1) / (k_l * k_l + k_n * k_n)) * dgdy * dgdy * k_l * Lx * cos(k_l * x / Lx + a_12 * t / Lt) * sin(k_l * x / Lx + a_12 * t / Lt) / (k_l * k_l + k_n * k_n);
  Q_dissipation = (Lz * Lz * k_l * k_l + Lx * Lx * k_n * k_n) * alpha * nu * dgdy * sin(k_l * x / Lx + a_12 * t / Lt) * cos(k_n * z / Lz + b_12 * t / Lt) * pow(Lx, -0.2e1) * pow(Lz, -0.2e1) - alpha * nu * d3gdy3 * sin(k_l * x / Lx + a_12 * t / Lt) * cos(k_n * z / Lz + b_12 * t / Lt);
  Qu = Q_time + Q_gradp + Q_convection + Q_dissipation;
  return(Qu);
}
