#include <math.h>

double SourceQ_v (
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
  double Qv;
  double Q_time;
  double Q_gradp;
  double Q_convection;
  double Q_dissipation;
  double g;
  double dgdy;
  double d2gdy2;
  double alpha;
  double beta;
  g = Ly * Ly * sin(PI * y / Ly) / PI - y * (Ly - y) + y * y * pow(Ly / 0.2e1 - y, 0.2e1) * pow(Ly - y, 0.2e1);
  dgdy = pow(Ly, 0.4e1) * y / 0.2e1 - 0.9e1 / 0.2e1 * pow(Ly, 0.3e1) * y * y + 0.13e2 * Ly * Ly * pow(y, 0.3e1) - 0.15e2 * Ly * pow(y, 0.4e1) + 0.6e1 * pow(y, 0.5e1) + Ly * cos(PI * y / Ly) - Ly + 0.2e1 * y;
  d2gdy2 = pow(Ly, 0.4e1) / 0.2e1 - 0.9e1 * pow(Ly, 0.3e1) * y + 0.39e2 * Ly * Ly * y * y - 0.60e2 * Ly * pow(y, 0.3e1) + 0.30e2 * pow(y, 0.4e1) - PI * sin(PI * y / Ly) + 0.2e1;
  alpha = -k_l / (k_l * k_l + k_n * k_n);
  beta = -k_n / (k_l * k_l + k_n * k_n);
  Q_time = -(a_12 * cos(b_12 * t + k_n * z) * sin(a_12 * t + k_l * x) + b_12 * cos(a_12 * t + k_l * x) * sin(b_12 * t + k_n * z)) * g;
  Q_gradp = (-alpha * k_l * pow(cos(b_12 * t + k_n * z), 0.2e1) * pow(sin(a_12 * t + k_l * x), 0.2e1) - beta * k_n * pow(cos(a_12 * t + k_l * x), 0.2e1) * pow(sin(b_12 * t + k_n * z), 0.2e1) + pow(cos(a_12 * t + k_l * x), 0.2e1) * pow(cos(b_12 * t + k_n * z), 0.2e1)) * dgdy * g;
  Q_convection = (alpha * k_l * pow(cos(b_12 * t + k_n * z), 0.2e1) * pow(sin(a_12 * t + k_l * x), 0.2e1) + beta * k_n * pow(cos(a_12 * t + k_l * x), 0.2e1) * pow(sin(b_12 * t + k_n * z), 0.2e1) - pow(cos(a_12 * t + k_l * x), 0.2e1) * pow(cos(b_12 * t + k_n * z), 0.2e1)) * dgdy * g;
  Q_dissipation = (-d2gdy2 + (k_l * k_l + k_n * k_n) * g) * nu * cos(a_12 * t + k_l * x) * cos(b_12 * t + k_n * z);
  Qv = Q_time + Q_gradp + Q_convection + Q_dissipation;
  return(Qv);
}
