#include <math.h>

double SourceQ_w (
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
  double Qw;
  double Q_time;
  double Q_gradp;
  double Q_convection;
  double Q_dissipation;
  double g;
  double dgdy;
  double d2gdy2;
  double d3gdy3;
  double alpha;
  double beta;
  g = Ly * Ly * sin(PI * y / Ly) / PI - y * (Ly - y) + y * y * pow(Ly / 0.2e1 - y, 0.2e1) * pow(Ly - y, 0.2e1);
  dgdy = pow(Ly, 0.4e1) * y / 0.2e1 - 0.9e1 / 0.2e1 * pow(Ly, 0.3e1) * y * y + 0.13e2 * Ly * Ly * pow(y, 0.3e1) - 0.15e2 * Ly * pow(y, 0.4e1) + 0.6e1 * pow(y, 0.5e1) + Ly * cos(PI * y / Ly) - Ly + 0.2e1 * y;
  d2gdy2 = pow(Ly, 0.4e1) / 0.2e1 - 0.9e1 * pow(Ly, 0.3e1) * y + 0.39e2 * Ly * Ly * y * y - 0.60e2 * Ly * pow(y, 0.3e1) + 0.30e2 * pow(y, 0.4e1) - PI * sin(PI * y / Ly) + 0.2e1;
  d3gdy3 = -0.9e1 * pow(Ly, 0.3e1) + 0.78e2 * Ly * Ly * y - 0.180e3 * Ly * y * y + 0.120e3 * pow(y, 0.3e1) - PI * PI * cos(PI * y / Ly) / Ly;
  alpha = -k_l / (k_l * k_l + k_n * k_n);
  beta = -k_n / (k_l * k_l + k_n * k_n);
  Q_time = -(a_12 * sin(a_12 * t + k_l * x) * sin(b_12 * t + k_n * z) - b_12 * cos(a_12 * t + k_l * x) * cos(b_12 * t + k_n * z)) * beta * dgdy;
  Q_gradp = beta * d2gdy2 * pow(cos(a_12 * t + k_l * x), 0.2e1) * cos(b_12 * t + k_n * z) * g * sin(b_12 * t + k_n * z) - (alpha * k_l * pow(sin(a_12 * t + k_l * x), 0.2e1) - beta * k_n * pow(cos(a_12 * t + k_l * x), 0.2e1)) * beta * dgdy * dgdy * cos(b_12 * t + k_n * z) * sin(b_12 * t + k_n * z);
  Q_convection = -beta * d2gdy2 * pow(cos(a_12 * t + k_l * x), 0.2e1) * cos(b_12 * t + k_n * z) * g * sin(b_12 * t + k_n * z) + (alpha * k_l * pow(sin(a_12 * t + k_l * x), 0.2e1) - beta * k_n * pow(cos(a_12 * t + k_l * x), 0.2e1)) * beta * dgdy * dgdy * cos(b_12 * t + k_n * z) * sin(b_12 * t + k_n * z);
  Q_dissipation = (-d3gdy3 + (k_l * k_l + k_n * k_n) * dgdy) * beta * nu * cos(a_12 * t + k_l * x) * sin(b_12 * t + k_n * z);
  Qw = Q_time + Q_gradp + Q_convection + Q_dissipation;
  return(Qw);
}
