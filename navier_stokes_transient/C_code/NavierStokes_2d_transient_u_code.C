#include <math.h>

double SourceQ_u (
  double x,
  double y,
  double t,
  double rho_0,
  double p_0,
  double u_0,
  double v_0,
  double omega,
  double epsilon,
  double mu)
{
  double Q_u;
  double Phi;
  Phi = x * x + y * y + omega * t;
  Q_u = (0.3e1 * sin(Phi) * sin(0.2e1 * Phi) + (0.4e1 * epsilon + 0.3e1) * sin(0.2e1 * Phi) + 0.2e1 * epsilon * (0.3e1 + epsilon) * cos(Phi)) * rho_0 * u_0 * u_0 * x + (-0.6e1 * pow(sin(Phi), 0.3e1) + (-0.4e1 * epsilon - 0.6e1) * pow(sin(Phi), 0.2e1) + (-0.3e1 * epsilon + 0.4e1) * sin(Phi) + 0.2e1 * epsilon * sin(0.2e1 * Phi) + epsilon * (0.2e1 * epsilon + 0.3e1) * cos(Phi) + 0.2e1 * epsilon + 0.3e1) * rho_0 * u_0 * v_0 * y + (0.2e1 * sin(Phi) + epsilon + 0.3e1 / 0.2e1) * rho_0 * u_0 * omega * cos(Phi) + 0.4e1 / 0.3e1 * mu * v_0 * x * y * cos(Phi) - 0.2e1 * p_0 * x * sin(Phi) - 0.2e1 / 0.3e1 * mu * ((-0.8e1 * x * x - 0.6e1 * y * y) * sin(Phi) + 0.7e1 * cos(Phi)) * u_0;
  return(Q_u);
}
#include <math.h>

double SourceQ_u (
  double x,
  double y,
  double t,
  double rho_0,
  double p_0,
  double u_0,
  double v_0,
  double omega,
  double epsilon,
  double mu)
{
  double Phi;
  double t1;
  double t14;
  double t19;
  double t2;
  double t22;
  double t32;
  double t4;
  double t6;
  double t9;
  t1 = x * x;
  t2 = y * y;
  Phi = t1 + t2 + omega * t;
  t4 = sin(Phi);
  t6 = sin(0.2e1 * Phi);
  t9 = 0.4e1 * epsilon;
  t14 = cos(Phi);
  t19 = u_0 * u_0;
  t22 = t4 * t4;
  t32 = 0.2e1 * epsilon;
  return((0.3e1 * t4 * t6 + (t9 + 0.3e1) * t6 + 0.2e1 * epsilon * (0.3e1 + epsilon) * t14) * rho_0 * t19 * x + (-0.6e1 * t22 * t4 + (-t9 - 0.6e1) * t22 + (-0.3e1 * epsilon + 0.4e1) * t4 + 0.2e1 * epsilon * t6 + epsilon * (t32 + 0.3e1) * t14 + t32 + 0.3e1) * rho_0 * u_0 * v_0 * y + (0.3e1 / 0.2e1 + 0.2e1 * t4 + epsilon) * rho_0 * u_0 * omega * t14 + 0.4e1 / 0.3e1 * mu * v_0 * x * y * t14 - 0.2e1 * p_0 * x * t4 - 0.2e1 / 0.3e1 * mu * ((-0.8e1 * t1 - 0.6e1 * t2) * t4 + 0.7e1 * t14) * u_0);
}
