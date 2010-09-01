#include <math.h>

double SourceQ_v (
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
  double Q_v;
  double Phi;
  Phi = x * x + y * y + omega * t;
  Q_v = (-0.6e1 * pow(sin(Phi), 0.3e1) + (-0.4e1 * epsilon - 0.6e1) * pow(sin(Phi), 0.2e1) + (-0.3e1 * epsilon + 0.4e1) * sin(Phi) + 0.2e1 * epsilon * sin(0.2e1 * Phi) + epsilon * (0.2e1 * epsilon + 0.3e1) * cos(Phi) + 0.2e1 * epsilon + 0.3e1) * rho_0 * u_0 * v_0 * x + (-0.8e1 * epsilon * pow(sin(Phi), 0.2e1) + (-0.3e1 * sin(0.2e1 * Phi) - 0.6e1 * epsilon) * sin(Phi) - 0.3e1 * sin(0.2e1 * Phi) + (0.2e1 * epsilon * epsilon + 0.2e1) * cos(Phi) + 0.4e1 * epsilon) * rho_0 * v_0 * v_0 * y + 0.4e1 / 0.3e1 * mu * u_0 * x * y * sin(Phi) + (-0.2e1 * pow(sin(Phi), 0.2e1) - 0.3e1 / 0.2e1 * sin(Phi) + epsilon * cos(Phi) + 0.1e1) * rho_0 * v_0 * omega - 0.2e1 * p_0 * y * sin(Phi) + 0.2e1 / 0.3e1 * mu * (0.7e1 * sin(Phi) + (0.6e1 * x * x + 0.8e1 * y * y) * cos(Phi)) * v_0;
  return(Q_v);
}
#include <math.h>

double SourceQ_v (
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
  double t15;
  double t18;
  double t2;
  double t21;
  double t30;
  double t34;
  double t4;
  double t40;
  double t5;
  double t8;
  t1 = x * x;
  t2 = y * y;
  Phi = t1 + t2 + omega * t;
  t4 = sin(Phi);
  t5 = t4 * t4;
  t8 = 0.4e1 * epsilon;
  t15 = sin(0.2e1 * Phi);
  t18 = 0.2e1 * epsilon;
  t21 = cos(Phi);
  t30 = 0.3e1 * t15;
  t34 = epsilon * epsilon;
  t40 = v_0 * v_0;
  return((-0.6e1 * t5 * t4 + (-t8 - 0.6e1) * t5 + (-0.3e1 * epsilon + 0.4e1) * t4 + 0.2e1 * epsilon * t15 + epsilon * (t18 + 0.3e1) * t21 + t18 + 0.3e1) * rho_0 * u_0 * v_0 * x + (-0.8e1 * epsilon * t5 + (-t30 - 0.6e1 * epsilon) * t4 - t30 + (0.2e1 * t34 + 0.2e1) * t21 + t8) * rho_0 * t40 * y + 0.4e1 / 0.3e1 * mu * u_0 * x * y * t4 + (0.1e1 - 0.2e1 * t5 - 0.3e1 / 0.2e1 * t4 + epsilon * t21) * rho_0 * v_0 * omega - 0.2e1 * p_0 * y * t4 + 0.2e1 / 0.3e1 * mu * (0.7e1 * t4 + (0.6e1 * t1 + 0.8e1 * t2) * t21) * v_0);
}
