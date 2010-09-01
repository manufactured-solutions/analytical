#include <math.h>

double SourceQ_rho (
  double x,
  double y,
  double t,
  double rho_0,
  double p_0,
  double u_0,
  double v_0,
  double omega,
  double epsilon)
{
  double Q_rho;
  Q_rho = cos(x * x + y * y + omega * t) * x * (0.2e1 * epsilon + 0.4e1 * sin(x * x + y * y + omega * t) + 0.3e1) * u_0 * rho_0 + (-0.4e1 * pow(sin(x * x + y * y + omega * t), 0.2e1) - 0.3e1 * sin(x * x + y * y + omega * t) + 0.2e1 * cos(x * x + y * y + omega * t) * epsilon + 0.2e1) * rho_0 * v_0 * y + rho_0 * cos(x * x + y * y + omega * t) * omega;
  return(Q_rho);
}
#include <math.h>

double SourceQ_rho (
  double x,
  double y,
  double t,
  double rho_0,
  double p_0,
  double u_0,
  double v_0,
  double omega,
  double epsilon)
{
  double t1;
  double t14;
  double t2;
  double t4;
  double t5;
  double t8;
  t1 = x * x;
  t2 = y * y;
  t4 = t1 + t2 + omega * t;
  t5 = cos(t4);
  t8 = sin(t4);
  t14 = t8 * t8;
  return(t5 * x * (0.2e1 * epsilon + 0.4e1 * t8 + 0.3e1) * u_0 * rho_0 + (-0.4e1 * t14 - 0.3e1 * t8 + 0.2e1 * t5 * epsilon + 0.2e1) * rho_0 * v_0 * y + rho_0 * t5 * omega);
}
