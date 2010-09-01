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
  double epsilon)
{
  double Q_u;
  Q_u = (0.3e1 * sin(x * x + y * y + omega * t) * sin(0.2e1 * x * x + 0.2e1 * y * y + 0.2e1 * omega * t) + 0.2e1 * epsilon * (0.3e1 + epsilon) * cos(x * x + y * y + omega * t) + (0.3e1 + 0.4e1 * epsilon) * sin(0.2e1 * x * x + 0.2e1 * y * y + 0.2e1 * omega * t)) * rho_0 * u_0 * u_0 * x - (0.6e1 * pow(sin(x * x + y * y + omega * t), 0.3e1) + (0.4e1 * epsilon + 0.6e1) * pow(sin(x * x + y * y + omega * t), 0.2e1) + (-0.4e1 + 0.3e1 * epsilon) * sin(x * x + y * y + omega * t) - epsilon * (0.2e1 * epsilon + 0.3e1) * cos(x * x + y * y + omega * t) - 0.2e1 * epsilon * sin(0.2e1 * x * x + 0.2e1 * y * y + 0.2e1 * omega * t) - 0.2e1 * epsilon - 0.3e1) * rho_0 * u_0 * v_0 * y + cos(x * x + y * y + omega * t) * omega * (0.2e1 * epsilon + 0.4e1 * sin(x * x + y * y + omega * t) + 0.3e1) * u_0 * rho_0 / 0.2e1 - 0.2e1 * p_0 * sin(x * x + y * y + omega * t) * x;
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
  double epsilon)
{
  double t1;
  double t12;
  double t15;
  double t2;
  double t20;
  double t23;
  double t31;
  double t4;
  double t5;
  double t7;
  t1 = x * x;
  t2 = y * y;
  t4 = t1 + t2 + omega * t;
  t5 = sin(t4);
  t7 = sin(0.2e1 * t4);
  t12 = cos(t4);
  t15 = 0.4e1 * epsilon;
  t20 = u_0 * u_0;
  t23 = t5 * t5;
  t31 = 0.2e1 * epsilon;
  return((0.3e1 * t5 * t7 + 0.2e1 * epsilon * (0.3e1 + epsilon) * t12 + (0.3e1 + t15) * t7) * rho_0 * t20 * x - (0.6e1 * t23 * t5 + (t15 + 0.6e1) * t23 + (-0.4e1 + 0.3e1 * epsilon) * t5 - epsilon * (t31 + 0.3e1) * t12 - 0.2e1 * epsilon * t7 - t31 - 0.3e1) * rho_0 * u_0 * v_0 * y + t12 * omega * (t31 + 0.4e1 * t5 + 0.3e1) * u_0 * rho_0 / 0.2e1 - 0.2e1 * p_0 * t5 * x);
}
