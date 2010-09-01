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
  double epsilon)
{
  double Q_v;
  Q_v = -(0.6e1 * pow(sin(x * x + y * y + omega * t), 0.3e1) + (0.6e1 + 0.4e1 * epsilon) * pow(sin(x * x + y * y + omega * t), 0.2e1) + (-0.4e1 + 0.3e1 * epsilon) * sin(x * x + y * y + omega * t) - epsilon * (0.3e1 + 0.2e1 * epsilon) * cos(x * x + y * y + omega * t) - 0.2e1 * epsilon * sin(0.2e1 * x * x + 0.2e1 * y * y + 0.2e1 * omega * t) - 0.3e1 - 0.2e1 * epsilon) * rho_0 * u_0 * v_0 * x - (0.8e1 * epsilon * pow(sin(x * x + y * y + omega * t), 0.2e1) + (0.3e1 * sin(0.2e1 * x * x + 0.2e1 * y * y + 0.2e1 * omega * t) + 0.6e1 * epsilon) * sin(x * x + y * y + omega * t) + (-0.2e1 * epsilon * epsilon - 0.2e1) * cos(x * x + y * y + omega * t) + 0.3e1 * sin(0.2e1 * x * x + 0.2e1 * y * y + 0.2e1 * omega * t) - 0.4e1 * epsilon) * rho_0 * v_0 * v_0 * y - (0.2e1 * pow(sin(x * x + y * y + omega * t), 0.2e1) + 0.3e1 / 0.2e1 * sin(x * x + y * y + omega * t) - epsilon * cos(x * x + y * y + omega * t) - 0.1e1) * rho_0 * v_0 * omega - 0.2e1 * p_0 * y * sin(x * x + y * y + omega * t);
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
  double epsilon)
{
  double t1;
  double t15;
  double t18;
  double t2;
  double t21;
  double t31;
  double t35;
  double t4;
  double t41;
  double t5;
  double t6;
  double t9;
  t1 = x * x;
  t2 = y * y;
  t4 = t1 + t2 + omega * t;
  t5 = sin(t4);
  t6 = t5 * t5;
  t9 = 0.4e1 * epsilon;
  t15 = 0.2e1 * epsilon;
  t18 = cos(t4);
  t21 = sin(0.2e1 * t4);
  t31 = 0.3e1 * t21;
  t35 = epsilon * epsilon;
  t41 = v_0 * v_0;
  return(-(0.6e1 * t6 * t5 + (0.6e1 + t9) * t6 + (-0.4e1 + 0.3e1 * epsilon) * t5 - epsilon * (0.3e1 + t15) * t18 - 0.2e1 * epsilon * t21 - 0.3e1 - t15) * rho_0 * u_0 * v_0 * x - (0.8e1 * epsilon * t6 + (t31 + 0.6e1 * epsilon) * t5 + (-0.2e1 * t35 - 0.2e1) * t18 + t31 - t9) * rho_0 * t41 * y - (-0.1e1 + 0.2e1 * t6 + 0.3e1 / 0.2e1 * t5 - epsilon * t18) * rho_0 * v_0 * omega - 0.2e1 * p_0 * y * t5);
}
