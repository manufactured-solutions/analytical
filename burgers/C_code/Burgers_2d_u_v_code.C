#include <math.h>

double SourceQ_u (
  double x,
  double y,
  double t,
  double u_0,
  double v_0,
  double omega,
  double epsilon,
  double nu)
{
  double Q_u;
  Q_u = (0.2e1 * sin(0.2e1 * x * x + 0.2e1 * y * y + 0.2e1 * omega * t) + 0.4e1 * epsilon * cos(x * x + y * y + omega * t)) * u_0 * u_0 * x + (-0.4e1 * pow(sin(x * x + y * y + omega * t), 0.2e1) - 0.2e1 * epsilon * sin(x * x + y * y + omega * t) + 0.2e1 * epsilon * cos(x * x + y * y + omega * t) + 0.2e1) * u_0 * v_0 * y + u_0 * cos(x * x + y * y + omega * t) * omega + (0.4e1 * x * x * sin(x * x + y * y + omega * t) + 0.4e1 * y * y * sin(x * x + y * y + omega * t) - 0.4e1 * cos(x * x + y * y + omega * t)) * u_0 * nu;
  return(Q_u);
}
#include <math.h>

double SourceQ_v (
  double x,
  double y,
  double t,
  double u_0,
  double v_0,
  double omega,
  double epsilon,
  double nu)
{
  double Q_v;
  Q_v = -(0.4e1 * epsilon * sin(x * x + y * y + omega * t) + 0.2e1 * sin(0.2e1 * x * x + 0.2e1 * y * y + 0.2e1 * omega * t)) * v_0 * v_0 * y + (-0.4e1 * pow(sin(x * x + y * y + omega * t), 0.2e1) - 0.2e1 * epsilon * sin(x * x + y * y + omega * t) + 0.2e1 * epsilon * cos(x * x + y * y + omega * t) + 0.2e1) * u_0 * v_0 * x - v_0 * sin(x * x + y * y + omega * t) * omega + (0.4e1 * x * x * cos(x * x + y * y + omega * t) + 0.4e1 * y * y * cos(x * x + y * y + omega * t) + 0.4e1 * sin(x * x + y * y + omega * t)) * v_0 * nu;
  return(Q_v);
}
