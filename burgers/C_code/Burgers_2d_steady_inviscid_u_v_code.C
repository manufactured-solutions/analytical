#include <math.h>

double SourceQ_u (
  double x,
  double y,
  double u_0,
  double v_0,
  double epsilon)
{
  double Q_u;
  Q_u = (0.2e1 * sin(0.2e1 * x * x + 0.2e1 * y * y) + 0.4e1 * epsilon * cos(x * x + y * y)) * u_0 * u_0 * x + (-0.4e1 * pow(sin(x * x + y * y), 0.2e1) - 0.2e1 * epsilon * sin(x * x + y * y) + 0.2e1 * epsilon * cos(x * x + y * y) + 0.2e1) * u_0 * v_0 * y;
  return(Q_u);
}
#include <math.h>

double SourceQ_v (
  double x,
  double y,
  double u_0,
  double v_0,
  double epsilon)
{
  double Q_v;
  Q_v = -(0.4e1 * epsilon * sin(x * x + y * y) + 0.2e1 * sin(0.2e1 * x * x + 0.2e1 * y * y)) * v_0 * v_0 * y + (-0.4e1 * pow(sin(x * x + y * y), 0.2e1) - 0.2e1 * epsilon * sin(x * x + y * y) + 0.2e1 * epsilon * cos(x * x + y * y) + 0.2e1) * u_0 * v_0 * x;
  return(Q_v);
}
