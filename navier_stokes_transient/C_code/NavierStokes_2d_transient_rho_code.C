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
  double Phi;
  Phi = x * x + y * y + omega * t;
  Q_rho = (0.4e1 * sin(Phi) + 0.2e1 * epsilon + 0.3e1) * rho_0 * u_0 * x * cos(Phi) - (0.4e1 * pow(sin(Phi), 0.2e1) + 0.3e1 * sin(Phi) - 0.2e1 * epsilon * cos(Phi) - 0.2e1) * rho_0 * v_0 * y + rho_0 * omega * cos(Phi);
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
  double Phi;
  double t1;
  double t10;
  double t13;
  double t2;
  double t4;
  t1 = x * x;
  t2 = y * y;
  Phi = t1 + t2 + omega * t;
  t4 = sin(Phi);
  t10 = cos(Phi);
  t13 = t4 * t4;
  return((0.4e1 * t4 + 0.2e1 * epsilon + 0.3e1) * rho_0 * u_0 * x * t10 - (0.4e1 * t13 + 0.3e1 * t4 - 0.2e1 * epsilon * t10 - 0.2e1) * rho_0 * v_0 * y + rho_0 * omega * t10);
}
