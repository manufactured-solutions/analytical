#include <math.h>

double SourceQ_e (
  double x,
  double y,
  double t,
  double rho_0,
  double p_0,
  double u_0,
  double v_0,
  double omega,
  double epsilon,
  double Gamma)
{
  double Q_e;
  Q_e = (0.4e1 * sin(omega * t + x * x + y * y) + epsilon + 0.9e1 / 0.2e1) * pow(sin(omega * t + x * x + y * y) + epsilon, 0.2e1) * x * cos(omega * t + x * x + y * y) * pow(u_0, 0.3e1) * rho_0 - y * v_0 * rho_0 * u_0 * u_0 * (sin(omega * t + x * x + y * y) + epsilon) * ((-0.2e1 * epsilon * epsilon - 0.6e1 * epsilon) * cos(omega * t + x * x + y * y) + 0.8e1 * pow(sin(omega * t + x * x + y * y), 0.3e1) + (0.4e1 * epsilon + 0.9e1) * pow(sin(omega * t + x * x + y * y), 0.2e1) - 0.3e1 * sin(0.2e1 * omega * t + 0.2e1 * x * x + 0.2e1 * y * y) * epsilon - 0.6e1 - 0.2e1 * epsilon + (0.3e1 * epsilon - 0.6e1) * sin(omega * t + x * x + y * y)) / 0.2e1 + (sin(omega * t + x * x + y * y) / 0.2e1 + epsilon / 0.2e1) * (0.3e1 * sin(omega * t + x * x + y * y) + 0.3e1 + epsilon) * cos(omega * t + x * x + y * y) * omega * rho_0 * u_0 * u_0 - x * u_0 * rho_0 * v_0 * v_0 * (cos(omega * t + x * x + y * y) + epsilon) * ((-0.2e1 * epsilon * epsilon - 0.3e1 * epsilon) * cos(omega * t + x * x + y * y) + 0.8e1 * pow(sin(omega * t + x * x + y * y), 0.3e1) + (0.6e1 * epsilon + 0.9e1) * pow(sin(omega * t + x * x + y * y), 0.2e1) - 0.2e1 * sin(0.2e1 * omega * t + 0.2e1 * x * x + 0.2e1 * y * y) * epsilon - 0.3e1 - 0.2e1 * epsilon + (0.6e1 * epsilon - 0.4e1) * sin(omega * t + x * x + y * y)) / 0.2e1 - y * pow(v_0, 0.3e1) * rho_0 * pow(cos(omega * t + x * x + y * y) + epsilon, 0.2e1) * (0.8e1 * pow(sin(omega * t + x * x + y * y), 0.2e1) - 0.2e1 - 0.2e1 * epsilon * cos(omega * t + x * x + y * y) + 0.9e1 * sin(omega * t + x * x + y * y)) / 0.2e1 - omega * rho_0 * v_0 * v_0 * (cos(omega * t + x * x + y * y) + epsilon) * (0.3e1 * pow(sin(omega * t + x * x + y * y), 0.2e1) - 0.1e1 - epsilon * cos(omega * t + x * x + y * y) + 0.3e1 * sin(omega * t + x * x + y * y)) / 0.2e1 - 0.2e1 * p_0 * x * Gamma * u_0 * (-0.2e1 * cos(omega * t + x * x + y * y) + 0.2e1 * pow(sin(omega * t + x * x + y * y), 0.2e1) - 0.1e1 + sin(omega * t + x * x + y * y) * epsilon) / (Gamma - 0.1e1) - 0.2e1 * y * p_0 * Gamma * sin(omega * t + x * x + y * y) * (0.2e1 * cos(omega * t + x * x + y * y) + epsilon + 0.2e1) * v_0 / (Gamma - 0.1e1) - omega * p_0 * sin(omega * t + x * x + y * y) / (Gamma - 0.1e1);
  return(Q_e);
}
#include <math.h>

double SourceQ_e (
  double x,
  double y,
  double t,
  double rho_0,
  double p_0,
  double u_0,
  double v_0,
  double omega,
  double epsilon,
  double Gamma)
{
  double t12;
  double t13;
  double t2;
  double t21;
  double t22;
  double t23;
  double t26;
  double t28;
  double t3;
  double t33;
  double t34;
  double t36;
  double t37;
  double t4;
  double t45;
  double t49;
  double t5;
  double t54;
  double t55;
  double t56;
  double t70;
  double t73;
  double t8;
  double t87;
  double t9;
  double t93;
  t2 = x * x;
  t3 = y * y;
  t4 = omega * t + t2 + t3;
  t5 = sin(t4);
  t8 = t5 + epsilon;
  t9 = t8 * t8;
  t12 = cos(t4);
  t13 = u_0 * u_0;
  t21 = epsilon * epsilon;
  t22 = 0.2e1 * t21;
  t23 = 0.6e1 * epsilon;
  t26 = t5 * t5;
  t28 = 0.8e1 * t26 * t5;
  t33 = sin(0.2e1 * t4);
  t34 = t33 * epsilon;
  t36 = 0.2e1 * epsilon;
  t37 = 0.3e1 * epsilon;
  t45 = 0.3e1 * t5;
  t49 = omega * rho_0;
  t54 = v_0 * v_0;
  t55 = t12 + epsilon;
  t56 = t54 * t55;
  t70 = t55 * t55;
  t73 = epsilon * t12;
  t87 = 0.2e1 * t12;
  t93 = 0.1e1 / (Gamma - 0.1e1);
  return((0.9e1 / 0.2e1 + 0.4e1 * t5 + epsilon) * t9 * x * t12 * t13 * u_0 * rho_0 - y * v_0 * rho_0 * t13 * t8 * ((-t22 - t23) * t12 + t28 + (0.4e1 * epsilon + 0.9e1) * t26 - 0.3e1 * t34 - 0.6e1 - t36 + (t37 - 0.6e1) * t5) / 0.2e1 + t8 * (t45 + 0.3e1 + epsilon) * t12 * t49 * t13 / 0.2e1 - x * u_0 * rho_0 * t56 * ((-t22 - t37) * t12 + t28 + (t23 + 0.9e1) * t26 - 0.2e1 * t34 - 0.3e1 - t36 + (t23 - 0.4e1) * t5) / 0.2e1 - y * t54 * v_0 * rho_0 * t70 * (0.8e1 * t26 - 0.2e1 - 0.2e1 * t73 + 0.9e1 * t5) / 0.2e1 - t49 * t56 * (0.3e1 * t26 - 0.1e1 - t73 + t45) / 0.2e1 - 0.2e1 * p_0 * x * Gamma * u_0 * (-t87 + 0.2e1 * t26 - 0.1e1 + t5 * epsilon) * t93 - 0.2e1 * y * p_0 * Gamma * t5 * (t87 + epsilon + 0.2e1) * v_0 * t93 - omega * p_0 * t5 * t93);
}
