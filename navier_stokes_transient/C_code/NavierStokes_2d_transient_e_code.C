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
  double mu,
  double Gamma)
{
  double Q_e;
  double Phi;
  Phi = x * x + y * y + omega * t;
  Q_e = (epsilon + 0.4e1 * sin(Phi) + 0.9e1 / 0.2e1) * pow(epsilon + sin(Phi), 0.2e1) * rho_0 * pow(u_0, 0.3e1) * x * cos(Phi) + (epsilon / 0.2e1 + sin(Phi) / 0.2e1) * (-0.8e1 * pow(sin(Phi), 0.3e1) + 0.3e1 * epsilon * sin(0.2e1 * Phi) + (0.2e1 * epsilon + 0.6e1) * epsilon * cos(Phi) + (-0.4e1 * epsilon - 0.9e1) * pow(sin(Phi), 0.2e1) + 0.2e1 * epsilon + (-0.3e1 * epsilon + 0.6e1) * sin(Phi) + 0.6e1) * rho_0 * u_0 * u_0 * v_0 * y + (epsilon / 0.2e1 + sin(Phi) / 0.2e1) * (epsilon + 0.3e1 * sin(Phi) + 0.3e1) * rho_0 * u_0 * u_0 * omega * cos(Phi) + (epsilon / 0.2e1 + cos(Phi) / 0.2e1) * (-0.8e1 * pow(sin(Phi), 0.3e1) + 0.2e1 * epsilon * sin(0.2e1 * Phi) + (0.2e1 * epsilon + 0.3e1) * epsilon * cos(Phi) + (-0.6e1 * epsilon - 0.9e1) * pow(sin(Phi), 0.2e1) + 0.2e1 * epsilon + (-0.6e1 * epsilon + 0.4e1) * sin(Phi) + 0.3e1) * rho_0 * u_0 * v_0 * v_0 * x + (epsilon * cos(Phi) - 0.4e1 * pow(sin(Phi), 0.2e1) - 0.9e1 / 0.2e1 * sin(Phi) + 0.1e1) * pow(epsilon + cos(Phi), 0.2e1) * rho_0 * pow(v_0, 0.3e1) * y + (epsilon / 0.2e1 + cos(Phi) / 0.2e1) * (epsilon * cos(Phi) - 0.3e1 * pow(sin(Phi), 0.2e1) - 0.3e1 * sin(Phi) + 0.1e1) * rho_0 * v_0 * v_0 * omega - (0.2e1 * epsilon + 0.4e1 * cos(Phi) + 0.4e1) * Gamma * p_0 * v_0 * y * sin(Phi) / (Gamma - 0.1e1) + 0.4e1 / 0.3e1 * mu * (epsilon * sin(Phi) + epsilon * cos(Phi) + 0.2e1 * sin(0.2e1 * Phi)) * u_0 * v_0 * x * y + Gamma * (-0.4e1 * epsilon * sin(Phi) - 0.8e1 * pow(sin(Phi), 0.2e1) + 0.8e1 * cos(Phi) + 0.4e1) * p_0 * u_0 * x / (Gamma - 0.1e1) / 0.2e1 - p_0 * omega * sin(Phi) / (Gamma - 0.1e1) + 0.2e1 / 0.3e1 * mu * (-0.7e1 * epsilon * cos(Phi) - 0.8e1 * x * x - 0.6e1 * y * y + (0.16e2 * x * x + 0.12e2 * y * y) * pow(sin(Phi), 0.2e1) + (0.8e1 * epsilon * x * x + 0.6e1 * epsilon * y * y) * sin(Phi) - 0.7e1 / 0.2e1 * sin(0.2e1 * Phi)) * u_0 * u_0 - 0.2e1 / 0.3e1 * mu * (-0.7e1 * epsilon * sin(Phi) - 0.6e1 * x * x - 0.8e1 * y * y + (0.12e2 * x * x + 0.16e2 * y * y) * pow(sin(Phi), 0.2e1) - 0.7e1 / 0.2e1 * sin(0.2e1 * Phi) + (-0.6e1 * epsilon * x * x - 0.8e1 * epsilon * y * y) * cos(Phi)) * v_0 * v_0 + 0.8e1 * k * (-0.16e2 * x * x - 0.16e2 * y * y + (0.8e1 * x * x + 0.8e1 * y * y + 0.6e1) * pow(sin(Phi), 0.2e1) + (-0.12e2 * x * x - 0.12e2 * y * y + 0.13e2) * sin(Phi) + (-0.3e1 * x * x - 0.3e1 * y * y + 0.4e1) * sin(0.2e1 * Phi) + (x * x + y * y + 0.12e2) * cos(Phi) + 0.6e1) * p_0 / R * pow(0.2e1 * sin(Phi) + 0.3e1, -0.3e1) / rho_0;
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
  double mu,
  double Gamma)
{
  double Phi;
  double t1;
  double t11;
  double t114;
  double t116;
  double t117;
  double t120;
  double t122;
  double t126;
  double t133;
  double t134;
  double t135;
  double t14;
  double t161;
  double t162;
  double t17;
  double t18;
  double t2;
  double t20;
  double t22;
  double t23;
  double t25;
  double t4;
  double t41;
  double t48;
  double t49;
  double t54;
  double t62;
  double t66;
  double t7;
  double t70;
  double t8;
  double t88;
  double t89;
  double t92;
  t1 = x * x;
  t2 = y * y;
  Phi = t1 + t2 + omega * t;
  t4 = sin(Phi);
  t7 = epsilon + t4;
  t8 = t7 * t7;
  t11 = u_0 * u_0;
  t14 = cos(Phi);
  t17 = t7 / 0.2e1;
  t18 = t4 * t4;
  t20 = 0.8e1 * t18 * t4;
  t22 = sin(0.2e1 * Phi);
  t23 = epsilon * t22;
  t25 = 0.2e1 * epsilon;
  t41 = 0.3e1 * t4;
  t48 = epsilon + t14;
  t49 = t48 / 0.2e1;
  t54 = 0.6e1 * epsilon;
  t62 = v_0 * v_0;
  t66 = epsilon * t14;
  t70 = t48 * t48;
  t88 = 0.1e1 / (Gamma - 0.1e1);
  t89 = t4 * t88;
  t92 = epsilon * t4;
  t114 = 0.8e1 * t1;
  t116 = 0.16e2 * t1;
  t117 = 0.12e2 * t2;
  t120 = epsilon * t1;
  t122 = epsilon * t2;
  t126 = 0.7e1 / 0.2e1 * t22;
  t133 = 0.8e1 * t2;
  t134 = 0.12e2 * t1;
  t135 = 0.16e2 * t2;
  t161 = 0.2e1 * t4 + 0.3e1;
  t162 = t161 * t161;
  return((0.9e1 / 0.2e1 + epsilon + 0.4e1 * t4) * t8 * rho_0 * t11 * u_0 * x * t14 + t17 * (-t20 + 0.3e1 * t23 + (t25 + 0.6e1) * epsilon * t14 + (-0.4e1 * epsilon - 0.9e1) * t18 + t25 + (-0.3e1 * epsilon + 0.6e1) * t4 + 0.6e1) * rho_0 * t11 * v_0 * y + t17 * (epsilon + t41 + 0.3e1) * rho_0 * t11 * omega * t14 + t49 * (-t20 + 0.2e1 * t23 + (t25 + 0.3e1) * epsilon * t14 + (-t54 - 0.9e1) * t18 + t25 + (-t54 + 0.4e1) * t4 + 0.3e1) * rho_0 * u_0 * t62 * x + (0.1e1 + t66 - 0.4e1 * t18 - 0.9e1 / 0.2e1 * t4) * t70 * rho_0 * t62 * v_0 * y + t49 * (t66 - 0.3e1 * t18 - t41 + 0.1e1) * rho_0 * t62 * omega - (0.4e1 + t25 + 0.4e1 * t14) * Gamma * p_0 * v_0 * y * t89 + 0.4e1 / 0.3e1 * mu * (t92 + t66 + 0.2e1 * t22) * u_0 * v_0 * x * y + Gamma * (-0.4e1 * t92 - 0.8e1 * t18 + 0.8e1 * t14 + 0.4e1) * p_0 * u_0 * x * t88 / 0.2e1 - p_0 * omega * t89 + 0.2e1 / 0.3e1 * mu * (-0.7e1 * t66 - t114 - 0.6e1 * t2 + (t116 + t117) * t18 + (0.8e1 * t120 + 0.6e1 * t122) * t4 - t126) * t11 - 0.2e1 / 0.3e1 * mu * (-0.7e1 * t92 - 0.6e1 * t1 - t133 + (t134 + t135) * t18 - t126 + (-0.6e1 * t120 - 0.8e1 * t122) * t14) * t62 + 0.8e1 * k * (-t116 - t135 + (t114 + t133 + 0.6e1) * t18 + (-t134 - t117 + 0.13e2) * t4 + (-0.3e1 * t1 - 0.3e1 * t2 + 0.4e1) * t22 + (t1 + t2 + 0.12e2) * t14 + 0.6e1) * p_0 / R / t162 / t161 / rho_0);
}
