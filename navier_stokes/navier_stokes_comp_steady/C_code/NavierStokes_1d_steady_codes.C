#include <math.h>

double SourceQ_rho (double x)
{
  double Q_rho;
  double RHO;
  double U;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  Q_rho = cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x / L + cos(a_rhox * PI * x / L) * U * a_rhox * PI * rho_x / L;
  return(Q_rho);
}
#include <math.h>

double SourceQ_u (double x)
{
  double Q_u;
  double RHO;
  double U;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  Q_u = 0.4e1 / 0.3e1 * a_ux * a_ux * PI * PI * u_x * mu * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L + 0.2e1 * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L;
  return(Q_u);
}
#include <math.h>

double SourceQ_e (double x)
{
  double Q_e;
  double RHO;
  double U;
  double P;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  P = p_0 + p_x * cos(a_px * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  Q_e = -0.4e1 / 0.3e1 * a_ux * a_ux * PI * PI * u_x * u_x * mu * pow(cos(a_ux * PI * x / L), 0.2e1) * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * a_ux * a_ux * PI * PI * u_x * mu * U * sin(a_ux * PI * x / L) * pow(L, -0.2e1) - 0.2e1 * a_rhox * a_rhox * PI * PI * rho_x * rho_x * k * P * pow(cos(a_rhox * PI * x / L), 0.2e1) * pow(L, -0.2e1) / R * pow(RHO, -0.3e1) - 0.2e1 * a_rhox * a_px * PI * PI * rho_x * p_x * k * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) + a_rhox * PI * rho_x * pow(U, 0.3e1) * cos(a_rhox * PI * x / L) / L / 0.2e1 + 0.3e1 / 0.2e1 * a_ux * PI * u_x * RHO * U * U * cos(a_ux * PI * x / L) / L - a_rhox * a_rhox * PI * PI * rho_x * k * P * sin(a_rhox * PI * x / L) * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) + a_px * a_px * PI * PI * p_x * k * cos(a_px * PI * x / L) * pow(L, -0.2e1) / R / RHO - a_px * PI * p_x * Gamma * U * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + a_ux * PI * u_x * Gamma * P * cos(a_ux * PI * x / L) / (Gamma - 0.1e1) / L;
  return(Q_e);
}
