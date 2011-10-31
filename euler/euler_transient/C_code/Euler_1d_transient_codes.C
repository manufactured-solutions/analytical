#include <math.h>

double SourceQ_rho (double x, double t)
{
  double Q_rho_t;
  double RHO;
  double U;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);
  Q_rho_t = cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t / L + cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x / L + cos(a_rhox * PI * x / L) * U * a_rhox * PI * rho_x / L;
  return(Q_rho_t);
}
#include <math.h>

double SourceQ_u (double x, double t)
{
  double Q_u_t;
  double RHO;
  double U;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);
  Q_u_t = cos(a_rhox * PI * x / L) * a_rhox * PI * rho_x * U * U / L + 0.2e1 * cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x * U / L + cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t * U / L - sin(a_ut * PI * t / L) * a_ut * PI * u_t * RHO / L - sin(a_px * PI * x / L) * a_px * PI * p_x / L;
  return(Q_u_t);
}
#include <math.h>

double SourceQ_e (double x, double t)
{
  double Q_e_t;
  double RHO;
  double U;
  double P;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  P = p_0 + p_x * cos(a_px * pi * x / L) + p_t * cos(a_pt * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);
  Q_e_t = cos(a_rhot * PI * t / L) * a_rhot * PI * rho_t * U * U / L / 0.2e1 - sin(a_ut * PI * t / L) * a_ut * PI * u_t * RHO * U / L - sin(a_pt * PI * t / L) * a_pt * PI * p_t / (Gamma - 0.1e1) / L + cos(a_rhox * PI * x / L) * pow(U, 0.3e1) * a_rhox * PI * rho_x / L / 0.2e1 + cos(a_ux * PI * x / L) * P * a_ux * PI * u_x * Gamma / (Gamma - 0.1e1) / L + 0.3e1 / 0.2e1 * cos(a_ux * PI * x / L) * RHO * U * U * a_ux * PI * u_x / L - sin(a_px * PI * x / L) * U * a_px * PI * p_x * Gamma / (Gamma - 0.1e1) / L;
  return(Q_e_t);
}
