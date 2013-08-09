#include <math.h>

double SourceQ_u (double x, double t)
{
  double RHO;
  double P;
  double U;
  double T;
  double MU;
  double DMu_Dx;
  double Q_u;
  double Q_u_convection;
  double Q_u_pressure;
  double Q_u_viscous;
  double Q_u_time;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_t * sin(a_rhot * PI * t / Lt);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / Lt);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_t * cos(a_pt * PI * t / Lt);
  T = P / RHO / R;
  MU = A_mu * pow(T, 0.3e1 / 0.2e1) / (T + B_mu);
  DMu_Dx = a_rhox * PI * rho_x * MU * MU * cos(a_rhox * PI * x / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhox * PI * rho_x * MU * cos(a_rhox * PI * x / L) / L / RHO + a_px * PI * p_x * MU * MU * sin(a_px * PI * x / L) / A_mu / L / R / RHO * pow(T, -0.3e1 / 0.2e1) - 0.3e1 / 0.2e1 * a_px * PI * p_x * MU * sin(a_px * PI * x / L) / L / R / RHO / T;
  Q_u_convection = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L + 0.2e1 * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L;
  Q_u_pressure = -a_px * PI * p_x * sin(a_px * PI * x / L) / L;
  Q_u_viscous = 0.4e1 / 0.3e1 * MU * a_ux * a_ux * PI * PI * u_x * sin(a_ux * PI * x / L) * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * DMu_Dx * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L;
  Q_u_time = a_rhot * PI * rho_t * U * cos(a_rhot * PI * t / Lt) / Lt - a_ut * PI * u_t * RHO * sin(a_ut * PI * t / Lt) / Lt;
  Q_u = Q_u_convection + Q_u_pressure + Q_u_viscous + Q_u_time;
  return(Q_u);
}
