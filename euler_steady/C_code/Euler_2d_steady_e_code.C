#include <math.h>

double SourceQ_e (double x, double y)
{
  double Qe;
  double RHO;
  double U;
  double V;
  double P;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L);
  Qe = -a_px * PI * p_x * Gamma * U * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + a_py * PI * p_y * Gamma * V * cos(a_py * PI * y / L) / (Gamma - 0.1e1) / L + (U * U + V * V) * a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L / 0.2e1 - (U * U + V * V) * a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L / 0.2e1 + (0.3e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * PI * RHO * U * V / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V * V / L / 0.2e1 + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * Gamma * P / (Gamma - 0.1e1) / L;
  return(Qe);
}
