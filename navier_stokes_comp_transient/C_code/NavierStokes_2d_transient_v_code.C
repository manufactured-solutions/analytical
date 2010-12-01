#include <math.h>

double SourceQ_v (double x, double y, double t)
{
  double Q_v;
  double RHO;
  double U;
  double V;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Q_v = a_vx * a_vx * PI * PI * v_x * mu * cos(a_vx * PI * x / L) * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * a_vy * a_vy * PI * PI * v_y * mu * sin(a_vy * PI * y / L) * pow(L, -0.2e1) + a_rhox * PI * rho_x * U * V * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * sin(a_rhoy * PI * y / L) / L - a_vx * PI * v_x * RHO * U * sin(a_vx * PI * x / L) / L + a_rhot * PI * rho_t * V * cos(a_rhot * PI * t / L) / L + a_vt * PI * v_t * RHO * cos(a_vt * PI * t / L) / L + a_py * PI * p_y * cos(a_py * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V / L;
  return(Q_v);
}
