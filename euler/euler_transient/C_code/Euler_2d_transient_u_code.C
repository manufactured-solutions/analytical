#include <math.h>

double SourceQ_u (double x, double y, double t)
{
  double Q_u_t;
  double RHO;
  double U;
  double V;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Q_u_t = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U / L + a_rhot * PI * rho_t * U * cos(a_rhot * PI * t / L) / L - a_ut * PI * u_t * RHO * sin(a_ut * PI * t / L) / L;
  return(Q_u_t);
}
