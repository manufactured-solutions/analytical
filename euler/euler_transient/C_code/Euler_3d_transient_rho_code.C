#include <math.h>

double SourceQ_rho (
  double x,
  double y,
  double z,
  double t)
{
  double Q_rho_t;
  double RHO;
  double U;
  double V;
  double W;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_t * sin(a_vt * PI * t / L);
  W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / L);
  Q_rho_t = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO / L + a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L;
  return(Q_rho_t);
}
