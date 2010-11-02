#include <math.h>

double SourceQ_u (double x, double y, double z)
{
  double Q_u;
  double RHO;
  double P;
  double U;
  double V;
  double W;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L);
  W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L);
  Q_u = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * U * W * cos(a_rhoz * PI * z / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L - a_uz * PI * u_z * RHO * W * sin(a_uz * PI * z / L) / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * U / L;
  return(Q_u);
}
