#include <math.h>

double SourceQ_rho (double x, double y, double t)
{
  double RHO;
  double U;
  double V;
  double Q_rho;
  double Q_rho_convection;
  double Q_rho_time;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / Lt);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / Lt);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / Lt);
  Q_rho_convection = cos(a_rhox * PI * x / L) * a_rhox * PI * rho_x * U / L - sin(a_rhoy * PI * y / L) * a_rhoy * PI * rho_y * V / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO / L;
  Q_rho_time = a_rhot * PI * rho_t * cos(a_rhot * PI * t / Lt) / Lt;
  Q_rho = Q_rho_convection + Q_rho_time;
  return(Q_rho);
}
