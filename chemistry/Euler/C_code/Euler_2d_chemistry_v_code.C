#include <math.h>

double SourceQ_v (
  double x,
  double y,
  double p_0,
  double p_x,
  double p_y,
  double rho_N_0,
  double rho_N_x,
  double rho_N_y,
  double rho_N2_0,
  double rho_N2_x,
  double rho_N2_y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double a_px,
  double a_py,
  double a_rho_N_x,
  double a_rho_N_y,
  double a_rho_N2_x,
  double a_rho_N2_y,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double L)
{
  double Q_v;
  Q_v = p_y * cos(a_py * PI * y / L) * a_py * PI / L - rho_N2_x * sin(a_rho_N2_x * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rho_N2_x * PI / L + rho_N2_y * cos(a_rho_N2_y * PI * y / L) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L), 0.2e1) * a_rho_N2_y * PI / L + rho_N_x * cos(a_rho_N_x * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rho_N_x * PI / L - rho_N_y * sin(a_rho_N_y * PI * y / L) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L), 0.2e1) * a_rho_N_y * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_y * cos(a_rho_N_y * PI * y / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_y * sin(a_rho_N2_y * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_ux * PI / L - v_x * sin(a_vx * PI * x / L) * (rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_y * cos(a_rho_N_y * PI * y / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_y * sin(a_rho_N2_y * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vx * PI / L + 0.2e1 * v_y * cos(a_vy * PI * y / L) * (rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_y * cos(a_rho_N_y * PI * y / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_y * sin(a_rho_N2_y * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_vy * PI / L;
  return(Q_v);
}
