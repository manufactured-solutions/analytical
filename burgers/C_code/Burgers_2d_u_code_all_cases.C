#include <math.h>

double SourceQ_u_transient_viscous (
  double x,
  double y,
  double t,
  double nu)
{
  double Qu_tv;
  double U;
  double V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qu_tv = -a_ut * PI * u_t * sin(a_ut * PI * t / L) / L - a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L + a_ux * a_ux * PI * PI * u_x * nu * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_uy * a_uy * PI * PI * u_y * nu * cos(a_uy * PI * y / L) * pow(L, -0.2e1);
  return(Qu_tv);
}
#include <math.h>

double SourceQ_u_steady_viscous (double x, double y, double nu)
{
  double Qu_sv;
  double U;
  double V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qu_sv = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L + a_ux * a_ux * PI * PI * u_x * nu * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_uy * a_uy * PI * PI * u_y * nu * cos(a_uy * PI * y / L) * pow(L, -0.2e1);
  return(Qu_sv);
}
#include <math.h>

double SourceQ_u_transient_inviscid (double x, double y, double t)
{
  double Qu_t;
  double U;
  double V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qu_t = -a_ut * PI * u_t * sin(a_ut * PI * t / L) / L - a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;
  return(Qu_t);
}
#include <math.h>

double SourceQ_u_steady_inviscid (double x, double y)
{
  double Qu_s;
  double U;
  double V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qu_s = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;
  return(Qu_s);
}
