#include <math.h>

double SourceQ_v_transient_viscous (
  double x,
  double y,
  double t,
  double nu)
{
  double Qv_tv;
  double U;
  double V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qv_tv = a_vt * PI * v_t * cos(a_vt * PI * t / L) / L - a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L + a_vx * a_vx * PI * PI * v_x * nu * cos(a_vx * PI * x / L) * pow(L, -0.2e1) + a_vy * a_vy * PI * PI * v_y * nu * sin(a_vy * PI * y / L) * pow(L, -0.2e1);
  return(Qv_tv);
}
#include <math.h>

double SourceQ_v_steady_viscous (double x, double y, double nu)
{
  double Qv_sv;
  double U;
  double V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qv_sv = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L + a_vx * a_vx * PI * PI * v_x * nu * cos(a_vx * PI * x / L) * pow(L, -0.2e1) + a_vy * a_vy * PI * PI * v_y * nu * sin(a_vy * PI * y / L) * pow(L, -0.2e1);
  return(Qv_sv);
}
#include <math.h>

double SourceQ_v_transient_inviscid (double x, double y, double t)
{
  double Qv_t;
  double U;
  double V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qv_t = a_vt * PI * v_t * cos(a_vt * PI * t / L) / L - a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;
  return(Qv_t);
}
#include <math.h>

double SourceQ_v_steady_inviscid (double x, double y)
{
  double Qv_s;
  double U;
  double V;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Qv_s = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;
  return(Qv_s);
}
