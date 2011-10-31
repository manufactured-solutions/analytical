#include <math.h>

double SourceQ_v (
  double x,
  double y,
  double t,
  double mu,
  double c_v1)
{
  double Q_v;
  double RHO;
  double U;
  double V;
  double NU_SA;
  double f_v1;
  double chi;
  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  Q_v = a_rhox * PI * rho_x * U * V * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * sin(a_rhoy * PI * y / L) / L - a_vx * PI * v_x * RHO * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V / L + a_py * PI * p_y * cos(a_py * PI * y / L) / L + (a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 / 0.3e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * f_v1 * PI * PI * RHO * NU_SA * pow(L, -0.2e1) + (-0.2e1 / 0.3e1 * a_ux * a_nusay * u_x * nu_sa_y * cos(a_ux * PI * x / L) * sin(a_nusay * PI * y / L) - a_uy * a_nusax * u_y * nu_sa_x * sin(a_uy * PI * y / L) * sin(a_nusax * PI * x / L) - a_vx * a_nusax * v_x * nu_sa_x * sin(a_vx * PI * x / L) * sin(a_nusax * PI * x / L) + 0.4e1 / 0.3e1 * a_vy * a_nusay * v_y * nu_sa_y * cos(a_vy * PI * y / L) * sin(a_nusay * PI * y / L)) * f_v1 * PI * PI * RHO * pow(L, -0.2e1) + (a_rhox * a_uy * rho_x * u_y * cos(a_rhox * PI * x / L) * sin(a_uy * PI * y / L) + a_rhox * a_vx * rho_x * v_x * cos(a_rhox * PI * x / L) * sin(a_vx * PI * x / L) - 0.2e1 / 0.3e1 * a_rhoy * a_ux * rho_y * u_x * sin(a_rhoy * PI * y / L) * cos(a_ux * PI * x / L) + 0.4e1 / 0.3e1 * a_rhoy * a_vy * rho_y * v_y * sin(a_rhoy * PI * y / L) * cos(a_vy * PI * y / L)) * f_v1 * PI * PI * NU_SA * pow(L, -0.2e1) + (pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1)) + f_v1) * (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * PI * PI * mu * pow(L, -0.2e1) / 0.3e1 + a_rhot * PI * rho_t * V * cos(a_rhot * PI * t / L) / L + a_vt * PI * v_t * RHO * cos(a_vt * PI * t / L) / L;
  return(Q_v);
}
