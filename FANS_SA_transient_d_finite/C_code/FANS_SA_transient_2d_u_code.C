#include <math.h>

double SourceQ_u (
  double x,
  double y,
  double t,
  double mu,
  double c_v1)
{
  double Q_u;
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
  Q_u = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L + (0.4e1 / 0.3e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * f_v1 * PI * PI * RHO * NU_SA * pow(L, -0.2e1) + (0.4e1 / 0.3e1 * a_ux * a_nusax * u_x * nu_sa_x * cos(a_ux * PI * x / L) * sin(a_nusax * PI * x / L) - a_uy * a_nusay * u_y * nu_sa_y * sin(a_uy * PI * y / L) * sin(a_nusay * PI * y / L) - a_vx * a_nusay * v_x * nu_sa_y * sin(a_vx * PI * x / L) * sin(a_nusay * PI * y / L) - 0.2e1 / 0.3e1 * a_vy * a_nusax * v_y * nu_sa_x * cos(a_vy * PI * y / L) * sin(a_nusax * PI * x / L)) * f_v1 * PI * PI * RHO * pow(L, -0.2e1) + (-0.4e1 / 0.3e1 * a_rhox * a_ux * rho_x * u_x * cos(a_rhox * PI * x / L) * cos(a_ux * PI * x / L) + 0.2e1 / 0.3e1 * a_rhox * a_vy * rho_x * v_y * cos(a_rhox * PI * x / L) * cos(a_vy * PI * y / L) - a_rhoy * a_uy * rho_y * u_y * sin(a_rhoy * PI * y / L) * sin(a_uy * PI * y / L) - a_rhoy * a_vx * rho_y * v_x * sin(a_rhoy * PI * y / L) * sin(a_vx * PI * x / L)) * f_v1 * PI * PI * NU_SA * pow(L, -0.2e1) + (pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1)) + f_v1) * (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * PI * PI * mu * pow(L, -0.2e1) / 0.3e1 + a_rhot * PI * rho_t * U * cos(a_rhot * PI * t / L) / L - a_ut * PI * u_t * RHO * sin(a_ut * PI * t / L) / L;
  return(Q_u);
}
