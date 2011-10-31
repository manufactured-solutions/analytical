#include <math.h>

double SourceQ_nu (
  double x,
  double y,
  double t,
  double mu,
  double c_b1,
  double c_b2,
  double sigma)
{
  double Q_nu;
  double RHO;
  double U;
  double V;
  double NU_SA;
  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  Q_nu = a_rhox * PI * rho_x * U * NU_SA * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * NU_SA * sin(a_rhoy * PI * y / L) / L - a_nusax * PI * nu_sa_x * RHO * U * sin(a_nusax * PI * x / L) / L - a_nusay * PI * nu_sa_y * RHO * V * sin(a_nusay * PI * y / L) / L + a_rhot * PI * rho_t * NU_SA * cos(a_rhot * PI * t / L) / L - a_nusat * PI * nu_sa_t * RHO * sin(a_nusat * PI * t / L) / L - (a_nusax * a_nusax * nu_sa_x * nu_sa_x * pow(sin(a_nusax * PI * x / L), 0.2e1) + a_nusay * a_nusay * nu_sa_y * nu_sa_y * pow(sin(a_nusay * PI * y / L), 0.2e1)) * (0.1e1 + c_b2) * PI * PI * RHO / sigma * pow(L, -0.2e1) + (a_rhox * a_nusax * rho_x * nu_sa_x * cos(a_rhox * PI * x / L) * sin(a_nusax * PI * x / L) - a_rhoy * a_nusay * rho_y * nu_sa_y * sin(a_rhoy * PI * y / L) * sin(a_nusay * PI * y / L)) * PI * PI * NU_SA / sigma * pow(L, -0.2e1) - c_b1 * sqrt(pow(-a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L), 0.2e1) * pow(L, -0.2e1)) * PI * RHO * NU_SA + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * NU_SA / L + (a_nusax * a_nusax * nu_sa_x * cos(a_nusax * PI * x / L) + a_nusay * a_nusay * nu_sa_y * cos(a_nusay * PI * y / L)) * (RHO * NU_SA + mu) * PI * PI / sigma * pow(L, -0.2e1);
  return(Q_nu);
}
