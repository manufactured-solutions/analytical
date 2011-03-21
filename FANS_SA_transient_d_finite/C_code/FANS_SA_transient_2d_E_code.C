#include <math.h>

double SourceQ_e (
  double x,
  double y,
  double t,
  double mu,
  double c_v1,
  double cp,
  double cv,
  double Pr_t,
  double Pr)
{
  double Q_E;
  double RHO;
  double U;
  double V;
  double P;
  double NU_SA;
  double chi;
  double f_v1;
  double R;
  double mu_t;
  NU_SA = nu_sa_0 + nu_sa_x * cos(a_nusax * PI * x / L) + nu_sa_y * cos(a_nusay * PI * y / L) + nu_sa_t * cos(a_nusat * PI * t / L);
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_t * cos(a_pt * PI * t / L);
  chi = RHO * NU_SA / mu;
  f_v1 = pow(chi, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1));
  mu_t = RHO * NU_SA * f_v1;
  R = cp - cv;
  Q_E = (0.4e1 / 0.3e1 * a_ux * a_nusax * u_x * nu_sa_x * cos(a_ux * PI * x / L) * sin(a_nusax * PI * x / L) - a_uy * a_nusay * u_y * nu_sa_y * sin(a_uy * PI * y / L) * sin(a_nusay * PI * y / L) - a_vx * a_nusay * v_x * nu_sa_y * sin(a_vx * PI * x / L) * sin(a_nusay * PI * y / L) - 0.2e1 / 0.3e1 * a_vy * a_nusax * v_y * nu_sa_x * cos(a_vy * PI * y / L) * sin(a_nusax * PI * x / L)) * PI * PI * f_v1 * RHO * U * pow(L, -0.2e1) + (-0.2e1 / 0.3e1 * a_ux * a_nusay * u_x * nu_sa_y * cos(a_ux * PI * x / L) * sin(a_nusay * PI * y / L) - a_uy * a_nusax * u_y * nu_sa_x * sin(a_uy * PI * y / L) * sin(a_nusax * PI * x / L) - a_vx * a_nusax * v_x * nu_sa_x * sin(a_vx * PI * x / L) * sin(a_nusax * PI * x / L) + 0.4e1 / 0.3e1 * a_vy * a_nusay * v_y * nu_sa_y * cos(a_vy * PI * y / L) * sin(a_nusay * PI * y / L)) * PI * PI * f_v1 * RHO * V * pow(L, -0.2e1) + (-0.4e1 / 0.3e1 * a_rhox * a_ux * rho_x * u_x * cos(a_rhox * PI * x / L) * cos(a_ux * PI * x / L) + 0.2e1 / 0.3e1 * a_rhox * a_vy * rho_x * v_y * cos(a_rhox * PI * x / L) * cos(a_vy * PI * y / L) - a_rhoy * a_uy * rho_y * u_y * sin(a_rhoy * PI * y / L) * sin(a_uy * PI * y / L) - a_rhoy * a_vx * rho_y * v_x * sin(a_rhoy * PI * y / L) * sin(a_vx * PI * x / L)) * PI * PI * f_v1 * U * NU_SA * pow(L, -0.2e1) + (a_rhox * a_uy * rho_x * u_y * cos(a_rhox * PI * x / L) * sin(a_uy * PI * y / L) + a_rhox * a_vx * rho_x * v_x * cos(a_rhox * PI * x / L) * sin(a_vx * PI * x / L) - 0.2e1 / 0.3e1 * a_rhoy * a_ux * rho_y * u_x * sin(a_rhoy * PI * y / L) * cos(a_ux * PI * x / L) + 0.4e1 / 0.3e1 * a_rhoy * a_vy * rho_y * v_y * sin(a_rhoy * PI * y / L) * cos(a_vy * PI * y / L)) * PI * PI * f_v1 * V * NU_SA * pow(L, -0.2e1) + (U * U + V * V) * a_rhot * PI * rho_t * cos(a_rhot * PI * t / L) / L / 0.2e1 + (f_v1 + pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1))) * (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * PI * PI * mu * U * pow(L, -0.2e1) / 0.3e1 + (f_v1 + pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1))) * (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * PI * PI * mu * V * pow(L, -0.2e1) / 0.3e1 - (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * PI * RHO * U * V / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * cp * PI * P / L / R - (mu_t / Pr_t + mu / Pr) * (-(a_px * a_px * p_x * cos(a_px * PI * x / L) + a_py * a_py * p_y * sin(a_py * PI * y / L)) * cp * PI * PI * pow(L, -0.2e1) / R / RHO + (a_rhox * a_rhox * rho_x * sin(a_rhox * PI * x / L) + a_rhoy * a_rhoy * rho_y * cos(a_rhoy * PI * y / L)) * cp * PI * PI * P * pow(L, -0.2e1) / R * pow(RHO, -0.2e1)) + (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * PI * PI * mu_t * U * pow(L, -0.2e1) / 0.3e1 + (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L)) * PI * PI * mu_t * V * pow(L, -0.2e1) / 0.3e1 + (0.3e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U * U / L / 0.2e1 + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V * V / L / 0.2e1 - (f_v1 + pow(c_v1, 0.3e1) / (pow(chi, 0.3e1) + pow(c_v1, 0.3e1))) * (0.4e1 * a_ux * a_ux * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) - 0.4e1 * a_ux * a_vy * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) + 0.3e1 * a_uy * a_uy * u_y * u_y * pow(sin(a_uy * PI * y / L), 0.2e1) + 0.6e1 * a_uy * a_vx * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) + 0.3e1 * a_vx * a_vx * v_x * v_x * pow(sin(a_vx * PI * x / L), 0.2e1) + 0.4e1 * a_vy * a_vy * v_y * v_y * pow(cos(a_vy * PI * y / L), 0.2e1)) * PI * PI * mu * pow(L, -0.2e1) / 0.3e1 + cp * a_py * PI * p_y * V * cos(a_py * PI * y / L) / L / R - cp * a_px * PI * p_x * U * sin(a_px * PI * x / L) / L / R - (0.4e1 * a_ux * a_ux * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) - 0.4e1 * a_ux * a_vy * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) + 0.3e1 * a_uy * a_uy * u_y * u_y * pow(sin(a_uy * PI * y / L), 0.2e1) + 0.6e1 * a_uy * a_vx * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) + 0.3e1 * a_vx * a_vx * v_x * v_x * pow(sin(a_vx * PI * x / L), 0.2e1) + 0.4e1 * a_vy * a_vy * v_y * v_y * pow(cos(a_vy * PI * y / L), 0.2e1)) * PI * PI * mu_t * pow(L, -0.2e1) / 0.3e1 + (U * U + V * V) * a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L / 0.2e1 + a_vt * PI * v_t * RHO * V * cos(a_vt * PI * t / L) / L - a_ut * PI * u_t * RHO * U * sin(a_ut * PI * t / L) / L + cv * a_rhot * PI * rho_t * P * cos(a_rhot * PI * t / L) / L / R / RHO - (a_px * a_nusax * p_x * nu_sa_x * sin(a_px * PI * x / L) * sin(a_nusax * PI * x / L) - a_py * a_nusay * p_y * nu_sa_y * cos(a_py * PI * y / L) * sin(a_nusay * PI * y / L)) * cp * PI * PI * mu_t * pow(L, -0.2e1) / Pr_t / R / RHO / NU_SA - (a_rhox * a_nusax * rho_x * nu_sa_x * cos(a_rhox * PI * x / L) * sin(a_nusax * PI * x / L) - a_rhoy * a_nusay * rho_y * nu_sa_y * sin(a_rhoy * PI * y / L) * sin(a_nusay * PI * y / L)) * cp * PI * PI * mu_t * P * pow(L, -0.2e1) / Pr_t / R * pow(RHO, -0.2e1) / NU_SA - (a_rhox * a_px * rho_x * p_x * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) + a_rhoy * a_py * rho_y * p_y * sin(a_rhoy * PI * y / L) * cos(a_py * PI * y / L)) * (Pr * mu_t + 0.2e1 * Pr_t * mu) * cp * PI * PI / Pr / Pr_t * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) - (a_rhox * a_rhox * rho_x * rho_x * pow(cos(a_rhox * PI * x / L), 0.2e1) + a_rhoy * a_rhoy * rho_y * rho_y * pow(sin(a_rhoy * PI * y / L), 0.2e1)) * (Pr * mu_t + 0.2e1 * Pr_t * mu) * cp * PI * PI * P / Pr / Pr_t * pow(L, -0.2e1) / R * pow(RHO, -0.3e1) - (U * U + V * V) * a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L / 0.2e1;
  return(Q_E);
}
