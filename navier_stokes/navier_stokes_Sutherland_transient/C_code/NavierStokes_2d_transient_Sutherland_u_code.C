#include <math.h>

double SourceQ_u (double x, double y, double t)
{
  double Q_u;
  double RHO;
  double U;
  double V;
  double P;
  double M1;
  double M2;
  double MU;
  RHO = r * rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_t * sin(a_rhot * PI * t / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / L);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_t * cos(a_pt * PI * t / L);
  M1 = A_mu * pow(P / R / RHO, 0.3e1 / 0.2e1);
  M2 = P / R / RHO + B_mu;
  MU = M1 / M2;
  Q_u = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L + a_rhot * PI * rho_t * U * cos(a_rhot * PI * t / L) / L - a_ut * PI * u_t * RHO * sin(a_ut * PI * t / L) / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L + (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L)) * PI * PI * MU * pow(L, -0.2e1) / 0.3e1 + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U / L + (0.4e1 * a_rhox * a_ux * rho_x * u_x * cos(a_rhox * PI * x / L) * cos(a_ux * PI * x / L) - 0.2e1 * a_rhox * a_vy * rho_x * v_y * cos(a_rhox * PI * x / L) * cos(a_vy * PI * y / L) + 0.3e1 * a_rhoy * a_uy * rho_y * u_y * sin(a_rhoy * PI * y / L) * sin(a_uy * PI * y / L) + 0.3e1 * a_rhoy * a_vx * rho_y * v_x * sin(a_rhoy * PI * y / L) * sin(a_vx * PI * x / L)) * (0.3e1 * B_mu * R * RHO + P) * PI * PI * MU / (B_mu * R * RHO + P) * pow(L, -0.2e1) / RHO / 0.6e1 + (0.4e1 * a_px * a_ux * p_x * u_x * sin(a_px * PI * x / L) * cos(a_ux * PI * x / L) - 0.2e1 * a_px * a_vy * p_x * v_y * sin(a_px * PI * x / L) * cos(a_vy * PI * y / L) + 0.3e1 * a_py * a_uy * p_y * u_y * cos(a_py * PI * y / L) * sin(a_uy * PI * y / L) + 0.3e1 * a_py * a_vx * p_y * v_x * cos(a_py * PI * y / L) * sin(a_vx * PI * x / L)) * (0.3e1 * B_mu * R * RHO + P) * PI * PI * MU / (B_mu * R * RHO + P) * pow(L, -0.2e1) / P / 0.6e1;
  return(Q_u);
}
