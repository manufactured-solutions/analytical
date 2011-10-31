#include <math.h>

double SourceQ_v (double x, double y)
{
  double Q_v;
  double RHO;
  double U;
  double V;
  double P;
  double M1;
  double M2;
  double MU;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L);
  M1 = A_mu * pow(P / R / RHO, 0.3e1 / 0.2e1);
  M2 = P / R / RHO + B_mu;
  MU = M1 / M2;
  Q_v = a_rhox * PI * rho_x * U * V * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * sin(a_rhoy * PI * y / L) / L - v_x * a_vx * PI * RHO * U * sin(a_vx * PI * x / L) / L + a_py * PI * p_y * cos(a_py * PI * y / L) / L + (0.3e1 * v_x * a_vx * a_vx * cos(a_vx * PI * x / L) + 0.4e1 * v_y * a_vy * a_vy * sin(a_vy * PI * y / L)) * PI * PI * MU * pow(L, -0.2e1) / 0.3e1 + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * v_y * a_vy * cos(a_vy * PI * y / L)) * PI * RHO * V / L + (-0.3e1 * a_rhox * a_uy * rho_x * u_y * cos(a_rhox * PI * x / L) * sin(a_uy * PI * y / L) + 0.2e1 * a_rhoy * a_ux * rho_y * u_x * sin(a_rhoy * PI * y / L) * cos(a_ux * PI * x / L) - 0.3e1 * v_x * a_rhox * a_vx * rho_x * cos(a_rhox * PI * x / L) * sin(a_vx * PI * x / L) - 0.4e1 * v_y * a_rhoy * a_vy * rho_y * sin(a_rhoy * PI * y / L) * cos(a_vy * PI * y / L)) * (0.3e1 * B_mu * R * RHO + P) * PI * PI * MU / (B_mu * R * RHO + P) * pow(L, -0.2e1) / RHO / 0.6e1 - (0.3e1 * a_px * a_uy * p_x * u_y * sin(a_px * PI * x / L) * sin(a_uy * PI * y / L) - 0.2e1 * a_py * a_ux * p_y * u_x * cos(a_py * PI * y / L) * cos(a_ux * PI * x / L) + 0.3e1 * v_x * a_px * a_vx * p_x * sin(a_px * PI * x / L) * sin(a_vx * PI * x / L) + 0.4e1 * v_y * a_py * a_vy * p_y * cos(a_py * PI * y / L) * cos(a_vy * PI * y / L)) * (0.3e1 * B_mu * R * RHO + P) * PI * PI * MU / (B_mu * R * RHO + P) * pow(L, -0.2e1) / P / 0.6e1;
  return(Q_v);
}
