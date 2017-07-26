#include <math.h>

double SourceQ_v (
  double x,
  double y,
  double z,
  double t)
{
  double RHO;
  double P;
  double U;
  double V;
  double W;
  double T;
  double MU;
  double DMu_Dx;
  double DMu_Dy;
  double DMu_Dz;
  double Q_v;
  double Q_v_convection;
  double Q_v_pressure;
  double Q_v_viscous;
  double Q_v_time;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / Lt);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / Lt);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_t * sin(a_vt * PI * t / Lt);
  W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / Lt);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_t * cos(a_pt * PI * t / Lt);
  T = P / RHO / R;
  MU = A_mu * pow(T, 0.3e1 / 0.2e1) / (T + B_mu);
  DMu_Dx = a_rhox * PI * rho_x * MU * MU * cos(a_rhox * PI * x / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhox * PI * rho_x * MU * cos(a_rhox * PI * x / L) / L / RHO + a_px * PI * p_x * MU * MU * sin(a_px * PI * x / L) / A_mu / L / R / RHO * pow(T, -0.3e1 / 0.2e1) - 0.3e1 / 0.2e1 * a_px * PI * p_x * MU * sin(a_px * PI * x / L) / L / R / RHO / T;
  DMu_Dy = -a_rhoy * PI * rho_y * MU * MU * sin(a_rhoy * PI * y / L) / A_mu / L / RHO / sqrt(T) + 0.3e1 / 0.2e1 * a_rhoy * PI * rho_y * MU * sin(a_rhoy * PI * y / L) / L / RHO - a_py * PI * p_y * MU * MU * cos(a_py * PI * y / L) / A_mu / L / R / RHO * pow(T, -0.3e1 / 0.2e1) + 0.3e1 / 0.2e1 * a_py * PI * p_y * MU * cos(a_py * PI * y / L) / L / R / RHO / T;
  DMu_Dz = a_rhoz * PI * rho_z * MU * MU * cos(a_rhoz * PI * z / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhoz * PI * rho_z * MU * cos(a_rhoz * PI * z / L) / L / RHO + a_pz * p_z * PI * MU * MU * sin(a_pz * PI * z / L) / A_mu / L / R / RHO * pow(T, -0.3e1 / 0.2e1) - 0.3e1 / 0.2e1 * a_pz * p_z * PI * MU * sin(a_pz * PI * z / L) / L / R / RHO / T;
  Q_v_convection = a_rhox * PI * rho_x * U * V * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * V * W * cos(a_rhoz * PI * z / L) / L - a_vx * PI * v_x * RHO * U * sin(a_vx * PI * x / L) / L + a_vz * PI * v_z * RHO * W * cos(a_vz * PI * z / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * V / L;
  Q_v_pressure = a_py * PI * p_y * cos(a_py * PI * y / L) / L;
  Q_v_viscous = (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L) + 0.3e1 * a_vz * a_vz * v_z * sin(a_vz * PI * z / L)) * PI * PI * MU * pow(L, -0.2e1) / 0.3e1 + (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * DMu_Dx * PI / L + 0.2e1 / 0.3e1 * (a_ux * u_x * cos(a_ux * PI * x / L) - 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dy * PI / L - (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * DMu_Dz * PI / L;
  Q_v_time = a_vt * v_t * PI * RHO * cos(a_vt * PI * t / Lt) / Lt + a_rhot * rho_t * PI * V * cos(a_rhot * PI * t / Lt) / Lt;
  Q_v = Q_v_convection + Q_v_pressure + Q_v_viscous + Q_v_time;
  return(Q_v);
}
