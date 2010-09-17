#include <math.h>

double SourceQ_e_CtrlC (
  double x,
  double y,
  double z,
  double u_0,
  double u_x,
  double u_y,
  double u_z,
  double v_0,
  double v_x,
  double v_y,
  double v_z,
  double w_0,
  double w_x,
  double w_y,
  double w_z,
  double rho_0,
  double rho_x,
  double rho_y,
  double rho_z,
  double p_0,
  double p_x,
  double p_y,
  double p_z,
  double a_px,
  double a_py,
  double a_pz,
  double a_rhox,
  double a_rhoy,
  double a_rhoz,
  double a_ux,
  double a_uy,
  double a_uz,
  double a_vx,
  double a_vy,
  double a_vz,
  double a_wx,
  double a_wy,
  double a_wz,
  double mu,
  double Gamma,
  double L)
{
  double Qe_divided;
  double RHO;
  double P;
  double U;
  double V;
  double W;
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L);
  RHO = rho_0 + rho_x * sin(a_rhox * x * PI / L) + rho_y * cos(a_rhoy * y * PI / L) + rho_z * sin(a_rhoz * z * PI / L);
  U = u_0 + u_x * sin(a_ux * x * PI / L) + u_y * cos(a_uy * y * PI / L) + u_z * cos(a_uz * z * PI / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L);
  W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L);
  Qe_divided = cos(a_rhox * x * PI / L) * (U * U + V * V + W * W) * a_rhox * PI * rho_x * U / L / 0.2e1 + cos(a_rhoz * z * PI / L) * (U * U + V * V + W * W) * a_rhoz * PI * rho_z * W / L / 0.2e1 - sin(a_rhoy * y * PI / L) * (U * U + V * V + W * W) * a_rhoy * PI * rho_y * V / L / 0.2e1 + ((U * U + 0.3e1 * V * V + W * W) * RHO / L / 0.2e1 + Gamma * P / L / (Gamma - 0.1e1)) * cos(a_vy * PI * y / L) * a_vy * PI * v_y + (-(U * U + V * V + 0.3e1 * W * W) * RHO / L / 0.2e1 - Gamma * P / L / (Gamma - 0.1e1)) * sin(a_wz * PI * z / L) * a_wz * PI * w_z + 0.2e1 * mu * cos(a_wx * PI * x / L) * sin(a_uz * z * PI / L) * a_uz * a_wx * PI * PI * u_z * w_x * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * cos(a_vy * PI * y / L) * sin(a_wz * PI * z / L) * a_vy * a_wz * PI * PI * v_y * w_z * pow(L, -0.2e1) - 0.2e1 * mu * cos(a_vz * PI * z / L) * cos(a_wy * PI * y / L) * a_vz * a_wy * PI * PI * v_z * w_y * pow(L, -0.2e1) - (0.2e1 * pow(sin(a_rhoy * y * PI / L), 0.2e1) * rho_y + cos(a_rhoy * y * PI / L) * RHO) * k * a_rhoy * a_rhoy * PI * PI * rho_y * P * pow(L, -0.2e1) / R * pow(RHO, -0.3e1) - 0.2e1 * cos(a_rhox * x * PI / L) * sin(a_px * PI * x / L) * k * a_rhox * a_px * PI * PI * rho_x * p_x * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) - 0.2e1 * sin(a_rhoy * y * PI / L) * cos(a_py * PI * y / L) * k * a_rhoy * a_py * PI * PI * rho_y * p_y * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) - 0.2e1 * cos(a_rhoz * z * PI / L) * sin(a_pz * PI * z / L) * k * a_rhoz * a_pz * PI * PI * rho_z * p_z * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) - (0.2e1 * pow(cos(a_rhox * x * PI / L), 0.2e1) * rho_x + sin(a_rhox * x * PI / L) * RHO) * k * a_rhox * a_rhox * PI * PI * rho_x * P * pow(L, -0.2e1) / R * pow(RHO, -0.3e1) - Gamma * sin(a_pz * PI * z / L) * a_pz * PI * p_z * W / L / (Gamma - 0.1e1) + Gamma * cos(a_py * PI * y / L) * a_py * PI * p_y * V / L / (Gamma - 0.1e1) - Gamma * sin(a_px * PI * x / L) * a_px * PI * p_x * U / L / (Gamma - 0.1e1) + ((0.3e1 * U * U + V * V + W * W) * RHO / L / 0.2e1 + Gamma * P / L / (Gamma - 0.1e1)) * cos(a_ux * x * PI / L) * a_ux * PI * u_x - (0.2e1 * pow(cos(a_rhoz * z * PI / L), 0.2e1) * rho_z + sin(a_rhoz * z * PI / L) * RHO) * k * a_rhoz * a_rhoz * PI * PI * rho_z * P * pow(L, -0.2e1) / R * pow(RHO, -0.3e1) + 0.4e1 / 0.3e1 * mu * cos(a_ux * x * PI / L) * cos(a_vy * PI * y / L) * a_ux * a_vy * PI * PI * u_x * v_y * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * cos(a_ux * x * PI / L) * sin(a_wz * PI * z / L) * a_ux * a_wz * PI * PI * u_x * w_z * pow(L, -0.2e1) - 0.2e1 * mu * sin(a_uy * y * PI / L) * sin(a_vx * PI * x / L) * a_uy * a_vx * PI * PI * u_y * v_x * pow(L, -0.2e1) + cos(a_vz * PI * z / L) * a_vz * PI * v_z * RHO * V * W / L - sin(a_uz * z * PI / L) * a_uz * PI * u_z * RHO * U * W / L + sin(a_py * PI * y / L) * k * a_py * a_py * PI * PI * p_y * pow(L, -0.2e1) / R / RHO - sin(a_vx * PI * x / L) * a_vx * PI * v_x * RHO * U * V / L + cos(a_pz * PI * z / L) * k * a_pz * a_pz * PI * PI * p_z * pow(L, -0.2e1) / R / RHO + cos(a_px * PI * x / L) * k * a_px * a_px * PI * PI * p_x * pow(L, -0.2e1) / R / RHO + cos(a_wy * PI * y / L) * a_wy * PI * w_y * RHO * V * W / L - sin(a_uy * y * PI / L) * a_uy * PI * u_y * RHO * U * V / L + (-pow(sin(a_uy * y * PI / L), 0.2e1) * u_y + cos(a_uy * y * PI / L) * U) * mu * a_uy * a_uy * PI * PI * u_y * pow(L, -0.2e1) + (-pow(sin(a_uz * z * PI / L), 0.2e1) * u_z + cos(a_uz * z * PI / L) * U) * mu * a_uz * a_uz * PI * PI * u_z * pow(L, -0.2e1) - (pow(sin(a_vx * PI * x / L), 0.2e1) * v_x - cos(a_vx * PI * x / L) * V) * mu * a_vx * a_vx * PI * PI * v_x * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * (pow(cos(a_vy * PI * y / L), 0.2e1) * v_y - sin(a_vy * PI * y / L) * V) * mu * a_vy * a_vy * PI * PI * v_y * pow(L, -0.2e1) - (pow(cos(a_vz * PI * z / L), 0.2e1) * v_z - sin(a_vz * PI * z / L) * V) * mu * a_vz * a_vz * PI * PI * v_z * pow(L, -0.2e1) + (-pow(cos(a_wx * PI * x / L), 0.2e1) * w_x + sin(a_wx * PI * x / L) * W) * mu * a_wx * a_wx * PI * PI * w_x * pow(L, -0.2e1) + (-pow(cos(a_wy * PI * y / L), 0.2e1) * w_y + sin(a_wy * PI * y / L) * W) * mu * a_wy * a_wy * PI * PI * w_y * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * (-pow(sin(a_wz * PI * z / L), 0.2e1) * w_z + cos(a_wz * PI * z / L) * W) * mu * a_wz * a_wz * PI * PI * w_z * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * (-pow(cos(a_ux * x * PI / L), 0.2e1) * u_x + sin(a_ux * x * PI / L) * U) * mu * a_ux * a_ux * PI * PI * u_x * pow(L, -0.2e1) + cos(a_wx * PI * x / L) * a_wx * PI * w_x * RHO * U * W / L;
  return(Qe_divided);
}
