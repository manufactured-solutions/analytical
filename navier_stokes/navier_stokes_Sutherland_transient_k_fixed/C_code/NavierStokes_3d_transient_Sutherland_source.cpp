
double eval_q_e(double x, double y, double z, double t)
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
 double kappa;
 double Q_e;
 double Q_e_convection;
 double Q_e_work_pressure;
 double Q_e_work_viscous;
 double Q_e_conduction;
 double Q_e_time;
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
 kappa = Gamma * R * MU / (Gamma - 0.1e1) / Pr;
 Q_e_convection = a_rhox * PI * rho_x * pow(U, 0.3e1) * cos(a_rhox * PI * x / L) / L / 0.2e1 + a_rhox * PI * rho_x * U * V * V * cos(a_rhox * PI * x / L) / L / 0.2e1 + a_rhox * PI * rho_x * U * W * W * cos(a_rhox * PI * x / L) / L / 0.2e1 - a_rhoy * PI * rho_y * U * U * V * sin(a_rhoy * PI * y / L) / L / 0.2e1 - a_rhoy * PI * rho_y * pow(V, 0.3e1) * sin(a_rhoy * PI * y / L) / L / 0.2e1 - a_rhoy * PI * rho_y * V * W * W * sin(a_rhoy * PI * y / L) / L / 0.2e1 + a_rhoz * PI * rho_z * U * U * W * cos(a_rhoz * PI * z / L) / L / 0.2e1 + a_rhoz * PI * rho_z * V * V * W * cos(a_rhoz * PI * z / L) / L / 0.2e1 + a_rhoz * PI * rho_z * pow(W, 0.3e1) * cos(a_rhoz * PI * z / L) / L / 0.2e1 - a_px * PI * p_x * U * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + a_py * PI * p_y * V * cos(a_py * PI * y / L) / (Gamma - 0.1e1) / L - a_pz * p_z * PI * W * sin(a_pz * PI * z / L) / (Gamma - 0.1e1) / L + (0.3e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * PI * RHO * U * V / L + (-a_uz * u_z * sin(a_uz * PI * z / L) + a_wx * w_x * cos(a_wx * PI * x / L)) * PI * RHO * U * W / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * V * V / L / 0.2e1 + (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * PI * RHO * V * W / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - 0.3e1 * a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * W * W / L / 0.2e1 + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * P / (Gamma - 0.1e1) / L;
 Q_e_work_pressure = -a_px * PI * p_x * U * sin(a_px * PI * x / L) / L + a_py * PI * p_y * V * cos(a_py * PI * y / L) / L - a_pz * p_z * PI * W * sin(a_pz * PI * z / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * P / L;
 Q_e_conduction = DMu_Dx * a_rhox * PI * rho_x * kappa * P * cos(a_rhox * PI * x / L) / L / R / MU * pow(RHO, -0.2e1) + DMu_Dx * a_px * PI * p_x * kappa * sin(a_px * PI * x / L) / L / R / MU / RHO - DMu_Dy * a_rhoy * PI * rho_y * kappa * P * sin(a_rhoy * PI * y / L) / L / R / MU * pow(RHO, -0.2e1) - DMu_Dy * a_py * PI * p_y * kappa * cos(a_py * PI * y / L) / L / R / MU / RHO + DMu_Dz * a_rhoz * PI * rho_z * kappa * P * cos(a_rhoz * PI * z / L) / L / R / MU * pow(RHO, -0.2e1) + a_pz * DMu_Dz * p_z * PI * kappa * sin(a_pz * PI * z / L) / L / R / MU / RHO + (a_px * a_px * p_x * cos(a_px * PI * x / L) + a_py * a_py * p_y * sin(a_py * PI * y / L) + a_pz * a_pz * p_z * cos(a_pz * PI * z / L)) * PI * PI * kappa * pow(L, -0.2e1) / R / RHO - (a_rhox * a_rhox * rho_x * sin(a_rhox * PI * x / L) + a_rhoy * a_rhoy * rho_y * cos(a_rhoy * PI * y / L) + a_rhoz * a_rhoz * rho_z * sin(a_rhoz * PI * z / L)) * PI * PI * kappa * P * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) - (0.2e1 * a_rhox * a_px * rho_x * p_x * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) + 0.2e1 * a_rhoy * a_py * rho_y * p_y * sin(a_rhoy * PI * y / L) * cos(a_py * PI * y / L) + 0.2e1 * a_pz * a_rhoz * p_z * rho_z * cos(a_rhoz * PI * z / L) * sin(a_pz * PI * z / L)) * PI * PI * kappa * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) - (0.2e1 * a_rhox * a_rhox * rho_x * rho_x * pow(cos(a_rhox * PI * x / L), 0.2e1) + 0.2e1 * a_rhoy * a_rhoy * rho_y * rho_y * pow(sin(a_rhoy * PI * y / L), 0.2e1) + 0.2e1 * a_rhoz * a_rhoz * rho_z * rho_z * pow(cos(a_rhoz * PI * z / L), 0.2e1)) * PI * PI * kappa * P * pow(L, -0.2e1) / R * pow(RHO, -0.3e1);
 Q_e_work_viscous = (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L) + 0.3e1 * a_uz * a_uz * u_z * cos(a_uz * PI * z / L)) * PI * PI * MU * U * pow(L, -0.2e1) / 0.3e1 + (0.3e1 * a_vx * a_vx * v_x * cos(a_vx * PI * x / L) + 0.4e1 * a_vy * a_vy * v_y * sin(a_vy * PI * y / L) + 0.3e1 * a_vz * a_vz * v_z * sin(a_vz * PI * z / L)) * PI * PI * MU * V * pow(L, -0.2e1) / 0.3e1 + (0.3e1 * a_wx * a_wx * w_x * sin(a_wx * PI * x / L) + 0.3e1 * a_wy * a_wy * w_y * sin(a_wy * PI * y / L) + 0.4e1 * a_wz * a_wz * w_z * cos(a_wz * PI * z / L)) * PI * PI * MU * W * pow(L, -0.2e1) / 0.3e1 - 0.2e1 / 0.3e1 * (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) - a_vy * v_y * cos(a_vy * PI * y / L) + a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dx * PI * U / L + (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * DMu_Dx * PI * V / L - (-a_uz * u_z * sin(a_uz * PI * z / L) + a_wx * w_x * cos(a_wx * PI * x / L)) * DMu_Dx * PI * W / L + (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * DMu_Dy * PI * U / L + 0.2e1 / 0.3e1 * (a_ux * u_x * cos(a_ux * PI * x / L) - 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dy * PI * V / L - (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * DMu_Dy * PI * W / L - (-a_uz * u_z * sin(a_uz * PI * z / L) + a_wx * w_x * cos(a_wx * PI * x / L)) * DMu_Dz * PI * U / L - (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * DMu_Dz * PI * V / L + 0.2e1 / 0.3e1 * (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) + 0.2e1 * a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dz * PI * W / L - (0.4e1 * a_ux * a_ux * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) - 0.4e1 * a_ux * a_vy * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) + 0.4e1 * a_ux * a_wz * u_x * w_z * cos(a_ux * PI * x / L) * sin(a_wz * PI * z / L) + 0.3e1 * a_uy * a_uy * u_y * u_y * pow(sin(a_uy * PI * y / L), 0.2e1) + 0.6e1 * a_uy * a_vx * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) + 0.3e1 * a_uz * a_uz * u_z * u_z * pow(sin(a_uz * PI * z / L), 0.2e1) - 0.6e1 * a_uz * a_wx * u_z * w_x * sin(a_uz * PI * z / L) * cos(a_wx * PI * x / L) + 0.3e1 * a_vx * a_vx * v_x * v_x * pow(sin(a_vx * PI * x / L), 0.2e1) + 0.4e1 * a_vy * a_vy * v_y * v_y * pow(cos(a_vy * PI * y / L), 0.2e1) + 0.4e1 * a_vy * a_wz * v_y * w_z * cos(a_vy * PI * y / L) * sin(a_wz * PI * z / L) + 0.3e1 * a_vz * a_vz * v_z * v_z * pow(cos(a_vz * PI * z / L), 0.2e1) + 0.6e1 * a_vz * a_wy * v_z * w_y * cos(a_vz * PI * z / L) * cos(a_wy * PI * y / L) + 0.3e1 * a_wx * a_wx * w_x * w_x * pow(cos(a_wx * PI * x / L), 0.2e1) + 0.3e1 * a_wy * a_wy * w_y * w_y * pow(cos(a_wy * PI * y / L), 0.2e1) + 0.4e1 * a_wz * a_wz * w_z * w_z * pow(sin(a_wz * PI * z / L), 0.2e1)) * PI * PI * MU * pow(L, -0.2e1) / 0.3e1;
 Q_e_time = -a_ut * u_t * PI * RHO * U * sin(a_ut * PI * t / Lt) / Lt + a_vt * v_t * PI * RHO * V * cos(a_vt * PI * t / Lt) / Lt - a_wt * w_t * PI * RHO * W * sin(a_wt * PI * t / Lt) / Lt + a_rhot * rho_t * PI * U * U * cos(a_rhot * PI * t / Lt) / Lt / 0.2e1 + a_rhot * rho_t * PI * V * V * cos(a_rhot * PI * t / Lt) / Lt / 0.2e1 + a_rhot * rho_t * PI * W * W * cos(a_rhot * PI * t / Lt) / Lt / 0.2e1 - a_pt * p_t * PI * sin(a_pt * PI * t / Lt) / (Gamma - 0.1e1) / Lt;
 Q_e = Q_e_convection + Q_e_work_pressure + Q_e_work_viscous + Q_e_conduction + Q_e_time;
 return(Q_e);
}

double eval_q_rho(double x, double y, double z, double t)
{
 double RHO;
 double U;
 double V;
 double W;
 double Q_rho;
 double Q_rho_convection;
 double Q_rho_time;
 RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_t * sin(a_rhot * PI * t / Lt);
 U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L) + u_t * cos(a_ut * PI * t / Lt);
 V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_t * sin(a_vt * PI * t / Lt);
 W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L) + w_t * cos(a_wt * PI * t / Lt);
 Q_rho_convection = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO / L;
 Q_rho_time = a_rhot * PI * rho_t * cos(a_rhot * PI * t / Lt) / Lt;
 Q_rho = Q_rho_convection + Q_rho_time;
 return(Q_rho);
}

double eval_q_u(double x, double y, double z, double t)
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
 double Q_u;
 double Q_u_convection;
 double Q_u_pressure;
 double Q_u_viscous;
 double Q_u_time;
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
 Q_u_convection = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * U * W * cos(a_rhoz * PI * z / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L - a_uz * PI * u_z * RHO * W * sin(a_uz * PI * z / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * U / L;
 Q_u_pressure = -a_px * PI * p_x * sin(a_px * PI * x / L) / L;
 Q_u_viscous = (0.4e1 * a_ux * a_ux * u_x * sin(a_ux * PI * x / L) + 0.3e1 * a_uy * a_uy * u_y * cos(a_uy * PI * y / L) + 0.3e1 * a_uz * a_uz * u_z * cos(a_uz * PI * z / L)) * PI * PI * MU * pow(L, -0.2e1) / 0.3e1 - 0.2e1 / 0.3e1 * (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) - a_vy * v_y * cos(a_vy * PI * y / L) + a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dx * PI / L + (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * DMu_Dy * PI / L + (a_uz * u_z * sin(a_uz * PI * z / L) - a_wx * w_x * cos(a_wx * PI * x / L)) * DMu_Dz * PI / L;
 Q_u_time = -a_ut * u_t * PI * RHO * sin(a_ut * PI * t / Lt) / Lt + a_rhot * rho_t * PI * U * cos(a_rhot * PI * t / Lt) / Lt;
 Q_u = Q_u_convection + Q_u_pressure + Q_u_viscous + Q_u_time;
 return(Q_u);
}
double eval_q_v(double x, double y, double z, double t)
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

double eval_q_w(double x, double y, double z, double t)
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
 double Q_w;
 double Q_w_convection;
 double Q_w_pressure;
 double Q_w_viscous;
 double Q_w_time;
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
 Q_w_convection = a_rhox * PI * rho_x * U * W * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * W * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * W * W * cos(a_rhoz * PI * z / L) / L + a_wx * PI * w_x * RHO * U * cos(a_wx * PI * x / L) / L + a_wy * PI * w_y * RHO * V * cos(a_wy * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - 0.2e1 * a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * W / L;
 Q_w_pressure = -a_pz * p_z * PI * sin(a_pz * PI * z / L) / L;
 Q_w_viscous = (0.3e1 * a_wx * a_wx * w_x * sin(a_wx * PI * x / L) + 0.3e1 * a_wy * a_wy * w_y * sin(a_wy * PI * y / L) + 0.4e1 * a_wz * a_wz * w_z * cos(a_wz * PI * z / L)) * PI * PI * MU * pow(L, -0.2e1) / 0.3e1 - (-a_uz * u_z * sin(a_uz * PI * z / L) + a_wx * w_x * cos(a_wx * PI * x / L)) * DMu_Dx * PI / L - (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * DMu_Dy * PI / L + 0.2e1 / 0.3e1 * (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) + 0.2e1 * a_wz * w_z * sin(a_wz * PI * z / L)) * DMu_Dz * PI / L;
 Q_w_time = -a_wt * w_t * PI * RHO * sin(a_wt * PI * t / Lt) / Lt + a_rhot * rho_t * PI * W * cos(a_rhot * PI * t / Lt) / Lt;
 Q_w = Q_w_convection + Q_w_pressure + Q_w_viscous + Q_w_time;
 return(Q_w);
}
