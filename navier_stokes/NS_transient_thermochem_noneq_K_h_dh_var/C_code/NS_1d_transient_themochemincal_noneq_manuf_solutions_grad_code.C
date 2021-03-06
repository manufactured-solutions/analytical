rho_an_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_t * cos(a_rho_N_t * PI * t / Lt);
rho_an_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_t * sin(a_rho_N2_t * PI * t / Lt);
rho_an = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_t * cos(a_rho_N_t * PI * t / Lt) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_t * sin(a_rho_N2_t * PI * t / Lt);
u_an = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / Lt);
T_an = T_0 + T_x * cos(a_Tx * PI * x / L) + T_t * cos(a_Tt * PI * t / Lt);
Tv_an = Tv_0 + Tv_x * cos(a_Tvx * PI * x / L) + Tv_t * sin(a_Tvt * PI * t / Lt);
grad_rho_an_N[0] = rho_N_x * cos(a_rho_N_x * PI * x / L) * a_rho_N_x * PI / L;
grad_rho_an_N[1] = 0;
grad_rho_an_N[2] = 0;
grad_rho_an_N2[0] = -rho_N2_x * sin(a_rho_N2_x * PI * x / L) * a_rho_N2_x * PI / L;
grad_rho_an_N2[1] = 0;
grad_rho_an_N2[2] = 0;
grad_rho_an[0] = rho_N_x * cos(a_rho_N_x * PI * x / L) * a_rho_N_x * PI / L - rho_N2_x * sin(a_rho_N2_x * PI * x / L) * a_rho_N2_x * PI / L;
grad_rho_an[1] = 0;
grad_rho_an[2] = 0;
grad_u_an[0] = u_x * cos(a_ux * PI * x / L) * a_ux * PI / L;
grad_u_an[1] = 0;
grad_u_an[2] = 0;
grad_T_an[0] = -T_x * sin(a_Tx * PI * x / L) * a_Tx * PI / L;
grad_T_an[1] = 0;
grad_T_an[2] = 0;
grad_Tv_an[0] = -Tv_x * sin(a_Tvx * PI * x / L) * a_Tvx * PI / L;
grad_Tv_an[1] = 0;
grad_Tv_an[2] = 0;
