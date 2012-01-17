// Source term for the momentum equation - u. 1D transient Navier-Stokes in thermochemical non-equilibrium. 
double SourceQ_u (
  double x,
  double t,
  double R,
  double M_N,
  double A_N,
  double B_N,
  double C_N,
  double A_N2,
  double B_N2,
  double C_N2)
{
  double Q_u;
  double RHO;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double Mu_N;
  double Mu_N2;
  double Phi_N;
  double Phi_N2;
  double Mtot;
  double Mu_mix;
  double DMu_mix_Dx;
  double Q_u_time;
  double Q_u_convection;
  double Q_u_gradp;
  double Q_u_work;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_t * cos(a_rho_N_t * PI * t / Lt);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_t * sin(a_rho_N2_t * PI * t / Lt);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / Lt);
  T = T_0 + T_x * cos(a_Tx * PI * x / L) + T_t * cos(a_Tt * PI * t / Lt);
    // "Aux relationships-----------------------------------------------"
;
  Mu_N = exp((A_N * log(T) + B_N) * log(T) + C_N) / 0.10e2;
  Mu_N2 = exp((A_N2 * log(T) + B_N2) * log(T) + C_N2) / 0.10e2;
  Mtot = 0.1e1 / (RHO_N / RHO / M_N + RHO_N2 / RHO / M_N / 0.2e1);
  Phi_N = RHO_N * Mtot / RHO / M_N + RHO_N2 * Mtot * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * sqrt(0.3e1) / RHO / M_N / 0.12e2;
  Phi_N2 = RHO_N * Mtot * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * sqrt(0.6e1) / RHO / M_N / 0.12e2 + RHO_N2 * Mtot / RHO / M_N / 0.2e1;
  Mu_mix = Mu_N * RHO_N * Mtot / RHO / M_N / Phi_N + Mu_N2 * RHO_N2 * Mtot / RHO / M_N / Phi_N2 / 0.2e1;
    // "Aux derivatives-----------------------------------------------"
;
  DMu_mix_Dx = -(0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * Mu_N * Mtot * RHO_N * sin(a_Tx * PI * x / L) / L / M_N / Phi_N / RHO / T - sqrt(0.6e1) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * a_rho_N_x * PI * rho_N_x * Mu_N2 * Mtot * Mtot * RHO_N2 * cos(a_rho_N_x * PI * x / L) / L * pow(M_N, -0.2e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) / 0.24e2 + sqrt(0.3e1) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * a_rho_N2_x * PI * rho_N2_x * Mu_N * Mtot * Mtot * RHO_N * sin(a_rho_N2_x * PI * x / L) / L * pow(M_N, -0.2e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) / 0.12e2 + sqrt(0.3e1) * pow(0.2e1, 0.1e1 / 0.4e1) * (0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1)) * (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * Mu_N * Mu_N * Mtot * Mtot * RHO_N * RHO_N2 * sin(a_Tx * PI * x / L) / sqrt(Mu_N / Mu_N2) / L * pow(M_N, -0.2e1) / Mu_N2 * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) / T / 0.12e2 - sqrt(0.6e1) * pow(0.2e1, 0.3e1 / 0.4e1) * (0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1) * (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * Mu_N2 * Mu_N2 * Mtot * Mtot * RHO_N * RHO_N2 * sin(a_Tx * PI * x / L) / sqrt(Mu_N2 / Mu_N) / L * pow(M_N, -0.2e1) / Mu_N * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) / T / 0.48e2 + (RHO * M_N * Phi_N - RHO_N * Mtot) * a_rho_N_x * PI * rho_N_x * Mu_N * Mtot * cos(a_rho_N_x * PI * x / L) / L * pow(M_N, -0.2e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) - (0.2e1 * RHO * M_N * Phi_N2 - RHO_N2 * Mtot) * a_rho_N2_x * PI * rho_N2_x * Mu_N2 * Mtot * sin(a_rho_N2_x * PI * x / L) / L * pow(M_N, -0.2e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) / 0.4e1 - (0.2e1 * Mu_N * Phi_N2 * RHO_N + Mu_N2 * Phi_N * RHO_N2) * (0.2e1 * A_N2 * log(T) + B_N2) * a_Tx * PI * T_x * Mtot * sin(a_Tx * PI * x / L) / L / M_N / Phi_N / Phi_N2 / RHO / T / 0.2e1 + sqrt(0.3e1) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * PI * Mu_N * pow(Mtot, 0.3e1) * RHO_N * RHO_N2 / L * pow(M_N, -0.3e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.3e1) / 0.24e2 + sqrt(0.6e1) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * PI * Mu_N2 * pow(Mtot, 0.3e1) * RHO_N * RHO_N2 / L * pow(M_N, -0.3e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.3e1) / 0.48e2 - (0.2e1 * Mu_N * Phi_N2 * RHO_N + Mu_N2 * Phi_N * RHO_N2) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * Mtot * Mtot / L * pow(M_N, -0.2e1) / Phi_N / Phi_N2 * pow(RHO, -0.2e1) / 0.4e1 + sqrt(0.3e1) * (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * PI * Mu_N * Mtot * Mtot * RHO_N * RHO_N2 / L * pow(M_N, -0.3e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.4e1) / 0.24e2 + sqrt(0.6e1) * (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * PI * Mu_N2 * Mtot * Mtot * RHO_N * RHO_N2 / L * pow(M_N, -0.3e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.4e1) / 0.48e2 + (0.4e1 * Mu_N * Phi_N2 * Phi_N2 * RHO_N * RHO_N + Mu_N2 * Phi_N * Phi_N * RHO_N2 * RHO_N2) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * pow(Mtot, 0.3e1) / L * pow(M_N, -0.3e1) * pow(Phi_N, -0.2e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.3e1) / 0.8e1 - (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (0.4e1 * M_N * Mu_N * Phi_N * Phi_N2 * Phi_N2 * RHO_N * RHO + 0.2e1 * M_N * Mu_N2 * Phi_N * Phi_N * Phi_N2 * RHO_N2 * RHO - 0.4e1 * Mu_N * Phi_N2 * Phi_N2 * Mtot * RHO_N * RHO_N - Mu_N2 * Phi_N * Phi_N * Mtot * RHO_N2 * RHO_N2) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * Mtot / L * pow(M_N, -0.3e1) * pow(Phi_N, -0.2e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.4e1) / 0.8e1;
    // "Contribution from the transient term to the total source term -----------------------------------------------"
;
  Q_u_time = -a_ut * PI * u_t * RHO * sin(a_ut * PI * t / Lt) / Lt - (a_rho_N_t * rho_N_t * sin(a_rho_N_t * PI * t / Lt) - a_rho_N2_t * rho_N2_t * cos(a_rho_N2_t * PI * t / Lt)) * PI * U / Lt;
    // "Contribution from the convective term to the total source term ---------------------------------------------"
;
  Q_u_convection = 0.2e1 * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L + (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * U * U / L;
    // "Contribution from the grad p term to the total source term --------------------------------------------"
;
  Q_u_gradp = -(0.2e1 * RHO_N + RHO_N2) * a_Tx * PI * T_x * R * sin(a_Tx * PI * x / L) / L / M_N / 0.2e1 + (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * T / L / M_N / 0.2e1;
    // "Contribution from the work done by viscous forces to the total source term -----------------------------------------------"
;
  Q_u_work = 0.4e1 / 0.3e1 * a_ux * a_ux * PI * PI * u_x * Mu_mix * sin(a_ux * PI * x / L) * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * DMu_mix_Dx * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L;
    // "Total source term -------------------------------------------------------------------------"
;
  Q_u = Q_u_time + Q_u_convection + Q_u_gradp + Q_u_work;
  return(Q_u);
}
