// Source term for the vibrational energy equation. 1D transient Navier-Stokes in thermochemical non-equilibrium. Equilibrium constant K must either be calculated previously or inside the procedure, as indicated bellow.
double SourceQ_eV (
  double x,
  double t,
  double q,
  double R,
  double M_N,
  double Cf1_N,
  double etaf1_N,
  double Ea_N,
  double Cf1_N2,
  double etaf1_N2,
  double Ea_N2,
  int energy_level_N,
  int energy_level_N2,
  double *heta_e_N,
  double *g_N,
  double *theta_e_N2,
  double *g_N2,
  double theta_v_N2,
  double A_N,
  double B_N,
  double C_N,
  double A_N2,
  double B_N2,
  double C_N2)
{
  double Q_eV;
  double RHO;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double TV;
  double kf1_N;
  double kf1_N2;
  double K;
  double R1;
  double w_dot_N;
  double w_dot_N2;
  double w_dot_V;
  double T_bar;
  double P;
  double tau_vib_N2_N2;
  double tau_vib_N2_N;
  double Mu_N;
  double Mu_N2;
  double Phi_N;
  double Phi_N2;
  double Mtot;
  double Cp;
  double Cv;
  double Kappa_mix;
  double Ds;
  double Kappa_ev_mix;
  double Kappa_tr_mix;
  double E_elec_N;
  double e_elec_N_num;
  double e_elec_N_den;
  double E_elec_N2;
  double e_elec_N2_num;
  double e_elec_N2_den;
  double E_vib_N2;
  double e_vib_eq_N2;
  int i;
  double eV;
  double eV_N;
  double eV_N2;
  double DKappa_mix_Dx;
  double DKappa_ev_Dx;
  double DKappa_tr_Dx;
  double DeV_Dx;
  double DeV_N_Dx;
  double DeV_N2_Dx;
  double DCp_Dx;
  double Sum_eN_thetae2_g_div_e;
  double Sum_eN_thetae3_g_div_e;
  double Sum_eN2_thetae2_g_div_e;
  double Sum_eN2_thetae3_g_div_e;
  double Q_eV_time;
  double Q_eV_convection;
  double Q_eV_production;
  double Q_eV_heatflux;
  double Q_eV_diffusion;
  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * PI * x / L) + rho_N_t * cos(a_rho_N_t * PI * t / Lt);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * PI * x / L) + rho_N2_t * sin(a_rho_N2_t * PI * t / Lt);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / Lt);
  T = T_0 + T_x * cos(a_Tx * PI * x / L) + T_t * cos(a_Tt * PI * t / Lt);
  TV = Tv_0 + Tv_x * cos(a_Tvx * PI * x / L) + Tv_t * sin(a_Tvt * PI * t / Lt);
  P = RHO_N * R * T / M_N + RHO_N2 * R * T / M_N / 0.2e1;
  T_bar = pow(T, q) * pow(TV, 0.1e1 - q);
    // " -----------------------------------------------"
;
  K = calculate_equilibrium_constant_K(T);
  kf1_N2 = Cf1_N2 * pow(T_bar, etaf1_N2) * exp(-Ea_N2 / R / T_bar);
  kf1_N = Cf1_N * pow(T_bar, etaf1_N) * exp(-Ea_N / R / T_bar);
  R1 = kf1_N * pow(RHO_N, 0.3e1) / K * pow(M_N, -0.3e1) - kf1_N * RHO_N2 * RHO_N * pow(M_N, -0.2e1) / 0.2e1 + kf1_N2 * RHO_N * RHO_N * RHO_N2 / K * pow(M_N, -0.3e1) / 0.2e1 - kf1_N2 * RHO_N2 * RHO_N2 * pow(M_N, -0.2e1) / 0.4e1;
  w_dot_N = -0.2e1 * M_N * R1;
  w_dot_N2 = 0.2e1 * M_N * R1;
  tau_vib_N2_N = exp(0.29e2 / 0.75000e5 * sqrt(0.6e1) * sqrt(M_N) * pow(theta_v_N2, 0.4e1 / 0.3e1) * (pow(T, -0.1e1 / 0.3e1) - pow(0.54e2, 0.1e1 / 0.4e1) * pow(M_N, 0.1e1 / 0.4e1) / 0.200e3) - 0.921e3 / 0.50e2) / P;
  tau_vib_N2_N2 = exp(0.29e2 / 0.25000e5 * sqrt(M_N) * pow(theta_v_N2, 0.4e1 / 0.3e1) * (pow(T, -0.1e1 / 0.3e1) - 0.3e1 / 0.200e3 * pow(M_N, 0.1e1 / 0.4e1)) - 0.921e3 / 0.50e2) / P;
  w_dot_V = E_elec_N * w_dot_N + (E_vib_N2 + E_elec_N2) * w_dot_N2 + (E_vib_eq_N2 - E_vib_N2) * RHO_N2 * (0.2e1 * RHO_N / (0.2e1 * RHO_N + RHO_N2) / tau_vib_N2_N + RHO_N2 / (0.2e1 * RHO_N + RHO_N2) / tau_vib_N2_N2);
    // "Vibrational, Electronic and Aux summations -----------------------------------------------"
;
  e_elec_N_num = 0.0e0;
  e_elec_N2_num = 0.0e0;
  e_elec_N_den = 0.0e0;
  e_elec_N2_den = 0.0e0;
  Sum_eN_thetae2_g_div_e = 0.0e0;
  Sum_eN_thetae3_g_div_e = 0.0e0;
  Sum_eN2_thetae2_g_div_e = 0.0e0;
  Sum_eN2_thetae3_g_div_e = 0.0e0;
  for (i = 1; i <= energy_level_N; i++)
  {
    e_elec_N_den = e_elec_N_den + g_N[i - 1] * exp(-theta_e_N[i - 1] / TV);
    e_elec_N_num = e_elec_N_num + theta_e_N[i - 1] * g_N[i - 1] * exp(-theta_e_N[i - 1] / TV);
    Sum_eN_thetae2_g_div_e = Sum_eN_thetae2_g_div_e + pow(theta_e_N[i - 1], 0.2e1) * g_N[i - 1] / exp(theta_e_N[i - 1] / TV);
    Sum_eN_thetae3_g_div_e = Sum_eN_thetae3_g_div_e + pow(theta_e_N[i - 1], 0.3e1) * g_N[i - 1] / exp(theta_e_N[i - 1] / TV);
  }
  for (i = 1; i <= energy_level_N2; i++)
  {
    e_elec_N2_den = e_elec_N2_den + g_N2[i - 1] * exp(-theta_e_N2[i - 1] / TV);
    e_elec_N2_num = e_elec_N2_num + theta_e_N2[i - 1] * g_N2[i - 1] * exp(-theta_e_N2[i - 1] / TV);
    Sum_eN2_thetae2_g_div_e = Sum_eN2_thetae2_g_div_e + pow(theta_e_N2[i - 1], 0.2e1) * g_N2[i - 1] / exp(theta_e_N2[i - 1] / TV);
    Sum_eN2_thetae3_g_div_e = Sum_eN2_thetae3_g_div_e + pow(theta_e_N2[i - 1], 0.3e1) * g_N2[i - 1] / exp(theta_e_N2[i - 1] / TV);
  }
  E_elec_N = e_elec_N_num * R / e_elec_N_den / M_N;
  E_elec_N2 = e_elec_N2_num * R / e_elec_N2_den / M_N / 0.2e1;
  E_vib_N2 = R * theta_v_N2 / M_N / (exp(theta_v_N2 / TV) - 0.1e1) / 0.2e1;
  eV_N = RHO_N * E_elec_N / RHO;
  eV_N2 = RHO_N2 * (E_vib_N2 + E_elec_N2) / RHO;
  eV = eV_N + eV_N2;
    // "Aux relationships-----------------------------------------------"
;
  Mu_N = exp((A_N * log(T) + B_N) * log(T) + C_N) / 0.10e2;
  Mu_N2 = exp((A_N2 * log(T) + B_N2) * log(T) + C_N2) / 0.10e2;
  Mtot = 0.1e1 / (RHO_N / RHO / M_N + RHO_N2 / RHO / M_N / 0.2e1);
  Phi_N = RHO_N * Mtot / RHO / M_N + RHO_N2 * Mtot * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * sqrt(0.3e1) / RHO / M_N / 0.12e2;
  Phi_N2 = RHO_N * Mtot * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * sqrt(0.6e1) / RHO / M_N / 0.12e2 + RHO_N2 * Mtot / RHO / M_N / 0.2e1;
  Cv = 0.3e1 / 0.2e1 * RHO_N * R / RHO / M_N + 0.5e1 / 0.4e1 * RHO_N2 * R / RHO / M_N + RHO_N * R * Sum_eN_thetae2_g_div_e / RHO * pow(TV, -0.2e1) / M_N / e_elec_N_den - RHO_N * E_elec_N * E_elec_N * M_N / RHO / R * pow(TV, -0.2e1) + 0.2e1 * RHO_N2 * exp(theta_v_N2 / TV) * M_N * E_vib_N2 * E_vib_N2 / RHO / R * pow(TV, -0.2e1) + RHO_N2 * R * Sum_eN2_thetae2_g_div_e / RHO * pow(TV, -0.2e1) / M_N / e_elec_N2_den / 0.2e1 - 0.2e1 * RHO_N2 * E_elec_N2 * E_elec_N2 * M_N / RHO / R * pow(TV, -0.2e1);
  Cp = Cv + R;
  Kappa_ev_mix = Mu_N * RHO_N * Mtot * R * Sum_eN_thetae2_g_div_e / RHO * pow(M_N, -0.2e1) / Phi_N * pow(TV, -0.2e1) / e_elec_N_den - Mu_N * RHO_N * Mtot * E_elec_N * E_elec_N / RHO / Phi_N / R * pow(TV, -0.2e1) + Mu_N2 * RHO_N2 * Mtot * R * Sum_eN2_thetae2_g_div_e / RHO * pow(M_N, -0.2e1) / Phi_N2 * pow(TV, -0.2e1) / e_elec_N2_den / 0.4e1 - Mu_N2 * RHO_N2 * Mtot * E_elec_N2 * E_elec_N2 / RHO / Phi_N2 / R * pow(TV, -0.2e1) + Mu_N2 * exp(theta_v_N2 / TV) * RHO_N2 * E_vib_N2 * E_vib_N2 * Mtot / RHO / Phi_N2 / R * pow(TV, -0.2e1);
  Kappa_tr_mix = 0.15e2 / 0.4e1 * Mu_N * RHO_N * Mtot * R / RHO * pow(M_N, -0.2e1) / Phi_N + 0.19e2 / 0.16e2 * Mu_N2 * RHO_N2 * Mtot * R / RHO * pow(M_N, -0.2e1) / Phi_N2;
  Kappa_mix = Kappa_ev_mix + Kappa_tr_mix;
  Ds = Le * Kappa_mix / RHO / Cp;
    // "Aux relationships-----------------------------------------------"
;
  DKappa_ev_Dx = -(0.2e1 * A_N2 * log(T) + B_N2) * a_Tx * PI * T_x * R * Mu_N * Mtot * RHO_N * sin(a_Tx * PI * x / L) * Sum_eN_thetae2_g_div_e / L * pow(M_N, -0.2e1) / Phi_N / RHO / T * pow(TV, -0.2e1) / e_elec_N_den - (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * R * Mu_N * Mtot * RHO_N * sin(a_Tx * PI * x / L) * Sum_eN_thetae2_g_div_e / L * pow(M_N, -0.2e1) / Phi_N / RHO / T * pow(TV, -0.2e1) / e_elec_N_den + sqrt(0.3e1) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * a_rho_N2_x * PI * rho_N2_x * R * Mu_N * Mtot * Mtot * RHO_N * sin(a_rho_N2_x * PI * x / L) * Sum_eN_thetae2_g_div_e / L * pow(M_N, -0.3e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) * pow(TV, -0.2e1) / e_elec_N_den / 0.12e2 + sqrt(0.3e1) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * PI * R * Mu_N * pow(Mtot, 0.3e1) * RHO_N * RHO_N2 * Sum_eN_thetae2_g_div_e / L * pow(M_N, -0.4e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.3e1) * pow(TV, -0.2e1) / e_elec_N_den / 0.24e2 + sqrt(0.3e1) * (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * PI * R * Mu_N * Mtot * Mtot * RHO_N * RHO_N2 * Sum_eN_thetae2_g_div_e / L * pow(M_N, -0.4e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.4e1) * pow(TV, -0.2e1) / e_elec_N_den / 0.24e2 - (M_N * Phi_N * RHO - RHO_N * Mtot) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * Mu_N * Mtot * Mtot * RHO_N * Sum_eN_thetae2_g_div_e / L * pow(M_N, -0.4e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.3e1) * pow(TV, -0.2e1) / e_elec_N_den / 0.2e1 - (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (M_N * Phi_N * RHO - RHO_N * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * Mu_N * Mtot * RHO_N * Sum_eN_thetae2_g_div_e / L * pow(M_N, -0.4e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.4e1) * pow(TV, -0.2e1) / e_elec_N_den / 0.2e1 + 0.2e1 * exp(theta_v_N2 / TV) * a_Tvx * PI * Tv_x * Mu_N2 * Mtot * E_vib_N2 * E_vib_N2 * RHO_N2 * sin(a_Tvx * PI * x / L) / L / R / Phi_N2 / RHO * pow(TV, -0.3e1) - (M_N * Phi_N * RHO - RHO_N * Mtot) * a_rho_N_x * PI * rho_N_x * Mu_N * Mtot * E_elec_N * E_elec_N * cos(a_rho_N_x * PI * x / L) / L / M_N / R * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) * pow(TV, -0.2e1) + sqrt(0.6e1) * (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (-E_elec_N2 * E_elec_N2 + exp(theta_v_N2 / TV) * E_vib_N2 * E_vib_N2) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * PI * Mu_N2 * Mtot * Mtot * RHO_N * RHO_N2 / L * pow(M_N, -0.2e1) / R * pow(Phi_N2, -0.2e1) * pow(RHO, -0.4e1) * pow(TV, -0.2e1) / 0.24e2 - sqrt(0.3e1) * (-0.4e1 * sqrt(0.3e1) * M_N * Mu_N * Phi_N * Phi_N2 * Phi_N2 * E_elec_N * E_elec_N * RHO_N * RHO - 0.4e1 * sqrt(0.3e1) * M_N * Mu_N2 * Phi_N * Phi_N * Phi_N2 * E_elec_N2 * E_elec_N2 * RHO_N2 * RHO + 0.4e1 * sqrt(0.3e1) * Mu_N * Phi_N2 * Phi_N2 * Mtot * E_elec_N * E_elec_N * RHO_N * RHO_N + pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * Mu_N * Phi_N2 * Phi_N2 * Mtot * E_elec_N * E_elec_N * RHO_N * RHO_N2 + 0.2e1 * sqrt(0.3e1) * Mu_N2 * Phi_N * Phi_N * Mtot * E_elec_N2 * E_elec_N2 * RHO_N2 * RHO_N2) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * Mtot * Mtot / L * pow(M_N, -0.2e1) / R * pow(Phi_N, -0.2e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.3e1) * pow(TV, -0.2e1) / 0.24e2 + sqrt(0.3e1) * pow(0.2e1, 0.1e1 / 0.4e1) * (0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1)) * (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * R * Mu_N * Mu_N * Mtot * Mtot * RHO_N * RHO_N2 * sin(a_Tx * PI * x / L) * Sum_eN_thetae2_g_div_e / sqrt(Mu_N / Mu_N2) / L * pow(M_N, -0.3e1) / Mu_N2 * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) / T * pow(TV, -0.2e1) / e_elec_N_den / 0.12e2 + (0.3e1 * M_N * E_elec_N + 0.2e1 * R * TV) * a_Tvx * PI * Tv_x * Mu_N * Mtot * RHO_N * sin(a_Tvx * PI * x / L) * Sum_eN_thetae2_g_div_e / L * pow(M_N, -0.2e1) / Phi_N / RHO * pow(TV, -0.4e1) / e_elec_N_den + (M_N * Phi_N * RHO - RHO_N * Mtot) * a_rho_N_x * PI * rho_N_x * R * Mu_N * Mtot * cos(a_rho_N_x * PI * x / L) * Sum_eN_thetae2_g_div_e / L * pow(M_N, -0.3e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) * pow(TV, -0.2e1) / e_elec_N_den - Sum_eN_thetae3_g_div_e * a_Tvx * PI * Tv_x * R * Mu_N * Mtot * RHO_N * sin(a_Tvx * PI * x / L) / L * pow(M_N, -0.2e1) / Phi_N / RHO * pow(TV, -0.4e1) / e_elec_N_den - Sum_eN2_thetae3_g_div_e * a_Tvx * PI * Tv_x * R * Mu_N2 * Mtot * RHO_N2 * sin(a_Tvx * PI * x / L) / L * pow(M_N, -0.2e1) / Phi_N2 / RHO * pow(TV, -0.4e1) / e_elec_N2_den / 0.4e1 - (0.2e1 * A_N2 * log(T) + B_N2) * a_Tx * PI * T_x * R * Mu_N2 * Mtot * RHO_N2 * sin(a_Tx * PI * x / L) * Sum_eN2_thetae2_g_div_e / L * pow(M_N, -0.2e1) / Phi_N2 / RHO / T * pow(TV, -0.2e1) / e_elec_N2_den / 0.4e1 - sqrt(0.6e1) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * a_rho_N_x * PI * rho_N_x * R * Mu_N2 * Mtot * Mtot * RHO_N2 * cos(a_rho_N_x * PI * x / L) * Sum_eN2_thetae2_g_div_e / L * pow(M_N, -0.3e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) * pow(TV, -0.2e1) / e_elec_N2_den / 0.48e2 - sqrt(0.6e1) * pow(0.2e1, 0.3e1 / 0.4e1) * (0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1) * (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * R * Mu_N2 * Mu_N2 * Mtot * Mtot * RHO_N * RHO_N2 * sin(a_Tx * PI * x / L) * Sum_eN2_thetae2_g_div_e / sqrt(Mu_N2 / Mu_N) / L * pow(M_N, -0.3e1) / Mu_N * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) / T * pow(TV, -0.2e1) / e_elec_N2_den / 0.96e2 + (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (0.2e1 * M_N * Mu_N * Phi_N * Phi_N2 * Phi_N2 * E_elec_N * E_elec_N * RHO_N * RHO + 0.2e1 * M_N * Mu_N2 * Phi_N * Phi_N * Phi_N2 * E_elec_N2 * E_elec_N2 * RHO_N2 * RHO - 0.2e1 * Mu_N * Phi_N2 * Phi_N2 * Mtot * E_elec_N * E_elec_N * RHO_N * RHO_N - Mu_N2 * Phi_N * Phi_N * Mtot * E_elec_N2 * E_elec_N2 * RHO_N2 * RHO_N2) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * Mtot / L * pow(M_N, -0.2e1) / R * pow(Phi_N, -0.2e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.4e1) * pow(TV, -0.2e1) / 0.4e1 + (0.3e1 * M_N * E_elec_N2 + R * TV) * a_Tvx * PI * Tv_x * Mu_N2 * Mtot * RHO_N2 * sin(a_Tvx * PI * x / L) * Sum_eN2_thetae2_g_div_e / L * pow(M_N, -0.2e1) / Phi_N2 / RHO * pow(TV, -0.4e1) / e_elec_N2_den / 0.2e1 - (0.2e1 * M_N * Phi_N2 * RHO - RHO_N2 * Mtot) * a_rho_N2_x * PI * rho_N2_x * R * Mu_N2 * Mtot * sin(a_rho_N2_x * PI * x / L) * Sum_eN2_thetae2_g_div_e / L * pow(M_N, -0.3e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) * pow(TV, -0.2e1) / e_elec_N2_den / 0.8e1 + sqrt(0.6e1) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * PI * R * Mu_N2 * pow(Mtot, 0.3e1) * RHO_N * RHO_N2 * Sum_eN2_thetae2_g_div_e / L * pow(M_N, -0.4e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.3e1) * pow(TV, -0.2e1) / e_elec_N2_den / 0.96e2 + sqrt(0.6e1) * (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * PI * R * Mu_N2 * Mtot * Mtot * RHO_N * RHO_N2 * Sum_eN2_thetae2_g_div_e / L * pow(M_N, -0.4e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.4e1) * pow(TV, -0.2e1) / e_elec_N2_den / 0.96e2 - (0.2e1 * M_N * Phi_N2 * RHO - RHO_N2 * Mtot) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * Mu_N2 * Mtot * Mtot * RHO_N2 * Sum_eN2_thetae2_g_div_e / L * pow(M_N, -0.4e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.3e1) * pow(TV, -0.2e1) / e_elec_N2_den / 0.16e2 - (0.2e1 * M_N * Phi_N2 * RHO - RHO_N2 * Mtot) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * exp(theta_v_N2 / TV) * PI * Mu_N2 * Mtot * Mtot * E_vib_N2 * E_vib_N2 * RHO_N2 / L * pow(M_N, -0.2e1) / R * pow(Phi_N2, -0.2e1) * pow(RHO, -0.3e1) * pow(TV, -0.2e1) / 0.4e1 + (0.2e1 * M_N * Phi_N2 * RHO - RHO_N2 * Mtot) * a_rho_N2_x * PI * rho_N2_x * Mu_N2 * Mtot * E_elec_N2 * E_elec_N2 * sin(a_rho_N2_x * PI * x / L) / L / M_N / R * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) * pow(TV, -0.2e1) / 0.2e1 - (0.2e1 * Mu_N * Phi_N2 * E_elec_N * E_elec_N * RHO_N + 0.2e1 * Mu_N2 * Phi_N * E_elec_N2 * E_elec_N2 * RHO_N2) * a_Tvx * PI * Tv_x * Mtot * sin(a_Tvx * PI * x / L) / L / R / Phi_N / Phi_N2 / RHO * pow(TV, -0.3e1) - 0.4e1 * pow(exp(theta_v_N2 / TV), 0.2e1) * a_Tvx * PI * Tv_x * M_N * Mu_N2 * Mtot * pow(E_vib_N2, 0.3e1) * RHO_N2 * sin(a_Tvx * PI * x / L) / L * pow(R, -0.2e1) / Phi_N2 / RHO * pow(TV, -0.4e1) + theta_v_N2 * exp(theta_v_N2 / TV) * a_Tvx * PI * Tv_x * Mu_N2 * Mtot * E_vib_N2 * E_vib_N2 * RHO_N2 * sin(a_Tvx * PI * x / L) / L / R / Phi_N2 / RHO * pow(TV, -0.4e1) - (0.2e1 * M_N * Phi_N2 * RHO - RHO_N2 * Mtot) * exp(theta_v_N2 / TV) * a_rho_N2_x * PI * rho_N2_x * Mu_N2 * Mtot * E_vib_N2 * E_vib_N2 * sin(a_rho_N2_x * PI * x / L) / L / M_N / R * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) * pow(TV, -0.2e1) / 0.2e1 - sqrt(0.6e1) * pow(0.2e1, 0.3e1 / 0.4e1) * (-E_elec_N2 * E_elec_N2 + exp(theta_v_N2 / TV) * E_vib_N2 * E_vib_N2) * (0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1) * (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * Mu_N2 * Mu_N2 * Mtot * Mtot * RHO_N * RHO_N2 * sin(a_Tx * PI * x / L) / sqrt(Mu_N2 / Mu_N) / L / M_N / R / Mu_N * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) / T * pow(TV, -0.2e1) / 0.24e2 + (Mu_N * Phi_N2 * E_elec_N * E_elec_N * RHO_N + Mu_N2 * Phi_N * E_elec_N2 * E_elec_N2 * RHO_N2 - Mu_N2 * exp(theta_v_N2 / TV) * Phi_N * E_vib_N2 * E_vib_N2 * RHO_N2) * (0.2e1 * A_N2 * log(T) + B_N2) * a_Tx * PI * T_x * Mtot * sin(a_Tx * PI * x / L) / L / R / Phi_N / Phi_N2 / RHO / T * pow(TV, -0.2e1) - sqrt(0.3e1) * (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * PI * Mu_N * Mtot * Mtot * E_elec_N * E_elec_N * RHO_N * RHO_N2 / L * pow(M_N, -0.2e1) / R * pow(Phi_N, -0.2e1) * pow(RHO, -0.4e1) * pow(TV, -0.2e1) / 0.24e2 + sqrt(0.6e1) * (-E_elec_N2 * E_elec_N2 + exp(theta_v_N2 / TV) * E_vib_N2 * E_vib_N2) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * PI * Mu_N2 * pow(Mtot, 0.3e1) * RHO_N * RHO_N2 / L * pow(M_N, -0.2e1) / R * pow(Phi_N2, -0.2e1) * pow(RHO, -0.3e1) * pow(TV, -0.2e1) / 0.24e2 - (0.2e1 * Mu_N * Phi_N2 * pow(E_elec_N, 0.3e1) * RHO_N + 0.4e1 * Mu_N2 * Phi_N * pow(E_elec_N2, 0.3e1) * RHO_N2) * a_Tvx * PI * Tv_x * M_N * Mtot * sin(a_Tvx * PI * x / L) / L * pow(R, -0.2e1) / Phi_N / Phi_N2 / RHO * pow(TV, -0.4e1) - (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (0.2e1 * M_N * Phi_N2 * RHO - RHO_N2 * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * exp(theta_v_N2 / TV) * PI * Mu_N2 * Mtot * E_vib_N2 * E_vib_N2 * RHO_N2 / L * pow(M_N, -0.2e1) / R * pow(Phi_N2, -0.2e1) * pow(RHO, -0.4e1) * pow(TV, -0.2e1) / 0.4e1 + (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * Mu_N * Mtot * E_elec_N * E_elec_N * RHO_N * sin(a_Tx * PI * x / L) / L / R / Phi_N / RHO / T * pow(TV, -0.2e1) - sqrt(0.3e1) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * a_rho_N2_x * PI * rho_N2_x * Mu_N * Mtot * Mtot * E_elec_N * E_elec_N * RHO_N * sin(a_rho_N2_x * PI * x / L) / L / M_N / R * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) * pow(TV, -0.2e1) / 0.12e2 - sqrt(0.3e1) * pow(0.2e1, 0.1e1 / 0.4e1) * (0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1)) * (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * Mu_N * Mu_N * Mtot * Mtot * E_elec_N * E_elec_N * RHO_N * RHO_N2 * sin(a_Tx * PI * x / L) / sqrt(Mu_N / Mu_N2) / L / M_N / R / Mu_N2 * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) / T * pow(TV, -0.2e1) / 0.12e2 - sqrt(0.6e1) * (-E_elec_N2 * E_elec_N2 + exp(theta_v_N2 / TV) * E_vib_N2 * E_vib_N2) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * a_rho_N_x * PI * rho_N_x * Mu_N2 * Mtot * Mtot * RHO_N2 * cos(a_rho_N_x * PI * x / L) / L / M_N / R * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) * pow(TV, -0.2e1) / 0.12e2 - (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (0.2e1 * M_N * Phi_N2 * RHO - RHO_N2 * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * Mu_N2 * Mtot * RHO_N2 * Sum_eN2_thetae2_g_div_e / L * pow(M_N, -0.4e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.4e1) * pow(TV, -0.2e1) / e_elec_N2_den / 0.16e2;
  DKappa_tr_Dx = -0.15e2 / 0.4e1 * (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * R * Mu_N * Mtot * RHO_N * sin(a_Tx * PI * x / L) / L * pow(M_N, -0.2e1) / Phi_N / RHO / T - 0.19e2 / 0.192e3 * sqrt(0.6e1) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * a_rho_N_x * PI * rho_N_x * R * Mu_N2 * Mtot * Mtot * RHO_N2 * cos(a_rho_N_x * PI * x / L) / L * pow(M_N, -0.3e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) + 0.5e1 / 0.16e2 * sqrt(0.3e1) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * a_rho_N2_x * PI * rho_N2_x * R * Mu_N * Mtot * Mtot * RHO_N * sin(a_rho_N2_x * PI * x / L) / L * pow(M_N, -0.3e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) + 0.5e1 / 0.16e2 * sqrt(0.3e1) * pow(0.2e1, 0.1e1 / 0.4e1) * (0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1)) * (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * R * Mu_N * Mu_N * Mtot * Mtot * RHO_N * RHO_N2 * sin(a_Tx * PI * x / L) / sqrt(Mu_N / Mu_N2) / L * pow(M_N, -0.3e1) / Mu_N2 * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) / T - 0.19e2 / 0.384e3 * sqrt(0.6e1) * pow(0.2e1, 0.3e1 / 0.4e1) * (0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1) * (0.2e1 * A_N * log(T) - 0.2e1 * A_N2 * log(T) + B_N - B_N2) * a_Tx * PI * T_x * R * Mu_N2 * Mu_N2 * Mtot * Mtot * RHO_N * RHO_N2 * sin(a_Tx * PI * x / L) / sqrt(Mu_N2 / Mu_N) / L * pow(M_N, -0.3e1) / Mu_N * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) / T + 0.15e2 / 0.4e1 * (M_N * Phi_N * RHO - RHO_N * Mtot) * a_rho_N_x * PI * rho_N_x * R * Mu_N * Mtot * cos(a_rho_N_x * PI * x / L) / L * pow(M_N, -0.3e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.2e1) - 0.19e2 / 0.32e2 * (0.2e1 * M_N * Phi_N2 * RHO - RHO_N2 * Mtot) * a_rho_N2_x * PI * rho_N2_x * R * Mu_N2 * Mtot * sin(a_rho_N2_x * PI * x / L) / L * pow(M_N, -0.3e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.2e1) - (0.60e2 * Mu_N * Phi_N2 * RHO_N + 0.19e2 * Mu_N2 * Phi_N * RHO_N2) * (0.2e1 * A_N2 * log(T) + B_N2) * a_Tx * PI * T_x * R * Mtot * sin(a_Tx * PI * x / L) / L * pow(M_N, -0.2e1) / Phi_N / Phi_N2 / RHO / T / 0.16e2 + 0.5e1 / 0.32e2 * sqrt(0.3e1) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * PI * R * Mu_N * pow(Mtot, 0.3e1) * RHO_N * RHO_N2 / L * pow(M_N, -0.4e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.3e1) + 0.19e2 / 0.384e3 * sqrt(0.6e1) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * PI * R * Mu_N2 * pow(Mtot, 0.3e1) * RHO_N * RHO_N2 / L * pow(M_N, -0.4e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.3e1) + 0.5e1 / 0.32e2 * sqrt(0.3e1) * (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N / Mu_N2) * pow(0.2e1, 0.1e1 / 0.4e1), 0.2e1) * PI * R * Mu_N * Mtot * Mtot * RHO_N * RHO_N2 / L * pow(M_N, -0.4e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.4e1) + 0.19e2 / 0.384e3 * sqrt(0.6e1) * (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * pow(0.1e1 + sqrt(Mu_N2 / Mu_N) * pow(0.2e1, 0.3e1 / 0.4e1) / 0.2e1, 0.2e1) * PI * R * Mu_N2 * Mtot * Mtot * RHO_N * RHO_N2 / L * pow(M_N, -0.4e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.4e1) - 0.15e2 / 0.8e1 * (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (M_N * Phi_N * RHO - RHO_N * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * Mu_N * Mtot * RHO_N / L * pow(M_N, -0.4e1) * pow(Phi_N, -0.2e1) * pow(RHO, -0.4e1) - 0.19e2 / 0.64e2 * (0.2e1 * RHO * M_N - 0.2e1 * RHO_N * Mtot - RHO_N2 * Mtot) * (0.2e1 * M_N * Phi_N2 * RHO - RHO_N2 * Mtot) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * Mu_N2 * Mtot * RHO_N2 / L * pow(M_N, -0.4e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.4e1) - (0.120e3 * M_N * Mu_N * Phi_N * Phi_N2 * Phi_N2 * RHO_N * RHO + 0.38e2 * M_N * Mu_N2 * Phi_N * Phi_N * Phi_N2 * RHO_N2 * RHO - 0.120e3 * Mu_N * Phi_N2 * Phi_N2 * Mtot * RHO_N * RHO_N - 0.19e2 * Mu_N2 * Phi_N * Phi_N * Mtot * RHO_N2 * RHO_N2) * (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * Mtot * Mtot / L * pow(M_N, -0.4e1) * pow(Phi_N, -0.2e1) * pow(Phi_N2, -0.2e1) * pow(RHO, -0.3e1) / 0.64e2;
  DKappa_mix_Dx = DKappa_ev_Dx + DKappa_tr_Dx;
  DeV_N_Dx = a_rho_N_x * PI * rho_N_x * E_elec_N * cos(a_rho_N_x * PI * x / L) / L / RHO + a_Tvx * PI * Tv_x * M_N * E_elec_N * E_elec_N * RHO_N * sin(a_Tvx * PI * x / L) / L / R / RHO * pow(TV, -0.2e1) - a_Tvx * PI * Tv_x * R * RHO_N * sin(a_Tvx * PI * x / L) * Sum_eN_thetae2_g_div_e / L / M_N / RHO * pow(TV, -0.2e1) / e_elec_N_den - (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * E_elec_N * RHO_N / L * pow(RHO, -0.2e1);
  DeV_N2_Dx = -(E_vib_N2 + E_elec_N2) * a_rho_N2_x * PI * rho_N2_x * sin(a_rho_N2_x * PI * x / L) / L / RHO - (-0.2e1 * E_elec_N2 * E_elec_N2 + 0.2e1 * exp(theta_v_N2 / TV) * E_vib_N2 * E_vib_N2) * a_Tvx * PI * Tv_x * M_N * RHO_N2 * sin(a_Tvx * PI * x / L) / L / R / RHO * pow(TV, -0.2e1) - a_Tvx * PI * Tv_x * R * RHO_N2 * sin(a_Tvx * PI * x / L) * Sum_eN2_thetae2_g_div_e / L / M_N / RHO * pow(TV, -0.2e1) / e_elec_N2_den / 0.2e1 - (E_vib_N2 + E_elec_N2) * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * RHO_N2 / L * pow(RHO, -0.2e1);
  DeV_Dx = DeV_N_Dx + DeV_N2_Dx;
  DCp_Dx = 0.2e1 * a_Tvx * PI * Tv_x * M_N * theta_v_N2 * exp(theta_v_N2 / TV) * E_vib_N2 * E_vib_N2 * RHO_N2 * sin(a_Tvx * PI * x / L) / L / R / RHO * pow(TV, -0.4e1) - a_rho_N_x * PI * rho_N_x * M_N * E_elec_N * E_elec_N * cos(a_rho_N_x * PI * x / L) / L / R / RHO * pow(TV, -0.2e1) + a_rho_N_x * PI * rho_N_x * R * cos(a_rho_N_x * PI * x / L) * Sum_eN_thetae2_g_div_e / L / M_N / RHO * pow(TV, -0.2e1) / e_elec_N_den - (-0.2e1 * E_elec_N2 * E_elec_N2 + 0.2e1 * exp(theta_v_N2 / TV) * E_vib_N2 * E_vib_N2) * a_rho_N2_x * PI * rho_N2_x * M_N * sin(a_rho_N2_x * PI * x / L) / L / R / RHO * pow(TV, -0.2e1) - a_rho_N2_x * PI * rho_N2_x * R * sin(a_rho_N2_x * PI * x / L) * Sum_eN2_thetae2_g_div_e / L / M_N / RHO * pow(TV, -0.2e1) / e_elec_N2_den / 0.2e1 + 0.2e1 * a_Tvx * PI * Tv_x * R * RHO_N * sin(a_Tvx * PI * x / L) * Sum_eN_thetae2_g_div_e / L / M_N / RHO * pow(TV, -0.3e1) / e_elec_N_den + a_Tvx * PI * Tv_x * R * RHO_N2 * sin(a_Tvx * PI * x / L) * Sum_eN2_thetae2_g_div_e / L / M_N / RHO * pow(TV, -0.3e1) / e_elec_N2_den + 0.3e1 * a_Tvx * PI * Tv_x * E_elec_N * RHO_N * sin(a_Tvx * PI * x / L) * Sum_eN_thetae2_g_div_e / L / RHO * pow(TV, -0.4e1) / e_elec_N_den + 0.3e1 * a_Tvx * PI * Tv_x * E_elec_N2 * RHO_N2 * sin(a_Tvx * PI * x / L) * Sum_eN2_thetae2_g_div_e / L / RHO * pow(TV, -0.4e1) / e_elec_N2_den + (0.4e1 * exp(theta_v_N2 / TV) * E_vib_N2 * E_vib_N2 * RHO_N2 - 0.2e1 * RHO_N * E_elec_N * E_elec_N - 0.4e1 * RHO_N2 * E_elec_N2 * E_elec_N2) * a_Tvx * PI * Tv_x * M_N * sin(a_Tvx * PI * x / L) / L / R / RHO * pow(TV, -0.3e1) - a_Tvx * PI * Tv_x * R * RHO_N * sin(a_Tvx * PI * x / L) * Sum_eN_thetae3_g_div_e / L / M_N / RHO * pow(TV, -0.4e1) / e_elec_N_den - a_Tvx * PI * Tv_x * R * RHO_N2 * sin(a_Tvx * PI * x / L) * Sum_eN2_thetae3_g_div_e / L / M_N / RHO * pow(TV, -0.4e1) / e_elec_N2_den / 0.2e1 + (0.6e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - 0.5e1 * a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R / L / M_N / RHO / 0.4e1 - (0.8e1 * pow(exp(theta_v_N2 / TV), 0.2e1) * pow(E_vib_N2, 0.3e1) * RHO_N2 + 0.2e1 * pow(E_elec_N, 0.3e1) * RHO_N + 0.8e1 * pow(E_elec_N2, 0.3e1) * RHO_N2) * a_Tvx * PI * Tv_x * M_N * M_N * sin(a_Tvx * PI * x / L) / L * pow(R, -0.2e1) / RHO * pow(TV, -0.4e1) - (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * (0.6e1 * RHO_N + 0.5e1 * RHO_N2) * PI * R / L / M_N * pow(RHO, -0.2e1) / 0.4e1 - (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * RHO_N * Sum_eN_thetae2_g_div_e / L / M_N * pow(RHO, -0.2e1) * pow(TV, -0.2e1) / e_elec_N_den - (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * R * RHO_N2 * Sum_eN2_thetae2_g_div_e / L / M_N * pow(RHO, -0.2e1) * pow(TV, -0.2e1) / e_elec_N2_den / 0.2e1 - (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * (0.2e1 * exp(theta_v_N2 / TV) * E_vib_N2 * E_vib_N2 * RHO_N2 - RHO_N * E_elec_N * E_elec_N - 0.2e1 * RHO_N2 * E_elec_N2 * E_elec_N2) * PI * M_N / L / R * pow(RHO, -0.2e1) * pow(TV, -0.2e1);
    // "Contribution from the transient term to the total source term ----------------------------------------------"
;
  Q_eV_time = -E_elec_N * a_rho_N_t * PI * rho_N_t * sin(a_rho_N_t * PI * t / Lt) / Lt + (E_vib_N2 + E_elec_N2) * a_rho_N2_t * PI * rho_N2_t * cos(a_rho_N2_t * PI * t / Lt) / Lt + a_Tvt * PI * Tv_t * R * RHO_N * cos(a_Tvt * PI * t / Lt) * Sum_eN_thetae2_g_div_e / Lt / M_N * pow(TV, -0.2e1) / e_elec_N_den + a_Tvt * PI * Tv_t * R * RHO_N2 * cos(a_Tvt * PI * t / Lt) * Sum_eN2_thetae2_g_div_e / Lt / M_N * pow(TV, -0.2e1) / e_elec_N2_den / 0.2e1 + (0.2e1 * exp(theta_v_N2 / TV) * E_vib_N2 * E_vib_N2 * RHO_N2 - RHO_N * E_elec_N * E_elec_N - 0.2e1 * RHO_N2 * E_elec_N2 * E_elec_N2) * a_Tvt * PI * Tv_t * M_N * cos(a_Tvt * PI * t / Lt) / Lt / R * pow(TV, -0.2e1);
    // "Contribution from the convective term to the total source term --------------------------------------------"
;
  Q_eV_convection = eV * a_ux * PI * u_x * RHO * cos(a_ux * PI * x / L) / L + DeV_Dx * RHO * U + (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * eV * PI * U / L;
    // "Contribution from the heat flux  term to the total source term --------------------------------------------"
;
  Q_eV_heatflux = a_Tvx * a_Tvx * PI * PI * Tv_x * Kappa_ev_mix * cos(a_Tvx * PI * x / L) * pow(L, -0.2e1) + DKappa_ev_Dx * a_Tvx * PI * Tv_x * sin(a_Tvx * PI * x / L) / L;
    // "Contribution from production  term to the total source term --------------------------------------------"
;
  Q_eV_production = -w_dot_V;
    // "Contribution from the diffusion term to the total source term ----------------------------------------------"
;
  Q_eV_diffusion = eV_N * a_rho_N_x * a_rho_N_x * PI * PI * rho_N_x * Ds * sin(a_rho_N_x * PI * x / L) * pow(L, -0.2e1) + eV_N2 * a_rho_N2_x * a_rho_N2_x * PI * PI * rho_N2_x * Ds * cos(a_rho_N2_x * PI * x / L) * pow(L, -0.2e1) + eV_N * DCp_Dx * a_rho_N_x * PI * rho_N_x * Ds * cos(a_rho_N_x * PI * x / L) / L / Cp - eV_N2 * DCp_Dx * a_rho_N2_x * PI * rho_N2_x * Ds * sin(a_rho_N2_x * PI * x / L) / L / Cp - DeV_N_Dx * a_rho_N_x * PI * rho_N_x * Ds * cos(a_rho_N_x * PI * x / L) / L + DeV_N2_Dx * a_rho_N2_x * PI * rho_N2_x * Ds * sin(a_rho_N2_x * PI * x / L) / L + (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - 0.2e1 * a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * eV_N * a_rho_N_x * PI * PI * rho_N_x * Ds * cos(a_rho_N_x * PI * x / L) * pow(L, -0.2e1) / RHO - (0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - 0.2e1 * a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * eV_N2 * a_rho_N2_x * PI * PI * rho_N2_x * Ds * sin(a_rho_N2_x * PI * x / L) * pow(L, -0.2e1) / RHO - DKappa_mix_Dx * eV_N * a_rho_N_x * PI * rho_N_x * Le * cos(a_rho_N_x * PI * x / L) / L / Cp / RHO + eV_N2 * DKappa_mix_Dx * a_rho_N2_x * PI * rho_N2_x * Le * sin(a_rho_N2_x * PI * x / L) / L / Cp / RHO - (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * eV_N * DCp_Dx * PI * Ds * RHO_N / L / Cp / RHO - (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * eV_N2 * DCp_Dx * PI * Ds * RHO_N2 / L / Cp / RHO + (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * DeV_N_Dx * PI * Ds * RHO_N / L / RHO + (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * DeV_N2_Dx * PI * Ds * RHO_N2 / L / RHO - (eV_N * RHO_N + eV_N2 * RHO_N2) * (a_rho_N_x * a_rho_N_x * rho_N_x * sin(a_rho_N_x * PI * x / L) + a_rho_N2_x * a_rho_N2_x * rho_N2_x * cos(a_rho_N2_x * PI * x / L)) * PI * PI * Ds * pow(L, -0.2e1) / RHO - (0.2e1 * eV_N * RHO_N + 0.2e1 * eV_N2 * RHO_N2) * pow(a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L), 0.2e1) * PI * PI * Ds * pow(L, -0.2e1) * pow(RHO, -0.2e1) + (eV_N * RHO_N + eV_N2 * RHO_N2) * DKappa_mix_Dx * (a_rho_N_x * rho_N_x * cos(a_rho_N_x * PI * x / L) - a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * PI * x / L)) * PI * Le / L / Cp * pow(RHO, -0.2e1);
    // "Total source term -------------------------------------------------------------------------"
;
  Q_eV = Q_eV_time + Q_eV_convection + Q_eV_production + Q_eV_heatflux + Q_eV_diffusion;
  return(Q_eV);
}
