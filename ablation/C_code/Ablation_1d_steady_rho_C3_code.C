#include <math.h>

double SourceQ_rho_C3 (
  double x,
  double W_C,
  double W_C3,
  double A_C3Enc,
  double E_aC3nc,
  double k_B,
  double beta_C3,
  double m_C3)
{
  double Q_rho_C3;
  double RHO;
  double RHO_C;
  double RHO_C3;
  double T;
  double P;
  double MF_C3;
  double MF_C3E;
  double Mdot_C3C;
  RHO_C3 = rho_C3_0 + rho_C3_x * cos(a_rho_C3_x * PI * x / L);
  RHO_C = rho_C_0 + rho_C_x * sin(a_rho_C_x * PI * x / L);
  RHO = RHO_C + RHO_C3;
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  P = R * T * (RHO_C / W_C + RHO_C3 / W_C3);
  MF_C3 = RHO_C3 / RHO;
  MF_C3E = A_C3Enc * exp(-E_aC3nc / T) / P;
  Mdot_C3C = sqrt(T * k_B / PI / m_C3) * sqrt(0.2e1) * (-MF_C3 + MF_C3E) * RHO * beta_C3 / 0.2e1;
  Q_rho_C3 = a_rho_C_x * PI * rho_C_x * D_C3 * RHO_C3 * cos(a_rho_C_x * PI * x / L) / L / RHO + a_rho_C3_x * PI * rho_C3_x * D_C3 * RHO_C * sin(a_rho_C3_x * PI * x / L) / L / RHO + a_rho_C_x * (D_C - D_C3) * PI * rho_C_x * RHO_C3 * RHO_C3 * cos(a_rho_C_x * PI * x / L) / L * pow(RHO, -0.2e1) + a_rho_C3_x * (D_C - D_C3) * PI * rho_C3_x * RHO_C * RHO_C3 * sin(a_rho_C3_x * PI * x / L) / L * pow(RHO, -0.2e1) + Mdot_C3C * (MF_C3 - 0.1e1);
  return(Q_rho_C3);
}
