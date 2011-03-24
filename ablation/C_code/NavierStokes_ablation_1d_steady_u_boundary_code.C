#include <math.h>

double DiscrepancyQ_u_vw (
  double x,
  double W_C,
  double W_C3,
  double A_C3Enc,
  double E_aC3nc,
  double k_B,
  double beta_C3,
  double m_C3)
{
  double Q_u_boundary;
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
  Q_u_boundary = -Mdot_C3C / RHO + U;
  return(Q_u_boundary);
}
