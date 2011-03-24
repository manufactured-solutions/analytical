#include <math.h>

double SourceQ_E (
  double x,
  double W_C,
  double W_C3,
  double A_C3Enc,
  double E_aC3nc,
  double k_B,
  double beta_C3,
  double m_C3,
  double qr,
  double alpha,
  double sigma,
  double epsilon,
  double Function_to_Calculate_h_C3)
{
  double Q_e;
  double P;
  double T;
  double RHO;
  double RHO_C;
  double RHO_C3;
  double MF_C3;
  double MF_C3E;
  double Mdot_C3C;
  double H_C3;
  RHO_C3 = rho_C3_0 + rho_C3_x * cos(a_rho_C3_x * PI * x / L);
  RHO_C = rho_C_0 + rho_C_x * sin(a_rho_C_x * PI * x / L);
  RHO = RHO_C + RHO_C3;
  T = T_0 + T_x * cos(a_Tx * PI * x / L);
  P = R * T * (RHO_C / W_C + RHO_C3 / W_C3);
  MF_C3 = RHO_C3 / RHO;
  MF_C3E = A_C3Enc * exp(-E_aC3nc / T) / P;
  Mdot_C3C = sqrt(T * k_B / PI / m_C3) * sqrt(0.2e1) * (-MF_C3 + MF_C3E) * RHO * beta_C3 / 0.2e1;
  H_C3 = Function_to_Calculate_h_C3;
  Q_e = k * a_Tx * PI * T_x * sin(a_Tx * PI * x / L) / L - sigma * epsilon * pow(T, 0.4e1) - Mdot_C3C * H_C3 + alpha * qr;
  return(Q_e);
}
