#include <math.h>

double SourceQ_u (double x)
{
  double Q_u;
  double RHO;
  double U;
  double P;
  double M1;
  double M2;
  double MU;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  P = p_0 + p_x * cos(a_px * pi * x / L);
  M1 = A_mu * pow(P / R / RHO, 0.3e1 / 0.2e1);
  M2 = P / R / RHO + B_mu;
  MU = M1 / M2;
  Q_u = 0.2e1 / 0.3e1 * (0.3e1 * B_mu * R * RHO + P) * a_rhox * a_ux * PI * PI * rho_x * u_x * MU * cos(a_rhox * PI * x / L) * cos(a_ux * PI * x / L) / (B_mu * R * RHO + P) * pow(L, -0.2e1) / RHO + 0.2e1 / 0.3e1 * (0.3e1 * B_mu * R * RHO + P) * a_px * a_ux * PI * PI * p_x * u_x * MU * sin(a_px * PI * x / L) * cos(a_ux * PI * x / L) / (B_mu * R * RHO + P) * pow(L, -0.2e1) / P + 0.4e1 / 0.3e1 * a_ux * a_ux * PI * PI * u_x * MU * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L + 0.2e1 * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L;
  return(Q_u);
}
