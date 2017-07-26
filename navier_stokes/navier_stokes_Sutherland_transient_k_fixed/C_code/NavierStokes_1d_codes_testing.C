#include <iostream>
#include <cmath>

using std::cout;
using std::endl;


// Global variables ================================
// Manufactured solutions --------------------------  
const double PI = M_PI;  
double L = 1.5, Lt=3.,
rho_0 =2., rho_x=3., a_rhox=5., rho_y=7., a_rhoy=11., rho_z=13., a_rhoz=17., rho_t=19, a_rhot=23,
p_0=29, p_x=31., a_px=37., p_y=41., a_py=43, p_z=47., a_pz=53, p_t=61., a_pt=67., 
u_0=2., u_x=13., a_ux=29., u_y=37., a_uy=-3.,  u_z=-17., a_uz= 83.,  u_t=-5., a_ut=19.,
v_0=3., v_x=11., a_vx=23., v_y=31., a_vy= 5.,  v_z= 19., a_vz=-47.,  v_t= 2., a_vt=13.,
w_0=5., w_x=7.,  a_wx=19., w_y=29., a_wy=-7.,  w_z= 23., a_wz= 11.,  w_t= 11., a_wt=17. ;

// double L = 12, Lt=12.,
// rho_0 =12., rho_x=12., a_rhox=12., rho_y=12., a_rhoy=12., rho_z=12., a_rhoz=12, rho_t=12, a_rhot=12.,
// p_0=12., p_x=12., a_px=12., p_y=12., a_py=12.,  p_z=12., a_pz=12.,  p_t=12., a_pt=12., 
// u_0=12., u_x=12., a_ux=12., u_y=12., a_uy=12.,  u_z=12., a_uz=12.,  u_t=12., a_ut=12.,
// v_0=12., v_x=12., a_vx=12., v_y=12., a_vy=12.,  v_z=12., a_vz=12.,  v_t=12., a_vt=12.,
// w_0=12., w_x=12., a_wx=12., w_y=12., a_wy=12.,  w_z=12., a_wz=12.,  w_t=12., a_wt=12. ;

// Model---------------------------------------------
double B_mu = 110.4,  A_mu=1.458e-6, Gamma  = 1.4 , R = 287, Pr  = 0.7;

// Headers =========================================
double SourceQ_e (double x, double t);
double SourceQ_u (double x, double t);
double SourceQ_rho (double x, double t);

// Main=============================================
int main(int argc, char **argv)
{
  double Qe=0., Qu=0., Qrho=0.;
  double x= 191., t=5.;//, y= 211., z= 163.;
  double p,u,T,e,rho;
  rho = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_t * sin(a_rhot * PI * t / Lt);
  p = p_0 + p_x * cos(a_px * PI * x / L) + p_t * cos(a_pt * PI * t / Lt);
  u = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / Lt);
  T=p/rho/R;
  e=R*T/(Gamma-1);
  
  cout << endl<<"rho_an= " << rho << endl <<
		"p_an= " << p << endl <<
		"u_an= " << u << endl <<
		"T=    " << T << endl <<
		"e   = " << e << endl;
		
  Qe= SourceQ_e ( x,  t);
  Qu= SourceQ_u ( x,  t);
  Qrho = SourceQ_rho ( x,  t);
  
  cout << endl<< "SourceQ_rho(x,t)= " << Qrho << endl <<
		"SourceQ_u(x,t)= "   << Qu << endl <<
		"SourceQ_e(x,t)= "   << Qe << endl;
   
  return 0;
}

//==================================================


double SourceQ_e (double x, double t)
{
  double RHO;
  double P;
  double U;
  double T;
  double MU;
  double DMu_Dx;
  double kappa;
  double Q_e;
  double Q_e_convection;
  double Q_e_work_pressure;
  double Q_e_work_viscous;
  double Q_e_conduction;
  double Q_e_time;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_t * sin(a_rhot * PI * t / Lt);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / Lt);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_t * cos(a_pt * PI * t / Lt);
  T = P / RHO / R;
  MU = A_mu * pow(T, 0.3e1 / 0.2e1) / (T + B_mu);
  DMu_Dx = a_rhox * PI * rho_x * MU * MU * cos(a_rhox * PI * x / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhox * PI * rho_x * MU * cos(a_rhox * PI * x / L) / L / RHO + a_px * PI * p_x * MU * MU * sin(a_px * PI * x / L) / A_mu / L / R / RHO * pow(T, -0.3e1 / 0.2e1) - 0.3e1 / 0.2e1 * a_px * PI * p_x * MU * sin(a_px * PI * x / L) / L / R / RHO / T;
  kappa = Gamma * R * MU / (Gamma - 0.1e1) / Pr;
  Q_e_convection = a_rhox * PI * rho_x * pow(U, 0.3e1) * cos(a_rhox * PI * x / L) / L / 0.2e1 + 0.3e1 / 0.2e1 * a_ux * PI * u_x * RHO * U * U * cos(a_ux * PI * x / L) / L - a_px * PI * p_x * U * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + a_ux * PI * u_x * P * cos(a_ux * PI * x / L) / (Gamma - 0.1e1) / L;
  Q_e_work_pressure = -a_px * PI * p_x * U * sin(a_px * PI * x / L) / L + a_ux * PI * u_x * P * cos(a_ux * PI * x / L) / L;
  Q_e_conduction = -0.2e1 * kappa * a_rhox * a_rhox * PI * PI * rho_x * rho_x * P * pow(cos(a_rhox * PI * x / L), 0.2e1) * pow(L, -0.2e1) / R * pow(RHO, -0.3e1) - 0.2e1 * kappa * a_rhox * a_px * PI * PI * rho_x * p_x * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) - kappa * a_rhox * a_rhox * PI * PI * rho_x * P * sin(a_rhox * PI * x / L) * pow(L, -0.2e1) / R * pow(RHO, -0.2e1) + kappa * a_px * a_px * PI * PI * p_x * cos(a_px * PI * x / L) * pow(L, -0.2e1) / R / RHO + kappa * DMu_Dx * a_rhox * PI * rho_x * P * cos(a_rhox * PI * x / L) / L / R / MU * pow(RHO, -0.2e1) + kappa * DMu_Dx * a_px * PI * p_x * sin(a_px * PI * x / L) / L / R / MU / RHO;
  Q_e_work_viscous = -0.4e1 / 0.3e1 * a_ux * a_ux * MU * PI * PI * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * a_ux * a_ux * MU * PI * PI * u_x * U * sin(a_ux * PI * x / L) * pow(L, -0.2e1);
  Q_e_time = -a_ut * u_t * PI * RHO * U * sin(a_ut * PI * t / Lt) / Lt + a_rhot * rho_t * PI * U * U * cos(a_rhot * PI * t / Lt) / Lt / 0.2e1 - a_pt * p_t * PI * sin(a_pt * PI * t / Lt) / (Gamma - 0.1e1) / Lt;
  Q_e = Q_e_convection + Q_e_work_pressure + Q_e_work_viscous + Q_e_conduction + Q_e_time;
  return(Q_e);
}
//==================================================
double SourceQ_u (double x, double t)
{
  double RHO;
  double P;
  double U;
  double T;
  double MU;
  double DMu_Dx;
  double Q_u;
  double Q_u_convection;
  double Q_u_pressure;
  double Q_u_viscous;
  double Q_u_time;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_t * sin(a_rhot * PI * t / Lt);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / Lt);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_t * cos(a_pt * PI * t / Lt);
  T = P / RHO / R;
  MU = A_mu * pow(T, 0.3e1 / 0.2e1) / (T + B_mu);
  DMu_Dx = a_rhox * PI * rho_x * MU * MU * cos(a_rhox * PI * x / L) / A_mu / L / RHO / sqrt(T) - 0.3e1 / 0.2e1 * a_rhox * PI * rho_x * MU * cos(a_rhox * PI * x / L) / L / RHO + a_px * PI * p_x * MU * MU * sin(a_px * PI * x / L) / A_mu / L / R / RHO * pow(T, -0.3e1 / 0.2e1) - 0.3e1 / 0.2e1 * a_px * PI * p_x * MU * sin(a_px * PI * x / L) / L / R / RHO / T;
  Q_u_convection = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L + 0.2e1 * a_ux * PI * u_x * RHO * U * cos(a_ux * PI * x / L) / L;
  Q_u_pressure = -a_px * PI * p_x * sin(a_px * PI * x / L) / L;
  Q_u_viscous = 0.4e1 / 0.3e1 * MU * a_ux * a_ux * PI * PI * u_x * sin(a_ux * PI * x / L) * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * DMu_Dx * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L;
  Q_u_time = a_rhot * PI * rho_t * U * cos(a_rhot * PI * t / Lt) / Lt - a_ut * PI * u_t * RHO * sin(a_ut * PI * t / Lt) / Lt;
  Q_u = Q_u_convection + Q_u_pressure + Q_u_viscous + Q_u_time;
  return(Q_u);
}
//==================================================


double SourceQ_rho (double x, double t)
{
  double RHO;
  double U;
  double Q_rho;
  double Q_rho_convection;
  double Q_rho_time;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_t * sin(a_rhot * PI * t / Lt);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_t * cos(a_ut * PI * t / Lt);
  Q_rho_convection = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L + a_ux * PI * u_x * RHO * cos(a_ux * PI * x / L) / L;
  Q_rho_time = a_rhot * rho_t * cos(a_rhot * PI * t / Lt) * PI / Lt;
  Q_rho = Q_rho_convection + Q_rho_time;
  return(Q_rho);
}
//==================================================
