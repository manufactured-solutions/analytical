/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#include "NS_Sutherland_scalar_transient_3d_v.h"
#include <math.h>

double Rho(double L, double Lt, double a_rhot, double a_rhox, double a_rhoy, double a_rhoz, double rho_0, double rho_t, double rho_x, double rho_y, double rho_z, double t, double x, double y, double z) {
  return rho_0 + rho_t*sin(M_PI*a_rhot*t/Lt) + rho_x*sin(M_PI*a_rhox*x/L) + rho_y*cos(M_PI*a_rhoy*y/L) + rho_z*sin(M_PI*a_rhoz*z/L);
}

double U(double L, double Lt, double a_ut, double a_ux, double a_uy, double a_uz, double t, double u_0, double u_t, double u_x, double u_y, double u_z, double x, double y, double z) {
  return u_0 + u_t*cos(M_PI*a_ut*t/Lt) + u_x*sin(M_PI*a_ux*x/L) + u_y*cos(M_PI*a_uy*y/L) + u_z*cos(M_PI*a_uz*z/L);
}

double V(double L, double Lt, double a_vt, double a_vx, double a_vy, double a_vz, double t, double v_0, double v_t, double v_x, double v_y, double v_z, double x, double y, double z) {
  return v_0 + v_t*sin(M_PI*a_vt*t/Lt) + v_x*cos(M_PI*a_vx*x/L) + v_y*sin(M_PI*a_vy*y/L) + v_z*sin(M_PI*a_vz*z/L);
}

double W(double L, double Lt, double a_wt, double a_wx, double a_wy, double a_wz, double t, double w_0, double w_t, double w_x, double w_y, double w_z, double x, double y, double z) {
  return w_0 + w_t*cos(M_PI*a_wt*t/Lt) + w_x*sin(M_PI*a_wx*x/L) + w_y*sin(M_PI*a_wy*y/L) + w_z*cos(M_PI*a_wz*z/L);
}

double P(double L, double Lt, double a_pt, double a_px, double a_py, double a_pz, double p_0, double p_t, double p_x, double p_y, double p_z, double t, double x, double y, double z) {
  return p_0 + p_t*cos(M_PI*a_pt*t/Lt) + p_x*cos(M_PI*a_px*x/L) + p_y*sin(M_PI*a_py*y/L) + p_z*cos(M_PI*a_pz*z/L);
}

double T(double P, double R, double Rho) {
  return P/(R*Rho);
}

double A_mu(double B_mu, double T_ref, double mu_ref) {
  return mu_ref*pow(T_ref,-1.5)*(B_mu + T_ref);
}

double Mu(double A_mu, double B_mu, double T) {
  return A_mu*pow(T,1.5)/(B_mu + T);
}

double DMu_Dx(double A_mu, double B_mu, double L, double R, double Rho, double T, double a_px, double a_rhox, double p_x, double rho_x, double x) {
  return A_mu*pow(T,0.5)*(-1.5*M_PI*T*a_rhox*rho_x*cos(M_PI*a_rhox*x/L)/(L*Rho) - 1.5*M_PI*a_px*p_x*sin(M_PI*a_px*x/L)/(L*R*Rho))/(B_mu + T) + A_mu*pow(T,1.5)*(M_PI*T*a_rhox*rho_x*cos(M_PI*a_rhox*x/L)/(L*Rho) + M_PI*a_px*p_x*sin(M_PI*a_px*x/L)/(L*R*Rho))/pow((B_mu + T),2);
}

double DMu_Dy(double A_mu, double B_mu, double L, double R, double Rho, double T, double a_py, double a_rhoy, double p_y, double rho_y, double y) {
  return A_mu*pow(T,0.5)*(1.5*M_PI*T*a_rhoy*rho_y*sin(M_PI*a_rhoy*y/L)/(L*Rho) + 1.5*M_PI*a_py*p_y*cos(M_PI*a_py*y/L)/(L*R*Rho))/(B_mu + T) + A_mu*pow(T,1.5)*(-M_PI*T*a_rhoy*rho_y*sin(M_PI*a_rhoy*y/L)/(L*Rho) - M_PI*a_py*p_y*cos(M_PI*a_py*y/L)/(L*R*Rho))/pow((B_mu + T),2);
}

double DMu_Dz(double A_mu, double B_mu, double L, double R, double Rho, double T, double a_pz, double a_rhoz, double p_z, double rho_z, double z) {
  return A_mu*pow(T,0.5)*(-1.5*M_PI*T*a_rhoz*rho_z*cos(M_PI*a_rhoz*z/L)/(L*Rho) - 1.5*M_PI*a_pz*p_z*sin(M_PI*a_pz*z/L)/(L*R*Rho))/(B_mu + T) + A_mu*pow(T,1.5)*(M_PI*T*a_rhoz*rho_z*cos(M_PI*a_rhoz*z/L)/(L*Rho) + M_PI*a_pz*p_z*sin(M_PI*a_pz*z/L)/(L*R*Rho))/pow((B_mu + T),2);
}

double Q_v_convection(double L, double Rho, double U, double V, double W, double a_rhox, double a_rhoy, double a_rhoz, double a_ux, double a_vx, double a_vy, double a_vz, double a_wz, double rho_x, double rho_y, double rho_z, double u_x, double v_x, double v_y, double v_z, double w_z, double x, double y, double z) {
  return M_PI*Rho*V*(a_ux*u_x*cos(M_PI*a_ux*x/L) - a_wz*w_z*sin(M_PI*a_wz*z/L) + 2*a_vy*v_y*cos(M_PI*a_vy*y/L))/L - M_PI*a_rhoy*rho_y*pow(V,2)*sin(M_PI*a_rhoy*y/L)/L + M_PI*Rho*W*a_vz*v_z*cos(M_PI*a_vz*z/L)/L + M_PI*U*V*a_rhox*rho_x*cos(M_PI*a_rhox*x/L)/L + M_PI*V*W*a_rhoz*rho_z*cos(M_PI*a_rhoz*z/L)/L - M_PI*Rho*U*a_vx*v_x*sin(M_PI*a_vx*x/L)/L;
}

double Q_v_gradp(double L, double a_py, double p_y, double y) {
  return M_PI*a_py*p_y*cos(M_PI*a_py*y/L)/L;
}

double Q_v_viscous(double DMu_Dx, double DMu_Dy, double DMu_Dz, double L, double Mu, double a_ux, double a_uy, double a_vx, double a_vy, double a_vz, double a_wy, double a_wz, double u_x, double u_y, double v_x, double v_y, double v_z, double w_y, double w_z, double x, double y, double z) {
  return M_PI*DMu_Dx*(a_uy*u_y*sin(M_PI*a_uy*y/L) + a_vx*v_x*sin(M_PI*a_vx*x/L))/L + M_PI*DMu_Dy*(-4*a_vy*v_y*cos(M_PI*a_vy*y/L)/3 - 2*a_wz*w_z*sin(M_PI*a_wz*z/L)/3 + 2*a_ux*u_x*cos(M_PI*a_ux*x/L)/3)/L + M_PI*DMu_Dz*(-a_vz*v_z*cos(M_PI*a_vz*z/L) - a_wy*w_y*cos(M_PI*a_wy*y/L))/L + Mu*pow(M_PI,2)*(v_x*pow(a_vx,2)*cos(M_PI*a_vx*x/L) + v_z*pow(a_vz,2)*sin(M_PI*a_vz*z/L) + 4*v_y*pow(a_vy,2)*sin(M_PI*a_vy*y/L)/3)/pow(L,2);
}

double Q_v_time(double Lt, double Rho, double V, double a_rhot, double a_vt, double rho_t, double t, double v_t) {
  return M_PI*Rho*a_vt*v_t*cos(M_PI*a_vt*t/Lt)/Lt + M_PI*V*a_rhot*rho_t*cos(M_PI*a_rhot*t/Lt)/Lt;
}

double Q_v(double Q_v_convection, double Q_v_gradp, double Q_v_time, double Q_v_viscous) {
  return Q_v_convection + Q_v_gradp + Q_v_time + Q_v_viscous;
}


