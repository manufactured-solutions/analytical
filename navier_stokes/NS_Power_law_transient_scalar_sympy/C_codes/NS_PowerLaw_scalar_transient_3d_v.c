/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#include "NS_PowerLaw_scalar_transient_3d_v.h"
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

double Mu(double T, double T_ref, double beta, double mu_ref) {
  return mu_ref*pow((T/T_ref),beta);
}

double DMu_Dx(double L, double R, double Rho, double T, double T_ref, double a_px, double a_rhox, double beta, double mu_ref, double p_x, double rho_x, double x) {
  return -M_PI*a_rhox*beta*mu_ref*rho_x*pow(T,beta)*pow((1.0/T_ref),beta)*cos(M_PI*a_rhox*x/L)/(L*Rho) - M_PI*a_px*beta*mu_ref*p_x*pow(T,beta)*pow((1.0/T_ref),beta)*sin(M_PI*a_px*x/L)/(L*R*Rho*T);
}

double DMu_Dy(double L, double R, double Rho, double T, double T_ref, double a_py, double a_rhoy, double beta, double mu_ref, double p_y, double rho_y, double y) {
  return M_PI*a_rhoy*beta*mu_ref*rho_y*pow(T,beta)*pow((1.0/T_ref),beta)*sin(M_PI*a_rhoy*y/L)/(L*Rho) + M_PI*a_py*beta*mu_ref*p_y*pow(T,beta)*pow((1.0/T_ref),beta)*cos(M_PI*a_py*y/L)/(L*R*Rho*T);
}

double DMu_Dz(double L, double R, double Rho, double T, double T_ref, double a_pz, double a_rhoz, double beta, double mu_ref, double p_z, double rho_z, double z) {
  return -M_PI*a_rhoz*beta*mu_ref*rho_z*pow(T,beta)*pow((1.0/T_ref),beta)*cos(M_PI*a_rhoz*z/L)/(L*Rho) - M_PI*a_pz*beta*mu_ref*p_z*pow(T,beta)*pow((1.0/T_ref),beta)*sin(M_PI*a_pz*z/L)/(L*R*Rho*T);
}

double Q_v_convection(double L, double Rho, double U, double V, double W, double a_rhox, double a_rhoy, double a_rhoz, double a_ux, double a_vx, double a_vy, double a_vz, double a_wz, double rho_x, double rho_y, double rho_z, double u_x, double v_x, double v_y, double v_z, double w_z, double x, double y, double z) {
  return M_PI*Rho*V*(a_ux*u_x*cos(M_PI*a_ux*x/L) - a_wz*w_z*sin(M_PI*a_wz*z/L) + 2*a_vy*v_y*cos(M_PI*a_vy*y/L))/L - M_PI*a_rhoy*rho_y*pow(V,2)*sin(M_PI*a_rhoy*y/L)/L + M_PI*Rho*W*a_vz*v_z*cos(M_PI*a_vz*z/L)/L + M_PI*U*V*a_rhox*rho_x*cos(M_PI*a_rhox*x/L)/L + M_PI*V*W*a_rhoz*rho_z*cos(M_PI*a_rhoz*z/L)/L - M_PI*Rho*U*a_vx*v_x*sin(M_PI*a_vx*x/L)/L;
}

double Q_v_gradp(double L, double a_py, double p_y, double y) {
  return M_PI*a_py*p_y*cos(M_PI*a_py*y/L)/L;
}

double Q_v_viscous(double DMu_Dx, double DMu_Dy, double DMu_Dz, double L, double LAMBDA, double Mu, double a_ux, double a_uy, double a_vx, double a_vy, double a_vz, double a_wy, double a_wz, double alpha, double u_x, double u_y, double v_x, double v_y, double v_z, double w_y, double w_z, double x, double y, double z) {
  return M_PI*DMu_Dx*(a_uy*u_y*sin(M_PI*a_uy*y/L) + a_vx*v_x*sin(M_PI*a_vx*x/L))/L + M_PI*DMu_Dz*(-a_vz*v_z*cos(M_PI*a_vz*z/L) - a_wy*w_y*cos(M_PI*a_wy*y/L))/L + Mu*pow(M_PI,2)*(v_x*pow(a_vx,2)*cos(M_PI*a_vx*x/L) + v_z*pow(a_vz,2)*sin(M_PI*a_vz*z/L) + 2*v_y*pow(a_vy,2)*sin(M_PI*a_vy*y/L))/pow(L,2) - M_PI*DMu_Dy*(2.0/3.0 - alpha)*(a_wz*w_z*sin(M_PI*a_wz*z/L) - a_ux*u_x*cos(M_PI*a_ux*x/L) - a_vy*v_y*cos(M_PI*a_vy*y/L))/L + LAMBDA*v_y*pow(M_PI,2)*pow(a_vy,2)*sin(M_PI*a_vy*y/L)/pow(L,2) - 2*M_PI*DMu_Dy*a_vy*v_y*cos(M_PI*a_vy*y/L)/L;
}

double Q_v_time(double Lt, double Rho, double V, double a_rhot, double a_vt, double rho_t, double t, double v_t) {
  return M_PI*Rho*a_vt*v_t*cos(M_PI*a_vt*t/Lt)/Lt + M_PI*V*a_rhot*rho_t*cos(M_PI*a_rhot*t/Lt)/Lt;
}

double Q_v(double Q_v_convection, double Q_v_gradp, double Q_v_time, double Q_v_viscous) {
  return Q_v_convection + Q_v_gradp + Q_v_time + Q_v_viscous;
}


