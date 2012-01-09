/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#include "NS_Sutherland_scalar_transient_3d_phi.h"
#include <math.h>

double U(double L, double Lt, double a_ut, double a_ux, double a_uy, double a_uz, double t, double u_0, double u_t, double u_x, double u_y, double u_z, double x, double y, double z) {
  return u_0 + u_t*cos(M_PI*a_ut*t/Lt) + u_x*sin(M_PI*a_ux*x/L) + u_y*cos(M_PI*a_uy*y/L) + u_z*cos(M_PI*a_uz*z/L);
}

double V(double L, double Lt, double a_vt, double a_vx, double a_vy, double a_vz, double t, double v_0, double v_t, double v_x, double v_y, double v_z, double x, double y, double z) {
  return v_0 + v_t*sin(M_PI*a_vt*t/Lt) + v_x*cos(M_PI*a_vx*x/L) + v_y*sin(M_PI*a_vy*y/L) + v_z*sin(M_PI*a_vz*z/L);
}

double W(double L, double Lt, double a_wt, double a_wx, double a_wy, double a_wz, double t, double w_0, double w_t, double w_x, double w_y, double w_z, double x, double y, double z) {
  return w_0 + w_t*cos(M_PI*a_wt*t/Lt) + w_x*sin(M_PI*a_wx*x/L) + w_y*sin(M_PI*a_wy*y/L) + w_z*cos(M_PI*a_wz*z/L);
}

double Phi(double L, double Lt, double a_phit, double a_phix, double a_phiy, double a_phiz, double phi_0, double phi_t, double phi_x, double phi_y, double phi_z, double t, double x, double y, double z) {
  return phi_0 + phi_t*cos(M_PI*a_phit*t/Lt) + phi_x*cos(M_PI*a_phix*x/L) + phi_y*cos(M_PI*a_phiy*y/L) + phi_z*sin(M_PI*a_phiz*z/L);
}

double Q_phi_time(double Lt, double Phi, double Rho, double a_phit, double a_rhot, double phi_t, double rho_t, double t) {
  return M_PI*Phi*a_rhot*rho_t*cos(M_PI*a_rhot*t/Lt)/Lt - M_PI*Rho*a_phit*phi_t*sin(M_PI*a_phit*t/Lt)/Lt;
}

double Q_phi_convection(double L, double Phi, double Rho, double U, double V, double W, double a_phix, double a_phiy, double a_phiz, double a_rhox, double a_rhoy, double a_rhoz, double a_ux, double a_vy, double a_wz, double phi_x, double phi_y, double phi_z, double rho_x, double rho_y, double rho_z, double u_x, double v_y, double w_z, double x, double y, double z) {
  return M_PI*Phi*Rho*(a_ux*u_x*cos(M_PI*a_ux*x/L) + a_vy*v_y*cos(M_PI*a_vy*y/L) - a_wz*w_z*sin(M_PI*a_wz*z/L))/L + M_PI*Phi*U*a_rhox*rho_x*cos(M_PI*a_rhox*x/L)/L + M_PI*Phi*W*a_rhoz*rho_z*cos(M_PI*a_rhoz*z/L)/L + M_PI*Rho*W*a_phiz*phi_z*cos(M_PI*a_phiz*z/L)/L - M_PI*Phi*V*a_rhoy*rho_y*sin(M_PI*a_rhoy*y/L)/L - M_PI*Rho*U*a_phix*phi_x*sin(M_PI*a_phix*x/L)/L - M_PI*Rho*V*a_phiy*phi_y*sin(M_PI*a_phiy*y/L)/L;
}

double Q_phi_diffusion(double Gamma, double L, double a_phix, double a_phiy, double a_phiz, double phi_x, double phi_y, double phi_z, double x, double y, double z) {
  return Gamma*pow(M_PI,2)*(phi_x*pow(a_phix,2)*cos(M_PI*a_phix*x/L) + phi_y*pow(a_phiy,2)*cos(M_PI*a_phiy*y/L) + phi_z*pow(a_phiz,2)*sin(M_PI*a_phiz*z/L))/pow(L,2);
}

double Q_phi(double Q_phi_convection, double Q_phi_diffusion, double Q_phi_time) {
  return Q_phi_convection + Q_phi_diffusion + Q_phi_time;
}


