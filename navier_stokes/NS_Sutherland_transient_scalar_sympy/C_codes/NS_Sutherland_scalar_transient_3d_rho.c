/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#include "NS_Sutherland_scalar_transient_3d_rho.h"
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

double P(double L, double Lt, double a_pt, double a_px, double a_py, double a_pz, double p_0, double p_t, double p_x, double p_y, double p_z, double t, double x, double y, double z) {
  return p_0 + p_t*cos(M_PI*a_pt*t/Lt) + p_x*cos(M_PI*a_px*x/L) + p_y*sin(M_PI*a_py*y/L) + p_z*cos(M_PI*a_pz*z/L);
}

double Q_rho_time(double Lt, double a_rhot, double rho_t, double t) {
  return M_PI*a_rhot*rho_t*cos(M_PI*a_rhot*t/Lt)/Lt;
}

double Q_rho_convection(double L, double Rho, double U, double V, double W, double a_rhox, double a_rhoy, double a_rhoz, double a_ux, double a_vy, double a_wz, double rho_x, double rho_y, double rho_z, double u_x, double v_y, double w_z, double x, double y, double z) {
  return M_PI*Rho*(a_ux*u_x*cos(M_PI*a_ux*x/L) + a_vy*v_y*cos(M_PI*a_vy*y/L) - a_wz*w_z*sin(M_PI*a_wz*z/L))/L + M_PI*U*a_rhox*rho_x*cos(M_PI*a_rhox*x/L)/L + M_PI*W*a_rhoz*rho_z*cos(M_PI*a_rhoz*z/L)/L - M_PI*V*a_rhoy*rho_y*sin(M_PI*a_rhoy*y/L)/L;
}

double Q_rho(double Q_rho_convection, double Q_rho_time) {
  return Q_rho_convection + Q_rho_time;
}


