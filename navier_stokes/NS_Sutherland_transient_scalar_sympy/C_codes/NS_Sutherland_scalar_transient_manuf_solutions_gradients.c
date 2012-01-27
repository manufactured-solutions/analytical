/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#include "NS_Sutherland_scalar_transient_manuf_solutions_gradients.h"
#include <math.h>

double u_an(double L, double Lt, double a_ut, double a_ux, double a_uy, double a_uz, double t, double u_0, double u_t, double u_x, double u_y, double u_z, double x, double y, double z) {
  return u_0 + u_t*cos(M_PI*a_ut*t/Lt) + u_x*sin(M_PI*a_ux*x/L) + u_y*cos(M_PI*a_uy*y/L) + u_z*cos(M_PI*a_uz*z/L);
}

double v_an(double L, double Lt, double a_vt, double a_vx, double a_vy, double a_vz, double t, double v_0, double v_t, double v_x, double v_y, double v_z, double x, double y, double z) {
  return v_0 + v_t*sin(M_PI*a_vt*t/Lt) + v_x*cos(M_PI*a_vx*x/L) + v_y*sin(M_PI*a_vy*y/L) + v_z*sin(M_PI*a_vz*z/L);
}

double w_an(double L, double Lt, double a_wt, double a_wx, double a_wy, double a_wz, double t, double w_0, double w_t, double w_x, double w_y, double w_z, double x, double y, double z) {
  return w_0 + w_t*cos(M_PI*a_wt*t/Lt) + w_x*sin(M_PI*a_wx*x/L) + w_y*sin(M_PI*a_wy*y/L) + w_z*cos(M_PI*a_wz*z/L);
}

double p_an(double L, double Lt, double a_pt, double a_px, double a_py, double a_pz, double p_0, double p_t, double p_x, double p_y, double p_z, double t, double x, double y, double z) {
  return p_0 + p_t*cos(M_PI*a_pt*t/Lt) + p_x*cos(M_PI*a_px*x/L) + p_y*sin(M_PI*a_py*y/L) + p_z*cos(M_PI*a_pz*z/L);
}

double rho_an(double L, double Lt, double a_rhot, double a_rhox, double a_rhoy, double a_rhoz, double rho_0, double rho_t, double rho_x, double rho_y, double rho_z, double t, double x, double y, double z) {
  return rho_0 + rho_t*sin(M_PI*a_rhot*t/Lt) + rho_x*sin(M_PI*a_rhox*x/L) + rho_y*cos(M_PI*a_rhoy*y/L) + rho_z*sin(M_PI*a_rhoz*z/L);
}

double phi_an(double L, double Lt, double a_phit, double a_phix, double a_phiy, double a_phiz, double phi_0, double phi_t, double phi_x, double phi_y, double phi_z, double t, double x, double y, double z) {
  return phi_0 + phi_t*cos(M_PI*a_phit*t/Lt) + phi_x*cos(M_PI*a_phix*x/L) + phi_y*cos(M_PI*a_phiy*y/L) + phi_z*sin(M_PI*a_phiz*z/L);
}

double du_dx(double L, double a_ux, double u_x, double x) {
  return M_PI*a_ux*u_x*cos(M_PI*a_ux*x/L)/L;
}

double du_dy(double L, double a_uy, double u_y, double y) {
  return -M_PI*a_uy*u_y*sin(M_PI*a_uy*y/L)/L;
}

double du_dz(double L, double a_uz, double u_z, double z) {
  return -M_PI*a_uz*u_z*sin(M_PI*a_uz*z/L)/L;
}

double dv_dx(double L, double a_vx, double v_x, double x) {
  return -M_PI*a_vx*v_x*sin(M_PI*a_vx*x/L)/L;
}

double dv_dy(double L, double a_vy, double v_y, double y) {
  return M_PI*a_vy*v_y*cos(M_PI*a_vy*y/L)/L;
}

double dv_dz(double L, double a_vz, double v_z, double z) {
  return M_PI*a_vz*v_z*cos(M_PI*a_vz*z/L)/L;
}

double dw_dx(double L, double a_wx, double w_x, double x) {
  return M_PI*a_wx*w_x*cos(M_PI*a_wx*x/L)/L;
}

double dw_dy(double L, double a_wy, double w_y, double y) {
  return M_PI*a_wy*w_y*cos(M_PI*a_wy*y/L)/L;
}

double dw_dz(double L, double a_wz, double w_z, double z) {
  return -M_PI*a_wz*w_z*sin(M_PI*a_wz*z/L)/L;
}

double dp_dx(double L, double a_px, double p_x, double x) {
  return -M_PI*a_px*p_x*sin(M_PI*a_px*x/L)/L;
}

double dp_dy(double L, double a_py, double p_y, double y) {
  return M_PI*a_py*p_y*cos(M_PI*a_py*y/L)/L;
}

double dp_dz(double L, double a_pz, double p_z, double z) {
  return -M_PI*a_pz*p_z*sin(M_PI*a_pz*z/L)/L;
}

double drho_dx(double L, double a_rhox, double rho_x, double x) {
  return M_PI*a_rhox*rho_x*cos(M_PI*a_rhox*x/L)/L;
}

double drho_dy(double L, double a_rhoy, double rho_y, double y) {
  return -M_PI*a_rhoy*rho_y*sin(M_PI*a_rhoy*y/L)/L;
}

double drho_dz(double L, double a_rhoz, double rho_z, double z) {
  return M_PI*a_rhoz*rho_z*cos(M_PI*a_rhoz*z/L)/L;
}

double dphi_dx(double L, double a_phix, double phi_x, double x) {
  return -M_PI*a_phix*phi_x*sin(M_PI*a_phix*x/L)/L;
}

double dphi_dy(double L, double a_phiy, double phi_y, double y) {
  return -M_PI*a_phiy*phi_y*sin(M_PI*a_phiy*y/L)/L;
}

double dphi_dz(double L, double a_phiz, double phi_z, double z) {
  return M_PI*a_phiz*phi_z*cos(M_PI*a_phiz*z/L)/L;
}


