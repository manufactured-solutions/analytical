/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#ifndef PROJECT__.._C_CODES_NS_SUTHERLAND_SCALAR_TRANSIENT_MANUF_SOLUTIONS_GRADIENTS__H
#define PROJECT__.._C_CODES_NS_SUTHERLAND_SCALAR_TRANSIENT_MANUF_SOLUTIONS_GRADIENTS__H

double u_an(double L, double Lt, double a_ut, double a_ux, double a_uy, double a_uz, double t, double u_0, double u_t, double u_x, double u_y, double u_z, double x, double y, double z);
double v_an(double L, double Lt, double a_vt, double a_vx, double a_vy, double a_vz, double t, double v_0, double v_t, double v_x, double v_y, double v_z, double x, double y, double z);
double w_an(double L, double Lt, double a_wt, double a_wx, double a_wy, double a_wz, double t, double w_0, double w_t, double w_x, double w_y, double w_z, double x, double y, double z);
double p_an(double L, double Lt, double a_pt, double a_px, double a_py, double a_pz, double p_0, double p_t, double p_x, double p_y, double p_z, double t, double x, double y, double z);
double rho_an(double L, double Lt, double a_rhot, double a_rhox, double a_rhoy, double a_rhoz, double rho_0, double rho_t, double rho_x, double rho_y, double rho_z, double t, double x, double y, double z);
double phi_an(double L, double Lt, double a_phit, double a_phix, double a_phiy, double a_phiz, double phi_0, double phi_t, double phi_x, double phi_y, double phi_z, double t, double x, double y, double z);
double du_dx(double L, double a_ux, double u_x, double x);
double du_dy(double L, double a_uy, double u_y, double y);
double du_dz(double L, double a_uz, double u_z, double z);
double dv_dx(double L, double a_vx, double v_x, double x);
double dv_dy(double L, double a_vy, double v_y, double y);
double dv_dz(double L, double a_vz, double v_z, double z);
double dw_dx(double L, double a_wx, double w_x, double x);
double dw_dy(double L, double a_wy, double w_y, double y);
double dw_dz(double L, double a_wz, double w_z, double z);
double dp_dx(double L, double a_px, double p_x, double x);
double dp_dy(double L, double a_py, double p_y, double y);
double dp_dz(double L, double a_pz, double p_z, double z);
double drho_dx(double L, double a_rhox, double rho_x, double x);
double drho_dy(double L, double a_rhoy, double rho_y, double y);
double drho_dz(double L, double a_rhoz, double rho_z, double z);
double dphi_dx(double L, double a_phix, double phi_x, double x);
double dphi_dy(double L, double a_phiy, double phi_y, double y);
double dphi_dz(double L, double a_phiz, double phi_z, double z);

#endif

