/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#ifndef PROJECT__.._C_CODES_NS_SUTHERLAND_SCALAR_TRANSIENT_3D_RHO__H
#define PROJECT__.._C_CODES_NS_SUTHERLAND_SCALAR_TRANSIENT_3D_RHO__H

double U(double L, double a_ut, double a_ux, double a_uy, double a_uz, double t, double u_0, double u_t, double u_x, double u_y, double u_z, double x, double y, double z);
double V(double L, double a_vt, double a_vx, double a_vy, double a_vz, double t, double v_0, double v_t, double v_x, double v_y, double v_z, double x, double y, double z);
double W(double L, double a_wt, double a_wx, double a_wy, double a_wz, double t, double w_0, double w_t, double w_x, double w_y, double w_z, double x, double y, double z);
double P(double L, double a_pt, double a_px, double a_py, double a_pz, double p_0, double p_t, double p_x, double p_y, double p_z, double t, double x, double y, double z);
double Q_rho_time(double Lt, double a_rhot, double rho_t, double t);
double Q_rho_convection(double L, double Rho, double U, double V, double W, double a_rhox, double a_rhoy, double a_rhoz, double a_ux, double a_vy, double a_wz, double rho_x, double rho_y, double rho_z, double u_x, double v_y, double w_z, double x, double y, double z);
double Q_rho(double Q_rho_convection, double Q_rho_time);

#endif

