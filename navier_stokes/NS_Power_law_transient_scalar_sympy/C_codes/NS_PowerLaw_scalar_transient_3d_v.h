/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#ifndef PROJECT__.._C_CODES_NS_POWERLAW_SCALAR_TRANSIENT_3D_V__H
#define PROJECT__.._C_CODES_NS_POWERLAW_SCALAR_TRANSIENT_3D_V__H

double Rho(double L, double Lt, double a_rhot, double a_rhox, double a_rhoy, double a_rhoz, double rho_0, double rho_t, double rho_x, double rho_y, double rho_z, double t, double x, double y, double z);
double U(double L, double Lt, double a_ut, double a_ux, double a_uy, double a_uz, double t, double u_0, double u_t, double u_x, double u_y, double u_z, double x, double y, double z);
double V(double L, double Lt, double a_vt, double a_vx, double a_vy, double a_vz, double t, double v_0, double v_t, double v_x, double v_y, double v_z, double x, double y, double z);
double W(double L, double Lt, double a_wt, double a_wx, double a_wy, double a_wz, double t, double w_0, double w_t, double w_x, double w_y, double w_z, double x, double y, double z);
double P(double L, double Lt, double a_pt, double a_px, double a_py, double a_pz, double p_0, double p_t, double p_x, double p_y, double p_z, double t, double x, double y, double z);
double T(double P, double R, double Rho);
double Mu(double T, double T_ref, double beta, double mu_ref);
double DMu_Dx(double L, double R, double Rho, double T, double T_ref, double a_px, double a_rhox, double beta, double mu_ref, double p_x, double rho_x, double x);
double DMu_Dy(double L, double R, double Rho, double T, double T_ref, double a_py, double a_rhoy, double beta, double mu_ref, double p_y, double rho_y, double y);
double DMu_Dz(double L, double R, double Rho, double T, double T_ref, double a_pz, double a_rhoz, double beta, double mu_ref, double p_z, double rho_z, double z);
double Q_v_convection(double L, double Rho, double U, double V, double W, double a_rhox, double a_rhoy, double a_rhoz, double a_ux, double a_vx, double a_vy, double a_vz, double a_wz, double rho_x, double rho_y, double rho_z, double u_x, double v_x, double v_y, double v_z, double w_z, double x, double y, double z);
double Q_v_gradp(double L, double a_py, double p_y, double y);
double Q_v_viscous(double DMu_Dx, double DMu_Dy, double DMu_Dz, double L, double LAMBDA, double Mu, double a_ux, double a_uy, double a_vx, double a_vy, double a_vz, double a_wy, double a_wz, double alpha, double u_x, double u_y, double v_x, double v_y, double v_z, double w_y, double w_z, double x, double y, double z);
double Q_v_time(double Lt, double Rho, double V, double a_rhot, double a_vt, double rho_t, double t, double v_t);
double Q_v(double Q_v_convection, double Q_v_gradp, double Q_v_time, double Q_v_viscous);

#endif

