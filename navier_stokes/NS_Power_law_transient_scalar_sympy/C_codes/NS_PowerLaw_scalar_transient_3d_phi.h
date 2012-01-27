/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#ifndef PROJECT__.._C_CODES_NS_POWERLAW_SCALAR_TRANSIENT_3D_PHI__H
#define PROJECT__.._C_CODES_NS_POWERLAW_SCALAR_TRANSIENT_3D_PHI__H

double U(double L, double Lt, double a_ut, double a_ux, double a_uy, double a_uz, double t, double u_0, double u_t, double u_x, double u_y, double u_z, double x, double y, double z);
double V(double L, double Lt, double a_vt, double a_vx, double a_vy, double a_vz, double t, double v_0, double v_t, double v_x, double v_y, double v_z, double x, double y, double z);
double W(double L, double Lt, double a_wt, double a_wx, double a_wy, double a_wz, double t, double w_0, double w_t, double w_x, double w_y, double w_z, double x, double y, double z);
double Phi(double L, double Lt, double a_phit, double a_phix, double a_phiy, double a_phiz, double phi_0, double phi_t, double phi_x, double phi_y, double phi_z, double t, double x, double y, double z);
double Q_phi_time(double Lt, double Phi, double Rho, double a_phit, double a_rhot, double phi_t, double rho_t, double t);
double Q_phi_convection(double L, double Phi, double Rho, double U, double V, double W, double a_phix, double a_phiy, double a_phiz, double a_rhox, double a_rhoy, double a_rhoz, double a_ux, double a_vy, double a_wz, double phi_x, double phi_y, double phi_z, double rho_x, double rho_y, double rho_z, double u_x, double v_y, double w_z, double x, double y, double z);
double Q_phi_diffusion(double Gamma, double L, double a_phix, double a_phiy, double a_phiz, double phi_x, double phi_y, double phi_z, double x, double y, double z);
double Q_phi(double Q_phi_convection, double Q_phi_diffusion, double Q_phi_time);

#endif

