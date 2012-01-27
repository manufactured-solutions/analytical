/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#ifndef PROJECT__.._C_CODES_NS_SUTHERLAND_SCALAR_TRANSIENT_3D_E__H
#define PROJECT__.._C_CODES_NS_SUTHERLAND_SCALAR_TRANSIENT_3D_E__H

double Rho(double L, double Lt, double a_rhot, double a_rhox, double a_rhoy, double a_rhoz, double rho_0, double rho_t, double rho_x, double rho_y, double rho_z, double t, double x, double y, double z);
double U(double L, double Lt, double a_ut, double a_ux, double a_uy, double a_uz, double t, double u_0, double u_t, double u_x, double u_y, double u_z, double x, double y, double z);
double V(double L, double Lt, double a_vt, double a_vx, double a_vy, double a_vz, double t, double v_0, double v_t, double v_x, double v_y, double v_z, double x, double y, double z);
double W(double L, double Lt, double a_wt, double a_wx, double a_wy, double a_wz, double t, double w_0, double w_t, double w_x, double w_y, double w_z, double x, double y, double z);
double P(double L, double Lt, double a_pt, double a_px, double a_py, double a_pz, double p_0, double p_t, double p_x, double p_y, double p_z, double t, double x, double y, double z);
double T(double P, double R, double Rho);
double E_t(double P, double Rho, double U, double V, double W, double gamma);
double A_mu(double B_mu, double T_ref, double mu_ref);
double Mu(double A_mu, double B_mu, double T);
double DMu_Dx(double A_mu, double B_mu, double L, double R, double Rho, double T, double a_px, double a_rhox, double p_x, double rho_x, double x);
double DMu_Dy(double A_mu, double B_mu, double L, double R, double Rho, double T, double a_py, double a_rhoy, double p_y, double rho_y, double y);
double DMu_Dz(double A_mu, double B_mu, double L, double R, double Rho, double T, double a_pz, double a_rhoz, double p_z, double rho_z, double z);
double kappa(double Mu, double Pr, double R, double gamma);
double Q_et_convection(double E_t, double L, double P, double Rho, double U, double V, double W, double a_px, double a_py, double a_pz, double a_rhox, double a_rhoy, double a_rhoz, double a_ux, double a_uy, double a_uz, double a_vx, double a_vy, double a_vz, double a_wx, double a_wy, double a_wz, double gamma, double p_x, double p_y, double p_z, double rho_x, double rho_y, double rho_z, double u_x, double u_y, double u_z, double v_x, double v_y, double v_z, double w_x, double w_y, double w_z, double x, double y, double z);
double Q_et_gradp(double L, double P, double U, double V, double W, double a_px, double a_py, double a_pz, double a_ux, double a_vy, double a_wz, double p_x, double p_y, double p_z, double u_x, double v_y, double w_z, double x, double y, double z);
double Q_et_viscous(double DMu_Dx, double DMu_Dy, double DMu_Dz, double L, double Mu, double U, double V, double W, double a_ux, double a_uy, double a_uz, double a_vx, double a_vy, double a_vz, double a_wx, double a_wy, double a_wz, double u_x, double u_y, double u_z, double v_x, double v_y, double v_z, double w_x, double w_y, double w_z, double x, double y, double z);
double Q_et_heatflux(double DMu_Dx, double DMu_Dy, double DMu_Dz, double L, double P, double Pr, double R, double Rho, double a_px, double a_py, double a_pz, double a_rhox, double a_rhoy, double a_rhoz, double gamma, double kappa, double p_x, double p_y, double p_z, double rho_x, double rho_y, double rho_z, double x, double y, double z);
double Q_et_time(double E_t, double Lt, double P, double Rho, double U, double V, double W, double a_pt, double a_rhot, double a_ut, double a_vt, double a_wt, double gamma, double p_t, double rho_t, double t, double u_t, double v_t, double w_t);
double Q_et(double Q_et_convection, double Q_et_gradp, double Q_et_heatflux, double Q_et_time, double Q_et_viscous);

#endif

