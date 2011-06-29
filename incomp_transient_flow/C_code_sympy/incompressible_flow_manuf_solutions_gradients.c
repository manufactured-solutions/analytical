/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#include "incompressible_flow_manuf_solutions_gradients.h"
#include <math.h>

double u_an(double a_u0, double a_ux, double a_uxy, double a_uxz, double a_uy, double a_uyz, double a_uz, double b_ux, double b_uxy, double b_uxz, double b_uy, double b_uyz, double b_uz, double c_ux, double c_uxy, double c_uxz, double c_uy, double c_uyz, double c_uz, double d_uxy, double d_uxz, double d_uyz, double e_uxy, double e_uxz, double e_uyz, double f_u0, double f_ux, double f_uxy, double f_uxz, double f_uy, double f_uyz, double f_uz, double g_u0, double g_ux, double g_uxy, double g_uxz, double g_uy, double g_uyz, double g_uz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return a_u0*cos(g_u0 + f_u0*t) + a_ux*cos(c_ux + b_ux*twopi_invLx*x)*cos(g_ux + f_ux*t) + a_uy*cos(g_uy + f_uy*t)*cos(c_uy + b_uy*twopi_invLy*y) + a_uz*cos(g_uz + f_uz*t)*cos(c_uz + b_uz*twopi_invLz*z) + a_uxy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t) + a_uxz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*cos(e_uxz + d_uxz*twopi_invLz*z) + a_uyz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z);
}

double v_an(double a_v0, double a_vx, double a_vxy, double a_vxz, double a_vy, double a_vyz, double a_vz, double b_vx, double b_vxy, double b_vxz, double b_vy, double b_vyz, double b_vz, double c_vx, double c_vxy, double c_vxz, double c_vy, double c_vyz, double c_vz, double d_vxy, double d_vxz, double d_vyz, double e_vxy, double e_vxz, double e_vyz, double f_v0, double f_vx, double f_vxy, double f_vxz, double f_vy, double f_vyz, double f_vz, double g_v0, double g_vx, double g_vxy, double g_vxz, double g_vy, double g_vyz, double g_vz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return a_v0*cos(g_v0 + f_v0*t) + a_vx*cos(c_vx + b_vx*twopi_invLx*x)*cos(g_vx + f_vx*t) + a_vy*cos(c_vy + b_vy*twopi_invLy*y)*cos(g_vy + f_vy*t) + a_vz*cos(g_vz + f_vz*t)*cos(c_vz + b_vz*twopi_invLz*z) + a_vxy*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x) + a_vxz*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x) + a_vyz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*cos(e_vyz + d_vyz*twopi_invLz*z);
}

double w_an(double a_w0, double a_wx, double a_wxy, double a_wxz, double a_wy, double a_wyz, double a_wz, double b_wx, double b_wxy, double b_wxz, double b_wy, double b_wyz, double b_wz, double c_wx, double c_wxy, double c_wxz, double c_wy, double c_wyz, double c_wz, double d_wxy, double d_wxz, double d_wyz, double e_wxy, double e_wxz, double e_wyz, double f_w0, double f_wx, double f_wxy, double f_wxz, double f_wy, double f_wyz, double f_wz, double g_w0, double g_wx, double g_wxy, double g_wxz, double g_wy, double g_wyz, double g_wz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return a_w0*cos(g_w0 + f_w0*t) + a_wx*cos(g_wx + f_wx*t)*cos(c_wx + b_wx*twopi_invLx*x) + a_wy*cos(c_wy + b_wy*twopi_invLy*y)*cos(g_wy + f_wy*t) + a_wz*cos(g_wz + f_wz*t)*cos(c_wz + b_wz*twopi_invLz*z) + a_wxy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*cos(e_wxy + d_wxy*twopi_invLy*y) + a_wxz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*cos(e_wxz + d_wxz*twopi_invLz*z) + a_wyz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*cos(e_wyz + d_wyz*twopi_invLz*z);
}

double du_dx(double a_ux, double a_uxy, double a_uxz, double b_ux, double b_uxy, double b_uxz, double c_ux, double c_uxy, double c_uxz, double d_uxy, double d_uxz, double e_uxy, double e_uxz, double f_ux, double f_uxy, double f_uxz, double g_ux, double g_uxy, double g_uxz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return -a_ux*b_ux*twopi_invLx*cos(g_ux + f_ux*t)*sin(c_ux + b_ux*twopi_invLx*x) - a_uxy*b_uxy*twopi_invLx*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t)*sin(c_uxy + b_uxy*twopi_invLx*x) - a_uxz*b_uxz*twopi_invLx*cos(g_uxz + f_uxz*t)*cos(e_uxz + d_uxz*twopi_invLz*z)*sin(c_uxz + b_uxz*twopi_invLx*x);
}

double du_dy(double a_uxy, double a_uy, double a_uyz, double b_uxy, double b_uy, double b_uyz, double c_uxy, double c_uy, double c_uyz, double d_uxy, double d_uyz, double e_uxy, double e_uyz, double f_uxy, double f_uy, double f_uyz, double g_uxy, double g_uy, double g_uyz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return -a_uy*b_uy*twopi_invLy*cos(g_uy + f_uy*t)*sin(c_uy + b_uy*twopi_invLy*y) - a_uxy*d_uxy*twopi_invLy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(g_uxy + f_uxy*t)*sin(e_uxy + d_uxy*twopi_invLy*y) - a_uyz*b_uyz*twopi_invLy*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z)*sin(c_uyz + b_uyz*twopi_invLy*y);
}

double du_dz(double a_uxz, double a_uyz, double a_uz, double b_uxz, double b_uyz, double b_uz, double c_uxz, double c_uyz, double c_uz, double d_uxz, double d_uyz, double e_uxz, double e_uyz, double f_uxz, double f_uyz, double f_uz, double g_uxz, double g_uyz, double g_uz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return -a_uz*b_uz*twopi_invLz*cos(g_uz + f_uz*t)*sin(c_uz + b_uz*twopi_invLz*z) - a_uxz*d_uxz*twopi_invLz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*sin(e_uxz + d_uxz*twopi_invLz*z) - a_uyz*d_uyz*twopi_invLz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*sin(e_uyz + d_uyz*twopi_invLz*z);
}

double dv_dx(double a_vx, double a_vxy, double a_vxz, double b_vx, double b_vxy, double b_vxz, double c_vx, double c_vxy, double c_vxz, double d_vxy, double d_vxz, double e_vxy, double e_vxz, double f_vx, double f_vxy, double f_vxz, double g_vx, double g_vxy, double g_vxz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return -a_vx*b_vx*twopi_invLx*cos(g_vx + f_vx*t)*sin(c_vx + b_vx*twopi_invLx*x) - a_vxy*b_vxy*twopi_invLx*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*sin(c_vxy + b_vxy*twopi_invLx*x) - a_vxz*b_vxz*twopi_invLx*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*sin(c_vxz + b_vxz*twopi_invLx*x);
}

double dv_dy(double a_vxy, double a_vy, double a_vyz, double b_vxy, double b_vy, double b_vyz, double c_vxy, double c_vy, double c_vyz, double d_vxy, double d_vyz, double e_vxy, double e_vyz, double f_vxy, double f_vy, double f_vyz, double g_vxy, double g_vy, double g_vyz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return -a_vy*b_vy*twopi_invLy*cos(g_vy + f_vy*t)*sin(c_vy + b_vy*twopi_invLy*y) - a_vxy*d_vxy*twopi_invLy*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x)*sin(e_vxy + d_vxy*twopi_invLy*y) - a_vyz*b_vyz*twopi_invLy*cos(g_vyz + f_vyz*t)*cos(e_vyz + d_vyz*twopi_invLz*z)*sin(c_vyz + b_vyz*twopi_invLy*y);
}

double dv_dz(double a_vxz, double a_vyz, double a_vz, double b_vxz, double b_vyz, double b_vz, double c_vxz, double c_vyz, double c_vz, double d_vxz, double d_vyz, double e_vxz, double e_vyz, double f_vxz, double f_vyz, double f_vz, double g_vxz, double g_vyz, double g_vz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return -a_vz*b_vz*twopi_invLz*cos(g_vz + f_vz*t)*sin(c_vz + b_vz*twopi_invLz*z) - a_vxz*d_vxz*twopi_invLz*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x)*sin(e_vxz + d_vxz*twopi_invLz*z) - a_vyz*d_vyz*twopi_invLz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*sin(e_vyz + d_vyz*twopi_invLz*z);
}

double dw_dx(double a_wx, double a_wxy, double a_wxz, double b_wx, double b_wxy, double b_wxz, double c_wx, double c_wxy, double c_wxz, double d_wxy, double d_wxz, double e_wxy, double e_wxz, double f_wx, double f_wxy, double f_wxz, double g_wx, double g_wxy, double g_wxz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return -a_wx*b_wx*twopi_invLx*cos(g_wx + f_wx*t)*sin(c_wx + b_wx*twopi_invLx*x) - a_wxy*b_wxy*twopi_invLx*cos(g_wxy + f_wxy*t)*cos(e_wxy + d_wxy*twopi_invLy*y)*sin(c_wxy + b_wxy*twopi_invLx*x) - a_wxz*b_wxz*twopi_invLx*cos(g_wxz + f_wxz*t)*cos(e_wxz + d_wxz*twopi_invLz*z)*sin(c_wxz + b_wxz*twopi_invLx*x);
}

double dw_dy(double a_wxy, double a_wy, double a_wyz, double b_wxy, double b_wy, double b_wyz, double c_wxy, double c_wy, double c_wyz, double d_wxy, double d_wyz, double e_wxy, double e_wyz, double f_wxy, double f_wy, double f_wyz, double g_wxy, double g_wy, double g_wyz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return -a_wy*b_wy*twopi_invLy*cos(g_wy + f_wy*t)*sin(c_wy + b_wy*twopi_invLy*y) - a_wxy*d_wxy*twopi_invLy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*sin(e_wxy + d_wxy*twopi_invLy*y) - a_wyz*b_wyz*twopi_invLy*cos(g_wyz + f_wyz*t)*cos(e_wyz + d_wyz*twopi_invLz*z)*sin(c_wyz + b_wyz*twopi_invLy*y);
}

double dw_dz(double a_wxz, double a_wyz, double a_wz, double b_wxz, double b_wyz, double b_wz, double c_wxz, double c_wyz, double c_wz, double d_wxz, double d_wyz, double e_wxz, double e_wyz, double f_wxz, double f_wyz, double f_wz, double g_wxz, double g_wyz, double g_wz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return -a_wz*b_wz*twopi_invLz*cos(g_wz + f_wz*t)*sin(c_wz + b_wz*twopi_invLz*z) - a_wxz*d_wxz*twopi_invLz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*sin(e_wxz + d_wxz*twopi_invLz*z) - a_wyz*d_wyz*twopi_invLz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*sin(e_wyz + d_wyz*twopi_invLz*z);
}


