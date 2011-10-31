/****************************************************************************** 
 *                      Code generated with sympy 0.6.7                       * 
 *                                                                            * 
 *              See http://www.sympy.org/ for more information.               * 
 *                                                                            * 
 *                       This file is part of 'project'                       * 
 ******************************************************************************/

#include "incompressible_flow_ource_Qcontinuity.h"
#include <math.h>

double Q_continuity(double a_ux, double a_uxy, double a_uxz, double a_vxy, double a_vy, double a_vyz, double a_wxz, double a_wyz, double a_wz, double b_ux, double b_uxy, double b_uxz, double b_vxy, double b_vy, double b_vyz, double b_wxz, double b_wyz, double b_wz, double c_ux, double c_uxy, double c_uxz, double c_vxy, double c_vy, double c_vyz, double c_wxz, double c_wyz, double c_wz, double d_uxy, double d_uxz, double d_vxy, double d_vyz, double d_wxz, double d_wyz, double e_uxy, double e_uxz, double e_vxy, double e_vyz, double e_wxz, double e_wyz, double f_ux, double f_uxy, double f_uxz, double f_vxy, double f_vy, double f_vyz, double f_wxz, double f_wyz, double f_wz, double g_ux, double g_uxy, double g_uxz, double g_vxy, double g_vy, double g_vyz, double g_wxz, double g_wyz, double g_wz, double t, double twopi_invLx, double twopi_invLy, double twopi_invLz, double x, double y, double z) {
  return -a_ux*b_ux*twopi_invLx*cos(g_ux + f_ux*t)*sin(c_ux + b_ux*twopi_invLx*x) - a_vy*b_vy*twopi_invLy*cos(g_vy + f_vy*t)*sin(c_vy + b_vy*twopi_invLy*y) - a_wz*b_wz*twopi_invLz*cos(g_wz + f_wz*t)*sin(c_wz + b_wz*twopi_invLz*z) - a_uxy*b_uxy*twopi_invLx*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t)*sin(c_uxy + b_uxy*twopi_invLx*x) - a_uxz*b_uxz*twopi_invLx*cos(g_uxz + f_uxz*t)*cos(e_uxz + d_uxz*twopi_invLz*z)*sin(c_uxz + b_uxz*twopi_invLx*x) - a_vxy*d_vxy*twopi_invLy*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x)*sin(e_vxy + d_vxy*twopi_invLy*y) - a_vyz*b_vyz*twopi_invLy*cos(g_vyz + f_vyz*t)*cos(e_vyz + d_vyz*twopi_invLz*z)*sin(c_vyz + b_vyz*twopi_invLy*y) - a_wxz*d_wxz*twopi_invLz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*sin(e_wxz + d_wxz*twopi_invLz*z) - a_wyz*d_wyz*twopi_invLz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*sin(e_wyz + d_wyz*twopi_invLz*z);
}


