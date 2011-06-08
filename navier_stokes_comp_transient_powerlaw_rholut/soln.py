#!/usr/bin/env python
# Using sympy, construct a desired analytical solution for some quantity \phi.
# Compute the analytical derivatives of \phi and output them as C code.  This
# logic will be used to compute the manufactured solutions rho, u, v, w, and T.

from sympy import *

# Coordinates
var('x y z t', real=True)

# Solution parameters used in the form of the analytical solution
var(""" a_0  a_x  a_xy  a_xz  a_y  a_yz  a_z
             b_x  b_xy  b_xz  b_y  b_yz  b_z
             c_x  c_xy  c_xz  c_y  c_yz  c_z
                  d_xy  d_xz       d_yz
                  e_xy  e_xz       e_yz
        f_0  f_x  f_xy  f_xz  f_y  f_yz  f_z
        g_0  g_x  g_xy  g_xz  g_y  g_yz  g_z """, real=True)

# Explicitly keep (2 * pi / L) terms together as indivisible tokens
var('twopi_invLx  twopi_invLy  twopi_invLz', real=True)

# Form the analytical solution and its derivatives
phi  = (
        a_0                                                                 *cos(f_0 *t + g_0 )
      + a_x  * cos(b_x *twopi_invLx*x + c_x )                               *cos(f_x *t + g_x )
      + a_xy * cos(b_xy*twopi_invLx*x + c_xy)*cos(d_xy*twopi_invLy*y + e_xy)*cos(f_xy*t + g_xy)
      + a_xz * cos(b_xz*twopi_invLx*x + c_xz)*cos(d_xz*twopi_invLz*z + e_xz)*cos(f_xz*t + g_xz)
      + a_y  * cos(b_y *twopi_invLy*y + c_y )                               *cos(f_y *t + g_y )
      + a_yz * cos(b_yz*twopi_invLy*y + c_yz)*cos(d_yz*twopi_invLz*z + e_yz)*cos(f_yz*t + g_yz)
      + a_z  * cos(b_z *twopi_invLz*z + c_z )                               *cos(f_z *t + g_z )
)
phi_t  = phi.diff(t)
phi_x  = phi.diff(x)
phi_y  = phi.diff(y)
phi_z  = phi.diff(z)
phi_xx = phi_x.diff(x)
phi_xy = phi_x.diff(y)
phi_xz = phi_x.diff(z)
phi_yy = phi_y.diff(y)
phi_yz = phi_y.diff(z)
phi_zz = phi_z.diff(z)

# Save these results to files 'soln.c' and 'soln.h' as C code
from sympy.utilities.codegen import codegen
codegen((
            ("phi",    phi   ),
            ("phi_t",  phi_t ),
            ("phi_x",  phi_x ),
            ("phi_xx", phi_xx),
            ("phi_xy", phi_xy),
            ("phi_xz", phi_xz),
            ("phi_y",  phi_y ),
            ("phi_yy", phi_yy),
            ("phi_yz", phi_yz),
            ("phi_z",  phi_z ),
            ("phi_zz", phi_zz),
        ),
        "C", "soln", header=False, to_files=True)
