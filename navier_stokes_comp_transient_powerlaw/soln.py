#!/usr/bin/env python
# Using sympy, construct a desired analytical solution for some quantity \phi.
# Compute the analytical derivatives of \phi and output them as C code.  This
# logic will be used to compute the manufactured solutions rho, u, v, w, and T.

from __future__ import division
import string
from sympy import *
from sympy.utilities.codegen import codegen

# Options used for all sympy.symbols calls
kwargs = { "real": True, "each_char": False }

# Coordinate definitions
x, y, z, t = symbols('x y z t', **kwargs)

# Solution parameters used in the form of the analytical solution
a_0, a_x, a_xy, a_xz, a_y, a_yz, a_z = symbols(string.join((
            'a_0', 'a_x', 'a_xy', 'a_xz', 'a_y', 'a_yz', 'a_z'
        )), **kwargs)
b_x, b_xy, b_xz, b_y, b_yz, b_z = symbols(string.join((
            'b_x', 'b_xy', 'b_xz', 'b_y', 'b_yz', 'b_z'
        )), **kwargs)
c_x, c_xy, c_xz, c_y, c_yz, c_z = symbols(string.join((
            'c_x', 'c_xy', 'c_xz', 'c_y', 'c_yz', 'c_z'
        )), **kwargs)
d_xy, d_xz, d_yz = symbols(string.join((
            'd_xy', 'd_xz', 'd_yz'
        )), **kwargs)
e_xy, e_xz, e_yz = symbols(string.join((
            'e_xy', 'e_xz', 'e_yz'
        )), **kwargs)
f_0, f_x, f_xy, f_xz, f_y, f_yz, f_z = symbols(string.join((
            'f_0', 'f_x', 'f_xy', 'f_xz', 'f_y', 'f_yz', 'f_z'
        )), **kwargs)
g_0, g_x, g_xy, g_xz, g_y, g_yz, g_z = symbols(string.join((
            'g_0', 'g_x', 'g_xy', 'g_xz', 'g_y', 'g_yz', 'g_z'
        )), **kwargs)

# Form the analytical solution and its derivatives
phi  = a_0                                         *sin(f_0 *t + g_0 ) \
     + a_x  * sin(b_x *x + c_x )                   *sin(f_x *t + g_x ) \
     + a_xy * sin(b_xy*x + c_xy)*sin(d_xy*y + e_xy)*sin(f_xy*t + g_xy) \
     + a_xz * sin(b_xz*x + c_xz)*sin(d_xz*z + e_xz)*sin(f_xz*t + g_xz) \
     + a_y  * sin(b_y *y + c_y )                   *sin(f_y *t + g_y ) \
     + a_yz * sin(b_yz*y + c_yz)*sin(d_yz*z + e_yz)*sin(f_yz*t + g_yz) \
     + a_z  * sin(b_z *z + c_z )                   *sin(f_z *t + g_z )
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

# Dump these results to files 'soln.c' and 'soln.h' as C code
codegen((                       \
            ("phi",    phi),    \
            ("phi_t",  phi_t),  \
            ("phi_x",  phi_x),  \
            ("phi_y",  phi_y),  \
            ("phi_z",  phi_z),  \
            ("phi_xx", phi_xx), \
            ("phi_xy", phi_xy), \
            ("phi_xz", phi_xz), \
            ("phi_yy", phi_yy), \
            ("phi_yz", phi_yz), \
            ("phi_zz", phi_zz), \
        ),                      \
        "C",                    \
        "soln",                 \
        header=False,           \
        to_files=True)

