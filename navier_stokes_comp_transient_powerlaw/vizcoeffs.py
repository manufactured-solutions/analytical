#!/usr/bin/env python
"""
SYNOPSIS

    vizcoeff.py [-h,--help]

DESCRIPTION

    This script is provides helpful definitions for interactively
    visualizing manufactured quantities within a Mayavi
    (http://code.enthought.com/projects/mayavi/) session.

VERSION

    $Id$
"""

import cPickle as pickle
import sys, os, traceback, optparse, string, textwrap, numpy
import nsctpl
from math import pi, sqrt
from enthought.mayavi import mlab

# Create a manufactured solution to which we'll bind nearly everything
m = nsctpl.manufactured_solution()

# Used to programmatically generate ufuncs of the appropriate types
ufuncs_to_define = [
    ("rho"    , m.rho.__call__, 4, 1),
    ("rho_t"  , m.rho._t      , 4, 1),
    ("rho_x"  , m.rho._x      , 4, 1),
    ("rho_xx" , m.rho._xx     , 4, 1),
    ("rho_xy" , m.rho._xy     , 4, 1),
    ("rho_xz" , m.rho._xz     , 4, 1),
    ("rho_y"  , m.rho._y      , 4, 1),
    ("rho_yy" , m.rho._yy     , 4, 1),
    ("rho_yz" , m.rho._yz     , 4, 1),
    ("rho_z"  , m.rho._z      , 4, 1),
    ("rho_zz" , m.rho._zz     , 4, 1),

    ("u"    , m.u.__call__, 4, 1),
    ("u_t"  , m.u._t      , 4, 1),
    ("u_x"  , m.u._x      , 4, 1),
    ("u_xx" , m.u._xx     , 4, 1),
    ("u_xy" , m.u._xy     , 4, 1),
    ("u_xz" , m.u._xz     , 4, 1),
    ("u_y"  , m.u._y      , 4, 1),
    ("u_yy" , m.u._yy     , 4, 1),
    ("u_yz" , m.u._yz     , 4, 1),
    ("u_z"  , m.u._z      , 4, 1),
    ("u_zz" , m.u._zz     , 4, 1),

    ("v"    , m.v.__call__, 4, 1),
    ("v_t"  , m.v._t      , 4, 1),
    ("v_x"  , m.v._x      , 4, 1),
    ("v_xx" , m.v._xx     , 4, 1),
    ("v_xy" , m.v._xy     , 4, 1),
    ("v_xz" , m.v._xz     , 4, 1),
    ("v_y"  , m.v._y      , 4, 1),
    ("v_yy" , m.v._yy     , 4, 1),
    ("v_yz" , m.v._yz     , 4, 1),
    ("v_z"  , m.v._z      , 4, 1),
    ("v_zz" , m.v._zz     , 4, 1),

    ("w"    , m.w.__call__, 4, 1),
    ("w_t"  , m.w._t      , 4, 1),
    ("w_x"  , m.w._x      , 4, 1),
    ("w_xx" , m.w._xx     , 4, 1),
    ("w_xy" , m.w._xy     , 4, 1),
    ("w_xz" , m.w._xz     , 4, 1),
    ("w_y"  , m.w._y      , 4, 1),
    ("w_yy" , m.w._yy     , 4, 1),
    ("w_yz" , m.w._yz     , 4, 1),
    ("w_z"  , m.w._z      , 4, 1),
    ("w_zz" , m.w._zz     , 4, 1),

    ("T"    , m.T.__call__, 4, 1),
    ("T_t"  , m.T._t      , 4, 1),
    ("T_x"  , m.T._x      , 4, 1),
    ("T_xx" , m.T._xx     , 4, 1),
    ("T_xy" , m.T._xy     , 4, 1),
    ("T_xz" , m.T._xz     , 4, 1),
    ("T_y"  , m.T._y      , 4, 1),
    ("T_yy" , m.T._yy     , 4, 1),
    ("T_yz" , m.T._yz     , 4, 1),
    ("T_z"  , m.T._z      , 4, 1),
    ("T_zz" , m.T._zz     , 4, 1),

    ("grad_rho" , m.grad_rho, 5, 1),
    ("grad_u"   , m.grad_u  , 5, 1),
    ("grad_v"   , m.grad_v  , 5, 1),
    ("grad_w"   , m.grad_w  , 5, 1),
    ("grad_T"   , m.grad_T  , 5, 1),
    ("e"        , m.e       , 4, 1),
    ("p"        , m.p       , 4, 1),
    ("mu"       , m.mu      , 4, 1),
    ("rhou"     , m.rhou    , 4, 1),
    ("rhov"     , m.rhov    , 4, 1),
    ("rhow"     , m.rhow    , 4, 1),
    ("rhoe"     , m.rhoe    , 4, 1),
    ("grad_e"   , m.grad_e  , 5, 1),
    ("grad_p"   , m.grad_p  , 5, 1),
    ("grad_mu"  , m.grad_mu , 5, 1),
    ("Q_rho"    , m.Q_rho   , 4, 1),
    ("Q_rhou"   , m.Q_rhou  , 4, 1),
    ("Q_rhov"   , m.Q_rhov  , 4, 1),
    ("Q_rhow"   , m.Q_rhow  , 4, 1),
    ("Q_rhoe"   , m.Q_rhoe  , 4, 1)
]

def ufunc_generator(name, base):
    def f(*args):
        b = numpy.broadcast(*args)
        out = numpy.empty(b.shape, dtype=numpy.float)
        out.flat = [base(*(map(float,arg))) for arg in b]
        return out
    return f

# Programmatically set a module variable for each generated ufunc
# See http://stackoverflow.com/questions/1429814
for (name, base, nin, nout) in ufuncs_to_define:
    setattr(sys.modules[__name__], name, ufunc_generator(name, base))

def grid(Nx, Ny, Nz):
    """Create a numpy.mgrid using the solution domain"""
    return numpy.mgrid[
        -m.Lx/2 : m.Lx/2 : Nx * 1j,
              0 : m.Ly   : Ny * 1j,
        -m.Lz/2 : m.Lz/2 : Nz * 1j
    ]


# Establish scenario parameters
m.gamma    = 1.4
m.R        = 287
m.beta     = 2./3.
m.mu_r     = 185.2e-7
m.T_r      = 300
m.k_r      = m.gamma*m.R*m.mu_r/(m.gamma - 1)/0.70
m.lambda_r = - 2./3. *m.mu_r
m.Lx       = 4*pi
m.Ly       = 2
m.Lz       = 4*pi/3

# Populate primitive solution parameters for an isothermal channel

# Y is the wall-normal direction
# Want u, v, w, and T constant at the y = {0, Ly} walls
# Accomplish using phase offset pi / 2 since cos(pi/2) = 0
for var in [m.u, m.v, m.w, m.T]:
    for phase in ['c_x', 'c_xy', 'c_xz', 'c_y', 'c_yz', 'c_z',
                         'e_xy', 'e_xz',        'e_yz',       ]:
        setattr(var, phase, pi/2)

# Set wavenumbers for all solution variables spatial variations
for var in [m.rho, m.u, m.v, m.w, m.T]:
    var.b_x = 3
    var.b_y = 1
    var.b_z = 2
    var.b_xy = var.b_xz = var.b_yz = 1
    var.d_x = 1
    var.d_y = 1
    var.d_z = 1
    var.d_xy = var.d_xz = var.d_yz = 1
    var.f_x  = var.f_y  = var.f_z  = 1
    var.f_xy = var.f_xz = var.f_yz = 1

m.rho.a_0  = 1
m.rho.a_x  = m.rho.a_0 / 11
m.rho.a_y  = m.rho.a_0 / 23
m.rho.a_z  = m.rho.a_0 / 31
m.rho.a_xy = m.rho.a_0 / 61
m.rho.a_xz = m.rho.a_0 / 67
m.rho.a_yz = m.rho.a_0 / 71

m.T.a_0   = m.T_r
m.T.a_x   = m.T.a_0 / 13
m.T.a_y   = m.T.a_0 / 29
m.T.a_z   = m.T.a_0 / 37
m.T.a_xy  = m.T.a_0 / 101
m.T.a_xz  = m.T.a_0 / 103
m.T.a_yz  = m.T.a_0 / 107

m.u.a_0  = 0
m.u.a_y  = 1.25 * sqrt(m.gamma * m.R * m.T.a_0) # Supersonic
m.u.a_x  = m.u.a_y / 503
m.u.a_z  = m.u.a_y / 509
m.u.a_xy = m.u.a_y / 983
m.u.a_xz = m.u.a_y / 991
m.u.a_yz = m.u.a_y / 997

m.v.a_0  = 0
m.v.a_x  = 2
m.v.a_y  = 3
m.v.a_z  = 7
m.v.a_xx = 2
m.v.a_xy = 3
m.v.a_yz = 5

m.w.a_0  = 0
m.w.a_x  = 2
m.w.a_y  = 3
m.w.a_z  = 7
m.w.a_xx = 2
m.w.a_xy = 3
m.w.a_yz = 5

# Compute a default grid based on a coarsening of Coleman et al JFM 1995
Nx = (144 / 1.5)
Ny = (119 / 1.5)
Nz = (80  / 1.5)
x, y, z = grid(Nx, Ny, Nz)


def plotfield(ufunc):
    mlab.clf()
    d = ufunc(x, y, z, t)
    f = mlab.contour3d(x, y, z, d, transparent = True)

## Plot up density field as a starting point
t = 0
plotfield(rho)

