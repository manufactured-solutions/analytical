#!/usr/bin/env python
"""
DESCRIPTION

    This script is provides helpful definitions for interactively visualizing
    manufactured quantities within a Mayavi
    (http://code.enthought.com/projects/mayavi/) or IPython session.  It is
    /not/ a standalone program and should be %run from within an 'python
    -wthread -pylab' session.

VERSION

    $Id$
"""

import cPickle as pickle
import sys, os, traceback, optparse, string, textwrap, numpy
import nsctpl
from math import pi, sqrt
from enthought.mayavi import mlab
from matplotlib import pyplot

# Create a manufactured solution to which we'll bind nearly everything
m = nsctpl.manufactured_solution()

# Establish isothermal channel scenario parameters
nsctpl.isothermal_channel(m)

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

# Extra mucking around necessary for returned ndarray to have float type
def ufunc_generator(name, base, nin, nout):
    ufunc = numpy.frompyfunc(base, nin, nout)
    def f(*args):
        b   = numpy.broadcast(*args)
        out = numpy.empty(b.shape, dtype=numpy.float)
        args = list(args)
        args.append(out)
        args = tuple(args)
        ufunc(*args)
        return out
    return f

# Programmatically set a module variable for each generated ufunc
# See http://stackoverflow.com/questions/1429814
for (name, base, nin, nout) in ufuncs_to_define:
    setattr(sys.modules[__name__], name, ufunc_generator(name, base, nin, nout))

def grid(Nx, Ny, Nz):
    """Create a numpy.mgrid using the solution domain"""
    return numpy.mgrid[
        -m.Lx/2 : m.Lx/2 : Nx * 1j,
              0 : m.Ly   : Ny * 1j,
        -m.Lz/2 : m.Lz/2 : Nz * 1j
    ]

# Compute a default grid based on a coarsened Coleman et al JFM 1995
Nx = int(144 / 2)
Ny = int(119 / 2)
Nz = int(80  / 2)
x, y, z = grid(Nx, Ny, Nz)
t = 0
tfinal = 1 / 10

# Plot a particular field's isocontours
def plotfield(ufunc):
    d = ufunc(x, y, z, t)
    print "field: min, max, mean, std = %g, %g, %g, %g" % (
            d.min(), d.max(), d.mean(), d.std() )
    print "lwall: min, max, mean, std = %g, %g, %g, %g" % (
            d[:,0,:].min(), d[:,0,:].max(), d[:,0,:].mean(), d[:,0,:].std() )
    print "uwall: min, max, mean, std = %g, %g, %g, %g" % (
            d[:,-1,:].min(), d[:,-1,:].max(), d[:,-1,:].mean(), d[:,-1,:].std() )
    mlab.clf()
    f = mlab.contour3d(x, y, z, d, transparent = True)
    return (d, f)

def fieldmin(ufunc, Nt = 100):
    times = numpy.linspace(0, tfinal, Nt)
    xlims = (times.min(), times.max())
    mins  = numpy.zeros_like(times)

    fig = pyplot.figure(num=None)
    pyplot.clf()
    for i in range(len(times)):
        d = ufunc(x, y, z, times[i])
        mins[i]  = d.min()

        if i < 2:
            continue

        pyplot.figure(fig.number)
        pyplot.hold(True);
        pyplot.plot(times[0:i-1], mins[0:i-1])
        pyplot.xlim(xlims)
        pyplot.draw()

    return (times, mins)

## Now, for example, one can plot up density contours using
# d, f = plotfield(rho)
## Or, for example, look at how minimum density changes with time using
# times, mins = fieldmin(rho)
