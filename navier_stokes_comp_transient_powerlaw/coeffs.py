#!/usr/bin/env python
# Provide a way to visualize various coefficient choices on the manufactured
# solution and associated forcing.  Note this is **not** strongly tied to
# soln.py or forcing.py and **will** need to change when those files change.

from math import pi, cos, sin
from enthought.traits.api import *
from enthought.traits.ui.api import Group, Include, Item, Tabbed, View
import numpy

# Traits shorthands for various coefficient types
Amplitude = Float(0.0)
Frequency = Range(0.0, None)
Phase     = Range(0.0, 2*pi, exclude_high=True)
Length    = Range(0.0, None)
Array3    = Array(numpy.float, shape=(None,None,None))

class Scenario(HasTraits):
    gamma    = Float(desc="Ratio of specific heats")
    R        = Float(desc="Gas constant (J/kg/K)")
    beta     = Float(desc="Viscosity power law exponent")
    mu_r     = Float(desc="Reference dynamic viscosity (Pa s)")
    T_r      = Float(desc="Reference temperature (K)")
    k_r      = Float(desc="Reference thermal conductivity (W/m/K)")
    lambda_r = Float(desc="Reference second viscosity (Pa s)")
    Lx       = Length(desc="Domain length in x direction (m)")
    Ly       = Length(desc="Domain length in y direction (m)")
    Lz       = Length(desc="Domain length in z direction (m)")

    def __init__(self):
        self.gamma    = 1.4
        self.R        = 287
        self.beta     = 2./3.
        self.mu_r     = 185.2e-7
        self.T_r      = 300
        self.k_r      = self.gamma*self.R*self.mu_r/(self.gamma - 1)/0.70
        self.lambda_r = - 2./3. *self.mu_r
        self.Lx       = 4*pi
        self.Ly       = 2
        self.Lz       = 4*pi/3

    traits_view = View(Group(Item(name = 'gamma'),
                             Item(name = 'R'),
                             Item(name = 'beta'),
                             Item(name = 'mu_r'),
                             Item(name = 'T_r'),
                             Item(name = 'k_r'),
                             Item(name = 'lambda_r'),
                             Item(name = 'Lx'),
                             Item(name = 'Ly'),
                             Item(name = 'Lz'),
                             label='Scenario parameters',
                             show_border = True))


class PrimitiveSolution(HasTraits):

    def __init__(self, scenario, name = ""):

        # Coefficients appearing in each primitive solution
        # from which we construct traits with instance-specific labels
        for pre, post, trait in [
            ("a_", "0" , Amplitude),
            ("f_", "0" , Frequency),
            ("g_", "0" , Phase),
            ("a_", "x" , Amplitude),
            ("b_", "x" , Frequency),
            ("c_", "x" , Phase),
            ("f_", "x" , Frequency),
            ("g_", "x" , Phase),
            ("a_", "xy", Amplitude),
            ("b_", "xy", Frequency),
            ("c_", "xy", Phase),
            ("d_", "xy", Frequency),
            ("e_", "xy", Phase),
            ("f_", "xy", Frequency),
            ("g_", "xy", Phase),
            ("a_", "xz", Amplitude),
            ("b_", "xz", Frequency),
            ("c_", "xz", Phase),
            ("d_", "xz", Frequency),
            ("e_", "xz", Phase),
            ("g_", "xz", Phase),
            ("f_", "xz", Frequency),
            ("a_", "y" , Amplitude),
            ("b_", "y" , Frequency),
            ("c_", "y" , Phase),
            ("f_", "y" , Frequency),
            ("g_", "y" , Phase),
            ("a_", "yz", Amplitude),
            ("b_", "yz", Frequency),
            ("c_", "yz", Phase),
            ("d_", "yz", Frequency),
            ("e_", "yz", Phase),
            ("f_", "yz", Frequency),
            ("g_", "yz", Phase),
            ("a_", "z" , Amplitude),
            ("b_", "z" , Frequency),
            ("c_", "z" , Phase),
            ("f_", "z" , Frequency),
            ("g_", "z" , Phase)
        ]:
            self.add_trait(pre + post, trait(label=pre + name + post))
            self.__dict__[pre + post] = self.trait(pre + post).default

        # Maintain live 2*pi/L variables for scenario.[Lx, Ly, Lz]
        self.update_length('Lx', scenario.Lx)
        self.update_length('Ly', scenario.Ly)
        self.update_length('Lz', scenario.Lz)
        scenario.on_trait_change(self.update_length, 'Lx, Ly, Lz')

        # Construct ufuncs based on scalar-only evaluation routines
        self.__ufunc = numpy.frompyfunc(self.__eval, 4, 1)
        self._t      = numpy.frompyfunc(self.__t   , 4, 1)
        self._x      = numpy.frompyfunc(self.__x   , 4, 1)
        self._xx     = numpy.frompyfunc(self.__xx  , 4, 1)
        self._xy     = numpy.frompyfunc(self.__xy  , 4, 1)
        self._xz     = numpy.frompyfunc(self.__xz  , 4, 1)
        self._y      = numpy.frompyfunc(self.__y   , 4, 1)
        self._yy     = numpy.frompyfunc(self.__yy  , 4, 1)
        self._yz     = numpy.frompyfunc(self.__yz  , 4, 1)
        self._z      = numpy.frompyfunc(self.__z   , 4, 1)
        self._zz     = numpy.frompyfunc(self.__zz  , 4, 1)

    def update_length(self, name, new):
        self.__dict__['twopi_inv' + name] = 2 * pi / new

    def __call__(self, x, y, z, t, out=None):
        """Return soln(x,y,z,t) given current parameters"""
        return self.__ufunc(x, y, z, t, out)

    def __eval(self, x, y, z, t):
        """Return soln(x,y,z,t) given current parameters"""
        return self.a_0*cos(self.g_0 + self.f_0*t) + self.a_x*cos(self.c_x + self.b_x*self.twopi_invLx*x)*cos(self.g_x + self.f_x*t) + self.a_y*cos(self.g_y + self.f_y*t)*cos(self.c_y + self.b_y*self.twopi_invLy*y) + self.a_z*cos(self.g_z + self.f_z*t)*cos(self.c_z + self.b_z*self.twopi_invLz*z) + self.a_xy*cos(self.g_xy + self.f_xy*t)*cos(self.c_xy + self.b_xy*self.twopi_invLx*x)*cos(self.e_xy + self.d_xy*self.twopi_invLy*y) + self.a_xz*cos(self.g_xz + self.f_xz*t)*cos(self.c_xz + self.b_xz*self.twopi_invLx*x)*cos(self.e_xz + self.d_xz*self.twopi_invLz*z) + self.a_yz*cos(self.c_yz + self.b_yz*self.twopi_invLy*y)*cos(self.e_yz + self.d_yz*self.twopi_invLz*z)*cos(self.g_yz + self.f_yz*t)

    def __t(self, x, y, z, t):
        """Return \partial_t soln(x,y,z,t) given current parameters"""
        return -self.a_0*self.f_0*sin(self.g_0 + self.f_0*t) - self.a_x*self.f_x*cos(self.c_x + self.b_x*self.twopi_invLx*x)*sin(self.g_x + self.f_x*t) - self.a_y*self.f_y*cos(self.c_y + self.b_y*self.twopi_invLy*y)*sin(self.g_y + self.f_y*t) - self.a_z*self.f_z*cos(self.c_z + self.b_z*self.twopi_invLz*z)*sin(self.g_z + self.f_z*t) - self.a_xy*self.f_xy*cos(self.c_xy + self.b_xy*self.twopi_invLx*x)*cos(self.e_xy + self.d_xy*self.twopi_invLy*y)*sin(self.g_xy + self.f_xy*t) - self.a_xz*self.f_xz*cos(self.c_xz + self.b_xz*self.twopi_invLx*x)*cos(self.e_xz + self.d_xz*self.twopi_invLz*z)*sin(self.g_xz + self.f_xz*t) - self.a_yz*self.f_yz*cos(self.c_yz + self.b_yz*self.twopi_invLy*y)*cos(self.e_yz + self.d_yz*self.twopi_invLz*z)*sin(self.g_yz + self.f_yz*t)

    def __x(self, x, y, z, t):
        """Return \partial_x soln(x,y,z,t) given current parameters"""
        return -self.a_x*self.b_x*self.twopi_invLx*cos(self.g_x + self.f_x*t)*sin(self.c_x + self.b_x*self.twopi_invLx*x) - self.a_xy*self.b_xy*self.twopi_invLx*cos(self.g_xy + self.f_xy*t)*cos(self.e_xy + self.d_xy*self.twopi_invLy*y)*sin(self.c_xy + self.b_xy*self.twopi_invLx*x) - self.a_xz*self.b_xz*self.twopi_invLx*cos(self.g_xz + self.f_xz*t)*cos(self.e_xz + self.d_xz*self.twopi_invLz*z)*sin(self.c_xz + self.b_xz*self.twopi_invLx*x)

    def __xx(self, x, y, z, t):
        """Return \partial_{xx} soln(x,y,z,t) given current parameters"""
        return -self.a_x*self.b_x**2*self.twopi_invLx**2*cos(self.c_x + self.b_x*self.twopi_invLx*x)*cos(self.g_x + self.f_x*t) - self.a_xy*self.b_xy**2*self.twopi_invLx**2*cos(self.g_xy + self.f_xy*t)*cos(self.c_xy + self.b_xy*self.twopi_invLx*x)*cos(self.e_xy + self.d_xy*self.twopi_invLy*y) - self.a_xz*self.b_xz**2*self.twopi_invLx**2*cos(self.g_xz + self.f_xz*t)*cos(self.c_xz + self.b_xz*self.twopi_invLx*x)*cos(self.e_xz + self.d_xz*self.twopi_invLz*z)

    def __xy(self, x, y, z, t):
        """Return \partial_{xy} soln(x,y,z,t) given current parameters"""
        return self.a_xy*self.b_xy*self.d_xy*self.twopi_invLx*self.twopi_invLy*cos(self.g_xy + self.f_xy*t)*sin(self.c_xy + self.b_xy*self.twopi_invLx*x)*sin(self.e_xy + self.d_xy*self.twopi_invLy*y)

    def __xz(self, x, y, z, t):
        """Return \partial_{xz} soln(x,y,z,t) given current parameters"""
        return self.a_xz*self.b_xz*self.d_xz*self.twopi_invLx*self.twopi_invLz*cos(self.g_xz + self.f_xz*t)*sin(self.c_xz + self.b_xz*self.twopi_invLx*x)*sin(self.e_xz + self.d_xz*self.twopi_invLz*z)

    def __y(self, x, y, z, t):
        """Return \partial_y soln(x,y,z,t) given current parameters"""
        return -self.a_y*self.b_y*self.twopi_invLy*cos(self.g_y + self.f_y*t)*sin(self.c_y + self.b_y*self.twopi_invLy*y) - self.a_xy*self.d_xy*self.twopi_invLy*cos(self.g_xy + self.f_xy*t)*cos(self.c_xy + self.b_xy*self.twopi_invLx*x)*sin(self.e_xy + self.d_xy*self.twopi_invLy*y) - self.a_yz*self.b_yz*self.twopi_invLy*cos(self.e_yz + self.d_yz*self.twopi_invLz*z)*cos(self.g_yz + self.f_yz*t)*sin(self.c_yz + self.b_yz*self.twopi_invLy*y)

    def __yy(self, x, y, z, t):
        """Return \partial_{yy} soln(x,y,z,t) given current parameters"""
        return -self.a_y*self.b_y**2*self.twopi_invLy**2*cos(self.g_y + self.f_y*t)*cos(self.c_y + self.b_y*self.twopi_invLy*y) - self.a_xy*self.d_xy**2*self.twopi_invLy**2*cos(self.g_xy + self.f_xy*t)*cos(self.c_xy + self.b_xy*self.twopi_invLx*x)*cos(self.e_xy + self.d_xy*self.twopi_invLy*y) - self.a_yz*self.b_yz**2*self.twopi_invLy**2*cos(self.c_yz + self.b_yz*self.twopi_invLy*y)*cos(self.e_yz + self.d_yz*self.twopi_invLz*z)*cos(self.g_yz + self.f_yz*t)

    def __yz(self, x, y, z, t):
        """Return \partial_{yz} soln(x,y,z,t) given current parameters"""
        return self.a_yz*self.b_yz*self.d_yz*self.twopi_invLy*self.twopi_invLz*cos(self.g_yz + self.f_yz*t)*sin(self.c_yz + self.b_yz*self.twopi_invLy*y)*sin(self.e_yz + self.d_yz*self.twopi_invLz*z)

    def __z(self, x, y, z, t):
        """Return \partial_z soln(x,y,z,t) given current parameters"""
        return -self.a_z*self.b_z*self.twopi_invLz*cos(self.g_z + self.f_z*t)*sin(self.c_z + self.b_z*self.twopi_invLz*z) - self.a_xz*self.d_xz*self.twopi_invLz*cos(self.g_xz + self.f_xz*t)*cos(self.c_xz + self.b_xz*self.twopi_invLx*x)*sin(self.e_xz + self.d_xz*self.twopi_invLz*z) - self.a_yz*self.d_yz*self.twopi_invLz*cos(self.c_yz + self.b_yz*self.twopi_invLy*y)*cos(self.g_yz + self.f_yz*t)*sin(self.e_yz + self.d_yz*self.twopi_invLz*z)

    def __zz(self, x, y, z, t):
        """Return \partial_{zz} soln(x,y,z,t) given current parameters"""
        -self.a_z*self.b_z**2*self.twopi_invLz**2*cos(self.g_z + self.f_z*t)*cos(self.c_z + self.b_z*self.twopi_invLz*z) - self.a_xz*self.d_xz**2*self.twopi_invLz**2*cos(self.g_xz + self.f_xz*t)*cos(self.c_xz + self.b_xz*self.twopi_invLx*x)*cos(self.e_xz + self.d_xz*self.twopi_invLz*z) - self.a_yz*self.d_yz**2*self.twopi_invLz**2*cos(self.c_yz + self.b_yz*self.twopi_invLy*y)*cos(self.e_yz + self.d_yz*self.twopi_invLz*z)*cos(self.g_yz + self.f_yz*t)

    traits_view = View(Group(Group(Group(Item(name = 'a_0'),
                                         Item(name = 'f_0'),
                                         Item(name = 'g_0'),
                                         show_border = True)),
                             Group(Group(Item(name = 'a_x'),
                                         Item(name = 'b_x'),
                                         Item(name = 'c_x'),
                                         Item(name = 'f_x'),
                                         Item(name = 'g_x'),
                                         show_border = True),
                                   Group(Item(name = 'a_xy'),
                                         Item(name = 'b_xy'),
                                         Item(name = 'c_xy'),
                                         Item(name = 'd_xy'),
                                         Item(name = 'e_xy'),
                                         Item(name = 'f_xy'),
                                         Item(name = 'g_xy'),
                                         show_border = True),
                                   Group(Item(name = 'a_xz'),
                                         Item(name = 'b_xz'),
                                         Item(name = 'c_xz'),
                                         Item(name = 'd_xz'),
                                         Item(name = 'e_xz'),
                                         Item(name = 'f_xz'),
                                         Item(name = 'g_xz'),
                                         show_border = True),
                                   show_border = False),
                             Group(Group(Item(name = 'a_y'),
                                         Item(name = 'b_y'),
                                         Item(name = 'c_y'),
                                         Item(name = 'f_y'),
                                         Item(name = 'g_y'),
                                         show_border = True),
                                   Group(Item(name = 'a_yz'),
                                         Item(name = 'b_yz'),
                                         Item(name = 'c_yz'),
                                         Item(name = 'd_yz'),
                                         Item(name = 'e_yz'),
                                         Item(name = 'f_yz'),
                                         Item(name = 'g_yz'),
                                         show_border = True),
                                   show_border = False),
                             Group(Group(Item(name = 'a_z'),
                                         Item(name = 'b_z'),
                                         Item(name = 'c_z'),
                                         Item(name = 'f_z'),
                                         Item(name = 'g_z'),
                                         show_border = True),
                                   show_border = False),
                             orientation = 'horizontal'))


class ManufacturedSolution(HasTraits):

    # Solution Parameters
    scenario    = Instance(Scenario)
    soln_rho    = Instance(PrimitiveSolution)
    soln_u      = Instance(PrimitiveSolution)
    soln_v      = Instance(PrimitiveSolution)
    soln_w      = Instance(PrimitiveSolution)
    soln_T      = Instance(PrimitiveSolution)

    # Manufactured solution computation details
    t = Float

    Nx, Ny, Nz = Int,    Int,    Int
    x,  y,  z  = Array3, Array3, Array3

    rho   , u   , v   , w   , T    = Array3, Array3, Array3, Array3, Array3
    rho_t , u_t , v_t , w_t , T_t  = Array3, Array3, Array3, Array3, Array3
    rho_x , u_x , v_x , w_x , T_x  = Array3, Array3, Array3, Array3, Array3
    rho_xx, u_xx, v_xx, w_xx, T_xx = Array3, Array3, Array3, Array3, Array3
    rho_xy, u_xy, v_xy, w_xy, T_xy = Array3, Array3, Array3, Array3, Array3
    rho_xz, u_xz, v_xz, w_xz, T_xz = Array3, Array3, Array3, Array3, Array3
    rho_y , u_y , v_y , w_y , T_y  = Array3, Array3, Array3, Array3, Array3
    rho_yy, u_yy, v_yy, w_yy, T_yy = Array3, Array3, Array3, Array3, Array3
    rho_yz, u_yz, v_yz, w_yz, T_yz = Array3, Array3, Array3, Array3, Array3
    rho_z , u_z , v_z , w_z , T_z  = Array3, Array3, Array3, Array3, Array3
    rho_zz, u_zz, v_zz, w_zz, T_zz = Array3, Array3, Array3, Array3, Array3

    def __init__(self):

        self.on_trait_change(self.grid_change, 'Nx,Ny,Nz,scenario.[Lx,Ly,Lz]')

        self.scenario = Scenario()
        self.soln_rho = PrimitiveSolution(self.scenario, 'rho')
        self.soln_u   = PrimitiveSolution(self.scenario, 'u'  )
        self.soln_v   = PrimitiveSolution(self.scenario, 'v'  )
        self.soln_w   = PrimitiveSolution(self.scenario, 'w'  )
        self.soln_T   = PrimitiveSolution(self.scenario, 'T'  )

        self.Nx = 4
        self.Ny = 4
        self.Nz = 4

    def grid_change(self):

        # Update the mesh grid
        self.x, self.y, self.z = numpy.mgrid[
                    -self.scenario.Lx/2 : self.scenario.Lx/2 : self.Nx * 1j,
                                      0 : self.scenario.Ly   : self.Ny * 1j,
                    -self.scenario.Lz/2 : self.scenario.Lz/2 : self.Nz * 1j
                ]

    def rho_change(self):
        self.rho = numpy.empty_like(self.x)
        self.soln_rho(self.x, self.y, self.z, self.t, out=self.rho)

# TEST
ms = ManufacturedSolution()
