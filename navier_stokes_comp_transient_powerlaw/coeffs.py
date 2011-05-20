#!/usr/bin/env python
# Provide a way to visualize various coefficient choices on the manufactured
# solution and associated forcing.

# Bring in symbolic solutions in phi, phi_t, phi_x, etc.
# Done prior to imports so some values will be overridden
execfile("soln.py")

from math import *
from enthought.traits.api import *
from enthought.traits.ui.api import View, Item, Group, Include, Tabbed

# Traits shorthands for various coefficient types
Amplitude = Float(0.0)
Frequency = Range(0.0, None)
Phase     = Range(0.0, 2*pi, exclude_high=True)
Length    = Range(0.0, None)

class Scenario(HasStrictTraits):
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
        self.lambda_r = self.mu_r / 3
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

        self._update_length('Lx', scenario.Lx)
        self._update_length('Ly', scenario.Ly)
        self._update_length('Lz', scenario.Lz)
        scenario.on_trait_change(self._update_length, 'Lx, Ly, Lz')

    def _update_length(self, name, new):
        self.__dict__[name] = new
        self.__dict__['twopi_inv' + name] = 2 * pi / new

    def _invoke_sympy_result(self, sympy_result, x, y, z, t):
        params = {'x': x, 'y': y, 'z': z, 't': t}
        params.update(self.__dict__)
        return sympy_result.subs(params).evalf(16,chop=True)

    def __call__(self, x, y, z, t):
        return self._invoke_sympy_result(phi, x, y, z, t)

    def _t(self, x, y, z, t):
        return self._invoke_sympy_result(phi_t, x, y, z, t)

    def _x(self, x, y, z, t):
        return self._invoke_sympy_result(phi_x, x, y, z, t)

    def _xx(self, x, y, z, t):
        return self._invoke_sympy_result(phi_xx, x, y, z, t)

    def _xy(self, x, y, z, t):
        return self._invoke_sympy_result(phi_xy, x, y, z, t)

    def _xz(self, x, y, z, t):
        return self._invoke_sympy_result(phi_xz, x, y, z, t)

    def _y(self, x, y, z, t):
        return self._invoke_sympy_result(phi_y, x, y, z, t)

    def _yy(self, x, y, z, t):
        return self._invoke_sympy_result(phi_yy, x, y, z, t)

    def _yz(self, x, y, z, t):
        return self._invoke_sympy_result(phi_yz, x, y, z, t)

    def _z(self, x, y, z, t):
        return self._invoke_sympy_result(phi_z, x, y, z, t)

    def _zz(self, x, y, z, t):
        return self._invoke_sympy_result(phi_zz, x, y, z, t)

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
    scenario    = Instance(Scenario, allow_none=False)
    rho         = Instance(PrimitiveSolution, allow_none=False)
    u           = Instance(PrimitiveSolution, allow_none=False)
    v           = Instance(PrimitiveSolution, allow_none=False)
    w           = Instance(PrimitiveSolution, allow_none=False)
    T           = Instance(PrimitiveSolution, allow_none=False)

    def __init__(self):
        self.scenario = Scenario()
        self.rho      = PrimitiveSolution(self.scenario, 'rho')
        self.u        = PrimitiveSolution(self.scenario, 'u'  )
        self.v        = PrimitiveSolution(self.scenario, 'v'  )
        self.w        = PrimitiveSolution(self.scenario, 'w'  )
        self.T        = PrimitiveSolution(self.scenario, 'T'  )

# TEST
ms = ManufacturedSolution()

