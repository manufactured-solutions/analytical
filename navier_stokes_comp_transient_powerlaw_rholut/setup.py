#!/usr/bin/env python

"""
setup.py file for nsctpl_rholut module (navier_stokes_comp_transient_powerlaw_rholut)

"""

from distutils.core import setup, Extension


nsctpl_rholut_module = Extension('_nsctpl_rholut',
                                 sources=['nsctpl_rholut_wrap.cxx'],
                                )

setup (name = 'nsctpl_rholut',
       version     = '0.1',
       author      = "Rhys Ulerich <rhys@ices.utexas.edu>",
       description = """Manufactured solution for navier_stokes_comp_transient_powerlaw""",
       ext_modules = [nsctpl_rholut_module],
       py_modules  = ["nsctpl_rholut"],
       )
