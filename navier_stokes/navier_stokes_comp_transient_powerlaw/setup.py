#!/usr/bin/env python

"""
setup.py file for nsctpl module (navier_stokes_comp_transient_powerlaw)
"""

from distutils.core import setup, Extension


nsctpl_module = Extension('_nsctpl',
                           sources=['nsctpl_wrap.cxx'],
                           )

setup (name = 'nsctpl',
       version = '0.1',
       author      = "Rhys Ulerich <rhys@ices.utexas.edu>",
       description = """Manufactured solution for navier_stokes_comp_transient_powerlaw""",
       ext_modules = [nsctpl_module],
       py_modules = ["nsctpl"],
       )
