#!/usr/bin/env python
# Compute manufactured solution and manufactured forcing at a particular point.
# Used as a test case to ensure C++ implementation consistent with Python, and
# to ensure floating point precision losses are acceptable.

# soln.py contains the nececssary sympy imports
execfile("soln.py")

# Work in (possibly slow) extended precision
PREC = 40;

# Use prime numbers for test data.  Slick generator for primes taken from
# http://stackoverflow.com/questions/567222/simple-prime-generator-in-python
def prime_generator():
    D = {}
    q = 2
    while True:
        if q not in D:
            yield q
            D[q * q] = [q]
        else:
            for p in D[q]:
                D.setdefault(p + q, []).append(p)
            del D[q]
        q += 1

# Test solution parameters.  These are not physically
# realizable but can be used to check the implementation.
p = prime_generator()
params = {
    'alpha'    : Real(p.next(), PREC),
    'beta'     : Real(p.next(), PREC),
    'gamma'    : Real(p.next(), PREC),
    'Ma'       : Real(p.next(), PREC),
    'Pr'       : Real(p.next(), PREC),
    'Re'       : Real(p.next(), PREC),
    'a_rho0'   : Real(p.next(), PREC),
    'a_rhox'   : Real(p.next(), PREC),
    'a_rhoxy'  : Real(p.next(), PREC),
    'a_rhoxz'  : Real(p.next(), PREC),
    'a_rhoy'   : Real(p.next(), PREC),
    'a_rhoyz'  : Real(p.next(), PREC),
    'a_rhoz'   : Real(p.next(), PREC),
    'b_rhox'   : Real(p.next(), PREC),
    'b_rhoxy'  : Real(p.next(), PREC),
    'b_rhoxz'  : Real(p.next(), PREC),
    'b_rhoy'   : Real(p.next(), PREC),
    'b_rhoyz'  : Real(p.next(), PREC),
    'b_rhoz'   : Real(p.next(), PREC),
    'c_rhox'   : Real(p.next(), PREC),
    'c_rhoxy'  : Real(p.next(), PREC),
    'c_rhoxz'  : Real(p.next(), PREC),
    'c_rhoy'   : Real(p.next(), PREC),
    'c_rhoyz'  : Real(p.next(), PREC),
    'c_rhoz'   : Real(p.next(), PREC),
    'd_rhoxy'  : Real(p.next(), PREC),
    'd_rhoxz'  : Real(p.next(), PREC),
    'd_rhoyz'  : Real(p.next(), PREC),
    'e_rhoxy'  : Real(p.next(), PREC),
    'e_rhoxz'  : Real(p.next(), PREC),
    'e_rhoyz'  : Real(p.next(), PREC),
    'f_rho0'   : Real(p.next(), PREC),
    'f_rhox'   : Real(p.next(), PREC),
    'f_rhoxy'  : Real(p.next(), PREC),
    'f_rhoxz'  : Real(p.next(), PREC),
    'f_rhoy'   : Real(p.next(), PREC),
    'f_rhoyz'  : Real(p.next(), PREC),
    'f_rhoz'   : Real(p.next(), PREC),
    'g_rho0'   : Real(p.next(), PREC),
    'g_rhox'   : Real(p.next(), PREC),
    'g_rhoxy'  : Real(p.next(), PREC),
    'g_rhoxz'  : Real(p.next(), PREC),
    'g_rhoy'   : Real(p.next(), PREC),
    'g_rhoyz'  : Real(p.next(), PREC),
    'g_rhoz'   : Real(p.next(), PREC),
    'a_u0'     : Real(p.next(), PREC),
    'a_ux'     : Real(p.next(), PREC),
    'a_uxy'    : Real(p.next(), PREC),
    'a_uxz'    : Real(p.next(), PREC),
    'a_uy'     : Real(p.next(), PREC),
    'a_uyz'    : Real(p.next(), PREC),
    'a_uz'     : Real(p.next(), PREC),
    'b_ux'     : Real(p.next(), PREC),
    'b_uxy'    : Real(p.next(), PREC),
    'b_uxz'    : Real(p.next(), PREC),
    'b_uy'     : Real(p.next(), PREC),
    'b_uyz'    : Real(p.next(), PREC),
    'b_uz'     : Real(p.next(), PREC),
    'c_ux'     : Real(p.next(), PREC),
    'c_uxy'    : Real(p.next(), PREC),
    'c_uxz'    : Real(p.next(), PREC),
    'c_uy'     : Real(p.next(), PREC),
    'c_uyz'    : Real(p.next(), PREC),
    'c_uz'     : Real(p.next(), PREC),
    'd_uxy'    : Real(p.next(), PREC),
    'd_uxz'    : Real(p.next(), PREC),
    'd_uyz'    : Real(p.next(), PREC),
    'e_uxy'    : Real(p.next(), PREC),
    'e_uxz'    : Real(p.next(), PREC),
    'e_uyz'    : Real(p.next(), PREC),
    'f_u0'     : Real(p.next(), PREC),
    'f_ux'     : Real(p.next(), PREC),
    'f_uxy'    : Real(p.next(), PREC),
    'f_uxz'    : Real(p.next(), PREC),
    'f_uy'     : Real(p.next(), PREC),
    'f_uyz'    : Real(p.next(), PREC),
    'f_uz'     : Real(p.next(), PREC),
    'g_u0'     : Real(p.next(), PREC),
    'g_ux'     : Real(p.next(), PREC),
    'g_uxy'    : Real(p.next(), PREC),
    'g_uxz'    : Real(p.next(), PREC),
    'g_uy'     : Real(p.next(), PREC),
    'g_uyz'    : Real(p.next(), PREC),
    'g_uz'     : Real(p.next(), PREC),
    'a_v0'     : Real(p.next(), PREC),
    'a_vx'     : Real(p.next(), PREC),
    'a_vxy'    : Real(p.next(), PREC),
    'a_vxz'    : Real(p.next(), PREC),
    'a_vy'     : Real(p.next(), PREC),
    'a_vyz'    : Real(p.next(), PREC),
    'a_vz'     : Real(p.next(), PREC),
    'b_vx'     : Real(p.next(), PREC),
    'b_vxy'    : Real(p.next(), PREC),
    'b_vxz'    : Real(p.next(), PREC),
    'b_vy'     : Real(p.next(), PREC),
    'b_vyz'    : Real(p.next(), PREC),
    'b_vz'     : Real(p.next(), PREC),
    'c_vx'     : Real(p.next(), PREC),
    'c_vxy'    : Real(p.next(), PREC),
    'c_vxz'    : Real(p.next(), PREC),
    'c_vy'     : Real(p.next(), PREC),
    'c_vyz'    : Real(p.next(), PREC),
    'c_vz'     : Real(p.next(), PREC),
    'd_vxy'    : Real(p.next(), PREC),
    'd_vxz'    : Real(p.next(), PREC),
    'd_vyz'    : Real(p.next(), PREC),
    'e_vxy'    : Real(p.next(), PREC),
    'e_vxz'    : Real(p.next(), PREC),
    'e_vyz'    : Real(p.next(), PREC),
    'f_v0'     : Real(p.next(), PREC),
    'f_vx'     : Real(p.next(), PREC),
    'f_vxy'    : Real(p.next(), PREC),
    'f_vxz'    : Real(p.next(), PREC),
    'f_vy'     : Real(p.next(), PREC),
    'f_vyz'    : Real(p.next(), PREC),
    'f_vz'     : Real(p.next(), PREC),
    'g_v0'     : Real(p.next(), PREC),
    'g_vx'     : Real(p.next(), PREC),
    'g_vxy'    : Real(p.next(), PREC),
    'g_vxz'    : Real(p.next(), PREC),
    'g_vy'     : Real(p.next(), PREC),
    'g_vyz'    : Real(p.next(), PREC),
    'g_vz'     : Real(p.next(), PREC),
    'a_w0'     : Real(p.next(), PREC),
    'a_wx'     : Real(p.next(), PREC),
    'a_wxy'    : Real(p.next(), PREC),
    'a_wxz'    : Real(p.next(), PREC),
    'a_wy'     : Real(p.next(), PREC),
    'a_wyz'    : Real(p.next(), PREC),
    'a_wz'     : Real(p.next(), PREC),
    'b_wx'     : Real(p.next(), PREC),
    'b_wxy'    : Real(p.next(), PREC),
    'b_wxz'    : Real(p.next(), PREC),
    'b_wy'     : Real(p.next(), PREC),
    'b_wyz'    : Real(p.next(), PREC),
    'b_wz'     : Real(p.next(), PREC),
    'c_wx'     : Real(p.next(), PREC),
    'c_wxy'    : Real(p.next(), PREC),
    'c_wxz'    : Real(p.next(), PREC),
    'c_wy'     : Real(p.next(), PREC),
    'c_wyz'    : Real(p.next(), PREC),
    'c_wz'     : Real(p.next(), PREC),
    'd_wxy'    : Real(p.next(), PREC),
    'd_wxz'    : Real(p.next(), PREC),
    'd_wyz'    : Real(p.next(), PREC),
    'e_wxy'    : Real(p.next(), PREC),
    'e_wxz'    : Real(p.next(), PREC),
    'e_wyz'    : Real(p.next(), PREC),
    'f_w0'     : Real(p.next(), PREC),
    'f_wx'     : Real(p.next(), PREC),
    'f_wxy'    : Real(p.next(), PREC),
    'f_wxz'    : Real(p.next(), PREC),
    'f_wy'     : Real(p.next(), PREC),
    'f_wyz'    : Real(p.next(), PREC),
    'f_wz'     : Real(p.next(), PREC),
    'g_w0'     : Real(p.next(), PREC),
    'g_wx'     : Real(p.next(), PREC),
    'g_wxy'    : Real(p.next(), PREC),
    'g_wxz'    : Real(p.next(), PREC),
    'g_wy'     : Real(p.next(), PREC),
    'g_wyz'    : Real(p.next(), PREC),
    'g_wz'     : Real(p.next(), PREC),
    'a_T0'     : Real(p.next(), PREC),
    'a_Tx'     : Real(p.next(), PREC),
    'a_Txy'    : Real(p.next(), PREC),
    'a_Txz'    : Real(p.next(), PREC),
    'a_Ty'     : Real(p.next(), PREC),
    'a_Tyz'    : Real(p.next(), PREC),
    'a_Tz'     : Real(p.next(), PREC),
    'b_Tx'     : Real(p.next(), PREC),
    'b_Txy'    : Real(p.next(), PREC),
    'b_Txz'    : Real(p.next(), PREC),
    'b_Ty'     : Real(p.next(), PREC),
    'b_Tyz'    : Real(p.next(), PREC),
    'b_Tz'     : Real(p.next(), PREC),
    'c_Tx'     : Real(p.next(), PREC),
    'c_Txy'    : Real(p.next(), PREC),
    'c_Txz'    : Real(p.next(), PREC),
    'c_Ty'     : Real(p.next(), PREC),
    'c_Tyz'    : Real(p.next(), PREC),
    'c_Tz'     : Real(p.next(), PREC),
    'd_Txy'    : Real(p.next(), PREC),
    'd_Txz'    : Real(p.next(), PREC),
    'd_Tyz'    : Real(p.next(), PREC),
    'e_Txy'    : Real(p.next(), PREC),
    'e_Txz'    : Real(p.next(), PREC),
    'e_Tyz'    : Real(p.next(), PREC),
    'f_T0'     : Real(p.next(), PREC),
    'f_Tx'     : Real(p.next(), PREC),
    'f_Txy'    : Real(p.next(), PREC),
    'f_Txz'    : Real(p.next(), PREC),
    'f_Ty'     : Real(p.next(), PREC),
    'f_Tyz'    : Real(p.next(), PREC),
    'f_Tz'     : Real(p.next(), PREC),
    'g_T0'     : Real(p.next(), PREC),
    'g_Tx'     : Real(p.next(), PREC),
    'g_Txy'    : Real(p.next(), PREC),
    'g_Txz'    : Real(p.next(), PREC),
    'g_Ty'     : Real(p.next(), PREC),
    'g_Tyz'    : Real(p.next(), PREC),
    'g_Tz'     : Real(p.next(), PREC),
    'Lx'       : Real(p.next(), PREC),
    'Ly'       : Real(p.next(), PREC),
    'Lz'       : Real(p.next(), PREC)
}

# Load test parameters as global names for evaluation purposes
globals().update(params)

# Prepare substitution dictionaries for each analytic solution.  Necessary
# since our analytic solution form is repeated for rho, u, v, w, and T.
params_rho = {
    'a_0'  : a_rho0,
    'a_x'  : a_rhox,
    'a_xy' : a_rhoxy,
    'a_xz' : a_rhoxz,
    'a_y'  : a_rhoy,
    'a_yz' : a_rhoyz,
    'a_z'  : a_rhoz,
    'b_x'  : b_rhox,
    'b_xy' : b_rhoxy,
    'b_xz' : b_rhoxz,
    'b_y'  : b_rhoy,
    'b_yz' : b_rhoyz,
    'b_z'  : b_rhoz,
    'c_x'  : c_rhox,
    'c_xy' : c_rhoxy,
    'c_xz' : c_rhoxz,
    'c_y'  : c_rhoy,
    'c_yz' : c_rhoyz,
    'c_z'  : c_rhoz,
    'd_xy' : d_rhoxy,
    'd_xz' : d_rhoxz,
    'd_yz' : d_rhoyz,
    'e_xy' : e_rhoxy,
    'e_xz' : e_rhoxz,
    'e_yz' : e_rhoyz,
    'f_0'  : f_rho0,
    'f_x'  : f_rhox,
    'f_xy' : f_rhoxy,
    'f_xz' : f_rhoxz,
    'f_y'  : f_rhoy,
    'f_yz' : f_rhoyz,
    'f_z'  : f_rhoz,
    'g_0'  : g_rho0,
    'g_x'  : g_rhox,
    'g_xy' : g_rhoxy,
    'g_xz' : g_rhoxz,
    'g_y'  : g_rhoy,
    'g_yz' : g_rhoyz,
    'g_z'  : g_rhoz,
}
params_u = {
    'a_0'  : a_u0,
    'a_x'  : a_ux,
    'a_xy' : a_uxy,
    'a_xz' : a_uxz,
    'a_y'  : a_uy,
    'a_yz' : a_uyz,
    'a_z'  : a_uz,
    'b_x'  : b_ux,
    'b_xy' : b_uxy,
    'b_xz' : b_uxz,
    'b_y'  : b_uy,
    'b_yz' : b_uyz,
    'b_z'  : b_uz,
    'c_x'  : c_ux,
    'c_xy' : c_uxy,
    'c_xz' : c_uxz,
    'c_y'  : c_uy,
    'c_yz' : c_uyz,
    'c_z'  : c_uz,
    'd_xy' : d_uxy,
    'd_xz' : d_uxz,
    'd_yz' : d_uyz,
    'e_xy' : e_uxy,
    'e_xz' : e_uxz,
    'e_yz' : e_uyz,
    'f_0'  : f_u0,
    'f_x'  : f_ux,
    'f_xy' : f_uxy,
    'f_xz' : f_uxz,
    'f_y'  : f_uy,
    'f_yz' : f_uyz,
    'f_z'  : f_uz,
    'g_0'  : g_u0,
    'g_x'  : g_ux,
    'g_xy' : g_uxy,
    'g_xz' : g_uxz,
    'g_y'  : g_uy,
    'g_yz' : g_uyz,
    'g_z'  : g_uz,
}
params_v = {
    'a_0'  : a_v0,
    'a_x'  : a_vx,
    'a_xy' : a_vxy,
    'a_xz' : a_vxz,
    'a_y'  : a_vy,
    'a_yz' : a_vyz,
    'a_z'  : a_vz,
    'b_x'  : b_vx,
    'b_xy' : b_vxy,
    'b_xz' : b_vxz,
    'b_y'  : b_vy,
    'b_yz' : b_vyz,
    'b_z'  : b_vz,
    'c_x'  : c_vx,
    'c_xy' : c_vxy,
    'c_xz' : c_vxz,
    'c_y'  : c_vy,
    'c_yz' : c_vyz,
    'c_z'  : c_vz,
    'd_xy' : d_vxy,
    'd_xz' : d_vxz,
    'd_yz' : d_vyz,
    'e_xy' : e_vxy,
    'e_xz' : e_vxz,
    'e_yz' : e_vyz,
    'f_0'  : f_v0,
    'f_x'  : f_vx,
    'f_xy' : f_vxy,
    'f_xz' : f_vxz,
    'f_y'  : f_vy,
    'f_yz' : f_vyz,
    'f_z'  : f_vz,
    'g_0'  : g_v0,
    'g_x'  : g_vx,
    'g_xy' : g_vxy,
    'g_xz' : g_vxz,
    'g_y'  : g_vy,
    'g_yz' : g_vyz,
    'g_z'  : g_vz,
}
params_w = {
    'a_0'  : a_w0,
    'a_x'  : a_wx,
    'a_xy' : a_wxy,
    'a_xz' : a_wxz,
    'a_y'  : a_wy,
    'a_yz' : a_wyz,
    'a_z'  : a_wz,
    'b_x'  : b_wx,
    'b_xy' : b_wxy,
    'b_xz' : b_wxz,
    'b_y'  : b_wy,
    'b_yz' : b_wyz,
    'b_z'  : b_wz,
    'c_x'  : c_wx,
    'c_xy' : c_wxy,
    'c_xz' : c_wxz,
    'c_y'  : c_wy,
    'c_yz' : c_wyz,
    'c_z'  : c_wz,
    'd_xy' : d_wxy,
    'd_xz' : d_wxz,
    'd_yz' : d_wyz,
    'e_xy' : e_wxy,
    'e_xz' : e_wxz,
    'e_yz' : e_wyz,
    'f_0'  : f_w0,
    'f_x'  : f_wx,
    'f_xy' : f_wxy,
    'f_xz' : f_wxz,
    'f_y'  : f_wy,
    'f_yz' : f_wyz,
    'f_z'  : f_wz,
    'g_0'  : g_w0,
    'g_x'  : g_wx,
    'g_xy' : g_wxy,
    'g_xz' : g_wxz,
    'g_y'  : g_wy,
    'g_yz' : g_wyz,
    'g_z'  : g_wz,
}
params_T = {
    'a_0'  : a_T0,
    'a_x'  : a_Tx,
    'a_xy' : a_Txy,
    'a_xz' : a_Txz,
    'a_y'  : a_Ty,
    'a_yz' : a_Tyz,
    'a_z'  : a_Tz,
    'b_x'  : b_Tx,
    'b_xy' : b_Txy,
    'b_xz' : b_Txz,
    'b_y'  : b_Ty,
    'b_yz' : b_Tyz,
    'b_z'  : b_Tz,
    'c_x'  : c_Tx,
    'c_xy' : c_Txy,
    'c_xz' : c_Txz,
    'c_y'  : c_Ty,
    'c_yz' : c_Tyz,
    'c_z'  : c_Tz,
    'd_xy' : d_Txy,
    'd_xz' : d_Txz,
    'd_yz' : d_Tyz,
    'e_xy' : e_Txy,
    'e_xz' : e_Txz,
    'e_yz' : e_Tyz,
    'f_0'  : f_T0,
    'f_x'  : f_Tx,
    'f_xy' : f_Txy,
    'f_xz' : f_Txz,
    'f_y'  : f_Ty,
    'f_yz' : f_Tyz,
    'f_z'  : f_Tz,
    'g_0'  : g_T0,
    'g_x'  : g_Tx,
    'g_xy' : g_Txy,
    'g_xz' : g_Txz,
    'g_y'  : g_Ty,
    'g_yz' : g_Tyz,
    'g_z'  : g_Tz,
}

# Choose where and when we compute the solution
# Note that we provide the (2*pi/L)-like terms expected by soln.py
params_xyzt = {
    'x': Rational(5, 10),
    'y': Rational(6, 10),
    'z': Rational(8, 10),
    't': Rational(9, 10),
    'Lx': Lx,
    'Ly': Ly,
    'Lz': Lz,
    'twopi_invLx': (2*pi / Lx),
    'twopi_invLy': (2*pi / Ly),
    'twopi_invLz': (2*pi / Lz)
}
params_rho.update(params_xyzt)
params_u.update(params_xyzt)
params_v.update(params_xyzt)
params_w.update(params_xyzt)
params_T.update(params_xyzt)

# Compute the analytic solution and derivatives for rho, u, v, w, and T
# These are the necessary inputs to forcing.py which we will source next
rho    = phi   .subs(params_rho)
rho_t  = phi_t .subs(params_rho)
rho_x  = phi_x .subs(params_rho)
rho_xx = phi_xx.subs(params_rho)
rho_xy = phi_xy.subs(params_rho)
rho_xz = phi_xz.subs(params_rho)
rho_y  = phi_y .subs(params_rho)
rho_yy = phi_yy.subs(params_rho)
rho_yz = phi_yz.subs(params_rho)
rho_z  = phi_z .subs(params_rho)
rho_zz = phi_zz.subs(params_rho)

u    = phi   .subs(params_u)
u_t  = phi_t .subs(params_u)
u_x  = phi_x .subs(params_u)
u_xx = phi_xx.subs(params_u)
u_xy = phi_xy.subs(params_u)
u_xz = phi_xz.subs(params_u)
u_y  = phi_y .subs(params_u)
u_yy = phi_yy.subs(params_u)
u_yz = phi_yz.subs(params_u)
u_z  = phi_z .subs(params_u)
u_zz = phi_zz.subs(params_u)

v    = phi   .subs(params_v)
v_t  = phi_t .subs(params_v)
v_x  = phi_x .subs(params_v)
v_xx = phi_xx.subs(params_v)
v_xy = phi_xy.subs(params_v)
v_xz = phi_xz.subs(params_v)
v_y  = phi_y .subs(params_v)
v_yy = phi_yy.subs(params_v)
v_yz = phi_yz.subs(params_v)
v_z  = phi_z .subs(params_v)
v_zz = phi_zz.subs(params_v)

w    = phi   .subs(params_w)
w_t  = phi_t .subs(params_w)
w_x  = phi_x .subs(params_w)
w_xx = phi_xx.subs(params_w)
w_xy = phi_xy.subs(params_w)
w_xz = phi_xz.subs(params_w)
w_y  = phi_y .subs(params_w)
w_yy = phi_yy.subs(params_w)
w_yz = phi_yz.subs(params_w)
w_z  = phi_z .subs(params_w)
w_zz = phi_zz.subs(params_w)

T    = phi   .subs(params_T)
T_t  = phi_t .subs(params_T)
T_x  = phi_x .subs(params_T)
T_xx = phi_xx.subs(params_T)
T_xy = phi_xy.subs(params_T)
T_xz = phi_xz.subs(params_T)
T_y  = phi_y .subs(params_T)
T_yy = phi_yy.subs(params_T)
T_yz = phi_yz.subs(params_T)
T_z  = phi_z .subs(params_T)
T_zz = phi_zz.subs(params_T)

# Compute the forcing using the same equations appearing in writeup.tex
# This also includes term-by-term intermediate values.
execfile("forcing.py")

# Prepare a list of the quantities that we should output.
# This makes assumptions about the names of temporaries in forcing.py,
# but that should not be a big problem.
qoi = """
    rho rho_t rho_x rho_xx rho_xy rho_xz rho_y rho_yy rho_yz rho_z rho_zz
    u   u_t   u_x   u_xx   u_xy   u_xz   u_y   u_yy   u_yz   u_z   u_zz
    v   v_t   v_x   v_xx   v_xy   v_xz   v_y   v_yy   v_yz   v_z   v_zz
    w   w_t   w_x   w_xx   w_xy   w_xz   w_y   w_yy   w_yz   w_z   w_zz
    T   T_t   T_x   T_xx   T_xy   T_xz   T_y   T_yy   T_yz   T_z   T_zz

    e e_x e_y e_z e_t
    p p_x p_y p_z mu
    mu_x mu_y mu_z
    lambda_ lambda_x lambda_y lambda_z
    qx qy qz
    qx_x qy_y qz_z

    rhou rhov rhow rhoe
    rhou_x rhov_y rhow_z
    rhou_t rhov_t rhow_t rhoe_t

    rhouu_x rhouv_y rhouw_z rhouv_x
    rhovv_y rhovw_z rhouw_x rhovw_y
    rhoww_z rhoue_x rhove_y rhowe_z

    tauxx tauyy tauzz tauxy tauxz tauyz
    tauxx_x tauyy_y tauzz_z tauxy_x tauxy_y tauxz_x tauxz_z tauyz_y tauyz_z

    pu_x pv_y pw_z
    utauxx_x vtauxy_x wtauxz_x
    utauxy_y vtauyy_y wtauyz_y
    utauxz_z vtauyz_z wtauzz_z

    Q_rho Q_rhou Q_rhov Q_rhow Q_rhoe
""".split()

# Output C++ namespace containing expected results in as long double constants
def dump(name, value):
    print "    const long double", name, "=",
    if isinstance(value, int):
        print value,
        sys.stdout.write(".0")
    else:
        print value.evalf(PREC),
    sys.stdout.write("L;\n");

print "namespace nsctpl_rholut {"
print "namespace test {"
dump('x', params_xyzt['x'])
dump('y', params_xyzt['y'])
dump('z', params_xyzt['z'])
dump('t', params_xyzt['t'])
for quantity in params:
    dump(quantity, globals()[quantity])
for quantity in qoi:
    dump(quantity, globals()[quantity])
    sys.stdout.flush()
print "} // end namespace test"
print "} // end namespace nsctpl_rholut"
