#!/usr/bin/env python
# Compute manufactured solution and manufactured forcing at a particular point.
# Used as a test case to ensure C++ implementation consistent with Python, and
# to ensure floating point precision losses are acceptable.

# soln.py contains the nececssary sympy imports
execfile("soln.py")

# Work in (possibly slow) extended precision
PREC = 30;

# Test solution parameters.  These are not physically
# realizable but can be used to check the implementation.
gamma    = Real(2, PREC)
R        = Real(3, PREC)
beta     = Real(5, PREC)
mu_r     = Real(7, PREC)
T_r      = Real(11, PREC)
k_r      = Real(13, PREC)
lambda_r = Real(17, PREC)
a_rho0   = Real(19, PREC)
a_rhox   = Real(23, PREC)
a_rhoxy  = Real(29, PREC)
a_rhoxz  = Real(31, PREC)
a_rhoy   = Real(37, PREC)
a_rhoyz  = Real(41, PREC)
a_rhoz   = Real(43, PREC)
b_rhox   = Real(47, PREC)
b_rhoxy  = Real(53, PREC)
b_rhoxz  = Real(59, PREC)
b_rhoy   = Real(61, PREC)
b_rhoyz  = Real(67, PREC)
b_rhoz   = Real(71, PREC)
c_rhox   = Real(73, PREC)
c_rhoxy  = Real(79, PREC)
c_rhoxz  = Real(83, PREC)
c_rhoy   = Real(89, PREC)
c_rhoyz  = Real(97, PREC)
c_rhoz   = Real(101, PREC)
d_rhoxy  = Real(103, PREC)
d_rhoxz  = Real(107, PREC)
d_rhoyz  = Real(109, PREC)
e_rhoxy  = Real(113, PREC)
e_rhoxz  = Real(127, PREC)
e_rhoyz  = Real(131, PREC)
f_rho0   = Real(137, PREC)
f_rhox   = Real(139, PREC)
f_rhoxy  = Real(149, PREC)
f_rhoxz  = Real(151, PREC)
f_rhoy   = Real(157, PREC)
f_rhoyz  = Real(163, PREC)
f_rhoz   = Real(167, PREC)
g_rho0   = Real(173, PREC)
g_rhox   = Real(179, PREC)
g_rhoxy  = Real(181, PREC)
g_rhoxz  = Real(191, PREC)
g_rhoy   = Real(193, PREC)
g_rhoyz  = Real(197, PREC)
g_rhoz   = Real(199, PREC)
a_u0     = Real(211, PREC)
a_ux     = Real(223, PREC)
a_uxy    = Real(227, PREC)
a_uxz    = Real(229, PREC)
a_uy     = Real(233, PREC)
a_uyz    = Real(239, PREC)
a_uz     = Real(241, PREC)
b_ux     = Real(251, PREC)
b_uxy    = Real(257, PREC)
b_uxz    = Real(263, PREC)
b_uy     = Real(269, PREC)
b_uyz    = Real(271, PREC)
b_uz     = Real(277, PREC)
c_ux     = Real(281, PREC)
c_uxy    = Real(283, PREC)
c_uxz    = Real(293, PREC)
c_uy     = Real(307, PREC)
c_uyz    = Real(311, PREC)
c_uz     = Real(313, PREC)
d_uxy    = Real(317, PREC)
d_uxz    = Real(331, PREC)
d_uyz    = Real(337, PREC)
e_uxy    = Real(347, PREC)
e_uxz    = Real(349, PREC)
e_uyz    = Real(353, PREC)
f_u0     = Real(359, PREC)
f_ux     = Real(367, PREC)
f_uxy    = Real(373, PREC)
f_uxz    = Real(379, PREC)
f_uy     = Real(383, PREC)
f_uyz    = Real(389, PREC)
f_uz     = Real(397, PREC)
g_u0     = Real(401, PREC)
g_ux     = Real(409, PREC)
g_uxy    = Real(419, PREC)
g_uxz    = Real(421, PREC)
g_uy     = Real(431, PREC)
g_uyz    = Real(433, PREC)
g_uz     = Real(439, PREC)
a_v0     = Real(443, PREC)
a_vx     = Real(449, PREC)
a_vxy    = Real(457, PREC)
a_vxz    = Real(461, PREC)
a_vy     = Real(463, PREC)
a_vyz    = Real(467, PREC)
a_vz     = Real(479, PREC)
b_vx     = Real(487, PREC)
b_vxy    = Real(491, PREC)
b_vxz    = Real(499, PREC)
b_vy     = Real(503, PREC)
b_vyz    = Real(509, PREC)
b_vz     = Real(521, PREC)
c_vx     = Real(523, PREC)
c_vxy    = Real(541, PREC)
c_vxz    = Real(547, PREC)
c_vy     = Real(557, PREC)
c_vyz    = Real(563, PREC)
c_vz     = Real(569, PREC)
d_vxy    = Real(571, PREC)
d_vxz    = Real(577, PREC)
d_vyz    = Real(587, PREC)
e_vxy    = Real(593, PREC)
e_vxz    = Real(599, PREC)
e_vyz    = Real(601, PREC)
f_v0     = Real(607, PREC)
f_vx     = Real(613, PREC)
f_vxy    = Real(617, PREC)
f_vxz    = Real(619, PREC)
f_vy     = Real(631, PREC)
f_vyz    = Real(641, PREC)
f_vz     = Real(643, PREC)
g_v0     = Real(647, PREC)
g_vx     = Real(653, PREC)
g_vxy    = Real(659, PREC)
g_vxz    = Real(661, PREC)
g_vy     = Real(673, PREC)
g_vyz    = Real(677, PREC)
g_vz     = Real(683, PREC)
a_w0     = Real(691, PREC)
a_wx     = Real(701, PREC)
a_wxy    = Real(709, PREC)
a_wxz    = Real(719, PREC)
a_wy     = Real(727, PREC)
a_wyz    = Real(733, PREC)
a_wz     = Real(739, PREC)
b_wx     = Real(743, PREC)
b_wxy    = Real(751, PREC)
b_wxz    = Real(757, PREC)
b_wy     = Real(761, PREC)
b_wyz    = Real(769, PREC)
b_wz     = Real(773, PREC)
c_wx     = Real(787, PREC)
c_wxy    = Real(797, PREC)
c_wxz    = Real(809, PREC)
c_wy     = Real(811, PREC)
c_wyz    = Real(821, PREC)
c_wz     = Real(823, PREC)
d_wxy    = Real(827, PREC)
d_wxz    = Real(829, PREC)
d_wyz    = Real(839, PREC)
e_wxy    = Real(853, PREC)
e_wxz    = Real(857, PREC)
e_wyz    = Real(859, PREC)
f_w0     = Real(863, PREC)
f_wx     = Real(877, PREC)
f_wxy    = Real(881, PREC)
f_wxz    = Real(883, PREC)
f_wy     = Real(887, PREC)
f_wyz    = Real(907, PREC)
f_wz     = Real(911, PREC)
g_w0     = Real(919, PREC)
g_wx     = Real(929, PREC)
g_wxy    = Real(937, PREC)
g_wxz    = Real(941, PREC)
g_wy     = Real(947, PREC)
g_wyz    = Real(953, PREC)
g_wz     = Real(967, PREC)
a_T0     = Real(971, PREC)
a_Tx     = Real(977, PREC)
a_Txy    = Real(983, PREC)
a_Txz    = Real(991, PREC)
a_Ty     = Real(997, PREC)
a_Tyz    = Real(1009, PREC)
a_Tz     = Real(1013, PREC)
b_Tx     = Real(1019, PREC)
b_Txy    = Real(1021, PREC)
b_Txz    = Real(1031, PREC)
b_Ty     = Real(1033, PREC)
b_Tyz    = Real(1039, PREC)
b_Tz     = Real(1049, PREC)
c_Tx     = Real(1051, PREC)
c_Txy    = Real(1061, PREC)
c_Txz    = Real(1063, PREC)
c_Ty     = Real(1069, PREC)
c_Tyz    = Real(1087, PREC)
c_Tz     = Real(1091, PREC)
d_Txy    = Real(1093, PREC)
d_Txz    = Real(1097, PREC)
d_Tyz    = Real(1103, PREC)
e_Txy    = Real(1109, PREC)
e_Txz    = Real(1117, PREC)
e_Tyz    = Real(1123, PREC)
f_T0     = Real(1129, PREC)
f_Tx     = Real(1151, PREC)
f_Txy    = Real(1153, PREC)
f_Txz    = Real(1163, PREC)
f_Ty     = Real(1171, PREC)
f_Tyz    = Real(1181, PREC)
f_Tz     = Real(1187, PREC)
g_T0     = Real(1193, PREC)
g_Tx     = Real(1201, PREC)
g_Txy    = Real(1213, PREC)
g_Txz    = Real(1217, PREC)
g_Ty     = Real(1223, PREC)
g_Tyz    = Real(1229, PREC)
g_Tz     = Real(1231, PREC)
Lx       = Real(1237, PREC)
Ly       = Real(1249, PREC)
Lz       = Real(1259, PREC)

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
    'z': Rational(7, 10),
    't': Rational(8, 10),
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

    gamma  R  beta  mu_r  T_r  k_r  lambda_r

    e e_x e_y e_z e_t
    p p_x p_y p_z mu
    mu_x mu_y mu_z
    lambda_ lambda_x lambda_y lambda_z
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

print "namespace nsctpl {"
print "namespace test {"
dump('x', params_xyzt['x'])
dump('y', params_xyzt['y'])
dump('z', params_xyzt['z'])
dump('t', params_xyzt['t'])
dump('Lx', params_xyzt['Lx'])
dump('Ly', params_xyzt['Ly'])
dump('Lz', params_xyzt['Lz'])
for quantity in qoi:
    dump(quantity, globals()[quantity])
    sys.stdout.flush()
print "} // end namespace test"
print "} // end namespace nsctpl"
